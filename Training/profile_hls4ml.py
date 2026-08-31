#!/usr/bin/env python3
import argparse
import importlib.util
import json
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1" # run on cpu for now for some weird version mismatch reason ;(

import numpy as np
import keras
from keras import layers
import hls4ml
from hls4ml.model.profiling import numerical
import matplotlib.pyplot as plt
import torch
from torch.utils.data import DataLoader


NMAX = 24
NFEAT = 18
HIDDEN = 16
PARALLELIZATION_FACTOR = 24


def load_training_module(path):
    spec = importlib.util.spec_from_file_location("mlpuppi_train", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pytorch_dense_kwargs(fan_in):
    limit = 1.0 / np.sqrt(fan_in)
    return {
        "kernel_initializer": keras.initializers.RandomUniform(minval=-limit, maxval=limit),
        "bias_initializer": keras.initializers.RandomUniform(minval=-limit, maxval=limit),
    }


def build_model():
    x = keras.Input(shape=(NMAX, NFEAT), name="x")
    mask = keras.Input(shape=(NMAX,), dtype="bool", name="mask")

    h = layers.Dense(HIDDEN, name="phi_dense1", **pytorch_dense_kwargs(NFEAT))(x)
    h = layers.BatchNormalization(momentum=0.9, epsilon=1e-5, name="phi_bn1")(h, mask=mask)
    h = layers.ReLU(name="phi_relu1")(h)

    h = layers.Dense(HIDDEN, name="phi_dense2", **pytorch_dense_kwargs(HIDDEN))(h)
    h = layers.BatchNormalization(momentum=0.9, epsilon=1e-5, name="phi_bn2")(h, mask=mask)
    h = layers.ReLU(name="phi_relu2")(h)

    h = layers.GlobalAveragePooling1D(name="masked_mean")(h, mask=mask)

    h = layers.Dense(HIDDEN, name="rho_dense1", **pytorch_dense_kwargs(HIDDEN))(h)
    h = layers.ReLU(name="rho_relu1")(h)
    h = layers.Dense(1, name="rho_output", **pytorch_dense_kwargs(HIDDEN))(h)
    y = layers.Activation("softplus", name="puppi_weight")(h)

    return keras.Model(inputs={"x": x, "mask": mask}, outputs=y, name="PuppiDeepSetMeanOnly")


def build_hls_keras_model():
    x = keras.Input(shape=(NMAX, NFEAT), name="x")
    mask = keras.Input(shape=(NMAX,), name="mask")
    scale = keras.Input(shape=(1,), name="scale")

    h = layers.Dense(HIDDEN, name="phi_dense1", **pytorch_dense_kwargs(NFEAT))(x)
    h = layers.BatchNormalization(momentum=0.9, epsilon=1e-5, name="phi_bn1")(h)
    h = layers.ReLU(name="phi_relu1")(h)

    h = layers.Dense(HIDDEN, name="phi_dense2", **pytorch_dense_kwargs(HIDDEN))(h)
    h = layers.BatchNormalization(momentum=0.9, epsilon=1e-5, name="phi_bn2")(h)
    h = layers.ReLU(name="phi_relu2")(h)

    mask_f = layers.Reshape((NMAX, 1), name="mask_reshape")(mask)
    h = layers.Multiply(name="apply_mask")([h, mask_f])

    h = layers.GlobalAveragePooling1D(name="global_average_pool")(h)
    h = layers.Multiply(name="masked_mean_scale")([h, scale])

    h = layers.Dense(HIDDEN, name="rho_dense1", **pytorch_dense_kwargs(HIDDEN))(h)
    h = layers.ReLU(name="rho_relu1")(h)
    h = layers.Dense(1, name="rho_output", **pytorch_dense_kwargs(HIDDEN))(h)
    y = layers.Activation("softplus", name="puppi_weight")(h)

    return keras.Model(inputs=[x, mask, scale], outputs=y, name="PuppiDeepSetHLS")


def plot_abs_ranges(names, arrays, title, outfile):
    p99 = []
    p999 = []
    xmax = []

    for arr in arrays:
        v = np.asarray(arr).reshape(-1)
        v = v[np.isfinite(v)]
        av = np.abs(v)
        p99.append(np.percentile(av, 99))
        p999.append(np.percentile(av, 99.9))
        xmax.append(np.max(av))

    y = np.arange(len(names))

    plt.figure(figsize=(10, max(6, 0.45 * len(names))))
    plt.scatter(p99, y, label="|x| p99")
    plt.scatter(p999, y, label="|x| p99.9")
    plt.scatter(xmax, y, label="|x| max")

    plt.yticks(y, names)
    plt.xscale("log")
    plt.xlabel("|x|")
    plt.title(title)
    plt.legend()
    plt.grid(axis="x", alpha=0.25)
    plt.tight_layout()
    plt.savefig(outfile, dpi=160)
    plt.close()


def save_fig(fig, path):
    if fig is not None:
        fig.savefig(path, bbox_inches="tight", dpi=160)
        plt.close(fig)


def print_ranges(name, x):
    x = np.asarray(x).reshape(-1)
    x = x[np.isfinite(x)]
    ax = np.abs(x)
    print(f"{name:24s} min={x.min(): .6g} max={x.max(): .6g} |x|p99={np.percentile(ax,99): .6g} |x|p99.9={np.percentile(ax,99.9): .6g} |x|max={ax.max(): .6g}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_script", required=True, help="Exact H16/NMAX24 Keras training script")
    parser.add_argument("--run_dir", required=True, help="Run directory containing model_best.weights.h5 and config.json")
    parser.add_argument("--data", required=True, help="Same dataset glob used for training")
    parser.add_argument("--batch_size", type=int, default=1024)
    parser.add_argument("--num_workers", type=int, default=4)
    parser.add_argument("--max_events", type=int, default=5000)
    parser.add_argument("--max_files", type=int, default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--plot", choices=["boxplot", "histogram"], default="boxplot")
    parser.add_argument("--output_dir", default=None)
    args = parser.parse_args()

    with open(os.path.join(args.run_dir, "config.json")) as f:
        run_cfg = json.load(f)

    nmax = int(run_cfg.get("nmax", NMAX))
    hidden_dim = int(run_cfg.get("hidden_dim", HIDDEN))
    if nmax != NMAX or hidden_dim != HIDDEN:
        raise RuntimeError(f"Expected frozen NMAX={NMAX}, HIDDEN={HIDDEN}; run config has NMAX={nmax}, HIDDEN={hidden_dim}")

    seed = int(run_cfg.get("seed", 42) if args.seed is None else args.seed)
    max_files = int(run_cfg.get("max_files", 10) if args.max_files is None else args.max_files)
    outdir = args.output_dir or os.path.join(args.run_dir, "hls4ml_profile")
    os.makedirs(outdir, exist_ok=True)

    np.random.seed(seed)
    torch.manual_seed(seed)

    train = load_training_module(args.train_script)
    train.NMAX = NMAX

    print("Loading dataset...")
    dataset = train.BestPuppiDataset(args.data, max_files=max_files)
    feature_idx = train.build_feature_indices(dataset.feature_names, train.INPUT_FEATURES)
    _, _, test_idx = train.make_splits(len(dataset), seed=seed)
    if args.max_events > 0:
        test_idx = test_idx[:args.max_events]

    loader = DataLoader(torch.utils.data.Subset(dataset, test_idx), batch_size=args.batch_size, shuffle=False, num_workers=args.num_workers, pin_memory=False)

    xs, masks = [], []
    n_seen = 0
    for data in loader:
        x, mask, _ = train.prepare_batch(data, feature_idx)
        xs.append(x.numpy().astype(np.float32, copy=False))
        masks.append(mask.numpy().astype(np.float32, copy=False))
        n_seen += x.shape[0]
        if args.max_events > 0 and n_seen >= args.max_events:
            break

    X = np.concatenate(xs, axis=0)
    M = np.concatenate(masks, axis=0)
    if args.max_events > 0:
        X = X[:args.max_events]
        M = M[:args.max_events]

    print("\n=== INPUT FEATURE RANGES ===")
    X_valid = X[M.astype(bool)]

    input_arrays = [X_valid[:, i] for i in range(X_valid.shape[1])]

    plot_abs_ranges(
        train.INPUT_FEATURES,
        input_arrays,
        "Input feature dynamic ranges",
        os.path.join(outdir, "input_feature_ranges.png"),
    )

    for i, name in enumerate(train.INPUT_FEATURES):
        print_ranges(name, X_valid[:, i])

    nvalid = np.sum(M, axis=1, keepdims=True)
    if np.any(nvalid <= 0):
        raise RuntimeError("Found event with zero valid candidates")
    S = (NMAX / nvalid).astype(np.float32)

    print(f"Profiling events: {len(X)}")
    print("x shape:", X.shape)
    print("mask shape:", M.shape)
    print("scale shape:", S.shape)

    # Exact trained standard model
    model = build_model()
    model.load_weights(os.path.join(args.run_dir, "model_best.weights.h5"))

    activation_names = ["phi_dense1", "phi_bn1", "phi_relu1", "phi_dense2", "phi_bn2", "phi_relu2", "masked_mean", "rho_dense1", "rho_relu1", "rho_output", "puppi_weight"]

    probe = keras.Model(inputs=model.inputs, outputs=[model.get_layer(name).output for name in activation_names])
    outputs = probe({"x": X, "mask": M.astype(bool)}, training=False)

    print("\n=== ACTIVATION RANGES ===")
    activation_arrays = []

    for name, out in zip(activation_names, outputs):
        arr = out.numpy()

        # Remove padded candidates for candidate-wise layers
        if arr.ndim == 3:
            arr = arr[M.astype(bool)]
        
        activation_arrays.append(arr)
        print_ranges(name, arr)

    plot_abs_ranges(
        activation_names,
        activation_arrays,
        "Activation dynamic ranges",
        os.path.join(outdir, "activation_ranges.png"),
    )

    # HLS-facing Keras model used for conversion; same trained weighted layers
    model_hls_keras = build_hls_keras_model()
    model_hls_keras.set_weights(model.get_weights())

    # Check exact masked-mean rewrite numerically before profiling
    y_ref = model({"x": X, "mask": M.astype(bool)}, training=False).numpy().reshape(-1)
    y_hls_keras = model_hls_keras([X, M, S], training=False).numpy().reshape(-1)
    delta = np.abs(y_ref - y_hls_keras)
    print(f"Keras rewrite check: max |delta| = {delta.max():.8g}, mean |delta| = {delta.mean():.8g}")

    config = hls4ml.utils.config_from_keras_model(model_hls_keras, granularity="name", backend="Vitis")
    config["LayerName"]["phi_dense1"]["ParallelizationFactor"] = PARALLELIZATION_FACTOR
    config["LayerName"]["phi_dense2"]["ParallelizationFactor"] = PARALLELIZATION_FACTOR

    hls_model = hls4ml.converters.convert_from_keras_model(
        model_hls_keras,
        hls_config=config,
        backend="Vitis",
        output_dir=os.path.join(outdir, "hls_project"),
        part="xcvu13p-flga2577-2-e",
        clock_period=2.778,
        io_type="io_parallel",
    )

    print("\nResolved HLS output precisions:")
    for layer in hls_model.get_layers():
        try:
            print(f"{layer.name:24s} {layer.get_output_variable().type.precision}")
        except Exception:
            pass

    # ------------------------------------------------------------------
    # hls4ml numerical profiling
    #
    # Current hls4ml numerical() internally converts X to one contiguous
    # ndarray for ModelGraph activation tracing. Our HLS-facing model has
    # three inputs with different shapes (x, mask, scale), so full
    # keras+hls activation profiling in ONE numerical() call is not
    # supported cleanly.
    #
    # We therefore do the two supported pieces separately:
    #   A) Keras-side weights + activations using hls4ml numerical()
    #   B) HLS ModelGraph weights + resolved fixed-point ranges
    #
    # This gives the float activation distributions needed for QKeras
    # precision choice, and the hls4ml fixed-point weight ranges.
    # ------------------------------------------------------------------

    print("\n=== hls4ml numerical: Keras weights + activations ===")
    wp_keras, _, ap_keras, _ = numerical(model=model_hls_keras, X=[X, M, S], plot=args.plot)
    save_fig(wp_keras, os.path.join(outdir, f"weights_keras_{args.plot}.png"))
    save_fig(ap_keras, os.path.join(outdir, f"activations_keras_{args.plot}.png"))

    print("\n=== hls4ml numerical: HLS weights + precision ranges ===")
    wp_before, wp_after, _, _ = numerical(hls_model=hls_model, plot=args.plot)
    save_fig(wp_before, os.path.join(outdir, f"weights_hls_before_{args.plot}.png"))
    save_fig(wp_after, os.path.join(outdir, f"weights_hls_after_{args.plot}.png"))

    print("\nSaved profiling plots:")
    print(os.path.join(outdir, f"weights_keras_{args.plot}.png"))
    print(os.path.join(outdir, f"activations_keras_{args.plot}.png"))
    print(os.path.join(outdir, f"weights_hls_before_{args.plot}.png"))
    print(os.path.join(outdir, f"weights_hls_after_{args.plot}.png"))
    print("\nNo hls_model.build() was run.")


if __name__ == "__main__":
    main()

