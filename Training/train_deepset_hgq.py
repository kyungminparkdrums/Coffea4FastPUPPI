#!/usr/bin/env python3

import argparse
import glob
import json
import math
import os
import random
import time
from datetime import datetime

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
import tensorflow as tf
import keras
from keras import layers
from hgq.layers import QDense, QBatchNormalization, QGlobalAveragePooling1D
from hgq.config import LayerConfigScope, QuantizerConfigScope
from tqdm import tqdm


NMAX = 24

INPUT_FEATURES = [
    "log_pt",
    "slog_px",
    "slog_py",
    "abs_eta",
    "charge",
    "nnvtx",
    "dxy",
    "z0",
    "puppiWeight",
    "pdg_abs_211",
    "pdg_abs_130",
    "pdg_abs_22",
    "pdg_abs_11",
    "pdg_abs_13",
    "deta",
    "dphi",
    "dr",
    "is_center",
]


class BestPuppiDataset(Dataset):
    def __init__(self, pattern: str, max_files=10):
        files = sorted(glob.glob(pattern))
        random.shuffle(files)
        self.files = files[:max_files]

        print("len(self.files)", len(self.files))
        print("Preloading dataset...")

        self.data = []
        self.metadata = None
        self.feature_names = None

        for f in self.files:
            t0 = time.time()
            payload = torch.load(f, map_location="cpu", weights_only=False)

            if isinstance(payload, dict) and "graphs" in payload:
                graphs = payload["graphs"]
                metadata = payload.get("metadata", {})
                feature_names = metadata.get("feature_names", None)

                if self.metadata is None:
                    self.metadata = metadata
                    self.feature_names = feature_names
                elif feature_names != self.feature_names:
                    raise RuntimeError(
                        f"Feature-name mismatch in {f}\n"
                        f"existing: {self.feature_names}\n"
                        f"new     : {feature_names}"
                    )
            else:
                graphs = payload
                metadata = {}
                feature_names = None

            self.data.extend(graphs)
            print(f"Loaded {f}: {len(graphs)} graphs in {time.time() - t0:.1f} s")

        print(f"Loaded {len(self.data)} graphs from {len(self.files)} files")

        if self.feature_names is None:
            in_dim = self.data[0].x.shape[1]
            self.feature_names = default_feature_names(in_dim)
            self.metadata = {
                "feature_names": self.feature_names,
                "note": "Inferred feature names from old-format dataset.",
            }

        print("Dataset feature names:")
        for i, name in enumerate(self.feature_names):
            print(f"  {i:2d}: {name}")

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        data = self.data[idx]
        x = data.x

        is_center_idx = self.feature_names.index("is_center")
        center_mask = x[:, is_center_idx] > 0.5

        if center_mask.sum().item() != 1:
            raise RuntimeError(f"Expected exactly one center, got {center_mask.sum().item()}")

        center_idx_local = torch.where(center_mask)[0].item()
        other_idx = torch.arange(x.size(0))
        other_idx = other_idx[other_idx != center_idx_local]
        order = torch.cat([torch.tensor([center_idx_local], dtype=torch.long), other_idx])
        x = x[order]

        n = x.size(0)

        if n > NMAX:
            x = x[:NMAX]
            n = NMAX

        x_fixed = torch.zeros(NMAX, x.size(1), dtype=x.dtype)
        mask = torch.zeros(NMAX, dtype=x.dtype)
        x_fixed[:n] = x
        mask[:n] = 1.0

        return {
            "x": x_fixed,
            "mask": mask,
            "y": data.y,
            "event_idx": data.event_idx,
            "center_idx": data.center_idx,
        }


def default_feature_names(in_dim):
    if in_dim == 15:
        return [
            "log_pt", "slog_px", "slog_py", "charge", "nnvtx", "puppiWeight",
            "pdg_abs_211", "pdg_abs_130", "pdg_abs_22", "pdg_abs_11", "pdg_abs_13",
            "deta", "dphi", "dr", "is_center",
        ]

    if in_dim == 16:
        return [
            "log_pt", "slog_px", "slog_py", "abs_eta", "charge", "nnvtx", "puppiWeight",
            "pdg_abs_211", "pdg_abs_130", "pdg_abs_22", "pdg_abs_11", "pdg_abs_13",
            "deta", "dphi", "dr", "is_center",
        ]

    return [f"feature_{i}" for i in range(in_dim)]


def build_feature_indices(dataset_feature_names, input_features):
    feature_to_idx = {name: i for i, name in enumerate(dataset_feature_names)}
    missing = [f for f in input_features if f not in feature_to_idx]

    if missing:
        raise RuntimeError(
            "Requested input features are missing from dataset metadata:\n"
            f"missing = {missing}\n"
            f"available = {dataset_feature_names}"
        )

    return [feature_to_idx[f] for f in input_features]


def pytorch_dense_kwargs(fan_in):
    bound = 1.0 / math.sqrt(fan_in)
    return {
        "kernel_initializer": keras.initializers.RandomUniform(minval=-bound, maxval=bound),
        "bias_initializer": keras.initializers.RandomUniform(minval=-bound, maxval=bound),
    }


def build_model(nmax=NMAX, in_dim=len(INPUT_FEATURES), hidden_dim=16, beta0=1e-5):
    # HGQ2 tutorial/default-style scopes:
    # - kbi/SAT_SYM for parameters
    # - kif/WRAP for datalane values
    # - EBOP regularization enabled with beta0
    with (
        QuantizerConfigScope(place="all", default_q_type="kbi", overflow_mode="SAT_SYM"),
        QuantizerConfigScope(place="datalane", default_q_type="kif", overflow_mode="WRAP"),
        LayerConfigScope(enable_ebops=True, beta0=beta0),
    ):
        x_input = keras.Input(shape=(nmax, in_dim), name="x")
        mask_input = keras.Input(shape=(nmax,), dtype="bool", name="mask")

        h = QDense(hidden_dim, name="phi_dense1", **pytorch_dense_kwargs(in_dim))(x_input)
        h = QBatchNormalization(momentum=0.9, epsilon=1e-5, name="phi_bn1")(h, mask=mask_input)
        h = layers.ReLU(name="phi_relu1")(h)

        h = QDense(hidden_dim, name="phi_dense2", **pytorch_dense_kwargs(hidden_dim))(h)
        h = QBatchNormalization(momentum=0.9, epsilon=1e-5, name="phi_bn2")(h, mask=mask_input)
        h = layers.ReLU(name="phi_relu2")(h)

        #h_mean = QGlobalAveragePooling1D(name="masked_mean")(h, mask=mask_input)
        h_mean = layers.GlobalAveragePooling1D(name="masked_mean")(h, mask=mask_input)

        out = QDense(hidden_dim, name="rho_dense1", **pytorch_dense_kwargs(hidden_dim))(h_mean)
        out = layers.ReLU(name="rho_relu1")(out)
        raw = QDense(1, name="rho_output", **pytorch_dense_kwargs(hidden_dim))(out)
        pred = layers.Activation("softplus", name="puppi_weight")(raw)

    return keras.Model(
        inputs={"x": x_input, "mask": mask_input},
        outputs=pred,
        name="PuppiDeepSetMeanOnlyHGQ",
    )


def puppi_loss(
    pred,
    target,
    loss_type="weighted_huber",
    alpha=1.0,
    eps=1e-3,
    zero_thr=0.05,
    pos_weight=9.0,
    pos_weight_threshold=0.05,
):
    pred = tf.reshape(pred, [-1])
    target = tf.reshape(target, [-1])

    error = pred - target
    abs_error = tf.abs(error)

    if loss_type in ("huber", "weighted_huber"):
        loss_abs = tf.where(abs_error < 1.0, 0.5 * tf.square(error), abs_error - 0.5)

        if loss_type == "weighted_huber":
            weight = tf.where(
                target > pos_weight_threshold,
                tf.cast(pos_weight, pred.dtype),
                tf.cast(1.0, pred.dtype),
            )
            loss_abs = loss_abs * weight

    elif loss_type == "mse":
        loss_abs = tf.square(error)

    elif loss_type == "mae":
        loss_abs = abs_error

    else:
        raise ValueError(f"Unknown loss_type: {loss_type}")

    loss = tf.reduce_mean(loss_abs)

    if alpha > 0.0:
        pos_mask = target > zero_thr
        pred_pos = tf.boolean_mask(pred, pos_mask)
        target_pos = tf.boolean_mask(target, pos_mask)

        def logratio_term():
            p = tf.maximum(pred_pos, tf.cast(eps, pred.dtype))
            t = tf.maximum(target_pos, tf.cast(eps, target.dtype))
            log_ratio = tf.math.log(p) - tf.math.log(t)
            return tf.reduce_mean(tf.square(log_ratio))

        loss_logratio = tf.cond(
            tf.size(pred_pos) > 0,
            logratio_term,
            lambda: tf.cast(0.0, pred.dtype),
        )
        loss = loss + alpha * loss_logratio

    return loss


def compute_diagnostic_metrics(pred, target, eps=1e-3, zero_thr=0.05, pred_thr=0.05):
    pred = np.asarray(pred).reshape(-1)
    target = np.asarray(target).reshape(-1)
    err = pred - target

    metrics = {}
    metrics["mse"] = float(np.mean(err ** 2))
    metrics["mae"] = float(np.mean(np.abs(err)))

    abs_err = np.abs(err)
    metrics["huber"] = float(
        np.mean(np.where(abs_err < 1.0, 0.5 * err ** 2, abs_err - 0.5))
    )

    zero_mask = target <= zero_thr
    pos_mask = target > zero_thr

    metrics["target_zero_frac"] = float(np.mean(zero_mask))
    metrics["target_pos_frac"] = float(np.mean(pos_mask))
    metrics["pred_mean"] = float(np.mean(pred))
    metrics["target_mean"] = float(np.mean(target))
    metrics["global_bias"] = float(np.mean(pred - target))

    if np.any(zero_mask):
        metrics["zero_mean_pred"] = float(np.mean(pred[zero_mask]))
        metrics["zero_median_pred"] = float(np.median(pred[zero_mask]))
        metrics["zero_nonzero_pred_frac"] = float(np.mean(pred[zero_mask] > pred_thr))
    else:
        metrics["zero_mean_pred"] = np.nan
        metrics["zero_median_pred"] = np.nan
        metrics["zero_nonzero_pred_frac"] = np.nan

    if np.any(pos_mask):
        log_ratio = np.log(pred[pos_mask] + eps) - np.log(target[pos_mask] + eps)
        ratio = (pred[pos_mask] + eps) / (target[pos_mask] + eps)

        metrics["pos_mae"] = float(np.mean(np.abs(pred[pos_mask] - target[pos_mask])))
        metrics["pos_bias"] = float(np.mean(pred[pos_mask] - target[pos_mask]))
        metrics["pos_logratio_mean"] = float(np.mean(log_ratio))
        metrics["pos_logratio_rms"] = float(np.sqrt(np.mean(log_ratio ** 2)))
        metrics["pos_ratio_median"] = float(np.median(ratio))
        metrics["pos_under_zero_frac"] = float(np.mean(pred[pos_mask] <= pred_thr))
    else:
        metrics["pos_mae"] = np.nan
        metrics["pos_bias"] = np.nan
        metrics["pos_logratio_mean"] = np.nan
        metrics["pos_logratio_rms"] = np.nan
        metrics["pos_ratio_median"] = np.nan
        metrics["pos_under_zero_frac"] = np.nan

    return metrics


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    tf.keras.utils.set_random_seed(seed)


def make_splits(n, seed=42):
    idx = list(range(n))
    random.Random(seed).shuffle(idx)

    n_train = int(0.7 * n)
    n_val = int(0.15 * n)

    return idx[:n_train], idx[n_train:n_train + n_val], idx[n_train + n_val:]


def make_run_dir():
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = f"runs_hgq/run_{ts}"
    os.makedirs(run_dir, exist_ok=True)
    return run_dir


def prepare_batch(data, feature_idx):
    x = data["x"][:, :, feature_idx].numpy().astype(np.float32, copy=False)
    mask = data["mask"].numpy().astype(bool, copy=False)
    target = data["y"].view(-1).numpy().astype(np.float32, copy=False)

    return (
        tf.convert_to_tensor(x),
        tf.convert_to_tensor(mask),
        tf.convert_to_tensor(target),
    )


def train_epoch(model, loader, optimizer, args, feature_idx):
    total_losses = []
    physics_losses = []
    ebops_losses = []

    for data in tqdm(loader, desc="Train", leave=False):
        x, mask, target = prepare_batch(data, feature_idx)

        with tf.GradientTape() as tape:
            pred = model({"x": x, "mask": mask}, training=True)

            physics_loss = puppi_loss(
                pred,
                target,
                loss_type=args.loss_type,
                alpha=args.penalty_alpha,
                eps=args.penalty_eps,
                zero_thr=args.penalty_zero_thr,
                pos_weight=args.pos_weight,
                pos_weight_threshold=args.pos_weight_threshold,
            )

            # HGQ2 attaches the EBOP regularization terms to model.losses.
            ebops_loss = (
                tf.add_n(model.losses)
                if model.losses
                else tf.cast(0.0, physics_loss.dtype)
            )

            total_loss = physics_loss + ebops_loss

        gradients = tape.gradient(total_loss, model.trainable_variables)
        optimizer.apply_gradients(zip(gradients, model.trainable_variables))

        total_losses.append(float(total_loss.numpy()))
        physics_losses.append(float(physics_loss.numpy()))
        ebops_losses.append(float(ebops_loss.numpy()))

    return {
        "loss": float(np.mean(total_losses)),
        "physics_loss": float(np.mean(physics_losses)),
        "ebops_loss": float(np.mean(ebops_losses)),
    }


def evaluate(model, loader, args, feature_idx, input_feature_name_to_idx):
    preds, targets, losses = [], [], []
    event_idx, center_idx, center_pt, center_eta = [], [], [], []

    for data in tqdm(loader, desc="Eval", leave=False):
        x, mask, target = prepare_batch(data, feature_idx)

        pred = model({"x": x, "mask": mask}, training=False)

        # Validation/test selection remains on the original physics loss only,
        # matching the float Keras baseline.
        loss = puppi_loss(
            pred,
            target,
            loss_type=args.loss_type,
            alpha=args.penalty_alpha,
            eps=args.penalty_eps,
            zero_thr=args.penalty_zero_thr,
            pos_weight=args.pos_weight,
            pos_weight_threshold=args.pos_weight_threshold,
        )

        pred_np = pred.numpy().reshape(-1)
        target_np = target.numpy().reshape(-1)

        preds.append(pred_np)
        targets.append(target_np)
        losses.append(float(loss.numpy()))

        event_idx.append(data["event_idx"].view(-1).numpy())
        center_idx.append(data["center_idx"].view(-1).numpy())

        x_np = x.numpy()
        center_pt.append(np.exp(x_np[:, 0, input_feature_name_to_idx["log_pt"]]))

        if "abs_eta" in input_feature_name_to_idx:
            center_eta.append(x_np[:, 0, input_feature_name_to_idx["abs_eta"]])

    out = {
        "loss": float(np.mean(losses)),
        "pred": np.concatenate(preds),
        "target": np.concatenate(targets),
        "event_idx": np.concatenate(event_idx),
        "center_idx": np.concatenate(center_idx),
        "center_pt": np.concatenate(center_pt),
    }

    if center_eta:
        out["center_eta"] = np.concatenate(center_eta)

    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch_size", type=int, default=1024)
    parser.add_argument("--num_workers", type=int, default=4)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max_files", type=int, default=10)
    parser.add_argument("--hidden_dim", type=int, default=16)

    # Only HGQ-specific training knob added in this first test.
    parser.add_argument("--hgq_beta0", type=float, default=1e-5)

    parser.add_argument(
        "--loss_type",
        choices=["huber", "mse", "mae", "weighted_huber"],
        default="weighted_huber",
    )
    parser.add_argument("--pos_weight", type=float, default=9.0)
    parser.add_argument("--pos_weight_threshold", type=float, default=0.05)
    parser.add_argument("--penalty_alpha", type=float, default=1.0)
    parser.add_argument("--penalty_eps", type=float, default=1e-3)
    parser.add_argument("--penalty_zero_thr", type=float, default=0.05)
    args = parser.parse_args()

    set_seed(args.seed)
    run_dir = make_run_dir()

    print("Run dir:", run_dir)
    print("TensorFlow:", tf.__version__)
    print("Keras:", keras.__version__)
    print("GPUs:", tf.config.list_physical_devices("GPU"))

    dataset = BestPuppiDataset(args.data, max_files=args.max_files)
    dataset_feature_names = dataset.feature_names
    feature_idx = build_feature_indices(dataset_feature_names, INPUT_FEATURES)
    input_feature_name_to_idx = {
        name: i for i, name in enumerate(INPUT_FEATURES)
    }

    print("Training input features:")
    for i, name in enumerate(INPUT_FEATURES):
        print(f"  {i:2d}: {name}  dataset_idx={feature_idx[i]}")

    train_idx, val_idx, test_idx = make_splits(len(dataset), seed=args.seed)

    train_loader = DataLoader(
        torch.utils.data.Subset(dataset, train_idx),
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=args.num_workers,
        pin_memory=False,
    )
    val_loader = DataLoader(
        torch.utils.data.Subset(dataset, val_idx),
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        pin_memory=False,
    )
    test_loader = DataLoader(
        torch.utils.data.Subset(dataset, test_idx),
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        pin_memory=False,
    )

    model = build_model(
        nmax=NMAX,
        in_dim=len(INPUT_FEATURES),
        hidden_dim=args.hidden_dim,
        beta0=args.hgq_beta0,
    )
    model.summary()

    n_trainable = int(
        np.sum([np.prod(v.shape) for v in model.trainable_variables])
    )
    print(f"Trainable parameters: {n_trainable:,}")

    optimizer = keras.optimizers.AdamW(
        learning_rate=args.lr,
        weight_decay=0.01,
        beta_1=0.9,
        beta_2=0.999,
        epsilon=1e-8,
    )

    config = vars(args).copy()
    config["dataset_metadata"] = dataset.metadata
    config["dataset_feature_names"] = dataset_feature_names
    config["input_features"] = INPUT_FEATURES
    config["feature_idx"] = feature_idx
    config["nmax"] = NMAX
    config["in_dim"] = len(INPUT_FEATURES)
    config["framework"] = "keras3_hgq2"
    config["model"] = "DeepSet_mean_only_HGQ2"
    config["aggregation"] = "masked_mean"
    config["center_branch"] = False
    config["sum_branch"] = False
    config["layernorm"] = False
    config["batchnorm"] = True
    config["batchnorm_momentum"] = 0.9
    config["batchnorm_epsilon"] = 1e-5
    config["batchnorm_padding_treatment"] = "valid_particles_only"
    config["prediction"] = "softplus(raw)"
    config["baseline"] = (
        "softplus + weighted_huber(pos_weight=9) + "
        "plain log-ratio penalty(alpha=1)"
    )
    config["quantization"] = "HGQ2 learned bitwidths"
    config["hgq_quantizer_all"] = {
        "default_q_type": "kbi",
        "overflow_mode": "SAT_SYM",
    }
    config["hgq_quantizer_datalane"] = {
        "default_q_type": "kif",
        "overflow_mode": "WRAP",
    }
    config["hgq_enable_ebops"] = True
    config["hgq_beta0"] = args.hgq_beta0
    config["train_size"] = len(train_idx)
    config["val_size"] = len(val_idx)
    config["test_size"] = len(test_idx)

    with open(os.path.join(run_dir, "config.json"), "w") as f:
        json.dump(config, f, indent=2)

    history = {
        "epoch": [],
        "train_loss": [],
        "train_physics_loss": [],
        "train_ebops_loss": [],
        "val_loss": [],
    }
    best_val_loss = float("inf")
    best_epoch = -1

    print("Start training!")

    for epoch in range(1, args.epochs + 1):
        train_out = train_epoch(
            model,
            train_loader,
            optimizer,
            args,
            feature_idx,
        )
        train_loss = train_out["loss"]

        val_out = evaluate(
            model,
            val_loader,
            args,
            feature_idx,
            input_feature_name_to_idx,
        )
        val_loss = val_out["loss"]

        diag = compute_diagnostic_metrics(
            val_out["pred"],
            val_out["target"],
            eps=args.penalty_eps,
            zero_thr=args.penalty_zero_thr,
            pred_thr=args.pos_weight_threshold,
        )

        history["epoch"].append(epoch)
        history["train_loss"].append(train_loss)
        history["train_physics_loss"].append(train_out["physics_loss"])
        history["train_ebops_loss"].append(train_out["ebops_loss"])
        history["val_loss"].append(val_loss)

        for k, v in diag.items():
            history.setdefault(f"val_{k}", []).append(v)

        print(
            f"Epoch {epoch:03d} | "
            f"train_total={train_loss:.5f} | "
            f"train_phys={train_out['physics_loss']:.5f} | "
            f"train_ebops={train_out['ebops_loss']:.5f} | "
            f"val={val_loss:.5f} | "
            f"val_mse={diag['mse']:.5f} | "
            f"val_huber={diag['huber']:.5f} | "
            f"zero_mean_pred={diag['zero_mean_pred']:.5f} | "
            f"zero_nonzero_frac={diag['zero_nonzero_pred_frac']:.4f} | "
            f"pos_logR_rms={diag['pos_logratio_rms']:.5f}"
        )

        np.savez(
            os.path.join(run_dir, f"val_outputs_epoch{epoch}.npz"),
            **val_out,
        )
        np.savez(
            os.path.join(run_dir, "loss_history.npz"),
            **{k: np.asarray(v) for k, v in history.items()},
        )
        model.save_weights(
            os.path.join(run_dir, f"model_epoch{epoch}.weights.h5")
        )

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_epoch = epoch
            model.save_weights(
                os.path.join(run_dir, "model_best.weights.h5")
            )
            np.savez(
                os.path.join(run_dir, "val_outputs_best.npz"),
                **val_out,
            )

    print(f"Best epoch: {best_epoch} with val_loss={best_val_loss:.6f}")

    model.load_weights(
        os.path.join(run_dir, "model_best.weights.h5")
    )

    test_out = evaluate(
        model,
        test_loader,
        args,
        feature_idx,
        input_feature_name_to_idx,
    )
    test_diag = compute_diagnostic_metrics(
        test_out["pred"],
        test_out["target"],
        eps=args.penalty_eps,
        zero_thr=args.penalty_zero_thr,
        pred_thr=args.pos_weight_threshold,
    )

    np.savez(
        os.path.join(run_dir, "test_outputs_best.npz"),
        **test_out,
    )

    with open(
        os.path.join(run_dir, "test_metrics_best.json"),
        "w",
    ) as f:
        json.dump(test_diag, f, indent=2)

    print("Test metrics from best model:")
    for k, v in test_diag.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()

