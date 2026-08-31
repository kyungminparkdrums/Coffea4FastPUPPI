import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

import numpy as np
import keras
from keras import layers
import hls4ml

NMAX = 24
NFEAT = 18
HIDDEN = 16
PARALLELIZATION_FACTOR = 12


def pytorch_dense_kwargs(fan_in):
    limit = 1.0 / np.sqrt(fan_in)
    return {
        "kernel_initializer": keras.initializers.RandomUniform(minval=-limit, maxval=limit),
        "bias_initializer": keras.initializers.RandomUniform(minval=-limit, maxval=limit),
    }


# ============================================================
# Latest standard Keras training/reference model
# ============================================================
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


# ============================================================
# HLS-facing inference model
#
# Same trained weighted layers as the standard model.
# The only structural rewrite is the masked mean:
#
#   mean_valid(h) = GAP(h * mask) * NMAX / Nvalid
#
# scale must therefore be NMAX / Nvalid for each event.
# ============================================================
def build_hls_model():
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


# ============================================================
# Load latest standard Keras weights
# ============================================================
model = build_model()
model.load_weights("models/model_best.weights.h5")

print("STANDARD KERAS MODEL LOADED")
model.summary()


# ============================================================
# Transfer trained weights to HLS-facing model
# ============================================================
model_hls = build_hls_model()
model_hls.set_weights(model.get_weights())

print("WEIGHTS TRANSFERRED TO HLS MODEL")
model_hls.summary()


# ============================================================
# hls4ml configuration
# ============================================================
config = hls4ml.utils.config_from_keras_model(model_hls, granularity="name", backend="Vitis")

print("\nHLS4ML LAYERS:")
for name in config["LayerName"]:
    print(name)

config["LayerName"]["phi_dense1"]["ParallelizationFactor"] = PARALLELIZATION_FACTOR
config["LayerName"]["phi_dense2"]["ParallelizationFactor"] = PARALLELIZATION_FACTOR

print("\nEMBEDDING PARALLELIZATION:")
print("phi_dense1:", config["LayerName"]["phi_dense1"])
print("phi_dense2:", config["LayerName"]["phi_dense2"])


# ============================================================
# Convert / compile / synthesize
# ============================================================
hls_model = hls4ml.converters.convert_from_keras_model(model_hls, hls_config=config, backend="Vitis", output_dir="hls4ml_float", part="xcvu13p-flga2577-2-e", clock_period=2.778, io_type="io_parallel")

hls_model.compile()
print("MODEL COMPILED")

hls_model.build()
print("MODEL BUILT")

