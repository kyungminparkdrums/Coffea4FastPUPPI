########################################
# --- PUPPI ML: Dataset validation --- #
########################################

import argparse
import glob
import json
import os
import random
import sys
import tempfile
import types
import warnings

import numpy as np
import torch
from tqdm import tqdm

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


LEGACY_FEATURE_NAMES = [
    "log_pt",
    "slog_px",
    "slog_py",
    "charge",
    "nnvtx",
    "puppi",
    "pdg_211",
    "pdg_130",
    "pdg_22",
    "pdg_11",
    "pdg_13",
    "deta",
    "dphi",
    "dr",
    "is_center",
]


def install_torch_geometric_stubs():
    if "torch_geometric.data.data" in sys.modules:
        return

    tg_mod = types.ModuleType("torch_geometric")
    tg_data_pkg = types.ModuleType("torch_geometric.data")
    tg_data_mod = types.ModuleType("torch_geometric.data.data")
    tg_storage_mod = types.ModuleType("torch_geometric.data.storage")

    class GlobalStorage(dict):
        pass

    class Data:
        def __getattr__(self, name):
            store = self.__dict__.get("_store", None)
            if store is None:
                raise AttributeError(name)

            mapping = getattr(store, "_mapping", {})
            if name in mapping:
                return mapping[name]

            raise AttributeError(name)

        def keys(self):
            store = self.__dict__.get("_store", None)
            if store is None:
                return []
            return list(getattr(store, "_mapping", {}).keys())

    class DataEdgeAttr:
        pass

    class DataTensorAttr:
        pass

    tg_data_mod.Data = Data
    tg_data_mod.DataEdgeAttr = DataEdgeAttr
    tg_data_mod.DataTensorAttr = DataTensorAttr
    tg_storage_mod.GlobalStorage = GlobalStorage

    sys.modules["torch_geometric"] = tg_mod
    sys.modules["torch_geometric.data"] = tg_data_pkg
    sys.modules["torch_geometric.data.data"] = tg_data_mod
    sys.modules["torch_geometric.data.storage"] = tg_storage_mod


def load_graphs(path):
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning, message="You are using `torch.load`")
            payload = torch.load(path, map_location="cpu")
    except ModuleNotFoundError as exc:
        if "torch_geometric" not in str(exc):
            raise
        install_torch_geometric_stubs()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning, message="You are using `torch.load`")
            payload = torch.load(path, map_location="cpu")

    if isinstance(payload, dict) and "graphs" in payload:
        return payload["graphs"], payload.get("metadata", {})

    return payload, {}


def get_graph_attr(graph, name):
    if hasattr(graph, name):
        value = getattr(graph, name)
        if torch.is_tensor(value):
            return value.detach().cpu().numpy()
        return value

    store = getattr(graph, "_store", None)
    mapping = getattr(store, "_mapping", {})
    value = mapping.get(name, None)
    if torch.is_tensor(value):
        return value.detach().cpu().numpy()
    return value


def infer_feature_names(x_dim):
    if x_dim == len(LEGACY_FEATURE_NAMES):
        return LEGACY_FEATURE_NAMES
    if x_dim == 16:
        return [
            "log_pt",
            "slog_px",
            "slog_py",
            "abs_eta",
            "charge",
            "nnvtx",
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
    return [f"feature_{i}" for i in range(x_dim)]


def pick_feature_index(feature_to_idx, *names):
    for name in names:
        if name in feature_to_idx:
            return feature_to_idx[name]
    return None


def finalize_axis(ax, xlabel, ylabel="Entries", title=None, logy=False):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, alpha=0.25)


def hist_or_empty(ax, values, bins, xlabel, title=None, logy=False, hist_range=None):
    if values.size:
        ax.hist(values, bins=bins, range=hist_range, histtype="stepfilled", alpha=0.75)
    finalize_axis(ax, xlabel, title=title, logy=logy)


def make_output_dir(path):
    os.makedirs(path, exist_ok=True)
    return os.path.abspath(path)


def sanitize_feature_filename(feature_name):
    return "".join(ch if ch.isalnum() or ch in ("_", "-", ".") else "_" for ch in feature_name)


def get_feature_plot_config(feature_name, stats):
    puppi_label = stats.get("puppi_label", "PUPPI weight")

    defaults = {
        "log_pt": {"xlabel": "log(pT)", "logy": True},
        "slog_px": {"xlabel": "signed log(px)", "logy": True},
        "slog_py": {"xlabel": "signed log(py)", "logy": True},
        "abs_eta": {"xlabel": "|eta|", "logy": True},
        "charge": {"xlabel": "Charge", "bins": np.arange(-1.5, 2.6, 1.0)},
        "nnvtx": {"xlabel": "nnVtxScore", "logy": True},
        "puppi": {"xlabel": puppi_label, "logy": True, "hist_range": (0.0, 1.05)},
        "puppiWeight": {"xlabel": puppi_label, "logy": True, "hist_range": (0.0, 1.05)},
        "deta": {"xlabel": "DeltaEta", "logy": True, "hist_range": (-0.31, 0.31)},
        "deta_to_center": {"xlabel": "DeltaEta", "logy": True, "hist_range": (-0.31, 0.31)},
        "dphi": {"xlabel": "DeltaPhi", "logy": True, "hist_range": (-0.31, 0.31)},
        "dphi_to_center": {"xlabel": "DeltaPhi", "logy": True, "hist_range": (-0.31, 0.31)},
        "dr": {"xlabel": "DeltaR to seed", "logy": True, "hist_range": (0.0, 0.31)},
        "dr_to_center": {"xlabel": "DeltaR to seed", "logy": True, "hist_range": (0.0, 0.31)},
        "is_seed": {"xlabel": "is_seed", "bins": np.arange(-0.5, 2.0, 1.0)},
        "is_center": {"xlabel": "is_center", "bins": np.arange(-0.5, 2.0, 1.0)},
        "is_chHad": {"xlabel": "is_chHad", "bins": np.arange(-0.5, 2.0, 1.0)},
        "is_nHad": {"xlabel": "is_nHad", "bins": np.arange(-0.5, 2.0, 1.0)},
        "is_neuHad": {"xlabel": "is_neuHad", "bins": np.arange(-0.5, 2.0, 1.0)},
        "is_gamma": {"xlabel": "is_gamma", "bins": np.arange(-0.5, 2.0, 1.0)},
        "is_ele": {"xlabel": "is_ele", "bins": np.arange(-0.5, 2.0, 1.0)},
        "is_mu": {"xlabel": "is_mu", "bins": np.arange(-0.5, 2.0, 1.0)},
        "pdg_211": {"xlabel": "pdg_211", "bins": np.arange(-0.5, 2.0, 1.0)},
        "pdg_130": {"xlabel": "pdg_130", "bins": np.arange(-0.5, 2.0, 1.0)},
        "pdg_22": {"xlabel": "pdg_22", "bins": np.arange(-0.5, 2.0, 1.0)},
        "pdg_11": {"xlabel": "pdg_11", "bins": np.arange(-0.5, 2.0, 1.0)},
        "pdg_13": {"xlabel": "pdg_13", "bins": np.arange(-0.5, 2.0, 1.0)},
        "pdg_abs_211": {"xlabel": "pdg_abs_211", "bins": np.arange(-0.5, 2.0, 1.0)},
        "pdg_abs_130": {"xlabel": "pdg_abs_130", "bins": np.arange(-0.5, 2.0, 1.0)},
        "pdg_abs_22": {"xlabel": "pdg_abs_22", "bins": np.arange(-0.5, 2.0, 1.0)},
        "pdg_abs_11": {"xlabel": "pdg_abs_11", "bins": np.arange(-0.5, 2.0, 1.0)},
        "pdg_abs_13": {"xlabel": "pdg_abs_13", "bins": np.arange(-0.5, 2.0, 1.0)},
    }
    config = {"xlabel": feature_name, "bins": 80, "logy": False, "hist_range": None}
    config.update(defaults.get(feature_name, {}))
    return config


def plot_input_variables(stats, outdir):
    feature_names = stats.get("dataset_feature_names", [])
    n_panels = len(feature_names) + 1
    ncols = 3
    nrows = int(np.ceil(n_panels / ncols))
    per_feature_dir = os.path.join(outdir, "InputVariables")
    os.makedirs(per_feature_dir, exist_ok=True)

    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 4.2 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for iax, feature_name in enumerate(feature_names):
        config = get_feature_plot_config(feature_name, stats)
        values = stats["feature_values"][feature_name]
        title = f"Input feature: {feature_name}"
        hist_or_empty(axes[iax], values, config["bins"], config["xlabel"], title=title, logy=config["logy"], hist_range=config["hist_range"])

        single_fig, single_ax = plt.subplots(figsize=(8, 6))
        hist_or_empty(single_ax, values, config["bins"], config["xlabel"], title=title, logy=config["logy"], hist_range=config["hist_range"])
        single_fig.tight_layout()
        single_fig.savefig(
            os.path.join(per_feature_dir, f"{sanitize_feature_filename(feature_name)}.png"),
            dpi=150,
        )
        plt.close(single_fig)

    particle_labels = ["ch had", "n had", "gamma", "ele", "mu"]
    axes[len(feature_names)].bar(particle_labels, stats["particle_counts"], alpha=0.8)
    finalize_axis(axes[len(feature_names)], "PF candidate type", ylabel="Candidates", title="PDG class composition")

    for ax in axes[n_panels:]:
        ax.set_visible(False)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "input_variables.png"), dpi=150)
    plt.close(fig)


def plot_target_vs_puppi(stats, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    hist_or_empty(
        axes[0],
        stats["seed_puppi"],
        80,
        stats.get("puppi_label", "Seed PUPPI weight"),
        title=f"Seed {stats.get('puppi_short_label', 'PUPPI')} distribution",
        logy=True,
        hist_range=(0.0, 1.05),
    )

    hist_or_empty(
        axes[1],
        stats["target"],
        100,
        "Target",
        title="Target distribution",
        logy=True,
    )

    residual = stats["target"] - stats["seed_puppi"]
    hist_or_empty(
        axes[2],
        residual,
        100,
        f"Target - seed {stats.get('puppi_short_label', 'PUPPI')}",
        title="Residual distribution",
        logy=True,
        hist_range=(-3.0, 3.0),
    )

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "target_vs_puppi.png"), dpi=150)
    plt.close(fig)


def plot_cone_multiplicity(stats, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    hist_or_empty(axes[0], stats["n_nodes"], 60, "PF candidates in cone", title="Total multiplicity", logy=True)
    hist_or_empty(axes[1], stats["n_neutral"], 60, "Neutral PF candidates in cone", title="Neutral multiplicity", logy=True)
    hist_or_empty(axes[2], stats["n_charged"], 60, "Charged PF candidates in cone", title="Charged multiplicity", logy=True)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "cone_multiplicity.png"), dpi=150)
    plt.close(fig)


def plot_pf_cloud_examples(example_graphs, outdir):
    examples_dir = os.path.join(outdir, "Clouds-around-NeutralSeed")
    os.makedirs(examples_dir, exist_ok=True)

    for iexample, example in enumerate(example_graphs):
        x = example["x"]
        cone_index = example["cone_index"]
        feature_to_idx = example["feature_to_idx"]

        deta = x[:, feature_to_idx["deta"]]
        dphi = x[:, feature_to_idx["dphi"]]
        pt = np.exp(x[:, feature_to_idx["log_pt"]])
        charge = x[:, feature_to_idx["charge"]]
        is_seed = x[:, feature_to_idx["is_seed"]] > 0.5
        is_charged = np.abs(charge) > 0.5
        is_neutral = ~is_charged
        title = f"PF cloud around neutral seed\ncone={cone_index}"

        fig, ax = plt.subplots(figsize=(9, 7))

        scatter = None
        if np.any(is_neutral & ~is_seed):
            scatter = ax.scatter(
                deta[is_neutral & ~is_seed],
                dphi[is_neutral & ~is_seed],
                c=pt[is_neutral & ~is_seed],
                s=40,
                cmap="viridis",
                marker="o",
                alpha=0.90,
                edgecolors="none",
                label="Neutral PF",
            )
        if np.any(is_charged):
            charged = ax.scatter(
                deta[is_charged],
                dphi[is_charged],
                c=pt[is_charged],
                s=55,
                cmap="viridis",
                marker="^",
                alpha=0.95,
                edgecolors="black",
                linewidths=0.3,
                label="Charged PF",
            )
            if scatter is None:
                scatter = charged
        if np.any(is_seed):
            seed = ax.scatter(
                deta[is_seed],
                dphi[is_seed],
                c=pt[is_seed],
                s=220,
                cmap="viridis",
                marker="*",
                alpha=1.0,
                edgecolors="red",
                linewidths=1.0,
                label="Neutral seed",
            )
            if scatter is None:
                scatter = seed

        ax.axhline(0.0, color="gray", linewidth=1.0, alpha=0.7)
        ax.axvline(0.0, color="gray", linewidth=1.0, alpha=0.7)
        finalize_axis(ax, "Δη", ylabel="Δφ", title=title)
        ax.set_xlim(-0.31, 0.31)
        ax.set_ylim(-0.31, 0.31)
        ax.legend(loc="upper right", fontsize=9)

        if scatter is not None:
            cbar = fig.colorbar(scatter, ax=ax, shrink=0.92, pad=0.02)
            cbar.set_label("PF candidate pT [GeV]")

        fig.tight_layout()
        fig.savefig(os.path.join(examples_dir, f"cone_{iexample:02d}_2d.png"), dpi=160)
        plt.close(fig)

        fig3d = plt.figure(figsize=(9, 7))
        ax3d = fig3d.add_subplot(1, 1, 1, projection="3d")

        scatter3d = None
        if np.any(is_neutral & ~is_seed):
            scatter3d = ax3d.scatter(
                deta[is_neutral & ~is_seed],
                dphi[is_neutral & ~is_seed],
                pt[is_neutral & ~is_seed],
                c=pt[is_neutral & ~is_seed],
                s=40,
                cmap="viridis",
                marker="o",
                alpha=0.90,
                edgecolors="none",
                label="Neutral PF",
            )
        if np.any(is_charged):
            charged3d = ax3d.scatter(
                deta[is_charged],
                dphi[is_charged],
                pt[is_charged],
                c=pt[is_charged],
                s=55,
                cmap="viridis",
                marker="^",
                alpha=0.95,
                edgecolors="black",
                linewidths=0.3,
                label="Charged PF",
            )
            if scatter3d is None:
                scatter3d = charged3d
        if np.any(is_seed):
            seed3d = ax3d.scatter(
                deta[is_seed],
                dphi[is_seed],
                pt[is_seed],
                c=pt[is_seed],
                s=220,
                cmap="viridis",
                marker="*",
                alpha=1.0,
                edgecolors="red",
                linewidths=1.0,
                label="Neutral seed",
            )
            if scatter3d is None:
                scatter3d = seed3d

        ax3d.set_xlabel("Δη")
        ax3d.set_ylabel("Δφ")
        ax3d.set_zlabel("PF pT [GeV]")
        ax3d.set_xlim(-0.31, 0.31)
        ax3d.set_ylim(-0.31, 0.31)
        ax3d.set_title(title, fontsize=11)
        ax3d.legend(loc="upper right", fontsize=9)

        if scatter3d is not None:
            cbar3d = fig3d.colorbar(scatter3d, ax=ax3d, shrink=0.82, pad=0.08)
            cbar3d.set_label("PF candidate pT [GeV]")

        fig3d.subplots_adjust(left=0.02, right=0.90, bottom=0.08, top=0.88)
        fig3d.savefig(os.path.join(examples_dir, f"cone_{iexample:02d}_3d.png"), dpi=160)
        plt.close(fig3d)


def write_summary(stats, outdir, files_used):
    summary_path = os.path.join(outdir, "summary.txt")
    with open(summary_path, "w", encoding="ascii") as handle:
        handle.write(f"Files used: {files_used}\n")
        handle.write(f"Cones processed: {len(stats['target'])}\n")
        handle.write(f"PF candidates processed: {len(stats['pt_all'])}\n")
        if len(stats["target"]):
            handle.write(f"Target mean: {np.mean(stats['target']):.6f}\n")
            handle.write(f"Target std: {np.std(stats['target']):.6f}\n")
        if len(stats["seed_puppi"]):
            handle.write(f"Seed {stats.get('puppi_short_label', 'PUPPI')} mean: {np.mean(stats['seed_puppi']):.6f}\n")
            handle.write(f"Seed {stats.get('puppi_short_label', 'PUPPI')} std: {np.std(stats['seed_puppi']):.6f}\n")


def main():
    parser = argparse.ArgumentParser(description="Validation plots for DEEPSets neutral-regression graph batches")
    parser.add_argument("--input", required=True, help="Glob pattern for the input .pt graph files")
    parser.add_argument("--output_dir", required=True, help="Directory where the plots and summary are written")
    parser.add_argument("--max_files", type=int, default=10, help="Maximum number of .pt files to load after shuffling, similar to training")
    parser.add_argument("--max_graphs", type=int, default=None, help="Optional cap on the number of cones")
    parser.add_argument("--n_examples", type=int, default=5, help="Number of explicit cone-cloud examples to save")
    parser.add_argument("--seed", type=int, default=42, help="Seed used when shuffling input files before applying --max_files")
    args = parser.parse_args()

    files = sorted(glob.glob(args.input))
    random.Random(args.seed).shuffle(files)
    if args.max_files is not None:
        files = files[: args.max_files]

    if not files:
        raise FileNotFoundError(f"No .pt files matched pattern: {args.input}")

    outdir = make_output_dir(args.output_dir)

    stats = {
        "pt_all": [],
        "log_pt_all": [],
        "puppi_all": [],
        "nnvtx_all": [],
        "charge_all": [],
        "deta_all": [],
        "dphi_all": [],
        "dr_all": [],
        "target": [],
        "seed_puppi": [],
        "n_nodes": [],
        "n_neutral": [],
        "n_charged": [],
        "particle_counts": np.zeros(5, dtype=np.int64),
        "puppi_label": "PUPPI weight",
        "puppi_short_label": "PUPPI",
        "dataset_feature_names": [],
        "feature_values": {},
    }
    example_graphs = []
    dataset_feature_names = None
    feature_to_idx = None
    skip_missing_xy = 0
    skip_invalid_seed = 0

    graphs_seen = 0

    for path in tqdm(files, desc="Loading .pt files"):
        graphs, metadata = load_graphs(path)
        for graph in graphs:
            x = get_graph_attr(graph, "x")
            y = get_graph_attr(graph, "y")

            if x is None or y is None:
                skip_missing_xy += 1
                continue

            x = np.asarray(x, dtype=np.float32)
            y = np.asarray(y, dtype=np.float32).reshape(-1)

            if x.ndim != 2:
                raise ValueError(f"Unexpected feature shape {x.shape} in {path}. Expected rank-2 node features.")

            current_feature_names = metadata.get("feature_names") or metadata.get("model_features")
            if isinstance(current_feature_names, str):
                try:
                    current_feature_names = json.loads(current_feature_names)
                except json.JSONDecodeError:
                    current_feature_names = [name.strip() for name in current_feature_names.split(",")]
            if current_feature_names is None:
                current_feature_names = infer_feature_names(x.shape[1])

            if dataset_feature_names is None:
                dataset_feature_names = current_feature_names
                feature_to_idx = {name: i for i, name in enumerate(dataset_feature_names)}
                stats["dataset_feature_names"] = list(dataset_feature_names)
                stats["feature_values"] = {name: [] for name in dataset_feature_names}
            elif current_feature_names != dataset_feature_names:
                raise ValueError(
                    f"Feature-name mismatch in {path}. "
                    f"Existing: {dataset_feature_names} "
                    f"New: {current_feature_names}"
                )

            log_pt_idx = pick_feature_index(feature_to_idx, "log_pt")
            charge_idx = pick_feature_index(feature_to_idx, "charge")
            nnvtx_idx = pick_feature_index(feature_to_idx, "nnvtx")
            puppi_idx = pick_feature_index(feature_to_idx, "puppiWeight", "puppi")
            deta_idx = pick_feature_index(feature_to_idx, "deta", "deta_to_center")
            dphi_idx = pick_feature_index(feature_to_idx, "dphi", "dphi_to_center")
            dr_idx = pick_feature_index(feature_to_idx, "dr", "dr_to_center")
            is_seed_idx = pick_feature_index(feature_to_idx, "is_seed", "is_center")

            pdg_indices = [
                pick_feature_index(feature_to_idx, "is_chHad", "pdg_211", "pdg_abs_211"),
                pick_feature_index(feature_to_idx, "is_nHad", "is_neuHad", "pdg_130", "pdg_abs_130"),
                pick_feature_index(feature_to_idx, "is_gamma", "pdg_22", "pdg_abs_22"),
                pick_feature_index(feature_to_idx, "is_ele", "pdg_11", "pdg_abs_11"),
                pick_feature_index(feature_to_idx, "is_mu", "pdg_13", "pdg_abs_13"),
            ]

            required = {
                "log_pt": log_pt_idx,
                "charge": charge_idx,
                "deta": deta_idx,
                "dphi": dphi_idx,
                "dr": dr_idx,
                "is_seed": is_seed_idx,
            }
            missing_required = [name for name, idx in required.items() if idx is None]
            if missing_required:
                raise ValueError(
                    f"Missing required features {missing_required} in {path}. "
                    f"Available: {dataset_feature_names}"
                )

            if puppi_idx is not None:
                stats["puppi_label"] = "PUPPI weight"
                stats["puppi_short_label"] = "PUPPI"

            log_pt = x[:, log_pt_idx]
            pt = np.exp(log_pt)
            charge = x[:, charge_idx]
            nnvtx = x[:, nnvtx_idx] if nnvtx_idx is not None else np.zeros_like(log_pt)
            puppi = x[:, puppi_idx] if puppi_idx is not None else np.zeros_like(log_pt)
            onehot = np.stack(
                [
                    x[:, idx] if idx is not None else np.zeros_like(log_pt)
                    for idx in pdg_indices
                ],
                axis=1,
            )
            deta = x[:, deta_idx]
            dphi = x[:, dphi_idx]
            dr = x[:, dr_idx]
            is_seed = x[:, is_seed_idx] > 0.5

            if is_seed.sum() != 1:
                skip_invalid_seed += 1
                continue

            seed_idx = np.argmax(is_seed)

            stats["pt_all"].append(pt)
            stats["log_pt_all"].append(log_pt)
            stats["puppi_all"].append(puppi)
            stats["nnvtx_all"].append(nnvtx)
            stats["charge_all"].append(charge)
            stats["deta_all"].append(deta)
            stats["dphi_all"].append(dphi)
            stats["dr_all"].append(dr)
            for feature_name, idx in feature_to_idx.items():
                stats["feature_values"][feature_name].append(x[:, idx])
            stats["target"].append(float(y[0]))
            stats["seed_puppi"].append(float(puppi[seed_idx]))
            stats["n_nodes"].append(x.shape[0])
            stats["n_charged"].append(int(np.sum(np.abs(charge) > 0.5)))
            stats["n_neutral"].append(int(np.sum(np.abs(charge) <= 0.5)))
            stats["particle_counts"] += np.rint(onehot.sum(axis=0)).astype(np.int64)
            if len(example_graphs) < args.n_examples:
                example_graphs.append(
                    {
                        "x": x.copy(),
                        "target": float(y[0]),
                        "seed_puppi": float(puppi[seed_idx]),
                        "cone_index": graphs_seen,
                        "feature_to_idx": {
                            "log_pt": log_pt_idx,
                            "charge": charge_idx,
                            "deta": deta_idx,
                            "dphi": dphi_idx,
                            "is_seed": is_seed_idx,
                        },
                    }
                )

            graphs_seen += 1
            if args.max_graphs is not None and graphs_seen >= args.max_graphs:
                break

        if args.max_graphs is not None and graphs_seen >= args.max_graphs:
            break

    # Keep label/config entries as plain Python objects; only numerical accumulators
    # should be converted to NumPy arrays here.
    array_fields = {
        "pt_all",
        "log_pt_all",
        "puppi_all",
        "nnvtx_all",
        "charge_all",
        "deta_all",
        "dphi_all",
        "dr_all",
    }
    scalar_fields = {"target", "seed_puppi", "n_nodes", "n_charged", "n_neutral"}

    for key in array_fields:
        value = stats[key]
        stats[key] = np.concatenate(value) if value else np.array([], dtype=np.float32)

    for key in scalar_fields:
        stats[key] = np.asarray(stats[key])

    for feature_name, values in list(stats["feature_values"].items()):
        stats["feature_values"][feature_name] = np.concatenate(values) if values else np.array([], dtype=np.float32)

    if graphs_seen == 0:
        raise RuntimeError(
            "No valid graphs were processed by validate_dataset.py. "
            f"skip_missing_xy={skip_missing_xy}, "
            f"skip_invalid_seed={skip_invalid_seed}, "
            f"dataset_feature_names={dataset_feature_names}"
        )

    plot_input_variables(stats, outdir)
    plot_target_vs_puppi(stats, outdir)
    plot_cone_multiplicity(stats, outdir)
    plot_pf_cloud_examples(example_graphs, outdir)
    write_summary(stats, outdir, len(files))

    print(f"Saved validation plots in: {outdir}")


if __name__ == "__main__":
    main()
