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
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm


# Upper bound on particle multiplicity near the neutral
NMAX = 32

# Dataset can contain a superset of features, but when we run the training some features can be commented out.
INPUT_FEATURES = [
    "log_pt",
    "slog_px",
    "slog_py",
    "abs_eta",
    "charge",
    "nnvtx",
    # "vz",
    "dxy",
    "z0",
    "puppiWeight",
    # "idProbPu",
    # "idProbEm",
    # "idProbPi",
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
        self.files = files[:max_files]  # each dataset currently contains ~50k neutrals

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

            self.data.extend(graphs)
            print(f"Loaded {f}: {len(graphs)} graphs in {time.time() - t0:.1f} s")

        print(f"Loaded {len(self.data)} graphs from {len(self.files)} files")

        if len(self.data) == 0:
            raise RuntimeError(f"No graphs were loaded from pattern: {pattern}")

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
        n_centers = center_mask.sum().item()

        if n_centers != 1:
            raise RuntimeError(f"Expected exactly one center, got {n_centers}")

        center_idx_local = torch.where(center_mask)[0].item()

        # Center first, then all other particles in original order
        other_idx = torch.arange(x.size(0), dtype=torch.long)
        other_idx = other_idx[other_idx != center_idx_local]
        order = torch.cat([torch.tensor([center_idx_local], dtype=torch.long), other_idx])
        x = x[order]

        n = x.size(0)

        # Deliberately do not truncate: fail if NMAX is insufficient.
        if n > NMAX:
            raise RuntimeError(f"Graph has {n} particles, exceeds NMAX={NMAX}")

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


class DeepSetRegressor(nn.Module):
    def __init__(self, in_dim, hidden_dim=256):
        super().__init__()

        self.phi = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
        )

        self.rho = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            #nn.Linear(hidden_dim * 3, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
        )


    def forward(self, x, mask):
        B, N, Fdim = x.shape

        valid = mask.bool()
        x_valid = x[valid]

        h_valid = self.phi(x_valid)

        h = torch.zeros(B, N, h_valid.size(1), device=x.device, dtype=h_valid.dtype)
        h[valid] = h_valid

        mask_f = mask.unsqueeze(-1).to(h.dtype)

        h_sum = h.sum(dim=1)
        counts = mask_f.sum(dim=1).clamp(min=1.0)
        h_mean = h_sum / counts
        h_center = h[:, 0, :]

        h_out = h_mean
        #h_out = torch.cat([h_sum, h_mean, h_center], dim=1)
        return self.rho(h_out).view(-1)

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def make_splits(n, seed=42):
    idx = list(range(n))
    random.Random(seed).shuffle(idx)

    n_train = int(0.7 * n)
    n_val = int(0.15 * n)

    return idx[:n_train], idx[n_train:n_train + n_val], idx[n_train + n_val:]


def make_run_dir():
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = f"runs/run_{ts}"
    os.makedirs(run_dir, exist_ok=True)
    return run_dir


def apply_output_activation(raw, activation="softplus"):
    if activation == "relu":
        return F.relu(raw)
    if activation == "softplus":
        return F.softplus(raw)
    if activation == "identity":
        return raw
    raise ValueError(f"Unknown output activation: {activation}")


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


def global_shift_penalty(pred):
    return pred.mean()


def puppi_loss(
    pred,
    target,
    loss_type="weighted_huber",
    alpha=1.0,
    mse_penalty_alpha=0.0,
    global_shift_alpha=0.0,
    eps=1e-3,
    zero_thr=0.05,
    raw=None,
    pos_weight=9.0,
    pos_weight_threshold=0.05,
):
    if loss_type == "huber":
        loss_abs = F.huber_loss(pred, target, reduction="none")

    elif loss_type == "weighted_huber":
        loss_abs = F.huber_loss(pred, target, reduction="none")

        weight = torch.ones_like(target)
        pos_mask = target > pos_weight_threshold
        weight[pos_mask] = pos_weight

        loss_abs = loss_abs * weight

    elif loss_type == "mse":
        loss_abs = F.mse_loss(pred, target, reduction="none")

    elif loss_type == "mae":
        loss_abs = F.l1_loss(pred, target, reduction="none")

    else:
        raise ValueError(f"Unknown loss_type: {loss_type}")

    loss = loss_abs.mean()

    if mse_penalty_alpha > 0.0:
        loss = loss + mse_penalty_alpha * torch.mean((pred - target) ** 2)

    # Plain symmetric log-ratio MSE on positive targets.
    if alpha > 0.0:
        pos_mask = target > zero_thr

        if torch.any(pos_mask):
            pred_pos = pred[pos_mask].clamp_min(eps)
            target_pos = target[pos_mask].clamp_min(eps)

            log_ratio = torch.log(pred_pos) - torch.log(target_pos)
            loss_logratio = torch.mean(log_ratio ** 2)
            loss = loss + alpha * loss_logratio

    if global_shift_alpha > 0.0:
        loss = loss + global_shift_alpha * global_shift_penalty(pred)

    return loss


def compute_diagnostic_metrics(pred, target, eps=1e-3, zero_thr=0.05, pred_thr=0.05):
    pred = np.asarray(pred)
    target = np.asarray(target)

    err = pred - target
    metrics = {}

    metrics["mse"] = float(np.mean(err ** 2))
    metrics["mae"] = float(np.mean(np.abs(err)))

    abs_err = np.abs(err)
    metrics["huber"] = float(np.mean(np.where(abs_err < 1.0, 0.5 * err**2, abs_err - 0.5)))

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


def train_epoch(model, loader, optimizer, device, args, feature_idx):
    model.train()
    total_loss = 0.0

    for data in tqdm(loader, desc="Train", leave=False):
        x = data["x"][:, :, feature_idx].to(device)
        mask = data["mask"].to(device)
        target = data["y"].view(-1).to(device)

        raw = model(x, mask)
        pred = apply_output_activation(raw, args.output_activation)

        loss = puppi_loss(
            pred,
            target,
            loss_type=args.loss_type,
            alpha=args.penalty_alpha,
            mse_penalty_alpha=args.mse_penalty_alpha,
            global_shift_alpha=args.global_shift_alpha,
            eps=args.penalty_eps,
            zero_thr=args.penalty_zero_thr,
            raw=raw,
            pos_weight=args.pos_weight,
            pos_weight_threshold=args.pos_weight_threshold,
        )

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    return total_loss / len(loader)


@torch.no_grad()
def evaluate(model, loader, device, args, feature_idx, input_feature_name_to_idx):
    model.eval()

    preds, raws, targets = [], [], []
    event_idx, center_idx, center_pt, center_eta = [], [], [], []
    losses = []

    for data in tqdm(loader, desc="Eval", leave=False):
        x = data["x"][:, :, feature_idx].to(device)
        mask = data["mask"].to(device)
        target = data["y"].view(-1).to(device)

        raw = model(x, mask)
        pred = apply_output_activation(raw, args.output_activation)
        pred_for_metrics = pred.clamp_min(0.0) if args.eval_clip_nonnegative else pred

        loss = puppi_loss(
            pred,
            target,
            loss_type=args.loss_type,
            alpha=args.penalty_alpha,
            mse_penalty_alpha=args.mse_penalty_alpha,
            global_shift_alpha=args.global_shift_alpha,
            eps=args.penalty_eps,
            zero_thr=args.penalty_zero_thr,
            raw=raw,
            pos_weight=args.pos_weight,
            pos_weight_threshold=args.pos_weight_threshold,
        )

        losses.append(loss.item())

        raws.append(raw.cpu())
        preds.append(pred_for_metrics.cpu())
        targets.append(target.cpu())

        event_idx.append(data["event_idx"].view(-1).cpu())
        center_idx.append(data["center_idx"].view(-1).cpu())

        center_pt.append(torch.exp(x[:, 0, input_feature_name_to_idx["log_pt"]]).cpu())

        if "abs_eta" in input_feature_name_to_idx:
            center_eta.append(x[:, 0, input_feature_name_to_idx["abs_eta"]].cpu())

    out = {
        "loss": float(np.mean(losses)),
        "raw": torch.cat(raws).numpy(),
        "pred": torch.cat(preds).numpy(),
        "target": torch.cat(targets).numpy(),
        "event_idx": torch.cat(event_idx).numpy(),
        "center_idx": torch.cat(center_idx).numpy(),
        "center_pt": torch.cat(center_pt).numpy(),
    }

    if len(center_eta) > 0:
        out["center_eta"] = torch.cat(center_eta).numpy()

    return out


@torch.no_grad()
def permutation_feature_importance(
    model,
    loader,
    device,
    args,
    feature_idx,
    input_feature_names,
    max_batches=20,
):
    model.eval()

    baseline_losses = []
    cached_batches = []

    for ibatch, data in enumerate(loader):
        if ibatch >= max_batches:
            break

        cached_batches.append({k: v.clone() for k, v in data.items()})

        x = data["x"][:, :, feature_idx].to(device)
        mask = data["mask"].to(device)
        target = data["y"].view(-1).to(device)

        raw = model(x, mask)
        pred = apply_output_activation(raw, args.output_activation)

        loss = puppi_loss(
            pred,
            target,
            loss_type=args.loss_type,
            alpha=args.penalty_alpha,
            mse_penalty_alpha=args.mse_penalty_alpha,
            global_shift_alpha=args.global_shift_alpha,
            eps=args.penalty_eps,
            zero_thr=args.penalty_zero_thr,
            raw=raw,
            pos_weight=args.pos_weight,
            pos_weight_threshold=args.pos_weight_threshold,
        )

        baseline_losses.append(loss.item())

    if len(baseline_losses) == 0:
        raise RuntimeError("No batches available for permutation feature importance.")

    baseline_loss = float(np.mean(baseline_losses))

    results = []
    protected_features = {"is_center"}

    for ifeat, name in enumerate(input_feature_names):
        if name in protected_features:
            results.append({
                "feature_index": ifeat,
                "feature": name,
                "baseline_loss": baseline_loss,
                "permuted_loss": np.nan,
                "delta_loss": np.nan,
                "note": "Skipped because permuting this feature breaks the fixed center definition.",
            })
            continue

        perm_losses = []

        for data_cpu in cached_batches:
            x = data_cpu["x"][:, :, feature_idx].clone()
            mask = data_cpu["mask"].clone()
            target = data_cpu["y"].view(-1).to(device)

            # Preserve padding: only permute this feature among valid particle rows.
            valid = mask.bool()
            feature_plane = x[:, :, ifeat]
            values = feature_plane[valid].clone()
            perm = torch.randperm(values.numel())
            feature_plane[valid] = values[perm]
            x[:, :, ifeat] = feature_plane

            x = x.to(device)
            mask = mask.to(device)

            raw = model(x, mask)
            pred = apply_output_activation(raw, args.output_activation)

            loss = puppi_loss(
                pred,
                target,
                loss_type=args.loss_type,
                alpha=args.penalty_alpha,
                mse_penalty_alpha=args.mse_penalty_alpha,
                global_shift_alpha=args.global_shift_alpha,
                eps=args.penalty_eps,
                zero_thr=args.penalty_zero_thr,
                raw=raw,
                pos_weight=args.pos_weight,
                pos_weight_threshold=args.pos_weight_threshold,
            )

            perm_losses.append(loss.item())

        perm_loss = float(np.mean(perm_losses))

        results.append({
            "feature_index": ifeat,
            "feature": name,
            "baseline_loss": baseline_loss,
            "permuted_loss": perm_loss,
            "delta_loss": perm_loss - baseline_loss,
        })

    finite_results = [x for x in results if np.isfinite(x["delta_loss"])]
    skipped_results = [x for x in results if not np.isfinite(x["delta_loss"])]
    finite_results = sorted(finite_results, key=lambda x: x["delta_loss"], reverse=True)
    return finite_results + skipped_results


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--data", required=True)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch_size", type=int, default=1024)
    parser.add_argument("--num_workers", type=int, default=4)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max_files", type=int, default=10)
    parser.add_argument("--hidden_dim", type=int, default=256)

    parser.add_argument(
        "--output_activation",
        choices=["relu", "softplus", "identity"],
        default="softplus",
    )
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

    parser.add_argument("--negative_penalty_beta", type=float, default=0.0)
    parser.add_argument("--mse_penalty_alpha", type=float, default=0.0)
    parser.add_argument("--global_shift_alpha", type=float, default=0.0)

    parser.add_argument("--feature_importance", action="store_true")
    parser.add_argument("--feature_importance_batches", type=int, default=20)
    parser.add_argument("--eval_clip_nonnegative", action="store_true")

    args = parser.parse_args()

    set_seed(args.seed)

    run_dir = make_run_dir()
    print("Run dir:", run_dir)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    dataset = BestPuppiDataset(args.data, max_files=args.max_files)

    dataset_feature_names = dataset.feature_names
    feature_idx = build_feature_indices(dataset_feature_names, INPUT_FEATURES)
    input_feature_name_to_idx = {name: i for i, name in enumerate(INPUT_FEATURES)}

    print("Training input features:")
    for i, name in enumerate(INPUT_FEATURES):
        print(f"  {i:2d}: {name}  dataset_idx={feature_idx[i]}")

    train_idx, val_idx, test_idx = make_splits(len(dataset), seed=args.seed)

    train_loader = DataLoader(
        torch.utils.data.Subset(dataset, train_idx),
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=args.num_workers,
        pin_memory=True,
    )

    val_loader = DataLoader(
        torch.utils.data.Subset(dataset, val_idx),
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=2,
    )

    test_loader = DataLoader(
        torch.utils.data.Subset(dataset, test_idx),
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=2,
    )

    model = DeepSetRegressor(
        in_dim=len(INPUT_FEATURES),
        hidden_dim=args.hidden_dim,
    ).to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr)

    config = vars(args).copy()
    config["dataset_metadata"] = dataset.metadata
    config["dataset_feature_names"] = dataset_feature_names
    config["input_features"] = INPUT_FEATURES
    config["feature_idx"] = feature_idx
    config["in_dim"] = len(INPUT_FEATURES)
    config["nmax"] = NMAX
    config["model"] = "DeepSetRegressor_fixedN_sum_mean_center"
    config["prediction"] = f"{args.output_activation}(raw)"
    config["baseline"] = (
        "softplus + weighted_huber(pos_weight=9) + "
        "plain log-ratio penalty(alpha=1)"
    )
    config["train_size"] = len(train_idx)
    config["val_size"] = len(val_idx)
    config["test_size"] = len(test_idx)

    with open(os.path.join(run_dir, "config.json"), "w") as f:
        json.dump(config, f, indent=2)

    history = {"epoch": [], "train_loss": [], "val_loss": []}
    best_val_loss = float("inf")
    best_epoch = -1

    print("Start training!")

    for epoch in range(1, args.epochs + 1):
        train_loss = train_epoch(
            model,
            train_loader,
            optimizer,
            device,
            args,
            feature_idx,
        )

        val_out = evaluate(
            model,
            val_loader,
            device,
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
        history["val_loss"].append(val_loss)

        for k, v in diag.items():
            history.setdefault(f"val_{k}", []).append(v)

        print(
            f"Epoch {epoch:03d} | "
            f"train={train_loss:.5f} | "
            f"val={val_loss:.5f} | "
            f"val_mse={diag['mse']:.5f} | "
            f"val_huber={diag['huber']:.5f} | "
            f"zero_mean_pred={diag['zero_mean_pred']:.5f} | "
            f"zero_nonzero_frac={diag['zero_nonzero_pred_frac']:.4f} | "
            f"pos_logR_rms={diag['pos_logratio_rms']:.5f}"
        )

        np.savez(os.path.join(run_dir, f"val_outputs_epoch{epoch}.npz"), **val_out)
        np.savez(
            os.path.join(run_dir, "loss_history.npz"),
            **{k: np.array(v) for k, v in history.items()},
        )

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_epoch = epoch

            torch.save(model.state_dict(), os.path.join(run_dir, "model_best.pt"))
            np.savez(os.path.join(run_dir, "val_outputs_best.npz"), **val_out)

        torch.save(
            model.state_dict(),
            os.path.join(run_dir, f"model_epoch{epoch}.pt"),
        )

    print(f"Best epoch: {best_epoch} with val_loss={best_val_loss:.6f}")

    best_path = os.path.join(run_dir, "model_best.pt")
    model.load_state_dict(
        torch.load(best_path, map_location=device, weights_only=True)
    )

    test_out = evaluate(
        model,
        test_loader,
        device,
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

    np.savez(os.path.join(run_dir, "test_outputs_best.npz"), **test_out)

    with open(os.path.join(run_dir, "test_metrics_best.json"), "w") as f:
        json.dump(test_diag, f, indent=2)

    print("Test metrics from best model:")
    for k, v in test_diag.items():
        print(f"  {k}: {v}")

    if args.feature_importance:
        print("Running permutation feature importance...")

        fi = permutation_feature_importance(
            model,
            val_loader,
            device,
            args,
            feature_idx=feature_idx,
            input_feature_names=INPUT_FEATURES,
            max_batches=args.feature_importance_batches,
        )

        with open(os.path.join(run_dir, "feature_importance.json"), "w") as f:
            json.dump(fi, f, indent=2)

        print("Top feature importances:")
        for item in fi[:15]:
            if np.isfinite(item["delta_loss"]):
                print(
                    f"  {item['feature_index']:2d} "
                    f"{item['feature']:15s} "
                    f"delta_loss={item['delta_loss']:.6f}"
                )
            else:
                print(
                    f"  {item['feature_index']:2d} "
                    f"{item['feature']:15s} skipped"
                )


if __name__ == "__main__":
    main()
