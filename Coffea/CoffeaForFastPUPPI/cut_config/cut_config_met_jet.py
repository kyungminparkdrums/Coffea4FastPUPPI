import os
import json
import awkward as ak
import numpy as np

# Select a matching checkpoint/config pair here (one model per analyzer run).
# The JSON controls input order, neighborhood dr, hidden_dim, activation and nmax.
# Set this directory to another training run to evaluate a different model.
_DEFAULT_RUN_DIR = "/home/phys/kyungmip/MLPUPPI/Coffea4FastPUPPI/Training/runs/run_20260916_085421"
ML_MODEL_PATH = os.environ.get(
    "ML_PUPPI_MODEL", os.path.join(_DEFAULT_RUN_DIR, "model_best.pt")
)
ML_CONFIG_PATH = os.environ.get(
    "ML_PUPPI_CONFIG", os.path.join(os.path.dirname(ML_MODEL_PATH), "config.json")
)


# Branch mapping
PUPPI_COLLECTION = "L1PuppiCands"
GEN_COLLECTION = "GenCands"
GENJET_COLLECTION = "GenJets"

BRANCH_MAP = {
    "pt": "L1PuppiCands_pt",
    "eta": "L1PuppiCands_eta",
    "phi": "L1PuppiCands_phi",
    "charge": "L1PuppiCands_charge",
    "pdgId": "L1PuppiCands_pdgId",
    "nnvtx": "L1PuppiCands_nnVtxScore",
    "vz": "L1PuppiCands_vz",
    "dxy": "L1PuppiCands_dxy",
    "z0": "L1PuppiCands_z0",
    "puppiWeight": "L1PuppiCands_puppiWeight",
    "idProbPu": "L1PuppiCands_hgcIdPu",
    "idProbEm": "L1PuppiCands_hgcIdEm",
    "idProbPi": "L1PuppiCands_hgcIdPi",
    "genPtSum0p2": "L1PuppiCands_genPtSum0p2",
    "recoDen_seedPlusAll0p2": "L1PuppiCands_recoDen_seedPlusAll0p2",
}

# NanoEvents exposes fields inside L1PuppiCands without the collection prefix.
COLLECTION_FIELD_MAP = {
    name: branch.removeprefix(f"{PUPPI_COLLECTION}_")
    for name, branch in BRANCH_MAP.items()
}

# Target branches
GENTARGET_GEN_SUM_FIELD = "genPtSum0p2"
GENTARGET_RECO_SUM_FIELD = "recoDen_seedPlusAll0p2"

# ML env setup
ML_DEVICE = os.environ.get("ML_PUPPI_DEVICE", "cuda")
ML_INFER_BATCH_GRAPHS = int(os.environ.get("ML_PUPPI_BATCH_GRAPHS", "8192"))

# some gen-level definitinos
DO_GEN_STATUS1 = True

GENMET_EDGES = np.linspace(0.0, 200.0, 20)
GENMET_BINS = [(float(lo), float(hi), f"bin{i:02d}")
               for i, (lo, hi) in enumerate(zip(GENMET_EDGES[:-1], GENMET_EDGES[1:]))]

STORED_GENMET_COLLECTION = os.environ.get("GENMET_COLLECTION", "genMetCentral")
STORED_GENMET_PT = os.environ.get("GENMET_PT_BRANCH", "genMetCentral_pt")
STORED_GENMET_PHI = os.environ.get("GENMET_PHI_BRANCH", "genMetCentral_phi")


def _stored_genmet(events):
    """Read scalar-per-event MET, either flat or grouped by NanoEvents."""
    if STORED_GENMET_PT in events.fields and STORED_GENMET_PHI in events.fields:
        pt, phi = events[STORED_GENMET_PT], events[STORED_GENMET_PHI]
    elif STORED_GENMET_COLLECTION in events.fields:
        collection = events[STORED_GENMET_COLLECTION]
        if not {"pt", "phi"}.issubset(collection.fields):
            raise ValueError("Stored GEN MET collection requires pt and phi")
        pt, phi = collection.pt, collection.phi
    else:
        raise KeyError(
            f"Stored GEN MET branches {STORED_GENMET_PT}, {STORED_GENMET_PHI} "
            f"or collection {STORED_GENMET_COLLECTION} not found. "
            "Set GENMET_PT_BRANCH/GENMET_PHI_BRANCH or GENMET_COLLECTION to the ntuple names."
        )
    if pt.ndim != 1 or phi.ndim != 1:
        raise ValueError("Stored GEN MET must contain one scalar per event")
    return pt * np.cos(phi), pt * np.sin(phi), pt


GENJET_PT_BINS = [
    (0, 20, "0_20"),
    (20, 40, "20_40"),
    (40, 60, "40_60"),
    (60, 80, "60_80"),
    (80, 100, "80_100"),
    (100, 150, "100_150"),
    (150, 200, "150_200"),
]

JET_RADIUS = 0.4
JET_MATCH_DR = 0.4
JET_MAX_BUILD = 16
JET_MAX_KEEP = 12
JET_MIN_PT = 5.0

DEFAULT_OUTPUT_ACTIVATION = "softplus"

_ML_MODEL = None
_ML_INFO = None


def _round_l1(x):
    return ak.values_astype(x / 0.25 + 0.5, np.int64) * 0.25


def _neutral_threshold(abs_eta):
    #return ak.where(abs_eta < 2.0, 1.0, 2.0)
    #return ak.where(abs_eta < 2.0, 1.0, ak.where(abs_eta < 2.4, 2.0, ak.where(abs_eta < 3.0, 4.0, 10.0)))
    return ak.where(abs_eta < 2.0, 1.0, ak.where(abs_eta < 2.5, 2.0, ak.where(abs_eta < 3.0, 4.0, 10.0)))

def _sum_px_py(pt, phi):
    px = ak.sum(pt * np.cos(phi), axis=1)
    py = ak.sum(pt * np.sin(phi), axis=1)
    return px, py


def _met_from_pt_phi(pt, phi):
    px, py = _sum_px_py(pt, phi)
    met_px = -px
    met_py = -py
    met_pt = np.sqrt(met_px**2 + met_py**2)
    return met_px, met_py, met_pt


def _central_met_from_pt_eta_phi(pt, eta, phi, eta_max=2.4):
    central_pt = ak.where(np.abs(eta) < eta_max, pt, 0.0)
    return _met_from_pt_phi(central_pt, phi)


def _met_resolution(reco_met_pt, gen_met_pt):
    return reco_met_pt - gen_met_pt


def _gen_met_from_gencands(gen, eta_max=2.4):
    mask = ak.ones_like(gen.pt, dtype=np.bool_)
    if DO_GEN_STATUS1 and hasattr(gen, "status"):
        mask = mask & (gen.status == 1)

    mask = mask & (np.abs(gen.eta) < eta_max)

    gen_sel = gen[mask]
    return _met_from_pt_phi(gen_sel.pt, gen_sel.phi)


def _effective_pt_from_weight(puppi, neutral_weight, eta_max=2.4):
    abs_eta = np.abs(puppi.eta)

    is_neutral = puppi.charge == 0
    is_charged = puppi.charge != 0

    neutral_pt = _round_l1(puppi.pt * neutral_weight)
    neutral_thr = _neutral_threshold(abs_eta)

    neutral_keep = is_neutral & (abs_eta < eta_max) & (neutral_pt >= neutral_thr)
    charged_keep = is_charged & (abs_eta < eta_max) & (puppi.nnVtxScore > 0.1)

    pt_for_met = ak.where(is_neutral, neutral_pt, puppi.pt)
    pt_for_met = ak.where(neutral_keep | charged_keep, pt_for_met, 0.0)
    return pt_for_met


def _current_puppi_effective_pt(puppi):
    return _effective_pt_from_weight(puppi, puppi.puppiWeight)


def _gentarget_weight(puppi):
    gen_sum = getattr(puppi, COLLECTION_FIELD_MAP[GENTARGET_GEN_SUM_FIELD])
    reco_sum = getattr(puppi, COLLECTION_FIELD_MAP[GENTARGET_RECO_SUM_FIELD])
    return ak.where(reco_sum > 0, gen_sum / reco_sum, 0.0)


def _gentarget_effective_pt(puppi):
    return _effective_pt_from_weight(puppi, _gentarget_weight(puppi))


def _compute_response_ratio(gen_met_pt, reco_met_pt):
    return ak.where(gen_met_pt > 0, reco_met_pt / gen_met_pt, 0.0)


def _event_singleton_record(mask, **fields):
    return ak.zip(
        {name: ak.singletons(ak.mask(value, mask)) for name, value in fields.items()},
        depth_limit=2,
    )


def _add_binned_met_objects(
    out, prefix,
    gen_met_px, gen_met_py, gen_met_pt,
    reco_met_px, reco_met_py, reco_met_pt,
):
    # One row per event, before histogram bins, clipping, or component masks.
    out[f"{prefix}_events"] = ak.zip({
        "gen_px": gen_met_px, "gen_py": gen_met_py, "gen_pt": gen_met_pt,
        "reco_px": reco_met_px, "reco_py": reco_met_py, "reco_pt": reco_met_pt,
    }, depth_limit=1)
    response_ratio = _compute_response_ratio(gen_met_pt, reco_met_pt)
    resolution = _met_resolution(reco_met_pt, gen_met_pt)
    safe_gen_pt = ak.where(gen_met_pt > 0, gen_met_pt, 1.0)
    gen_ux = gen_met_px / safe_gen_pt
    gen_uy = gen_met_py / safe_gen_pt

    reco_parallel = reco_met_px * gen_ux + reco_met_py * gen_uy
    reco_perpendicular = -reco_met_px * gen_uy + reco_met_py * gen_ux
    parallel_diff = reco_parallel - gen_met_pt

    valid_gen = (
        (gen_met_pt > 0) & np.isfinite(gen_met_pt)
        & np.isfinite(gen_met_px) & np.isfinite(gen_met_py)
        & np.isfinite(reco_met_px) & np.isfinite(reco_met_py)
    )
    reco_parallel = ak.where(valid_gen, reco_parallel, 0.0)
    reco_perpendicular = ak.where(valid_gen, reco_perpendicular, 0.0)
    parallel_diff = ak.where(valid_gen, parallel_diff, 0.0)

    for lo, hi, tag in GENMET_BINS:
        mask = (gen_met_pt >= lo) & (gen_met_pt < hi)

        out[f"{prefix}_genmet_{tag}"] = _event_singleton_record(
            mask,
            gen_met_pt=gen_met_pt,
            reco_met_pt=reco_met_pt,
            response_ratio=response_ratio,
            met_diff=resolution,
            reco_met_parallel=reco_parallel,
            reco_met_perpendicular=ak.mask(reco_perpendicular, valid_gen),
            met_parallel_diff=ak.mask(parallel_diff, valid_gen),
            # GEN perpendicular projection is zero: residual = reco perpendicular.
            # These are unscaled residuals in GeV; widths are computed per event
            # in the plotter with np.std(..., ddof=0).
        )

    return out


def _seeded_cone_jets(pt, eta, phi):
    """Following Federica's notebook: seeded-cone clustering, event by event."""
    pt_events = ak.to_list(pt)
    eta_events = ak.to_list(eta)
    phi_events = ak.to_list(phi)
    output = []

    for event_pt, event_eta, event_phi in zip(pt_events, eta_events, phi_events):
        event_pt = np.asarray(event_pt, dtype=np.float64)
        event_eta = np.asarray(event_eta, dtype=np.float64)
        event_phi = np.asarray(event_phi, dtype=np.float64)

        valid = (
            np.isfinite(event_pt)
            & np.isfinite(event_eta)
            & np.isfinite(event_phi)
            & (event_pt > 0.0)
        )
        event_pt = event_pt[valid]
        event_eta = event_eta[valid]
        event_phi = event_phi[valid]
        remaining = np.ones(len(event_pt), dtype=np.bool_)
        jets = []

        while remaining.any() and len(jets) < JET_MAX_BUILD:
            remaining_idx = np.flatnonzero(remaining)
            seed_idx = remaining_idx[np.argmax(event_pt[remaining_idx])]
            dphi = _wrap_phi(event_phi[seed_idx] - event_phi)
            dr2 = (event_eta[seed_idx] - event_eta) ** 2 + dphi ** 2
            constituents = np.flatnonzero((dr2 <= JET_RADIUS ** 2) & remaining)
            constituent_pt = event_pt[constituents]
            jet_pt = constituent_pt.sum()

            if jet_pt > 0.0:
                jet_eta = np.sum(constituent_pt * event_eta[constituents]) / jet_pt
                jet_phi = np.arctan2(
                    np.sum(constituent_pt * np.sin(event_phi[constituents])),
                    np.sum(constituent_pt * np.cos(event_phi[constituents])),
                )
                jets.append({
                    "pt": float(jet_pt),
                    "eta": float(jet_eta),
                    "phi": float(jet_phi),
                })

            remaining[constituents] = False

        jets.sort(key=lambda jet: jet["pt"], reverse=True)
        output.append([jet for jet in jets[:JET_MAX_KEEP] if jet["pt"] > JET_MIN_PT])

    return ak.zip({
        field: ak.Array([[jet[field] for jet in event] for event in output])
        for field in ("pt", "eta", "phi")
    })


def _match_reco_to_gen_jets(reco_jets, gen_jets):
    """Match each reco jet to its closest GEN jet, as in Federica's notebook."""
    reco_events = ak.to_list(reco_jets)
    gen_events = ak.to_list(gen_jets)
    output = []

    for reco_event, gen_event in zip(reco_events, gen_events):
        matched = []
        for reco in reco_event:
            if not gen_event:
                continue

            gen_eta = np.asarray([jet["eta"] for jet in gen_event], dtype=np.float64)
            gen_phi = np.asarray([jet["phi"] for jet in gen_event], dtype=np.float64)
            deta = float(reco["eta"]) - gen_eta
            dphi = _wrap_phi(float(reco["phi"]) - gen_phi)
            dr = np.sqrt(deta ** 2 + dphi ** 2)
            best = int(np.argmin(dr))

            if dr[best] < JET_MATCH_DR:
                gen = gen_event[best]
                gen_pt = float(gen["pt"])
                reco_pt = float(reco["pt"])
                matched.append({
                    "gen_jet_pt": gen_pt,
                    "reco_jet_pt": reco_pt,
                    "response_ratio": reco_pt / gen_pt if gen_pt > 0 else 0.0,
                    "pt_residual_relative": (
                        (reco_pt - gen_pt) / gen_pt if gen_pt > 0 else 0.0
                    ),
                    "deltaR": float(dr[best]),
                })
        output.append(matched)

    return ak.zip({
        field: ak.Array([
            [record[field] for record in event]
            for event in output
        ])
        for field in (
            "gen_jet_pt", "reco_jet_pt", "response_ratio",
            "pt_residual_relative", "deltaR",
        )
    })


def _add_binned_jet_objects(out, prefix, matches):
    for lo, hi, tag in GENJET_PT_BINS:
        mask = (matches.gen_jet_pt >= lo) & (matches.gen_jet_pt < hi)
        out[f"{prefix}_genjetpt_{tag}"] = matches[mask]
    return out


def _wrap_phi(dphi):
    return (dphi + np.pi) % (2.0 * np.pi) - np.pi


def _load_ml_config():
    # Do not silently evaluate a model with fallback features or radius.
    if not ML_CONFIG_PATH or not os.path.isfile(ML_CONFIG_PATH):
        raise FileNotFoundError(f"ML training config not found: {ML_CONFIG_PATH}")
    with open(ML_CONFIG_PATH) as f:
        cfg = json.load(f)
    metadata = cfg.get("dataset_metadata", {})
    dataset_feature_names = cfg.get("dataset_feature_names")
    if dataset_feature_names is None:
        dataset_feature_names = metadata.get("feature_names")
    input_features = cfg.get("input_features")
    for label, names in (("dataset_feature_names", dataset_feature_names),
                         ("input_features", input_features)):
        if (not isinstance(names, list) or not names
                or not all(isinstance(name, str) for name in names)
                or len(set(names)) != len(names)):
            raise ValueError(f"{label} must be a nonempty list of unique feature names")

    name_to_idx = {name: i for i, name in enumerate(dataset_feature_names)}
    missing = [name for name in input_features if name not in name_to_idx]
    if missing:
        raise ValueError(f"ML input features missing from dataset metadata: {missing}")
    # Build the feature selection from names, in the model's declared order.
    # Never treat the dataset's column order as the model's input order.
    feature_idx = [name_to_idx[name] for name in input_features]
    saved_idx = cfg.get("feature_idx")
    if saved_idx is not None and saved_idx != feature_idx:
        raise ValueError(
            "feature_idx disagrees with input_features/dataset_feature_names; "
            "use the unmodified config.json from the checkpoint's training run"
        )
    in_dim = int(cfg.get("in_dim", len(input_features)))
    if in_dim != len(input_features):
        raise ValueError(f"in_dim={in_dim} does not match {len(input_features)} input features")

    # Training JSON uses dataset_metadata.dr; also accept a top-level dr.
    # dr is needed even when it is not a model input feature.
    dr_value = cfg.get("dr", metadata.get("dr"))
    if dr_value is None:
        raise ValueError("Training config must specify dataset_metadata.dr or dr")
    dr = float(dr_value)
    if not np.isfinite(dr) or dr <= 0:
        raise ValueError(f"Neighborhood dr must be finite and positive, got {dr}")
    if ("dr" in cfg and "dr" in metadata
            and not np.isclose(dr, float(metadata["dr"]), rtol=0, atol=1e-12)):
        raise ValueError("Conflicting dr values in config and dataset_metadata")

    # Use training nmax, never dataset_metadata.max_nodes.
    for key in ("nmax", "hidden_dim"):
        value = cfg.get(key)
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"Training config must specify a positive integer {key}")
    nmax = cfg["nmax"]
    hidden_dim = cfg["hidden_dim"]
    return {
        "dataset_feature_names": dataset_feature_names,
        "input_features": input_features,
        "feature_idx": np.asarray(feature_idx, dtype=np.int64),
        "in_dim": in_dim,
        "hidden_dim": hidden_dim,
        "output_activation": cfg.get("output_activation", DEFAULT_OUTPUT_ACTIVATION),
        "dr": dr,
        "nmax": nmax,
    }



def _apply_output_activation(raw, activation):
    import torch.nn.functional as F

    if activation == "softplus":
        return F.softplus(raw)
    if activation == "relu":
        return F.relu(raw)
    if activation == "identity":
        return raw

    raise ValueError(f"Unknown output activation: {activation}")


def _load_ml_model():
    global _ML_MODEL, _ML_INFO

    if _ML_MODEL is not None:
        return _ML_MODEL, _ML_INFO

    if not ML_MODEL_PATH:
        raise RuntimeError(
            "ML_PUPPI_MODEL is not set. Example: export ML_PUPPI_MODEL=/path/to/model_best.pt"
        )

    import torch
    import torch.nn as nn

    class DeepSetRegressor(nn.Module):
        def __init__(self, in_dim, hidden_dim):
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
                nn.ReLU(),
                nn.Linear(hidden_dim, 1),
            )

        def forward(self, x, mask):
            batch_size, n_nodes, _ = x.shape
            valid = mask.bool()
            x_valid = x[valid]
            h_valid = self.phi(x_valid)

            h = torch.zeros(
                batch_size, n_nodes, h_valid.size(1),
                device=x.device, dtype=h_valid.dtype,
            )
            h[valid] = h_valid
            mask_f = mask.unsqueeze(-1).to(h.dtype)
            counts = mask_f.sum(dim=1).clamp(min=1.0)
            h_mean = h.sum(dim=1) / counts
            return self.rho(h_mean).view(-1)

    info = _load_ml_config()

    model = DeepSetRegressor(
        in_dim=info["in_dim"],
        hidden_dim=info["hidden_dim"],
    )

    try:
        state = torch.load(
            ML_MODEL_PATH, map_location=ML_DEVICE, weights_only=True
        )
    except TypeError:
        state = torch.load(ML_MODEL_PATH, map_location=ML_DEVICE)

    if isinstance(state, dict) and "state_dict" in state:
        state = state["state_dict"]
    elif isinstance(state, dict) and "model_state_dict" in state:
        state = state["model_state_dict"]

    model.load_state_dict(state)
    model.to(ML_DEVICE)
    model.eval()

    _ML_MODEL = model
    _ML_INFO = info

    return _ML_MODEL, _ML_INFO


def _build_base_feature_table(
    pt, eta, phi, charge, pdg, puppi_w, nnvtx, vz, dxy, z0,
    id_prob_pu, id_prob_em, id_prob_pi,
    dataset_feature_names
):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    abs_pdg = np.abs(pdg.astype(np.int32))
    is_neutral = charge == 0

    out = np.zeros(
        (len(pt), len(dataset_feature_names)),
        dtype=np.float32
    )
    name_to_col = {
        name: i for i, name in enumerate(dataset_feature_names)
    }

    def set_col(name, values):
        if name in name_to_col:
            out[:, name_to_col[name]] = values.astype(np.float32)

    set_col("log_pt", np.log1p(pt))
    set_col("slog_px", np.sign(px) * np.log1p(np.abs(px)))
    set_col("slog_py", np.sign(py) * np.log1p(np.abs(py)))
    set_col("abs_eta", np.abs(eta))
    set_col("eta", eta)
    set_col("phi", phi)
    set_col("charge", charge)
    set_col("nnvtx", np.where(is_neutral, -1.0, nnvtx))
    set_col("vz", np.where(is_neutral, -1.0, vz))
    set_col("dxy", np.where(is_neutral, -1.0, dxy))
    set_col(
        "log_dxy",
        np.where(
            is_neutral,
            -1.0,
            np.sign(dxy) * np.log1p(np.abs(dxy))
        )
    )
    set_col("z0", np.where(is_neutral, -1.0, z0))
    set_col("puppiWeight", puppi_w)
    set_col("idProbPu", id_prob_pu)
    set_col("idProbEm", id_prob_em)
    set_col("idProbPi", id_prob_pi)
    set_col("pdg_abs_211", (abs_pdg == 211).astype(np.float32))
    set_col("pdg_abs_130", (abs_pdg == 130).astype(np.float32))
    set_col("pdg_abs_22", (abs_pdg == 22).astype(np.float32))
    set_col("pdg_abs_11", (abs_pdg == 11).astype(np.float32))
    set_col("pdg_abs_13", (abs_pdg == 13).astype(np.float32))

    return out, name_to_col


def _predict_graphs_in_batches(model, info, graph_arrays):
    import torch

    preds = []
    batch_graphs = max(1, ML_INFER_BATCH_GRAPHS)

    for start in range(0, len(graph_arrays), batch_graphs):
        graphs = graph_arrays[start:start + batch_graphs]

        n_graphs = len(graphs)
        nmax = info["nmax"]
        in_dim = info["in_dim"]
        x_np = np.zeros((n_graphs, nmax, in_dim), dtype=np.float32)
        mask_np = np.zeros((n_graphs, nmax), dtype=np.float32)

        for i, graph in enumerate(graphs):
            n_nodes = min(len(graph), nmax)
            x_np[i, :n_nodes] = graph[:n_nodes]
            mask_np[i, :n_nodes] = 1.0

        x = torch.from_numpy(x_np).to(ML_DEVICE)
        mask = torch.from_numpy(mask_np).to(ML_DEVICE)

        with torch.no_grad():
            raw = model(x, mask)
            pred = _apply_output_activation(
                raw, info["output_activation"]
            )

        preds.append(
            pred.detach().cpu().numpy().reshape(-1)
        )

    if preds:
        return np.concatenate(preds, axis=0)

    return np.array([], dtype=np.float32)


def _ml_neutral_weight_array(puppi):
    model, info = _load_ml_model()

    dataset_feature_names = info["dataset_feature_names"]
    feature_idx = info["feature_idx"]
    dr_cut = float(info["dr"])
    nmax = int(info["nmax"])

    feature_branches = {
        name: COLLECTION_FIELD_MAP[name]
        for name in (
            "nnvtx", "vz", "dxy", "z0", "puppiWeight",
            "idProbPu", "idProbEm", "idProbPi",
        )
    }
    missing_branches = [
        branch
        for feature, branch in feature_branches.items()
        if feature in info["input_features"] and not hasattr(puppi, branch)
    ]
    if missing_branches:
        raise RuntimeError(
            "L1PuppiCands is missing branches required by this model: "
            f"{missing_branches}"
        )

    pt_list = ak.to_list(puppi.pt)
    eta_list = ak.to_list(puppi.eta)
    phi_list = ak.to_list(puppi.phi)
    charge_list = ak.to_list(puppi.charge)

    pdg_list = (
        ak.to_list(puppi.pdgId)
        if hasattr(puppi, "pdgId")
        else None
    )
    puppi_w_list = (
        ak.to_list(puppi.puppiWeight)
        if hasattr(puppi, "puppiWeight")
        else None
    )
    nnvtx_list = (
        ak.to_list(puppi.nnVtxScore)
        if hasattr(puppi, "nnVtxScore")
        else None
    )
    vz_list = (
        ak.to_list(puppi.vz)
        if hasattr(puppi, "vz")
        else None
    )
    dxy_list = (
        ak.to_list(puppi.dxy)
        if hasattr(puppi, "dxy")
        else None
    )
    z0_list = (
        ak.to_list(puppi.z0)
        if hasattr(puppi, "z0")
        else None
    )
    id_prob_pu_list = (
        ak.to_list(getattr(puppi, COLLECTION_FIELD_MAP["idProbPu"]))
        if hasattr(puppi, COLLECTION_FIELD_MAP["idProbPu"])
        else None
    )
    id_prob_em_list = (
        ak.to_list(getattr(puppi, COLLECTION_FIELD_MAP["idProbEm"]))
        if hasattr(puppi, COLLECTION_FIELD_MAP["idProbEm"])
        else None
    )
    id_prob_pi_list = (
        ak.to_list(getattr(puppi, COLLECTION_FIELD_MAP["idProbPi"]))
        if hasattr(puppi, COLLECTION_FIELD_MAP["idProbPi"])
        else None
    )

    graph_arrays = []
    refs = []
    weight_events = []

    for iev in range(len(pt_list)):
        n = len(pt_list[iev])

        if n == 0:
            weight_events.append(
                np.asarray([], dtype=np.float32)
            )
            continue

        pt = np.asarray(pt_list[iev], dtype=np.float32)
        eta = np.asarray(eta_list[iev], dtype=np.float32)
        phi = np.asarray(phi_list[iev], dtype=np.float32)
        charge = np.asarray(charge_list[iev], dtype=np.float32)

        pdg = (
            np.asarray(pdg_list[iev], dtype=np.int32)
            if pdg_list is not None
            else np.zeros(n, dtype=np.int32)
        )
        puppi_w = (
            np.asarray(puppi_w_list[iev], dtype=np.float32)
            if puppi_w_list is not None
            else np.ones(n, dtype=np.float32)
        )
        nnvtx = (
            np.asarray(nnvtx_list[iev], dtype=np.float32)
            if nnvtx_list is not None
            else np.full(n, -1.0, dtype=np.float32)
        )
        vz = (
            np.asarray(vz_list[iev], dtype=np.float32)
            if vz_list is not None
            else np.full(n, -1.0, dtype=np.float32)
        )
        dxy = (
            np.asarray(dxy_list[iev], dtype=np.float32)
            if dxy_list is not None
            else np.full(n, -1.0, dtype=np.float32)
        )
        z0 = (
            np.asarray(z0_list[iev], dtype=np.float32)
            if z0_list is not None
            else np.full(n, -1.0, dtype=np.float32)
        )
        id_prob_pu = (
            np.asarray(id_prob_pu_list[iev], dtype=np.float32)
            if id_prob_pu_list is not None
            else np.zeros(n, dtype=np.float32)
        )
        id_prob_em = (
            np.asarray(id_prob_em_list[iev], dtype=np.float32)
            if id_prob_em_list is not None
            else np.zeros(n, dtype=np.float32)
        )
        id_prob_pi = (
            np.asarray(id_prob_pi_list[iev], dtype=np.float32)
            if id_prob_pi_list is not None
            else np.zeros(n, dtype=np.float32)
        )

        event_weights = puppi_w.copy()

        base, name_to_col = _build_base_feature_table(
            pt, eta, phi, charge, pdg,
            puppi_w, nnvtx, vz, dxy, z0,
            id_prob_pu, id_prob_em, id_prob_pi,
            dataset_feature_names,
        )

        deta_all = eta[:, None] - eta[None, :]
        dphi_all = _wrap_phi(
            phi[:, None] - phi[None, :]
        )
        dr_all = np.sqrt(
            deta_all * deta_all + dphi_all * dphi_all
        )

        neutral_centers = np.where(
            (charge == 0) & (np.abs(eta) < 2.5)
        )[0]

        for center_idx in neutral_centers:
            dist = dr_all[:, center_idx]
            idx = np.nonzero(dist < dr_cut)[0]

            if len(idx) == 0:
                continue

            # Match BestPuppiDataset.__getitem__: center first, followed by
            # other particles in their original graph order, then truncate.
            other_idx = idx[idx != center_idx]
            idx = np.concatenate(([center_idx], other_idx))[:nmax]

            x_full = base[idx].copy()

            if "deta" in name_to_col:
                x_full[:, name_to_col["deta"]] = (
                    eta[idx] - eta[center_idx]
                )

            if "dphi" in name_to_col:
                x_full[:, name_to_col["dphi"]] = _wrap_phi(
                    phi[idx] - phi[center_idx]
                )

            if "dr" in name_to_col:
                x_full[:, name_to_col["dr"]] = dist[idx]

            if "is_center" in name_to_col:
                x_full[:, name_to_col["is_center"]] = 0.0
                center_pos = np.where(idx == center_idx)[0]

                if len(center_pos) == 0:
                    continue

                x_full[
                    center_pos[0],
                    name_to_col["is_center"]
                ] = 1.0

            graph_arrays.append(
                x_full[:, feature_idx]
            )
            refs.append((iev, int(center_idx)))

        weight_events.append(event_weights)

    if len(graph_arrays) > 0:
        pred = _predict_graphs_in_batches(
            model, info, graph_arrays
        )

        for value, (iev, icand) in zip(pred, refs):
            weight_events[iev][icand] = float(value)

    return ak.Array(weight_events)


def _ml_puppi_effective_pt(puppi):
    ml_weight = _ml_neutral_weight_array(puppi)
    return _effective_pt_from_weight(puppi, ml_weight)


def build_objects(events):
    return {
        "puppi": getattr(events, PUPPI_COLLECTION),
        "gen": getattr(events, GEN_COLLECTION),
        "gen_jets": getattr(events, GENJET_COLLECTION),
    }


def cut_met_response(e, o):
    out = dict(o)

    puppi = out["puppi"]
    gen = out["gen"]
    gen_jets = out["gen_jets"]
    stored_gen_px, stored_gen_py, stored_gen_pt = _stored_genmet(e)

    # GEN MET: central |eta| < 2.4, barrel |eta| < 1.4
    gen_central_met_px, gen_central_met_py, gen_central_met_pt = (
        _gen_met_from_gencands(gen, eta_max=2.4)
    )
    gen_barrel_met_px, gen_barrel_met_py, gen_barrel_met_pt = (
        _gen_met_from_gencands(gen, eta_max=1.4)
    )

    # Existing plain MET objects are aliases of central |eta| < 2.4.
    gen_met_px = gen_central_met_px
    gen_met_py = gen_central_met_py
    gen_met_pt = gen_central_met_pt

    # Handmade PUPPI
    puppi_eff_pt = _current_puppi_effective_pt(puppi)

    (
        puppi_central_met_px,
        puppi_central_met_py,
        puppi_central_met_pt,
    ) = _central_met_from_pt_eta_phi(
        puppi_eff_pt, puppi.eta, puppi.phi, eta_max=2.4
    )

    (
        puppi_barrel_met_px,
        puppi_barrel_met_py,
        puppi_barrel_met_pt,
    ) = _central_met_from_pt_eta_phi(
        puppi_eff_pt, puppi.eta, puppi.phi, eta_max=1.4
    )

    puppi_met_px = puppi_central_met_px
    puppi_met_py = puppi_central_met_py
    puppi_met_pt = puppi_central_met_pt

    # GenTarget PUPPI
    gentarget_weight = _gentarget_weight(puppi)
    gentarget_eff_pt = _effective_pt_from_weight(puppi, gentarget_weight)

    (
        gentarget_central_met_px,
        gentarget_central_met_py,
        gentarget_central_met_pt,
    ) = _central_met_from_pt_eta_phi(
        gentarget_eff_pt, puppi.eta, puppi.phi, eta_max=2.4
    )

    (
        gentarget_barrel_met_px,
        gentarget_barrel_met_py,
        gentarget_barrel_met_pt,
    ) = _central_met_from_pt_eta_phi(
        gentarget_eff_pt, puppi.eta, puppi.phi, eta_max=1.4
    )

    gentarget_met_px = gentarget_central_met_px
    gentarget_met_py = gentarget_central_met_py
    gentarget_met_pt = gentarget_central_met_pt

    # ML PUPPI
    ml_weight = _ml_neutral_weight_array(puppi)
    ml_eff_pt = _effective_pt_from_weight(puppi, ml_weight)

    (
        ml_central_met_px,
        ml_central_met_py,
        ml_central_met_pt,
    ) = _central_met_from_pt_eta_phi(
        ml_eff_pt, puppi.eta, puppi.phi, eta_max=2.4
    )

    (
        ml_barrel_met_px,
        ml_barrel_met_py,
        ml_barrel_met_pt,
    ) = _central_met_from_pt_eta_phi(
        ml_eff_pt, puppi.eta, puppi.phi, eta_max=1.4
    )

    ml_met_px = ml_central_met_px
    ml_met_py = ml_central_met_py
    ml_met_pt = ml_central_met_pt

    # Central |eta| < 2.4
    out = _add_binned_met_objects(
        out, "puppi_central_met_perf",
        gen_central_met_px, gen_central_met_py, gen_central_met_pt,
        puppi_central_met_px, puppi_central_met_py, puppi_central_met_pt,
    )
    out = _add_binned_met_objects(
        out, "gentarget_central_met_perf",
        gen_central_met_px, gen_central_met_py, gen_central_met_pt,
        gentarget_central_met_px, gentarget_central_met_py, gentarget_central_met_pt,
    )
    out = _add_binned_met_objects(
        out, "ml_central_met_perf",
        gen_central_met_px, gen_central_met_py, gen_central_met_pt,
        ml_central_met_px, ml_central_met_py, ml_central_met_pt,
    )

    # Barrel |eta| < 1.4
    out = _add_binned_met_objects(
        out, "puppi_barrel_met_perf",
        gen_barrel_met_px, gen_barrel_met_py, gen_barrel_met_pt,
        puppi_barrel_met_px, puppi_barrel_met_py, puppi_barrel_met_pt,
    )
    out = _add_binned_met_objects(
        out, "gentarget_barrel_met_perf",
        gen_barrel_met_px, gen_barrel_met_py, gen_barrel_met_pt,
        gentarget_barrel_met_px, gentarget_barrel_met_py, gentarget_barrel_met_pt,
    )
    out = _add_binned_met_objects(
        out, "ml_barrel_met_perf",
        gen_barrel_met_px, gen_barrel_met_py, gen_barrel_met_pt,
        ml_barrel_met_px, ml_barrel_met_py, ml_barrel_met_pt,
    )

    # Existing plain names: aliases of central |eta| < 2.4
    out = _add_binned_met_objects(
        out, "puppi_met_perf",
        gen_met_px, gen_met_py, gen_met_pt,
        puppi_met_px, puppi_met_py, puppi_met_pt,
    )
    out = _add_binned_met_objects(
        out, "gentarget_met_perf",
        gen_met_px, gen_met_py, gen_met_pt,
        gentarget_met_px, gentarget_met_py, gentarget_met_pt,
    )
    out = _add_binned_met_objects(
        out, "ml_met_perf",
        gen_met_px, gen_met_py, gen_met_pt,
        ml_met_px, ml_met_py, ml_met_pt,
    )

    # Keep the stored-reference comparison separate from GenCands-based MET.
    for method, reco_px, reco_py, reco_pt in (
        ("puppi", puppi_central_met_px, puppi_central_met_py, puppi_central_met_pt),
        ("gentarget", gentarget_central_met_px, gentarget_central_met_py, gentarget_central_met_pt),
        ("ml", ml_central_met_px, ml_central_met_py, ml_central_met_pt),
    ):
        _add_binned_met_objects(
            out, f"{method}_storedgen_met_perf",
            stored_gen_px, stored_gen_py, stored_gen_pt, reco_px, reco_py, reco_pt,
        )

    # Notebook-equivalent seeded-cone jets. MET remains central |eta| < 2.4;
    # jet inputs follow the notebook candidate acceptance through |eta| < 2.5.
    puppi_jet_pt = _effective_pt_from_weight(
        puppi, puppi.puppiWeight, eta_max=2.5
    )
    gentarget_jet_pt = _effective_pt_from_weight(
        puppi, gentarget_weight, eta_max=2.5
    )
    ml_jet_pt = _effective_pt_from_weight(
        puppi, ml_weight, eta_max=2.5
    )

    reco_jet_collections = {
        "puppi_jet_perf": _seeded_cone_jets(
            puppi_jet_pt, puppi.eta, puppi.phi
        ),
        "gentarget_jet_perf": _seeded_cone_jets(
            gentarget_jet_pt, puppi.eta, puppi.phi
        ),
        "ml_jet_perf": _seeded_cone_jets(
            ml_jet_pt, puppi.eta, puppi.phi
        ),
    }

    for prefix, reco_jets in reco_jet_collections.items():
        matches = _match_reco_to_gen_jets(reco_jets, gen_jets)
        out = _add_binned_jet_objects(out, prefix, matches)

    return out


CUTFLOW = [
    ("met_response", [cut_met_response]),
]
