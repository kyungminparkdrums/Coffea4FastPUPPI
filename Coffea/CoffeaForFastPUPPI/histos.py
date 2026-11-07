# Coffea4FP/histos.py
import awkward as ak
import numpy as np


def _has_field(coll, field):
    return (coll is not None) and hasattr(coll, field)


def _to_numpy_1d(x, *, default=0):
    if x is None:
        return np.array([])

    # fully flatten (safe for jagged / option / union)
    x = ak.flatten(x, axis=None)

    # replace None
    x = ak.fill_none(x, default)

    # convert to numpy
    arr = ak.to_numpy(x, allow_missing=False)

    # bool → int (boost_histogram hates bools)
    if arr.dtype == np.bool_:
        arr = arr.astype(np.int8)

    return arr


def _ratio_per_event(numer_coll, denom_coll, field="pt"):
    """
    Per-event ratio: ak.num(numer[field]) / ak.num(denom[field])
    Returns numpy 1D (nEvents,)
    """
    if numer_coll is None or denom_coll is None:
        return np.array([])
    if not _has_field(numer_coll, field) or not _has_field(denom_coll, field):
        return np.array([])

    denom = ak.num(getattr(denom_coll, field), axis=1)
    numer = ak.num(getattr(numer_coll, field), axis=1)

    denom = ak.where(denom == 0, 1, denom)
    return ak.to_numpy(numer / denom)


def _resolution(matched_reco, matched_true, field="pt"):
    """
    Per-object resolution: matched_reco.pt / matched_true.pt (flattened).
    Returns numpy 1D
    """
    if matched_reco is None or matched_true is None:
        return np.array([])
    if not _has_field(matched_reco, field) or not _has_field(matched_true, field):
        return np.array([])

    res = getattr(matched_reco, field) / getattr(matched_true, field)
    return _to_numpy_1d(res, default=0)


def fill_histo(selected_obj, histos, out):
    """
    selected_obj: dict[str, awkward array/records] (collections)
    histos: dict of expanded hist config entries
    out: dict-like accumulator of Hist objects
    """

    # --- Derived per-event quantities ---
    isMatched = {}        # purity proxy: matched_reco/reco
    isReconstructed = {}  # efficiency proxy: matched_gen/gen
    resolution = {}       # per-object pt ratio for matched reco

    # purity proxies (reco side)
    if "pf" in selected_obj and "matched_pf" in selected_obj:
        isMatched["pf"] = _ratio_per_event(selected_obj["matched_pf"], selected_obj["pf"], field="pt")
    if "puppi" in selected_obj and "matched_puppi" in selected_obj:
        isMatched["puppi"] = _ratio_per_event(selected_obj["matched_puppi"], selected_obj["puppi"], field="pt")
    if "tkele" in selected_obj and "matched_tkele" in selected_obj:
        isMatched["tkele"] = _ratio_per_event(selected_obj["matched_tkele"], selected_obj["tkele"], field="pt")

    # efficiency proxies (gen side)
    if "gen" in selected_obj and "matched_gen" in selected_obj:
        isReconstructed["gen"] = _ratio_per_event(selected_obj["matched_gen"], selected_obj["gen"], field="pt")
    if "genel" in selected_obj and "matched_genel" in selected_obj:
        isReconstructed["genel"] = _ratio_per_event(selected_obj["matched_genel"], selected_obj["genel"], field="pt")

    # resolutions (matched reco vs its matched truth)
    if "matched_pf" in selected_obj and "matched_pfTrue" in selected_obj:
        resolution["matched_pf"] = _resolution(selected_obj["matched_pf"], selected_obj["matched_pfTrue"], field="pt")
    if "matched_puppi" in selected_obj and "matched_puppiTrue" in selected_obj:
        resolution["matched_puppi"] = _resolution(selected_obj["matched_puppi"], selected_obj["matched_puppiTrue"], field="pt")
    if "matched_tkele" in selected_obj and "matched_tkeleTrue" in selected_obj:
        resolution["matched_tkele"] = _resolution(selected_obj["matched_tkele"], selected_obj["matched_tkeleTrue"], field="pt")

    # --- Fill loop ---
    for hname, hcfg in histos.items():
        variables = hcfg["variables"]

        obj_key = hcfg.get("object", None)
        if obj_key is None:
            # fallback heuristic (shouldn't happen if YAML expanded properly)
            parts = hname.split("_")
            obj_key = parts[0] if "matched" not in hname else "_".join(parts[:2])

        if obj_key not in selected_obj:
            continue

        coll = selected_obj[obj_key]

        data = {}
        for var in variables:
            # derived variables
            if var == "multiplicity":
                if _has_field(coll, "pt"):
                    data[var] = ak.to_numpy(ak.num(coll.pt, axis=1))
                else:
                    data[var] = np.array([])

            elif var == "pdgId":
                if _has_field(coll, "pdgId"):
                    data[var] = _to_numpy_1d(np.abs(coll.pdgId), default=0)
                else:
                    data[var] = np.array([])

            elif var == "isMatched":
                data[var] = isMatched.get(obj_key, np.array([]))

            elif var == "isReconstructed":
                data[var] = isReconstructed.get(obj_key, np.array([]))

            elif var == "resolution":
                data[var] = resolution.get(obj_key, np.array([]))

            else:
                # direct field
                if _has_field(coll, var):
                    data[var] = _to_numpy_1d(getattr(coll, var), default=0)
                else:
                    data[var] = np.array([])

        # skip empty
        if any(len(v) == 0 for v in data.values()):
            continue

        # sanity: if 2D, require same lengths
        if len(data.keys()) == 2:
            lens = [len(v) for v in data.values()]
            if len(set(lens)) != 1:
                continue

        out[hname].fill(**data)

    return out
