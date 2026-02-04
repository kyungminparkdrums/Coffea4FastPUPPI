import awkward as ak
import numpy as np

from utils import utils
from cut_config import ETA_RANGES, cut_range, apply_evt_mask

# ============================================================
# USER CONFIG (only edit here)
# ============================================================

MODE = "tkele"
ETA_REGION = "barrel"
PT_MIN = 1.0
DR_CUT = 0.1

APPLY_GEN_TWO_PROMPT_OS = True  # signal-only; safely skipped if no gen

PAIR_DVZ_MAX = 1.0
PAIR_MAX_RELPFISO = 1   # set None to disable
PAIR_MAX_RELPUISO = 0.5    # set None to disable
BESTPAIR_SCORE = "min_dvz"     # "pt" / "min_dvz" / "min_iso_pf"
BESTPAIR_PTMIN = 10     # e.g. 10.0


# ============================================================
# Helpers
# ============================================================

def has_gen(events):
    return hasattr(events, "GenEl")

def build_objects(events):
    obj = {
        "tkele": events.TkEleL2,
    }
    # Normalize pt to ptCorr for TkEleL2
    if hasattr(obj["tkele"], "ptCorr"):
        obj["tkele"] = ak.with_field(obj["tkele"], obj["tkele"].ptCorr, "pt")

    obj["tkele"] = add_rel_iso(obj["tkele"])

    if has_gen(events):
        obj["genel"] = events.GenEl

    return obj


def add_rel_iso(coll):
    if coll is None:
        return None
    if not hasattr(coll, "pt"):
        return coll

    pt = coll.pt
    safe_pt = ak.where(pt > 0, pt, 1.0)

    if hasattr(coll, "pfIso"):
        coll = ak.with_field(coll, coll.pfIso / safe_pt, "relPfIso")
    if hasattr(coll, "puppiIso"):
        coll = ak.with_field(coll, coll.puppiIso / safe_pt, "relPuppiIso")

    return coll


def _delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    return (dphi + np.pi) % (2 * np.pi) - np.pi


def build_pairs_from_electrons(ele):
    """
    Generic: build ee pairs from an electron-like collection (pt,eta,phi,(vz optional)).
    Output has: l1,l2 and pair vars pt,eta,phi,vz,mass,delta_*.
    """
    pairs = ak.combinations(ele, 2, fields=["l1", "l2"])
    l1, l2 = pairs.l1, pairs.l2

    pt1, eta1, phi1 = l1.pt, l1.eta, l1.phi
    pt2, eta2, phi2 = l2.pt, l2.eta, l2.phi

    d_eta = eta1 - eta2
    d_phi = _delta_phi(phi1, phi2)

    m2 = 2.0 * pt1 * pt2 * (np.cosh(d_eta) - np.cos(d_phi))
    m2 = ak.where(m2 < 0, 0, m2)
    mass = np.sqrt(m2)

    px = pt1 * np.cos(phi1) + pt2 * np.cos(phi2)
    py = pt1 * np.sin(phi1) + pt2 * np.sin(phi2)
    pt = np.sqrt(px * px + py * py)
    phi = np.arctan2(py, px)

    pz = pt1 * np.sinh(eta1) + pt2 * np.sinh(eta2)
    eta = np.arcsinh(ak.where(pt == 0, 0, pz / pt))

    if hasattr(l1, "vz") and hasattr(l2, "vz"):
        vz = 0.5 * (l1.vz + l2.vz)
        d_vz = l1.vz - l2.vz
    else:
        vz = ak.zeros_like(pt)
        d_vz = ak.zeros_like(pt)

    d_pt = pt1 - pt2

    pairs = ak.with_field(pairs, pt, "pt")
    pairs = ak.with_field(pairs, eta, "eta")
    pairs = ak.with_field(pairs, phi, "phi")
    pairs = ak.with_field(pairs, vz, "vz")
    pairs = ak.with_field(pairs, mass, "mass")

    pairs = ak.with_field(pairs, abs(d_eta), "delta_eta")
    pairs = ak.with_field(pairs, abs(d_phi), "delta_phi")
    pairs = ak.with_field(pairs, abs(d_vz), "delta_vz")
    pairs = ak.with_field(pairs, abs(d_pt), "delta_pt")

    # normalized variables
    eps = 1e-12
    pt_safe  = ak.where(pt == 0, eps, pt)
    eta_safe = ak.where(eta == 0, eps, abs(eta))
    phi_safe = ak.where(phi == 0, eps, abs(phi))

    pairs = ak.with_field(pairs, abs(d_pt)  / pt_safe,  "delta_pt_over_pt")
    pairs = ak.with_field(pairs, abs(d_eta) / eta_safe, "delta_eta_over_eta")
    pairs = ak.with_field(pairs, abs(d_phi) / phi_safe, "delta_phi_over_phi")

    # pt/m
    pt_over_mass = ak.where(mass > 0, pt / mass, 0.0)
    pairs = ak.with_field(pairs, pt_over_mass, "pt_over_mass")

    if hasattr(l1, "charge") and hasattr(l2, "charge"):
        charge_prod = l1.charge * l2.charge
        pairs = ak.with_field(pairs, charge_prod, "charge_prod")

    eta_prod = l1.eta * l2.eta
    pairs = ak.with_field(pairs, eta_prod, "eta_prod")

    # --- max relIso per pair ---
    if hasattr(l1, "relPfIso") and hasattr(l2, "relPfIso"):
        pairs = ak.with_field(pairs, np.maximum(l1.relPfIso, l2.relPfIso), "max_relPfIso")

    if hasattr(l1, "relPuppiIso") and hasattr(l2, "relPuppiIso"):
        pairs = ak.with_field(pairs, np.maximum(l1.relPuppiIso, l2.relPuppiIso), "max_relPuppiIso")

    d_r = np.sqrt(d_eta*d_eta + d_phi*d_phi)
    pairs = ak.with_field(pairs, d_r, "delta_r")

    return pairs


def build_lead_sub_from_pairs(pair_coll):
    """
    From pair_coll.l1/l2 → returns (lead, sub) jagged [events][pairs]
    """
    if pair_coll is None or (not hasattr(pair_coll, "l1")) or (not hasattr(pair_coll, "l2")):
        return None, None
    l1, l2 = pair_coll.l1, pair_coll.l2
    if not (hasattr(l1, "pt") and hasattr(l2, "pt")):
        return None, None
    lead_mask = (l1.pt >= l2.pt)
    lead = ak.where(lead_mask, l1, l2)
    sub  = ak.where(lead_mask, l2, l1)
    return lead, sub


def _safe_abs(x):
    return abs(x)

def add_pair_level_iso(pairs):
    """
    Requires electrons already have relPfIso and/or relPuppiIso.
    Adds:
      max_relPfIso, max_relPuppiIso per pair
    """
    if pairs is None or (not hasattr(pairs, "l1")) or (not hasattr(pairs, "l2")):
        return pairs

    l1, l2 = pairs.l1, pairs.l2

    if hasattr(l1, "relPfIso") and hasattr(l2, "relPfIso"):
        pairs = ak.with_field(pairs, np.maximum(l1.relPfIso, l2.relPfIso), "max_relPfIso")
    if hasattr(l1, "relPuppiIso") and hasattr(l2, "relPuppiIso"):
        pairs = ak.with_field(pairs, np.maximum(l1.relPuppiIso, l2.relPuppiIso), "max_relPuppiIso")

    return pairs


def _as_jagged01(per_event_option):
    # [evt] option-record -> [evt][0/1]
    return ak.fill_none(ak.singletons(per_event_option), [])



# ============================================================
# CUTS
# ============================================================

def cut_base(events, obj):
    return obj


def cut_evt_gen_two_prompt_os(events, obj):
    """
    Signal-only event selection. Safe for bkg (no GenEl).
    """
    if not APPLY_GEN_TWO_PROMPT_OS:
        return obj
    if not has_gen(events):
        return obj

    mask_genPrompt = (events.GenEl.prompt == 2)
    evt_mask = (ak.num(events.GenEl.pt[mask_genPrompt], axis=1) == 2)
    evt_mask = evt_mask & (ak.sum(events.GenEl.charge[mask_genPrompt], axis=1) == 0)
    evt_mask = evt_mask & (ak.num(events.GenEl.pt, axis=1) == 2)

    return apply_evt_mask(obj, evt_mask)


# ---- linked reco/true cutting helpers ----

MATCHED_RECO_TRUE_PAIRS = [("matched_tkele", "matched_tkeleTrue")]

def _cut_range_linked_reco_true(obj, reco_key, true_key, var, vmin=None, vmax=None, doAbs=False):
    if reco_key not in obj or true_key not in obj:
        return obj
    reco = obj[reco_key]
    true = obj[true_key]
    if reco is None or true is None or (not hasattr(reco, var)):
        return obj

    vals = getattr(reco, var)
    if doAbs:
        vals = abs(vals)

    mask = ak.ones_like(vals, dtype=bool)
    if vmin is not None:
        mask = mask & (vals >= vmin)
    if vmax is not None:
        mask = mask & (vals < vmax)

    out = dict(obj)
    out[reco_key] = reco[mask]
    out[true_key] = true[mask]
    return out


def _apply_base_cut_range(events, obj, *, var, vmin=None, vmax=None, doAbs=False):
    out = obj
    for key in ("genel", "tkele"):
        if key in out:
            out = cut_range(events, out, key, var, vmin=vmin, vmax=vmax, doAbs=doAbs)
    return out


def _apply_linked_matched_cuts(obj, *, var, vmin=None, vmax=None, doAbs=False):
    out = obj
    for reco_key, true_key in MATCHED_RECO_TRUE_PAIRS:
        out = _cut_range_linked_reco_true(out, reco_key, true_key, var, vmin=vmin, vmax=vmax, doAbs=doAbs)
    return out


def cut_pt(events, obj):
    out = _apply_base_cut_range(events, obj, var="pt", vmin=PT_MIN, vmax=None, doAbs=False)
    out = _apply_linked_matched_cuts(out, var="pt", vmin=PT_MIN, vmax=None, doAbs=False)

    if "tkele" in out:
        out["tkele"] = add_rel_iso(out["tkele"])
    if "matched_tkele" in out:
        out["matched_tkele"] = add_rel_iso(out["matched_tkele"])
    if "nonMatched_tkele" in out:
        out["nonMatched_tkele"] = add_rel_iso(out["nonMatched_tkele"])

    return out


def cut_eta(events, obj):
    eta_min, eta_max = ETA_RANGES[ETA_REGION]
    out = _apply_base_cut_range(events, obj, var="eta", vmin=eta_min, vmax=eta_max, doAbs=True)
    out = _apply_linked_matched_cuts(out, var="eta", vmin=eta_min, vmax=eta_max, doAbs=True)

    if "tkele" in out:
        out["tkele"] = add_rel_iso(out["tkele"])
    if "matched_tkele" in out:
        out["matched_tkele"] = add_rel_iso(out["matched_tkele"])
    if "nonMatched_tkele" in out:
        out["nonMatched_tkele"] = add_rel_iso(out["nonMatched_tkele"])

    return out


def cut_build_pairs(events, obj):
    """
    Build reco pairs + reco lead/sub always.
    Also build gen pairs + gen lead/sub when gen exists.
    """
    out = dict(obj)

    # reco
    if "tkele" in out:
        out["tkelePair"] = build_pairs_from_electrons(out["tkele"])
        out["tkelePair"] = add_pair_level_iso(out["tkelePair"])
        out["tkeleLead"], out["tkeleSub"] = build_lead_sub_from_pairs(out["tkelePair"])

    # gen (signal only)
    if "genel" in out:
        out["genelPair"] = build_pairs_from_electrons(out["genel"])
        out["genelLead"], out["genelSub"] = build_lead_sub_from_pairs(out["genelPair"])

    return out


def cut_add_matching(events, obj):
    """
    Signal-only gen matching. Safe no-op for bkg.
    """
    out = dict(obj)
    if not has_gen(events):
        return out
    if "genel" not in out or "tkele" not in out:
        return out

    # gen->reco efficiency
    out["matched_genel"], out["nonMatched_genel"] = utils.get_genMatched(
        out["genel"], out["tkele"], typ="Gen", dr_cut=DR_CUT
    )

    matched_mask, matched_gen_idx = utils.match_reco_to_gen_indices(
        out["genel"], out["tkele"], dr_cut=DR_CUT
    )

    # keep only highest-pt reco per gen
    best_mask = utils.keep_highest_pt_reco_per_gen(
        out["tkele"], matched_mask, matched_gen_idx, pt_field="pt"
    )
    matched_mask = matched_mask & best_mask

    # decorate reco electrons
    tkele = out["tkele"]
    tkele = ak.with_field(tkele, matched_mask, "isGenMatched")
    tkele = ak.with_field(tkele, matched_gen_idx, "matchedGenIdx")
    out["tkele"] = tkele

    # matched / nonMatched reco
    out["matched_tkele"] = tkele[matched_mask]
    out["nonMatched_tkele"] = tkele[~matched_mask]

    # matched truth for reco electrons
    safe_idx = ak.where(matched_gen_idx < 0, 0, matched_gen_idx)
    matched_true_all = out["genel"][safe_idx]          # aligned to reco shape
    out["matched_tkeleTrue"] = matched_true_all[matched_mask]

    # ------------------------------------------------------------
    # rebuild pairs from the UPDATED tkele
    # ------------------------------------------------------------
    pair_all = build_pairs_from_electrons(out["tkele"])
    out["tkelePair_genMatchAware"] = pair_all  # optional debug; can omit

    l1 = pair_all.l1
    l2 = pair_all.l2

    both_matched = (l1.isGenMatched & l2.isGenMatched)
    distinct = (l1.matchedGenIdx != l2.matchedGenIdx)
    pair_matched_mask = both_matched & distinct

    out["matched_tkelePair"] = pair_all[pair_matched_mask]
    out["nonMatched_tkelePair"] = pair_all[~pair_matched_mask]

    # lead/sub from matched/nonmatched pairs
    mlead, msub = build_lead_sub_from_pairs(out["matched_tkelePair"])
    if mlead is not None:
        out["matched_tkeleLead"] = mlead
    if msub is not None:
        out["matched_tkeleSub"] = msub

    nlead, nsub = build_lead_sub_from_pairs(out["nonMatched_tkelePair"])
    if nlead is not None:
        out["nonMatched_tkeleLead"] = nlead
    if nsub is not None:
        out["nonMatched_tkeleSub"] = nsub

    out["matched_tkele"] = add_rel_iso(out["matched_tkele"])
    out["nonMatched_tkele"] = add_rel_iso(out["nonMatched_tkele"])

    if "matched_tkeleLead" in out:
        out["matched_tkeleLead"] = add_rel_iso(out["matched_tkeleLead"])
    if "matched_tkeleSub" in out:
        out["matched_tkeleSub"] = add_rel_iso(out["matched_tkeleSub"])
    if "nonMatched_tkeleLead" in out:
        out["nonMatched_tkeleLead"] = add_rel_iso(out["nonMatched_tkeleLead"])
    if "nonMatched_tkeleSub" in out:
        out["nonMatched_tkeleSub"] = add_rel_iso(out["nonMatched_tkeleSub"])


    return out


def cut_pair_os(events, obj, pair_key="tkelePair"):
    out = dict(obj)
    if pair_key not in out:
        return out
    p = out[pair_key]
    if p is None or not hasattr(p, "charge_prod"):
        return out

    out[pair_key] = p[p.charge_prod < 0]
    out = cut_refresh_pair_match_views(out, pair_key=pair_key)
    out = cut_pick_best_pair(out, pair_key=pair_key, score=BESTPAIR_SCORE)
    return out


def cut_pair_dvz(events, obj, pair_key="tkelePair", dvz_max=1.0):
    out = dict(obj)
    if pair_key not in out:
        return out
    p = out[pair_key]
    if p is None or not hasattr(p, "delta_vz"):
        return out

    out[pair_key] = p[p.delta_vz < dvz_max]
    out = cut_refresh_pair_match_views(out, pair_key=pair_key)
    out = cut_pick_best_pair(out, pair_key=pair_key, score=BESTPAIR_SCORE)
    return out


def cut_pair_iso(events, obj, pair_key="tkelePair", max_relPfIso=None, max_relPuppiIso=None):
    out = dict(obj)
    if pair_key not in out:
        return out
    p = out[pair_key]
    if p is None:
        return out

    mask = ak.ones_like(p.pt, dtype=bool)
    if max_relPfIso is not None and hasattr(p, "max_relPfIso"):
        mask = mask & (p.max_relPfIso < max_relPfIso)
    if max_relPuppiIso is not None and hasattr(p, "max_relPuppiIso"):
        mask = mask & (p.max_relPuppiIso < max_relPuppiIso)

    out[pair_key] = p[mask]
    out = cut_refresh_pair_match_views(out, pair_key=pair_key)
    out = cut_pick_best_pair(out, pair_key=pair_key, score=BESTPAIR_SCORE)
    return out


def pick_best_pair(pair_coll, score="min_dvz"):
    if pair_coll is None:
        return None

    # choose idx [evt,1]
    if score == "pt":
        idx = ak.argmax(pair_coll.pt, axis=1, keepdims=True)
    elif score == "min_dvz":
        idx = ak.argmin(pair_coll.delta_vz, axis=1, keepdims=True)
    elif score == "min_iso_pf":
        idx = ak.argmin(pair_coll.max_relPfIso, axis=1, keepdims=True)
    else:
        raise ValueError(score)

    best = pair_coll[idx]      # [evt,1] (can be empty per evt)
    best = ak.firsts(best)     # [evt] option-record
    return best


def add_pair_genmatch_flag(pairs):
    """
    Adds pair.isGenMatched.
    If no gen info on legs, sets False.
    """
    if pairs is None or (not hasattr(pairs, "l1")):
        return pairs

    l1, l2 = pairs.l1, pairs.l2

    if not (hasattr(l1, "isGenMatched") and hasattr(l2, "isGenMatched")):
        return ak.with_field(pairs, ak.zeros_like(pairs.pt, dtype=bool), "isGenMatched")

    both = (l1.isGenMatched & l2.isGenMatched)

    if hasattr(l1, "matchedGenIdx") and hasattr(l2, "matchedGenIdx"):
        both = both & (l1.matchedGenIdx != l2.matchedGenIdx)

    return ak.with_field(pairs, both, "isGenMatched")


def cut_refresh_pair_match_views(obj, pair_key="tkelePair"):
    """
    Keeps matched/nonMatched pair collections in-sync with current (cut) tkelePair.
    """
    out = dict(obj)
    if pair_key not in out:
        return out
    p = out[pair_key]
    if p is None or (not hasattr(p, "isGenMatched")):
        return out

    out["matched_tkelePair"] = p[p.isGenMatched]
    out["nonMatched_tkelePair"] = p[~p.isGenMatched]
    return out


def cut_pick_best_pair(obj, pair_key="tkelePair", out_key="best_tkelePair", score="min_dvz"):
    """
    Produces best pair as [evt][0/1] so your counting/histos work.
    """
    out = dict(obj)
    if pair_key not in out:
        return out
    p = out[pair_key]
    if p is None:
        return out

    n = ak.num(p, axis=1)
    has_any = (n > 0)

    if score == "pt":
        idx = ak.argmax(p.pt, axis=1, keepdims=True)
    elif score == "min_dvz":
        idx = ak.argmin(p.delta_vz, axis=1, keepdims=True)
    elif score == "min_iso_pf":
        idx = ak.argmin(p.max_relPfIso, axis=1, keepdims=True)
    else:
        raise ValueError(score)

    best = p[idx]                      # [evt,1] (empty per evt possible)
    best_opt = ak.firsts(best)         # [evt] option-record
    best_j01 = _as_jagged01(best_opt)  # [evt][0/1]

    out[out_key] = best_j01
    return out


def cut_event_on_bestpair(events, obj, best_key="best_tkelePair", *, ptmin=None, massmin=None, massmax=None):
    out = dict(obj)
    if best_key not in out:
        return out

    # Always convert to option-per-event
    bp = ak.firsts(out[best_key])   # safe for [evt,1], [evt,N], or option

    # start with "has a best pair"
    mask = ~ak.is_none(bp)

    if ptmin is not None and hasattr(bp, "pt"):
        pt = ak.fill_none(bp.pt, -1.0)
        mask = mask & (pt >= ptmin)

    if (massmin is not None or massmax is not None) and hasattr(bp, "mass"):
        m = ak.fill_none(bp.mass, -1.0)
        if massmin is not None:
            mask = mask & (m >= massmin)
        if massmax is not None:
            mask = mask & (m < massmax)

    return apply_evt_mask(out, mask)



def add_best_lead_sub(obj, best_key="best_tkelePair"):
    out = dict(obj)
    if best_key not in out:
        return out
    bp = out[best_key]  # [evt][0/1]
    if bp is None or (not hasattr(bp, "l1")):
        return out

    l1, l2 = bp.l1, bp.l2
    lead_mask = (l1.pt >= l2.pt)
    out["best_tkeleLead"] = ak.where(lead_mask, l1, l2)  # [evt][0/1]
    out["best_tkeleSub"]  = ak.where(lead_mask, l2, l1)
    return out


def cut_veto_if_no_bestpair(events, obj, best_key="best_tkelePair"):
    """
    Veto events with no best pair.
    Expects best_key to be jagged [evt][0/1].
    """
    if best_key not in obj:
        return obj

    bp = obj[best_key]
    if bp is None:
        return obj

    evt_mask = ak.num(bp, axis=1) > 0
    return apply_evt_mask(obj, evt_mask)


def cut_genmatch_after_buildpairs(events, obj):
    """
    Run right after cut_build_pairs.
    Decorate tkele with gen-match flags, then rebuild tkelePair with isGenMatched.
    Also fills matched/nonMatched tkele and matched_tkeleTrue for later resolution plots.
    """
    out = dict(obj)

    if not has_gen(events):
        # bkg: ensure pair has isGenMatched=False for consistent plotting
        if "tkelePair" in out:
            out["tkelePair"] = add_pair_genmatch_flag(out["tkelePair"])
            out = cut_refresh_pair_match_views(out)
        return out

    if "genel" not in out or "tkele" not in out:
        return out

    # reco->gen indices
    matched_mask, matched_gen_idx = utils.match_reco_to_gen_indices(
        out["genel"], out["tkele"], dr_cut=DR_CUT
    )

    # best-pt pruning
    best_mask = utils.keep_highest_pt_reco_per_gen(
        out["tkele"], matched_mask, matched_gen_idx, pt_field="pt"
    )
    matched_mask = matched_mask & best_mask

    # decorate reco electrons
    tkele = out["tkele"]
    tkele = ak.with_field(tkele, matched_mask, "isGenMatched")
    tkele = ak.with_field(tkele, matched_gen_idx, "matchedGenIdx")
    out["tkele"] = tkele

    out["matched_tkele"] = add_rel_iso(tkele[matched_mask])
    out["nonMatched_tkele"] = add_rel_iso(tkele[~matched_mask])

    # matched truth (for resolution)
    safe_idx = ak.where(matched_gen_idx < 0, 0, matched_gen_idx)
    matched_true_all = out["genel"][safe_idx]
    out["matched_tkeleTrue"] = matched_true_all[matched_mask]

    # IMPORTANT: rebuild pairs from updated tkele (legs have isGenMatched/matchedGenIdx)
    pair_all = build_pairs_from_electrons(out["tkele"])
    pair_all = add_pair_genmatch_flag(pair_all)
    out["tkelePair"] = pair_all

    # keep matched/nonMatched pair views in sync
    out = cut_refresh_pair_match_views(out)

    if "matched_tkelePair" in out:
        mlead, msub = build_lead_sub_from_pairs(out["matched_tkelePair"])
        if mlead is not None:
            out["matched_tkeleLead"] = add_rel_iso(mlead)
        if msub is not None:
            out["matched_tkeleSub"] = add_rel_iso(msub)

    if "nonMatched_tkelePair" in out:
        nlead, nsub = build_lead_sub_from_pairs(out["nonMatched_tkelePair"])
        if nlead is not None:
            out["nonMatched_tkeleLead"] = add_rel_iso(nlead)
        if nsub is not None:
            out["nonMatched_tkeleSub"] = add_rel_iso(nsub)

    return out


def cut_bestpair_leg_pt(events, obj,
                        lead_key="best_tkeleLead",
                        sub_key="best_tkeleSub",
                        lead_pt_min=5.0,
                        sub_pt_min=4.0):
    out = dict(obj)
    if (lead_key not in out) or (sub_key not in out):
        return out

    lead = out[lead_key]
    sub  = out[sub_key]
    if lead is None or sub is None:
        return out

    lead1 = ak.firsts(lead)  # [evt] option-record
    sub1  = ak.firsts(sub)

    mask = (~ak.is_none(lead1)) & (~ak.is_none(sub1))

    if hasattr(lead1, "pt"):
        mask = mask & (ak.fill_none(lead1.pt, -1.0) > lead_pt_min)
    if hasattr(sub1, "pt"):
        mask = mask & (ak.fill_none(sub1.pt, -1.0) > sub_pt_min)

    return apply_evt_mask(out, mask)


def cut_bestpair_leg_idscore(events, obj,
                             lead_key="best_tkeleLead",
                             sub_key="best_tkeleSub",
                             id_min=-0.1):
    out = dict(obj)
    if (lead_key not in out) or (sub_key not in out):
        return out

    lead = out[lead_key]
    sub  = out[sub_key]
    if lead is None or sub is None:
        return out

    lead1 = ak.firsts(lead)
    sub1  = ak.firsts(sub)

    mask = (~ak.is_none(lead1)) & (~ak.is_none(sub1))

    if hasattr(lead1, "idScore"):
        mask = mask & (ak.fill_none(lead1.idScore, -999.0) > id_min)
    if hasattr(sub1, "idScore"):
        mask = mask & (ak.fill_none(sub1.idScore, -999.0) > id_min)

    return apply_evt_mask(out, mask)



# ============================================================
# CUTFLOW
# ============================================================


CUTFLOW = [
    ("cut0_None",            [cut_base]),
    ("cut1_gen_twoPromptOS", [cut_evt_gen_two_prompt_os]),

    (f"cut2_eta_{ETA_REGION}", [cut_eta]),
    (f"cut3_pt_{int(PT_MIN)}", [cut_pt]),

    ("cut4_buildPairs",      [cut_build_pairs]),
    ("cut4b_genMatchAfterPairs", [cut_genmatch_after_buildpairs]),
    ("cut4c_pickBestPairInit", [lambda e,o: cut_pick_best_pair(o, score=BESTPAIR_SCORE)]),

    # pair cuts keep matched/nonMatched pair views AND bestPair updated
    ("cut5_pairOS",          [lambda e,o: cut_pair_os(e, o, "tkelePair")]),
    ("cut5b_vetoNoBestPair", [cut_veto_if_no_bestpair]),
    
    (f"cut6_pairDVZ_{PAIR_DVZ_MAX}", [lambda e,o: cut_pair_dvz(e, o, "tkelePair", PAIR_DVZ_MAX)]),
    ("cut6b_vetoNoBestPair", [cut_veto_if_no_bestpair]),
    
    ("cut7_pairIsoPF",         [lambda e,o: cut_pair_iso(e, o, "tkelePair", max_relPfIso=PAIR_MAX_RELPFISO)]),
    ("cut7b_vetoNoBestPair", [cut_veto_if_no_bestpair]),

    ("cut8_pairIsoPUPPI",         [lambda e,o: cut_pair_iso(e, o, "tkelePair", max_relPuppiIso=PAIR_MAX_RELPUISO)]),
    ("cut8b_vetoNoBestPair", [cut_veto_if_no_bestpair]),
    ("cut8c_bestLeadSub",        [lambda e,o: add_best_lead_sub(o, best_key="best_tkelePair")]),

    ("cut9_bestPairLegPt",      [lambda e,o: cut_bestpair_leg_pt(e, o, lead_pt_min=5.0, sub_pt_min=4.0)]),
    ("cut10_bestPairLegIdScore", [lambda e,o: cut_bestpair_leg_idscore(e, o, id_min=-0.1)]),

]

