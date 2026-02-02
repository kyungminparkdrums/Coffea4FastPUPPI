import awkward as ak
import numpy as np

from cut_config.cut_config_tkele import add_best_lead_sub
from cut_config import ETA_RANGES, cut_range
from cut_config import ETA_RANGES, cut_range, apply_evt_mask

# ============================================================
# USER CONFIG (only edit here)
# ============================================================

MODE = "tkele_bkg"
ETA_REGION = "barrel"
PT_MIN = 1.0

PAIR_DVZ_MAX = 1.0
PAIR_MAX_RELPFISO = 1     # None to disable
PAIR_MAX_RELPUISO = 0.5    # None to disable

BESTPAIR_SCORE = "min_dvz"  # "pt" / "min_dvz" / "min_iso_pf"
BESTPAIR_PTMIN = 10.0       # None to disable


# ============================================================
# Helpers
# ============================================================

def add_rel_iso(coll):
    if coll is None or (not hasattr(coll, "pt")):
        return coll

    pt = coll.pt
    safe_pt = ak.where(pt > 0, pt, 1.0)

    if hasattr(coll, "pfIso"):
        coll = ak.with_field(coll, coll.pfIso / safe_pt, "relPfIso")
    if hasattr(coll, "puppiIso"):
        coll = ak.with_field(coll, coll.puppiIso / safe_pt, "relPuppiIso")

    return coll

def _as_jagged01(per_event_option):
    # [evt] option-record -> [evt][0/1] (never None)
    return ak.fill_none(ak.singletons(per_event_option), [])


def pick_best_pair(pair_coll, score="min_dvz"):
    """
    Returns [evt] option-record best pair.
    """
    if pair_coll is None:
        return None

    if score == "pt":
        idx = ak.argmax(pair_coll.pt, axis=1, keepdims=True)
    elif score == "min_dvz":
        idx = ak.argmin(pair_coll.delta_vz, axis=1, keepdims=True)
    elif score == "min_iso_pf":
        idx = ak.argmin(pair_coll.max_relPfIso, axis=1, keepdims=True)
    else:
        raise ValueError(score)

    best = pair_coll[idx]      # [evt,1] (empty where no pairs)
    return ak.firsts(best)     # [evt] option-record


def cut_pick_best_pair(events, obj, pair_key="tkelePair", out_key="best_tkelePair", score="min_dvz"):
    out = dict(obj)
    if pair_key not in out:
        return out
    best_opt = pick_best_pair(out[pair_key], score=score)  # [evt] option-record
    out[out_key] = _as_jagged01(best_opt)                  # [evt][0/1]
    return out


# ============================================================
# OBJECT BUILDING
# ============================================================

def build_objects(events):
    obj = {}

    if not hasattr(events, "TkEleL2"):
        return obj

    tkele = events.TkEleL2

    # Normalize pt to ptCorr
    if hasattr(tkele, "ptCorr"):
        tkele = ak.with_field(tkele, tkele.ptCorr, "pt")

    # ADD HERE
    tkele = add_rel_iso(tkele)

    obj["tkele"] = tkele
    return obj

def _delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    return (dphi + np.pi) % (2 * np.pi) - np.pi


# ============================================================
# Pair building helpers
# ============================================================
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

    # massless approx
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

    # charge product (only if charge exists)
    if hasattr(l1, "charge") and hasattr(l2, "charge"):
        charge_prod = l1.charge * l2.charge
        pairs = ak.with_field(pairs, charge_prod, "charge_prod")

    # eta product
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


# ============================================================
# CUTS
# ============================================================

def cut_base(events, obj):
    return obj

def cut_eta(events, obj):
    eta_min, eta_max = ETA_RANGES[ETA_REGION]
    out = obj
    if "tkele" in out:
        out = cut_range(events, out, "tkele", "eta", vmin=eta_min, vmax=eta_max, doAbs=True)
        out["tkele"] = add_rel_iso(out["tkele"])   # ADD
    return out

def cut_pt(events, obj):
    out = obj
    if "tkele" in out:
        out = cut_range(events, out, "tkele", "pt", vmin=PT_MIN, vmax=None, doAbs=False)
        out["tkele"] = add_rel_iso(out["tkele"])   # ADD
    return out

def cut_build_pairs(events, obj):
    out = dict(obj)
    if "tkele" not in out:
        return out

    out["tkelePair"] = build_pairs_from_electrons(out["tkele"])
    out["tkeleLead"], out["tkeleSub"] = build_lead_sub_from_pairs(out["tkelePair"])

    # ADD (lead/sub are electron records too)
    out["tkeleLead"] = add_rel_iso(out["tkeleLead"])
    out["tkeleSub"]  = add_rel_iso(out["tkeleSub"])

    return out

def cut_pair_os(obj, pair_key="tkelePair"):
    out = dict(obj)
    if pair_key not in out:
        return out
    p = out[pair_key]
    if p is None or (not hasattr(p, "charge_prod")):
        return out
    out[pair_key] = p[p.charge_prod < 0]
    return out


def cut_pair_dvz(obj, pair_key="tkelePair", dvz_max=1.0):
    out = dict(obj)
    if pair_key not in out:
        return out
    p = out[pair_key]
    if p is None or (not hasattr(p, "delta_vz")):
        return out
    out[pair_key] = p[p.delta_vz < dvz_max]
    return out


def cut_pair_max_reliso(obj, pair_key="tkelePair", max_relPfIso=None, max_relPuppiIso=None):
    out = dict(obj)
    if pair_key not in out:
        return out
    p = out[pair_key]
    if p is None:
        return out

    mask = ak.ones_like(p.pt, dtype=bool)
    if (max_relPfIso is not None) and hasattr(p, "max_relPfIso"):
        mask = mask & (p.max_relPfIso < max_relPfIso)
    if (max_relPuppiIso is not None) and hasattr(p, "max_relPuppiIso"):
        mask = mask & (p.max_relPuppiIso < max_relPuppiIso)

    out[pair_key] = p[mask]
    return out


def cut_veto_if_no_bestpair(events, obj, best_key="best_tkelePair"):
    if best_key not in obj:
        return obj
    bp = obj[best_key]
    if bp is None:
        return obj
    evt_mask = ak.num(bp, axis=1) > 0
    return apply_evt_mask(obj, evt_mask)


def cut_event_on_bestpair(events, obj, best_key="best_tkelePair", *, ptmin=None):
    out = dict(obj)
    if best_key not in out:
        return out

    bp = ak.firsts(out[best_key])  # [evt] option-record
    mask = ~ak.is_none(bp)

    if ptmin is not None and hasattr(bp, "pt"):
        mask = mask & (ak.fill_none(bp.pt, -1.0) >= ptmin)

    return apply_evt_mask(out, mask)


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
    ("cut0_None",              [cut_base]),
    (f"cut1_eta_{ETA_REGION}", [cut_eta]),
    (f"cut2_pt_{int(PT_MIN)}", [cut_pt]),

    ("cut3_buildPairs",        [cut_build_pairs,
                                lambda e,o: cut_pick_best_pair(e, o, score=BESTPAIR_SCORE)]),

    ("cut4_pairOS",            [lambda e,o: cut_pair_os(o, "tkelePair"),
                                lambda e,o: cut_pick_best_pair(e, o, score=BESTPAIR_SCORE)]),

    (f"cut5_pairDVZ_{PAIR_DVZ_MAX}", [lambda e,o: cut_pair_dvz(o, "tkelePair", PAIR_DVZ_MAX),
                                     lambda e,o: cut_pick_best_pair(e, o, score=BESTPAIR_SCORE)]),

    ("cut6_pairIsoPF",           [lambda e,o: cut_pair_max_reliso(o, "tkelePair",
                                                                 max_relPfIso=PAIR_MAX_RELPFISO),
                                lambda e,o: cut_pick_best_pair(e, o, score=BESTPAIR_SCORE)]),

    ("cut7_pairIsoPUPPI",           [lambda e,o: cut_pair_max_reliso(o, "tkelePair",
                                                                 max_relPuppiIso=PAIR_MAX_RELPUISO),
                                lambda e,o: cut_pick_best_pair(e, o, score=BESTPAIR_SCORE)]),


    ("cut8_vetoNoBestPair",    [cut_veto_if_no_bestpair]),
    ("cut8c_bestLeadSub",        [lambda e,o: add_best_lead_sub(o, best_key="best_tkelePair")]),

    ("cut9_bestPairLegPt",      [lambda e,o: cut_bestpair_leg_pt(e, o, lead_pt_min=5.0, sub_pt_min=4.0)]),
    ("cut10_bestPairLegIdScore", [lambda e,o: cut_bestpair_leg_idscore(e, o, id_min=-0.1)]),

]
