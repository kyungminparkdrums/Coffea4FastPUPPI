import awkward as ak
from utils import utils

from cut_config import ETA_RANGES, cut_range, cut_equal, apply_evt_mask


# ============================================================
# USER CONFIG (only edit here)
# ============================================================

# Choose one:
#  - "pfpuppi": use GenCands/L1PFCands/L1PuppiCands only
#  - "tkele"  : use GenEl/TkEleL2 only


#MODE = "pfpuppi"   # or "tkele"
MODE = "tkele"   # or "tkele"

#ETA_REGION = "endcap"   # "barrel" / "endcap" / "endcapNoTk" / "hf"
ETA_REGION = "barrel"   # "barrel" / "endcap" / "endcapNoTk" / "hf"

PT_MIN = 5.0
DR_CUT = 0.1

# For MODE="pfpuppi" only (GenCands has pdgId)
GEN_PDGID = 130      # set None to disable

# Optional reco pdgId (pf/puppi only)
RECO_PDGID = None    # set None to disable

# Optional event-level selection (GenEl-only example)
APPLY_GEN_TWO_PROMPT_OS = False

# ============================================================
# OBJECT BUILDING (mode-exclusive)
# ============================================================

def build_objects(events):
    """
    Build only the collections needed for the selected MODE.
    """
    if MODE == "pfpuppi":
        return {
            "pf": events.L1PFCands,
            "puppi": events.L1PuppiCands,
            "gen": events.GenCands,
        }

    if MODE == "tkele":
        obj = {
            "genel": events.GenEl,
            "tkele": events.TkEleL2,
        }
        # For TkEleL2, use ptCorr.
        # If TkEleL2 already has "pt", overwrite for consistency.
        obj["tkele"] = ak.with_field(obj["tkele"], obj["tkele"].ptCorr, "pt")

        return obj

    raise ValueError(f"Unknown MODE='{MODE}'. Use 'pfpuppi' or 'tkele'.")


# ============================================================
# CUTS
# ============================================================

def cut_base(events, obj):
    return obj


def cut_evt_gen_two_prompt_os(events, obj):
    """
    Event-level GenEl selection (only for MODE='tkele'):

      mask_genPrompt = GenEl_prompt == 2
      exactly 2 prompt electrons
      sum charge == 0
    """
    if MODE != "tkele":
        return obj

    mask_genPrompt = (events.GenEl_prompt == 2)
    evt_mask = (ak.num(events.GenEl_pt[mask_genPrompt], axis=1) == 2)
    evt_mask = evt_mask & (ak.sum(events.GenEl_charge, axis=1) == 0)
    return apply_evt_mask(obj, evt_mask)


def cut_pt(events, obj):
    out = obj
    for k in list(obj.keys()):
        out = cut_range(events, out, k, "pt", vmin=PT_MIN, vmax=None, doAbs=False)
    return out


def cut_eta(events, obj):
    eta_min, eta_max = ETA_RANGES[ETA_REGION]
    out = obj
    for k in list(obj.keys()):
        out = cut_range(events, out, k, "eta", vmin=eta_min, vmax=eta_max, doAbs=True)
    return out


def cut_gen_pdgid(events, obj):
    """
    Only applies in MODE='pfpuppi' because GenCands has pdgId.
    """
    if MODE != "pfpuppi":
        return obj
    if GEN_PDGID is None:
        return obj
    out = obj
    out = cut_equal(events, out, "gen", "pdgId", GEN_PDGID, doAbs=True)
    return out


def cut_reco_pdgid(events, obj):
    """
    Optional reco pdgId selection for pf/puppi only.
    """
    if MODE != "pfpuppi":
        return obj
    if RECO_PDGID is None:
        return obj
    out = obj
    out = cut_equal(events, out, "pf", "pdgId", RECO_PDGID, doAbs=True)
    out = cut_equal(events, out, "puppi", "pdgId", RECO_PDGID, doAbs=True)
    return out


def cut_add_matching(events, obj):
    """
    Mode-exclusive matching:
      - pfpuppi: gen<->pf/puppi matching
      - tkele  : genel<->tkele matching
    """
    out = dict(obj)

    if MODE == "pfpuppi":
        out["matched_gen"], out["nonMatched_gen"] = utils.get_genMatched(
            out["gen"], out["pf"], typ="Gen", dr_cut=DR_CUT
        )
        out["matched_pf"], out["nonMatched_pf"], out["matched_pfTrue"] = utils.get_genMatched(
            out["gen"], out["pf"], typ="Reco", dr_cut=DR_CUT
        )
        out["matched_puppi"], out["nonMatched_puppi"], out["matched_puppiTrue"] = utils.get_genMatched(
            out["gen"], out["puppi"], typ="Reco", dr_cut=DR_CUT
        )

    if MODE == "tkele":
        out["matched_genel"], out["nonMatched_genel"] = utils.get_genMatched(
            out["genel"], out["tkele"], typ="Gen", dr_cut=DR_CUT
        )
        out["matched_tkele"], out["nonMatched_tkele"], out["matched_tkeleTrue"] = utils.get_genMatched(
            out["genel"], out["tkele"], typ="Reco", dr_cut=DR_CUT
        )
    
    return out


# ============================================================
# CUTFLOW (mode-dependent)
# ============================================================

# No cut
CUTFLOW = [
    ("base", [cut_base]),
]

if MODE == "tkele" and APPLY_GEN_TWO_PROMPT_OS:
    CUTFLOW.append(("gen_twoPromptOS", [cut_evt_gen_two_prompt_os]))

CUTFLOW += [
    ("pt", [cut_pt]),
    ("eta", [cut_eta]),
]

if MODE == "pfpuppi":
    CUTFLOW += [
        ("genPdgId", [cut_gen_pdgid]),
    ]

CUTFLOW += [
    ("matching", [cut_add_matching]),
]
