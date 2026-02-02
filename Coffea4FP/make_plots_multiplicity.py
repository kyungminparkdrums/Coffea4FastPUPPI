import os
import shlex
import subprocess

# ----------------------------
# Samples
# ----------------------------
SIG_MASSES = [2, 5, 10, 15, 20, 30]  # GeV
SIG_FILE_TPL = "coffea/tkele_{m}.coffea"

#sig_cut  = "cut5_pairOS"
#sig_cut  = "cut6_pairDVZ_1.0"
#sig_cut  = "cut7_pairIsoPF"
sig_cut  = "cut8_pairIsoPUPPI"
sig_hist_default = None  # use default histogram name from file

signals = [
    {
        "mass": m,
        "infile": SIG_FILE_TPL.format(m=m),
        "label": f"Signal: {m} GeV",
    }
    for m in SIG_MASSES
]

bkg = {
    "infile": "coffea/test_bkg_tkele.coffea",
    "cut": "cut4_pairOS",
    #"cut": "cut5_pairDVZ_1.0",
    #"cut": "cut7_pairIsoPUPPI",
    "label": "Bkg ee",
}

# ----------------------------
# Plot settings (global defaults)
# ----------------------------
GLOBAL = dict(
    logy=True,
    density=True,
)

# ----------------------------
# Plot “jobs”
# Add new plots by appending entries here.
# ----------------------------
PLOTS = [
    dict(
        name="tkelePair_mass",
        sig_hist="matched_tkelePair_mass",
        bkg_hist="tkelePair_mass",
        xtitle="M(ee) [GeV]",
        xmin=0,
        xmax=60,
        out="plots/cutIsoPUPPI/tkelePair_mass.png",
        **GLOBAL,
    ),
    # another plot
    # dict(
    #     name="tkelePair_pt",
    #     sig_hist="matched_tkelePair_pt",
    #     bkg_hist="tkelePair_pt",
    #     xtitle="pT(ee)",
    #     xmin=0,
    #     xmax=200,
    #     out="plots/overlay_tkelePair_pt.png",
    #     logy=True,
    #     density=True,
    # ),
]

PLOTS.append(dict(
    name="tkelePair_multiplicity",
    sig_hist="tkelePair_multiplicity",
    bkg_hist="tkelePair_multiplicity",
    xtitle=rf"N(reco ee)",
    xmin=0,
    xmax=15,
    out="plots/cutIsoPUPPI/tkelePair_multiplicity.png",
    logy=False,
    density=True,
))

'''
PLOTS.append(dict(
    name="tkelePair_isGenMatched",
    sig_hist="tkelePair_isGenMatched",
    bkg_hist="tkelePair_isGenMatched",
    xtitle=rf"Is Gen Matched (reco ee)",
    xmin=0,
    xmax=1,
    out="plots/cutIsoPUPPI/tkelePair_isGenMatched.png",
    logy=False,
    density=True,
))
'''


# ----------------------------
# Helpers
# ----------------------------
def build_cmd(job):
    cmd = ["python", "plot_overlay.py"]

    # signals (sorted by mass)
    for s in sorted(signals, key=lambda d: d["mass"]):
        cmd += [s["infile"], sig_cut, job["sig_hist"], s["label"]]

    # background
    cmd += [bkg["infile"], bkg["cut"], job["bkg_hist"], bkg["label"]]

    # output + title
    cmd += ["--out", job["out"], "--xtitle", job.get("xtitle", job["name"])]

    # flags
    if job.get("logy", False):
        cmd += ["--logy"]
    if job.get("density", False):
        cmd += ["--density"]
    if job.get("xmin", None) is not None:
        cmd += ["--xmin", str(job["xmin"])]
    if job.get("xmax", None) is not None:
        cmd += ["--xmax", str(job["xmax"])]

    return cmd


def run(cmd):
    # pretty print
    print(" ".join(shlex.quote(x) for x in cmd))
    subprocess.run(cmd, check=True)


# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    os.makedirs("plots", exist_ok=True)

    for job in PLOTS:
        run(build_cmd(job))
