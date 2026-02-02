import os
import shlex
import subprocess

# ----------------------------
# Samples
# ----------------------------
SIG_MASSES = [2, 5, 10, 15, 20, 30]  # GeV
SIG_FILE_TPL = "coffea/tkele_{m}.coffea"

#sig_cut  = "cut10_bestPairLegIdScore"
sig_cut  = "cut8c_bestLeadSub"
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
    "cut": "cut10_bestPairLegIdScore",
    #"cut": "cut10_bestPairLegIdScore",
    "label": "Bkg",
}

# ----------------------------
# Plot settings (global defaults)
# ----------------------------
GLOBAL = dict(
    logy=True,
    density=False,
)

# ----------------------------
# Plot “jobs”
# Add new plots by appending entries here.
# ----------------------------
PLOTS = [
    dict(
        name="best_tkelePair_pt",
        sig_hist="best_tkelePair_pt",
        bkg_hist="best_tkelePair_pt",
        xtitle="pT(ee) [GeV]",
        xmin=0,
        xmax=60,
        out="plots/bestPair/tkelePair_pt.png",
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
    name="best_tkelePair_mass",
    sig_hist="best_tkelePair_mass",
    bkg_hist="best_tkelePair_mass",
    xtitle=rf"$M(ee)$ [GeV]",
    xmin=0,
    xmax=60,
    out="plots/bestPair/tkelePair_mass.png",
    logy=True,
    density=False,
))

PLOTS.append(dict(
    name="best_tkelePair_delta_pt",
    sig_hist="best_tkelePair_delta_pt",
    bkg_hist="best_tkelePair_delta_pt",
    xtitle=rf"$\Delta p_T(ee)$",
    xmin=0,
    xmax=60,
    out="plots/bestPair/tkelePair_delta_pt.png",
    logy=True,
    density=False,
))

PLOTS.append(dict(
    name="best_tkelePair_delta_phi",
    sig_hist="best_tkelePair_delta_phi",
    bkg_hist="best_tkelePair_delta_phi",
    xtitle=rf"$\Delta \phi(ee)$",
    xmin=0,
    xmax=3,
    out="plots/bestPair/tkelePair_delta_phi.png",
    logy=True,
    density=False,
))

PLOTS.append(dict(
    name="best_tkelePair_delta_r",
    sig_hist="best_tkelePair_delta_r",
    bkg_hist="best_tkelePair_delta_r",
    xtitle=rf"$\Delta r(ee)$",
    xmin=0,
    xmax=5,
    out="plots/bestPair/tkelePair_delta_r.png",
    logy=True,
    density=False,
))

PLOTS.append(dict(
    name="best_tkelePair_eta_prod",
    sig_hist="best_tkelePair_eta_prod",
    bkg_hist="best_tkelePair_eta_prod",
    xtitle=rf"$\eta_1 \times \eta_2$",
    xmin=-5,
    xmax=5,
    out="plots/bestPair/tkelePair_eta_prod.png",
    logy=True,
    density=False,
))

PLOTS.append(dict(
    name="best_tkeleLead_pt",
    sig_hist="best_tkeleLead_pt",
    bkg_hist="best_tkeleLead_pt",
    xtitle=rf"$p_T(lead)$ [GeV]",
    xmin=0,
    xmax=60,
    out="plots/bestPair/tkeleLead_pt.png",
    logy=True,
    density=False,
))

PLOTS.append(dict(
    name="best_tkeleSub_pt",
    sig_hist="best_tkeleSub_pt",
    bkg_hist="best_tkeleSub_pt",
    xtitle=rf"$p_T(sub)$ [GeV]",
    xmin=0,
    xmax=60,
    out="plots/bestPair/tkeleSub_pt.png",
    logy=True,
    density=False,
))

PLOTS.append(dict(
    name="best_tkeleLead_idScore",
    sig_hist="best_tkeleLead_idScore",
    bkg_hist="best_tkeleLead_idScore",
    xtitle=rf"ID Score (lead)$",
    xmin=-1,
    xmax=1,
    out="plots/bestPair/tkeleLead_idScore.png",
    logy=True,
    density=False,
))

PLOTS.append(dict(
    name="best_tkeleSub_idScore",
    sig_hist="best_tkeleSub_idScore",
    bkg_hist="best_tkeleSub_idScore",
    xtitle=rf"ID Score (sub)$",
    xmin=-1,
    xmax=1,
    out="plots/bestPair/tkeleSub_idScore.png",
    logy=True,
    density=False,
))





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

    cmd += ["--scale"]

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
