import os
import shlex
import subprocess

# ----------------------------
# Samples
# ----------------------------
SIG_MASSES = [2, 5, 10, 15, 20, 30]  # GeV
SIG_FILE_TPL = "coffea/tkele_{m}.coffea"

sig_cut  = "cut4b_genMatchAfterPairs"
sig_hist_default = None  # use default histogram name from file

signals = [
    {
        "mass": m,
        "infile": SIG_FILE_TPL.format(m=m),
        "label": f"Gen-matched ee: {m} GeV",
    }
    for m in SIG_MASSES
]


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
        xtitle="M(ee) [GeV]",
        xmin=0,
        xmax=60,
        out="plots/sig_bkg/tkelePair_mass.png",
        **GLOBAL,
    ),
    # another plot
    dict(
         name="tkelePair_pt",
         sig_hist="matched_tkelePair_pt",
         xtitle="pT(ee)",
         xmin=0,
         xmax=60,
         out="plots/sig/tkelePair_pt.png",
         logy=True,
         density=True,
     ),
]

PLOTS.append(dict(
    name="genelPair_pt",
    sig_hist="genelPair_pt",
    xtitle=rf"$p_T$(Gen ee)",
    xmin=0,
    xmax=60,
    out="plots/sig/genelPair_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_eta",
    sig_hist="genelPair_eta",
    xtitle=rf"$\eta$(Gen ee)",
    xmin=-5,
    xmax=5,
    out="plots/sig/genelPair_eta.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_vz",
    sig_hist="genelPair_vz",
    xtitle=rf"$v_z$(Gen ee)",
    xmin=-15,
    xmax=15,
    out="plots/sig/genelPair_vz.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_mass",
    sig_hist="genelPair_mass",
    xtitle=rf"$M$(Gen ee)",
    xmin=0,
    xmax=50,
    out="plots/sig/genelPair_mass.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_deltaR",
    sig_hist="genelPair_deltaR",
    xtitle=rf"$\Delta R$(Gen ee)",
    xmin=0,
    xmax=5,
    out="plots/sig/genelPair_deltaR.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_deltaR",
    sig_hist="genelPair_deltaR",
    xtitle=rf"$\Delta R$(Gen ee)",
    xmin=0,
    xmax=5,
    out="plots/sig/genelPair_deltaR.png",
    logy=False,
    density=True,
))




PLOTS.append(dict(
    name="genel_multiplicity",
    sig_hist="genel_multiplicity",
    xtitle=rf"N(Gen electrons)",
    xmin=0,
    xmax=10,
    out="plots/sig/genel_multiplicity.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="genelLead_pt",
    sig_hist="genelLead_pt",
    xtitle=rf"$p_T$(Gen leading electron)",
    xmin=0,
    xmax=60,
    out="plots/sig/genelLead_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelSub_pt",
    sig_hist="genelSub_pt",
    xtitle=rf"$p_T$(Gen subleading electron)",
    xmin=0,
    xmax=60,
    out="plots/sig/genelSub_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelLead_eta",
    sig_hist="genelLead_eta",
    xtitle=rf"$\eta$(Gen leading electron)",
    xmin=-5,
    xmax=5,
    out="plots/sig/genelLead_eta.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelSub_eta",
    sig_hist="genelSub_eta",
    xtitle=rf"$\eta$(Gen subleading electron)",
    xmin=-5,
    xmax=5,
    out="plots/sig/genelSub_eta.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelLead_phi",
    sig_hist="genelLead_phi",
    xtitle=rf"$\phi$(Gen leading electron)",
    xmin=-5,
    xmax=5,
    out="plots/sig/genelLead_phi.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelSub_phi",
    sig_hist="genelSub_phi",
    xtitle=rf"$\phi$(Gen subleading electron)",
    xmin=-5,
    xmax=5,
    out="plots/sig/genelSub_phi.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelLead_vz",
    sig_hist="genelLead_vz",
    xtitle=rf"$v_z$(Gen leading electron)",
    xmin=-100,
    xmax=100,
    out="plots/sig/genelLead_vz.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelSub_vz",
    sig_hist="genelSub_vz",
    xtitle=rf"$v_z$(Gen subleading electron)",
    xmin=-100,
    xmax=100,
    out="plots/sig/genelSub_vz.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelLead_charge",
    sig_hist="genelLead_charge",
    xtitle=rf"$q$(Gen leading electron)",
    xmin=-1,
    xmax=1,
    out="plots/sig/genelLead_charge.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelSub_charge",
    sig_hist="genelSub_charge",
    xtitle=rf"$q$(Gen subleading electron)",
    xmin=-1,
    xmax=1,
    out="plots/sig/genelSub_charge.png",
    logy=True,
    density=True,
))



# ----------------------------
# Helpers
# ----------------------------
def build_cmd(job):
    cmd = ["python", "plot_overlay.py"]

    # signals (sorted by mass)
    for s in sorted(signals, key=lambda d: d["mass"]):
        cmd += [s["infile"], sig_cut, job["sig_hist"], s["label"]]

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
