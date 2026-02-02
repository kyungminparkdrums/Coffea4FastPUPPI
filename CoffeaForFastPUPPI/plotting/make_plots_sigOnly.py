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
        "label": f"Signal: {m} GeV",
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
        xmax=40,
        out="plots/sig_bkg/tkelePair_mass.png",
        **GLOBAL,
    ),
    # another plot
    # dict(
    #     name="tkelePair_pt",
    #     sig_hist="matched_tkelePair_pt",
    #     xtitle="pT(ee)",
    #     xmin=0,
    #     xmax=200,
    #     out="plots/overlay_tkelePair_pt.png",
    #     logy=True,
    #     density=True,
    # ),
]

PLOTS.append(dict(
    name="genel_multiplicity",
    sig_hist="genel_multiplicity",
    xtitle=rf"N(Gen electrons)",
    xmin=0,
    xmax=10,
    out="plots/sig_genlevel/genel_multiplicity.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="genelLead_pt",
    sig_hist="genelLead_pt",
    xtitle=rf"$p_T$(Gen leading electron)",
    xmin=0,
    xmax=60,
    out="plots/sig_genlevel/genelLead_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelSub_pt",
    sig_hist="genelSub_pt",
    xtitle=rf"$p_T$(Gen subleading electron)",
    xmin=0,
    xmax=60,
    out="plots/sig_genlevel/genelSub_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelLead_eta",
    sig_hist="genelLead_eta",
    xtitle=rf"$\eta$(Gen leading electron)",
    xmin=-5,
    xmax=5,
    out="plots/sig_genlevel/genelLead_eta.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelSub_eta",
    sig_hist="genelSub_eta",
    xtitle=rf"$\eta$(Gen subleading electron)",
    xmin=-5,
    xmax=5,
    out="plots/sig_genlevel/genelSub_eta.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelLead_phi",
    sig_hist="genelLead_phi",
    xtitle=rf"$\phi$(Gen leading electron)",
    xmin=-5,
    xmax=5,
    out="plots/sig_genlevel/genelLead_phi.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelSub_phi",
    sig_hist="genelSub_phi",
    xtitle=rf"$\phi$(Gen subleading electron)",
    xmin=-5,
    xmax=5,
    out="plots/sig_genlevel/genelSub_phi.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelLead_vz",
    sig_hist="genelLead_vz",
    xtitle=rf"$v_z$(Gen leading electron)",
    xmin=-20,
    xmax=20,
    out="plots/sig_genlevel/genelLead_vz.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelSub_vz",
    sig_hist="genelSub_vz",
    xtitle=rf"$v_z$(Gen subleading electron)",
    xmin=-20,
    xmax=20,
    out="plots/sig_genlevel/genelSub_vz.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelLead_charge",
    sig_hist="genelLead_charge",
    xtitle=rf"$q$(Gen leading electron)",
    xmin=-1,
    xmax=1,
    out="plots/sig_genlevel/genelLead_charge.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelSub_charge",
    sig_hist="genelSub_charge",
    xtitle=rf"$q$(Gen subleading electron)",
    xmin=-1,
    xmax=1,
    out="plots/sig_genlevel/genelSub_charge.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_pt",
    sig_hist="genelPair_pt",
    xtitle=rf"$p_T$(Gen electron pair)",
    xmin=0,
    xmax=60,
    out="plots/sig_genlevel/genelPair_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_eta",
    sig_hist="genelPair_eta",
    xtitle=rf"$\eta$(Gen electron pair)",
    xmin=-5,
    xmax=5,
    out="plots/sig_genlevel/genelPair_eta.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_phi",
    sig_hist="genelPair_phi",
    xtitle=rf"$\phi$(Gen electron pair)",
    xmin=-5,
    xmax=5,
    out="plots/sig_genlevel/genelPair_phi.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_vz",
    sig_hist="genelPair_vz",
    xtitle=rf"$v_z$(Gen electron pair)",
    xmin=-20,
    xmax=20,
    out="plots/sig_genlevel/genelPair_vz.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_mass",
    sig_hist="genelPair_mass",
    xtitle=rf"$m$(Gen electron pair)",
    xmin=0,
    xmax=60,
    out="plots/sig_genlevel/genelPair_mass.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_delta_r",
    sig_hist="genelPair_delta_r",
    xtitle=rf"$\Delta r$(Gen electron pair)",
    xmin=0,
    xmax=5,
    out="plots/sig_genlevel/genelPair_delta_r.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_delta_pt",
    sig_hist="genelPair_delta_pt",
    xtitle=rf"$\Delta p_T$(Gen electron pair)",
    xmin=0,
    xmax=60,
    out="plots/sig_genlevel/genelPair_delta_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_delta_eta",
    sig_hist="genelPair_delta_eta",
    xtitle=rf"$\Delta \eta$(Gen electron pair)",
    xmin=0,
    xmax=5,
    out="plots/sig_genlevel/genelPair_delta_eta.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_delta_phi",
    sig_hist="genelPair_delta_phi",
    xtitle=rf"$\Delta \phi$(Gen electron pair)",
    xmin=0,
    xmax=3,
    out="plots/sig_genlevel/genelPair_delta_phi.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_delta_vz",
    sig_hist="genelPair_delta_vz",
    xtitle=rf"$\Delta v_z$(Gen electron pair)",
    xmin=-20,
    xmax=20,
    out="plots/sig_genlevel/genelPair_delta_vz.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_delta_pt_over_pt",
    sig_hist="genelPair_delta_pt_over_pt",
    xtitle=rf"$\Delta p_T / p_T$(Gen electron pair)",
    xmin=0,
    xmax=1,
    out="plots/sig_genlevel/genelPair_delta_pt_over_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_delta_eta_over_eta",
    sig_hist="genelPair_delta_eta_over_eta",
    xtitle=rf"$\Delta \eta / \eta$(Gen electron pair)",
    xmin=0,
    xmax=5,
    out="plots/sig_genlevel/genelPair_delta_eta_over_eta.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="genelPair_delta_phi_over_phi",
    sig_hist="genelPair_delta_phi_over_phi",
    xtitle=rf"$\Delta \phi / \phi$(Gen electron pair)",
    xmin=0,
    xmax=5,
    out="plots/sig_genlevel/genelPair_delta_phi_over_phi.png",
    logy=True,
    density=True,
))


PLOTS.append(dict(
    name="matched_tkelePair_multiplicity",
    sig_hist="matched_tkelePair_multiplicity",
    xtitle=rf"has matched pair",
    xmin=0,
    xmax=2,
    out="plots/sig_genlevel/matched_tkelePair_multiplicity.png",
    logy=False,
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
