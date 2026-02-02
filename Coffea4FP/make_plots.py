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

bkg = {
    "infile": "coffea/test_bkg_tkele.coffea",
    "cut": "cut3_buildPairs",
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
        out="plots/sig_bkg/tkelePair_mass.png",
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
    name="tkelePair_pt",
    sig_hist="matched_tkelePair_pt",
    bkg_hist="tkelePair_pt",
    xtitle=rf"$p_T(ee)$",
    xmin=0,
    xmax=60,
    out="plots/sig_bkg/tkelePair_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_eta",
    sig_hist="matched_tkelePair_eta",
    bkg_hist="tkelePair_eta",
    xtitle=rf"$\eta(ee)$",
    xmin=-5,
    xmax=5,
    out="plots/sig_bkg/tkelePair_eta.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_phi",
    sig_hist="matched_tkelePair_phi",
    bkg_hist="tkelePair_phi",
    xtitle=rf"$\phi(ee)$",
    xmin=-5,
    xmax=5,
    out="plots/sig_bkg/tkelePair_phi.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_vz",
    sig_hist="matched_tkelePair_vz",
    bkg_hist="tkelePair_vz",
    xtitle=rf"$v_z(ee)$",
    xmin=-100,
    xmax=100,
    out="plots/sig_bkg/tkelePair_vz.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_delta_r",
    sig_hist="matched_tkelePair_delta_r",
    bkg_hist="tkelePair_delta_r",
    xtitle=rf"$\Delta r(ee)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkelePair_delta_r.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_mass",
    sig_hist="matched_tkelePair_mass",
    bkg_hist="tkelePair_mass",
    xtitle=rf"$M(ee)$ [GeV]",
    xmin=0,
    xmax=60,
    out="plots/sig_bkg/tkelePair_mass.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_delta_pt",
    sig_hist="matched_tkelePair_delta_pt",
    bkg_hist="tkelePair_delta_pt",
    xtitle=rf"$\Delta p_T(ee)$",
    xmin=0,
    xmax=60,
    out="plots/sig_bkg/tkelePair_delta_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_delta_eta",
    sig_hist="matched_tkelePair_delta_eta",
    bkg_hist="tkelePair_delta_eta",
    xtitle=rf"$\Delta \eta(ee)$",
    xmin=0,
    xmax=4,
    out="plots/sig_bkg/tkelePair_delta_eta.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_delta_phi",
    sig_hist="matched_tkelePair_delta_phi",
    bkg_hist="tkelePair_delta_phi",
    xtitle=rf"$\Delta \phi(ee)$",
    xmin=0,
    xmax=3,
    out="plots/sig_bkg/tkelePair_delta_phi.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_delta_vz",
    sig_hist="matched_tkelePair_delta_vz",
    bkg_hist="tkelePair_delta_vz",
    xtitle=rf"$\Delta v_z(ee)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkelePair_delta_vz.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_delta_pt_over_pt",
    sig_hist="matched_tkelePair_delta_pt_over_pt",
    bkg_hist="tkelePair_delta_pt_over_pt",
    xtitle=rf"$\Delta p_T(ee)$",
    xmin=0,
    xmax=1,
    out="plots/sig_bkg/tkelePair_delta_pt_over_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_delta_eta_over_eta",
    sig_hist="matched_tkelePair_delta_eta_over_eta",
    bkg_hist="tkelePair_delta_eta_over_eta",
    xtitle=rf"$\Delta \eta(ee)$",
    xmin=0,
    xmax=4,
    out="plots/sig_bkg/tkelePair_delta_eta_over_eta.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_delta_phi_over_phi",
    sig_hist="matched_tkelePair_delta_phi_over_phi",
    bkg_hist="tkelePair_delta_phi_over_phi",
    xtitle=rf"$\Delta \phi(ee)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkelePair_delta_phi_over_phi.png",
    logy=True,
    density=True,
))


PLOTS.append(dict(
    name="tkelePair_pt_over_mass",
    sig_hist="matched_tkelePair_pt_over_mass",
    bkg_hist="tkelePair_pt_over_mass",
    xtitle=rf"$p_T(ee)/M(ee)$",
    xmin=0,
    xmax=10,
    out="plots/sig_bkg/tkelePair_pt_over_mass.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_multiplicity",
    sig_hist="tkelePair_multiplicity",
    bkg_hist="tkelePair_multiplicity",
    xtitle=rf"N(reco ee)",
    xmin=0,
    xmax=50,
    out="plots/sig_bkg/tkelePair_multiplicity.png",
    logy=False,
    density=True,
))

PLOTS.append(dict(
    name="tkele_multiplicity",
    sig_hist="tkele_multiplicity",
    bkg_hist="tkele_multiplicity",
    xtitle=rf"N(reco e)",
    xmin=0,
    xmax=15,
    out="plots/sig_bkg/tkele_multiplicity.png",
    logy=False,
    density=True,
))


PLOTS.append(dict(
    name="tkeleLead_pt",
    sig_hist="matched_tkeleLead_pt",
    bkg_hist="tkeleLead_pt",
    xtitle=rf"$p_T(Leading electron)$ [GeV]",
    xmin=0,
    xmax=60,
    out="plots/sig_bkg/tkeleLead_pt.png",
    logy=True,
    density=True,
))


PLOTS.append(dict(
    name="tkeleSub_pt",
    sig_hist="matched_tkeleSub_pt",
    bkg_hist="tkeleSub_pt",
    xtitle=rf"$p_T(Subleading electron)$ [GeV]",
    xmin=0,
    xmax=60,
    out="plots/sig_bkg/tkeleSub_pt.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleLead_eta",
    sig_hist="matched_tkeleLead_eta",
    bkg_hist="tkeleLead_eta",
    xtitle=rf"$\eta(Leading electron)$",
    xmin=-5,
    xmax=5,
    out="plots/sig_bkg/tkeleLead_eta.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleSub_eta",
    sig_hist="matched_tkeleSub_eta",
    bkg_hist="tkeleSub_eta",
    xtitle=rf"$\eta(Subleading electron)$",
    xmin=-5,
    xmax=5,
    out="plots/sig_bkg/tkeleSub_eta.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleLead_vz",
    sig_hist="matched_tkeleLead_vz",
    bkg_hist="tkeleLead_vz",
    xtitle=rf"$v_z(Leading electron)$",
    xmin=-5,
    xmax=5,
    out="plots/sig_bkg/tkeleLead_vz.png",
    logy=True,
    density=True,
))


PLOTS.append(dict(
    name="tkeleSub_vz",
    sig_hist="matched_tkeleSub_vz",
    bkg_hist="tkeleSub_vz",
    xtitle=rf"$v_z(Subleading electron)$",
    xmin=-5,
    xmax=5,
    out="plots/sig_bkg/tkeleSub_vz.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleLead_idScore",
    sig_hist="matched_tkeleLead_idScore",
    bkg_hist="tkeleLead_idScore",
    xtitle=rf"$idScore(Leading electron)$",
    xmin=-1,
    xmax=1,
    out="plots/sig_bkg/tkeleLead_idScore.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleSub_idScore",
    sig_hist="matched_tkeleSub_idScore",
    bkg_hist="tkeleSub_idScore",
    xtitle=rf"$idScore(Subleading electron)$",
    xmin=-1,
    xmax=1,
    out="plots/sig_bkg/tkeleSub_idScore.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleLead_pfIso",
    sig_hist="matched_tkeleLead_pfIso",
    bkg_hist="tkeleLead_pfIso",
    xtitle=rf"$pfIso(Leading electron)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkeleLead_pfIso.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleSub_pfIso",
    sig_hist="matched_tkeleSub_pfIso",
    bkg_hist="tkeleSub_pfIso",
    xtitle=rf"$pfIso(Subleading electron)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkeleSub_pfIso.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleLead_puppiIso",
    sig_hist="matched_tkeleLead_puppiIso",
    bkg_hist="tkeleLead_puppiIso",
    xtitle=rf"$puppiIso(Leading electron)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkeleLead_puppiIso.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleSub_puppiIso",
    sig_hist="matched_tkeleSub_puppiIso",
    bkg_hist="tkeleSub_puppiIso",
    xtitle=rf"$puppiIso(Subleading electron)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkeleSub_puppiIso.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleLead_relPfIso",
    sig_hist="matched_tkeleLead_relPfIso",
    bkg_hist="tkeleLead_relPfIso",
    xtitle=rf"$relPfIso(Leading electron)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkeleLead_relPfIso.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleSub_relPfIso",
    sig_hist="matched_tkeleSub_relPfIso",
    bkg_hist="tkeleSub_relPfIso",
    xtitle=rf"$relPfIso(Subleading electron)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkeleSub_relPfIso.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleLead_relPuppiIso",
    sig_hist="matched_tkeleLead_relPuppiIso",
    bkg_hist="tkeleLead_relPuppiIso",
    xtitle=rf"$relPuppiIso(Leading electron)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkeleLead_relPuppiIso.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkeleSub_relPuppiIso",
    sig_hist="matched_tkeleSub_relPuppiIso",
    bkg_hist="tkeleSub_relPuppiIso",
    xtitle=rf"$relPuppiIso(Subleading electron)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkeleSub_relPuppiIso.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_max_relPfIso",
    sig_hist="matched_tkelePair_max_relPfIso",
    bkg_hist="tkelePair_max_relPfIso",
    xtitle=rf"$max(relPfIso)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkelePair_max_relPfIso.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_max_relPuppiIso",
    sig_hist="matched_tkelePair_max_relPuppiIso",
    bkg_hist="tkelePair_max_relPuppiIso",
    xtitle=rf"$max(relPuppiIso)$",
    xmin=0,
    xmax=5,
    out="plots/sig_bkg/tkelePair_max_relPuppiIso.png",
    logy=True,
    density=True,
))


PLOTS.append(dict(
    name="tkelePair_charge_prod",
    sig_hist="matched_tkelePair_charge_prod",
    bkg_hist="tkelePair_charge_prod",
    xtitle=rf"$q1*q2$",
    xmin=-1,
    xmax=2,
    out="plots/sig_bkg/tkelePair_charge_prod.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_eta_prod",
    sig_hist="matched_tkelePair_eta_prod",
    bkg_hist="tkelePair_eta_prod",
    xtitle=rf"$\eta_1*\eta_2$",
    xmin=-5,
    xmax=5,
    out="plots/sig_bkg/tkelePair_eta_prod.png",
    logy=True,
    density=True,
))

PLOTS.append(dict(
    name="tkelePair_pt_over_mass",
    sig_hist="matched_tkelePair_pt_over_mass",
    bkg_hist="tkelePair_pt_over_mass",
    xtitle=rf"$pT(ee)/m(ee)$",
    xmin=0,
    xmax=20,
    out="plots/sig_bkg/tkelePair_pt_over_mass.png",
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
