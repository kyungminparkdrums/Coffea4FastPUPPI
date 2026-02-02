import os
import argparse
import warnings

# No GUI popup on macOS
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import coffea.util as util
import mplhep as hep

hep.style.use("CMS")  # 1) CMS style


warnings.filterwarnings(
    "ignore",
    message=r".*All sumw are zero!.*",
    category=RuntimeWarning,
)


def print_cutflow(cutflow_dict, doMultiplicity=False):
    stages = list(cutflow_dict.keys())
    keys = set()
    for s in stages:
        keys |= set(cutflow_dict[s].keys())

    if doMultiplicity:
        preferred = [
            "gen", "matched_gen",
            "pf", "matched_pf",
            "puppi", "matched_puppi",

            "genel", "matched_genel",
            "tkele", "matched_tkele",

            "tkelePair", "matched_tkelePair",
        ]  
        print("\n=== OBJECT MULTIPLICITY ===")
    else:
        preferred = [ "events" ]
        print("\n=== EVENT CUTFLOW ===")

    cols = [k for k in preferred if k in keys]
    #cols = [k for k in preferred if k in keys] + sorted([k for k in keys if k not in preferred])

    header = f"{'stage':<18}" + "".join([f"{c:>12}" for c in cols])
    print(header)
    for s in stages:
        row = f"{s:<18}"
        for c in cols:
            row += f"{cutflow_dict[s].get(c, 0):>12}"
        print(row)
    print("")


def _ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def _hist_is_empty(h):
    try:
        vals = np.asarray(h.values(flow=False), dtype=float)
        if vals.size == 0:
            return True
        s = np.nansum(vals)
        return (not np.isfinite(s)) or (s == 0.0)
    except Exception:
        return True


def plot_1d_hist(h, out_png, *, title=None, xlabel=None, ylabel="Entries", logy=False):
    fig, ax = plt.subplots()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r".*All sumw are zero!.*", category=RuntimeWarning)
        h.plot1d(ax=ax)

    # CMS label
    hep.cms.text("Preliminary", ax=ax)

    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if logy:
        ax.set_yscale("log")

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_2d_hist(h, out_png, *, title=None):
    fig, ax = plt.subplots()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r".*All sumw are zero!.*", category=RuntimeWarning)
        h.plot2d(ax=ax)

    hep.cms.text("Preliminary", ax=ax)

    if title:
        ax.set_title(title)

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_ratio_hist(num_h, den_h, out_png, *, title, ylabel="ratio", ylim=(0.0, 1.05)):
    num = np.asarray(num_h.values(flow=False), dtype=float)
    den = np.asarray(den_h.values(flow=False), dtype=float)

    ax0 = den_h.axes[0]
    edges = np.asarray(ax0.edges, dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])

    ratio = np.full_like(num, np.nan, dtype=float)
    m = den > 0
    ratio[m] = num[m] / den[m]

    xlabel = getattr(ax0, "label", ax0.name)

    fig, ax = plt.subplots()
    ax.step(centers, ratio, where="mid")

    hep.cms.text("Preliminary", ax=ax)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(*ylim)

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def add_efficiency_purity_plots(stage_hists, out_stage_dir):
    def _make_ratio(num_key, den_key, out_name, title, ylabel, ylim=(0.0, 1.05)):
        if num_key not in stage_hists or den_key not in stage_hists:
            return 0
        num = stage_hists[num_key]
        den = stage_hists[den_key]
        if _hist_is_empty(num) or _hist_is_empty(den):
            return 0

        plot_ratio_hist(
            num, den,
            os.path.join(out_stage_dir, out_name),
            title=title,
            ylabel=ylabel,
            ylim=ylim,
        )
        return 1

    recipes = [
        ("matched_gen_pt",   "gen_pt",   "efficiency_gen_pt.png",
         "Efficiency vs pT (matched_gen_pt / gen_pt)", "efficiency", (0.0, 1.05)),

        ("matched_pf_pt",    "pf_pt",    "purity_pf_pt.png",
         "Purity vs pT (matched_pf_pt / pf_pt)", "purity", (0.0, 1.05)),

        ("matched_puppi_pt", "puppi_pt", "purity_puppi_pt.png",
         "Purity vs pT (matched_puppi_pt / puppi_pt)", "purity", (0.0, 1.05)),
    
        ("matched_genel_pt",   "genel_pt",   "efficiency_genel_pt.png",
         "Efficiency vs pT (matched_genel_pt / genel_pt)", "efficiency", (0.0, 1.05)),

        ("matched_tkele_pt",    "tkele_pt",    "purity_tkele_pt.png",
         "Purity vs pT (matched_tkele_pt / tkele_pt)", "purity", (0.0, 1.05)),
    ]

    made = 0
    for num_key, den_key, out_name, title, ylabel, ylim in recipes:
        made += _make_ratio(num_key, den_key, out_name, title, ylabel, ylim)
    return made


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("coffea_file", help="Path to .coffea file")
    ap.add_argument("--stage", default=None, help="Plot only this stage (default: all)")
    ap.add_argument("--outdir", default="plots", help="Output directory for PNGs")
    ap.add_argument("--tag", default=None, help="Override folder name under outdir")
    args = ap.parse_args()

    res = util.load(args.coffea_file)

    hists_by_stage = res.get("hists", {})
    cutflow = res.get("cutflow", {})

    # Optional sidecar plot config (no Hist.metadata usage)
    # expected shape: plot_cfg[stage][histname] = {"logy": bool}
    plot_cfg = res.get("plot_cfg", {}) or {}

    if cutflow:
        print_cutflow(cutflow, doMultiplicity=False) # event cutflow
        print_cutflow(cutflow, doMultiplicity=True)  # object multiplicity

    stages = list(hists_by_stage.keys())
    if args.stage is not None:
        if args.stage not in stages:
            raise KeyError(f"Stage '{args.stage}' not found. Available: {stages}")
        stages = [args.stage]

    stem = args.tag
    if stem is None:
        stem = os.path.splitext(os.path.basename(args.coffea_file))[0]

    for stage in stages:
        stage_hists = hists_by_stage[stage]
        out_stage_dir = os.path.join(args.outdir, stem, stage)
        _ensure_dir(out_stage_dir)

        wrote = 0
        skipped_empty = 0
        skipped_nd = 0

        stage_cfg = plot_cfg.get(stage, {}) if isinstance(plot_cfg, dict) else {}

        for hname, h in stage_hists.items():
            ndim = len(h.axes)
            if ndim == 0 or ndim > 2:
                skipped_nd += 1
                continue

            if _hist_is_empty(h):
                skipped_empty += 1
                continue

            # logy from sidecar (no metadata touching)
            logy = False
            cfg = stage_cfg.get(hname, {}) if isinstance(stage_cfg, dict) else {}
            if isinstance(cfg, dict):
                logy = bool(cfg.get("logy", False))

            out_png = os.path.join(out_stage_dir, f"{hname}.png")

            if ndim == 1:
                xlabel = getattr(h.axes[0], "label", h.axes[0].name)
                plot_1d_hist(h, out_png, title=hname, xlabel=xlabel, logy=logy)
            else:
                plot_2d_hist(h, out_png, title=hname)

            wrote += 1

        made = add_efficiency_purity_plots(stage_hists, out_stage_dir)

        print(
            f"[plot.py] Stage '{stage}': wrote {wrote}, skipped empty {skipped_empty}, "
            f"skipped ND {skipped_nd}, derived ratio plots {made} -> {out_stage_dir}"
        )

    return 0


if __name__ == "__main__":
    main()
    