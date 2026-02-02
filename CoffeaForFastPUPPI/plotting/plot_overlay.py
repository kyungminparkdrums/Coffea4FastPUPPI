import os
import argparse
import warnings
import re

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import coffea.util as util
import mplhep as hep

hep.style.use("CMS")

warnings.filterwarnings(
    "ignore",
    message=r".*All sumw are zero!.*",
    category=RuntimeWarning,
)

# Fixed colors for signal masses
MASS_COLORS = {
    2:  "#377eb8",  # blue
    5:  "#4daf4a",  # green
    10: "#e41a1c",  # red
    15: "#984ea3",  # purple
    20: "#ff7f00",  # orange
    30: "#ef7aa9",  # pink
}
BKG_COLOR = "black"

# -------------------------
# Normalization
# -------------------------
LUMI = 400.0  # fb^-1

SIG_XSEC_PB = {
    2: 44.16,
    5: 11.64,
    10: 6.683,
    15: 26.97,
    20: 17.41,
    30: 8.462,
}
SIG_XSEC_FB = {m: xs * 1000.0 for m, xs in SIG_XSEC_PB.items()}  # 1 pb = 1000 fb

NORM_BKG = 31.5 * 60 * 24 * 60 * 60 * 1_000_000  # fixed convention


def _hist_is_empty(h):
    try:
        vals = np.asarray(h.values(flow=False), dtype=float)
        return (vals.size == 0) or (np.nansum(vals) == 0.0)
    except Exception:
        return True


def _get_edges_and_vals(h):
    ax0 = h.axes[0]
    edges = np.asarray(ax0.edges, dtype=float)
    vals = np.asarray(h.values(flow=False), dtype=float)
    return edges, vals, ax0


def _extract_mass(label: str):
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*GeV", label)
    return float(m.group(1)) if m else None


def _infer_scale(is_bkg: bool, label: str, infile: str):
    if is_bkg:
        return float(NORM_BKG)

    m = _extract_mass(label)
    if m is None:
        m2 = re.search(r"[_\-](\d+)(?:GeV)?\b", os.path.basename(infile))
        m = float(m2.group(1)) if m2 else None

    if m is None or int(m) not in SIG_XSEC_FB:
        raise ValueError(f"Cannot infer signal mass/xsec for label='{label}' infile='{infile}'")

    return SIG_XSEC_FB[int(m)] * LUMI


def _step_from_edges(ax, edges, vals, *, label, density=False, linestyle="-", linewidth=2.0, color=None, scale=1.0):
    vals = np.asarray(vals, dtype=float) * float(scale)
    edges = np.asarray(edges, dtype=float)

    if density:
        area = np.sum(vals * np.diff(edges))
        if area > 0:
            vals = vals / area

    # extend last bin so the last segment is drawn
    y = np.r_[vals, vals[-1]]
    ax.step(edges, y, where="post", label=label, linestyle=linestyle, linewidth=linewidth, color=color)


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Overlay 1D histograms from multiple .coffea files on one canvas.\n\n"
            "Usage:\n"
            "  python plot_overlay.py infile stage hist label [infile stage hist label ...] --out out.png\n\n"
            "Each sample is a 4-tuple: infile stage hist label"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    ap.add_argument("items", nargs="+", help="Repeated groups of: infile stage hist label")
    ap.add_argument("--out", default="overlay.png", help="Output PNG path")
    ap.add_argument("--logy", action="store_true", help="Log y axis")
    ap.add_argument("--density", action="store_true", help="Normalize each hist to unit area")
    ap.add_argument("--xmin", type=float, default=None, help="Lower x-axis limit")
    ap.add_argument("--xmax", type=float, default=None, help="Upper x-axis limit")
    ap.add_argument("--xtitle", default=None, help="Plot x-axis title override")
    ap.add_argument("--ytitle", default=None, help="Plot y-axis title override")
    ap.add_argument("--scale", action="store_true", help="Scale histograms by xsec*lumi (signal) and fixed norm (bkg)")

    args = ap.parse_args()

    if len(args.items) % 4 != 0:
        raise SystemExit(
            "ERROR: Expected items in groups of 4 (infile stage hist label).\n"
            "Example:\n"
            '  python plot_overlay.py a.coffea cut hist "label A" b.coffea cut hist "label B" --out out.png'
        )

    samples = []
    for i in range(0, len(args.items), 4):
        infile, stage, histname, label = args.items[i:i + 4]
        samples.append((infile, stage, histname, label))

    loaded = {}
    hlist = []

    for infile, stage, histname, label in samples:
        if infile not in loaded:
            loaded[infile] = util.load(infile)

        res = loaded[infile]
        stage_hists = (res.get("hists", {}) or {}).get(stage, {})
        if histname not in stage_hists:
            raise KeyError(
                f"Hist '{histname}' not found in file='{infile}', stage='{stage}'. "
                f"Available (first 30): {list(stage_hists.keys())[:30]}"
            )

        h = stage_hists[histname]
        if len(h.axes) != 1:
            raise ValueError(f"'{histname}' in '{infile}:{stage}' is not 1D.")
        if _hist_is_empty(h):
            raise ValueError(f"Histogram '{histname}' in '{infile}:{stage}' is empty.")

        edges, vals, ax0 = _get_edges_and_vals(h)
        hlist.append((edges, vals, ax0, label, infile))

    # Same binning check
    ref_edges = hlist[0][0]
    for edges, _, _, label, _ in hlist[1:]:
        if edges.shape != ref_edges.shape or not np.allclose(edges, ref_edges, rtol=0, atol=1e-12):
            raise ValueError("Binnings differ between overlaid histograms.")

    outdir = os.path.dirname(args.out)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    fig, ax = plt.subplots()

    for edges, vals, _, label, infile in hlist:
        is_bkg = ("bkg" in label.lower()) or ("background" in label.lower())
        m = _extract_mass(label)

        linestyle = "--" if is_bkg else "-"
        color = BKG_COLOR if is_bkg else MASS_COLORS.get(int(m) if m is not None else -1, None)
        linewidth = 2.2  # a bit thinner than before, but still visible

        scale = _infer_scale(is_bkg, label, infile) if args.scale else 1.0

        _step_from_edges(
            ax, edges, vals,
            label=label,
            density=args.density,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color,
            scale=scale,
        )

    hep.cms.text("Phase-2 Simulation", ax=ax)

    xlabel = getattr(hlist[0][2], "label", hlist[0][2].name)
    ax.set_xlabel(args.xtitle if args.xtitle else xlabel)

    if args.ytitle:
        ax.set_ylabel(args.ytitle)
    else:
        if args.density:
            ax.set_ylabel("A.U.")
        else:
            ax.set_ylabel("Events" if args.scale else "Entries")

    if args.xmin is not None or args.xmax is not None:
        ax.set_xlim(
            left=args.xmin if args.xmin is not None else ax.get_xlim()[0],
            right=args.xmax if args.xmax is not None else ax.get_xlim()[1],
        )

    if args.logy:
        ax.set_yscale("log")

    # legend under plot, but tighter
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.10),
        #bbox_to_anchor=(0.5, -0.14),
        ncol=2,
        frameon=False,
        handlelength=2.6,
        columnspacing=1.4,
    )

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.26)

    fig.savefig(args.out)
    plt.close(fig)

    print("\nSaved:", args.out)
    return 0


if __name__ == "__main__":
    main()
