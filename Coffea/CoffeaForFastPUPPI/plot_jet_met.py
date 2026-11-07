import argparse
import os

import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt


import mplhep as hep

hep.style.use("CMS")

# Must match GENMET_BINS in the cut config exactly.
GENMET_EDGES = np.linspace(0., 200., 20)
GENMET_ORDER = [(f"bin{i:02d}", float(lo), float(hi))
                for i, (lo, hi) in enumerate(zip(GENMET_EDGES[:-1], GENMET_EDGES[1:]))]


METHODS_CENTRAL = [
    ("puppi_central_met_perf", "Baseline PUPPI", "o"),
    ("gentarget_central_met_perf", "Gen-target PUPPI", "s"),
    ("ml_central_met_perf", "Neutral ML PUPPI", "^"),
]

METHODS_BARREL = [
    ("puppi_barrel_met_perf", "Baseline PUPPI", "o"),
    ("gentarget_barrel_met_perf", "Gen-target PUPPI", "s"),
    ("ml_barrel_met_perf", "Neutral ML PUPPI", "^"),
]


def cms_label(ax):
    hep.cms.label(
        llabel="Private Work (Phase-2 Simulation)",
        rlabel="PU 200 (14 TeV)", data=False, com=14,
        loc=0, ax=ax, fontsize=20,
        #loc=0, ax=ax
    )


def unwrap_hist(h):
    return h._h if hasattr(h, "_h") else h


def get_stage_hists(output, stage):
    hists = output["hists"]

    if stage is None:
        stages = list(hists.keys())
        if len(stages) != 1:
            print("Available stages:")
            for s in stages:
                print(" ", s)
            raise RuntimeError("Multiple stages found. Please pass --stage.")
        stage = stages[0]

    return hists[stage], stage


def values_edges(h):
    h = unwrap_hist(h)
    values = np.asarray(h.values(flow=False), dtype=float)
    edges = np.asarray(h.axes[0].edges, dtype=float)
    if values.ndim != 1:
        raise ValueError("Expected a 1D histogram; select dataset/category axes first")
    if np.any(values < 0):
        raise ValueError("These distribution statistics require nonnegative bin counts")
    return values, edges


def hist_mean(h):
    values, edges = values_edges(h)
    centers = 0.5 * (edges[:-1] + edges[1:])
    n = np.sum(values)

    if n <= 0:
        return np.nan, 0.0

    mean = np.sum(values * centers) / n
    return mean, n


def hist_std(h):
    values, edges = values_edges(h)
    centers = 0.5 * (edges[:-1] + edges[1:])
    n = np.sum(values)

    if n <= 0:
        return np.nan, 0.0

    mean = np.sum(values * centers) / n
    std = np.sqrt(np.sum(values * (centers - mean) ** 2) / n)
    return std, n


class PlotInputs(dict):
    """Histogram distributions plus unbinned MET matrices from the processor."""
    def __init__(self, hists, met_events):
        super().__init__(hists)
        self.met_events = met_events


def event_met_curve(inputs, prefix):
    if prefix not in inputs.met_events:
        raise RuntimeError(f"Missing event-level MET data for {prefix}; rerun the updated processor and cut config")
    saved = inputs.met_events[prefix]
    a = np.asarray(saved.value if hasattr(saved, "value") else saved, dtype=float)
    if a.ndim != 2 or a.shape[1] != 6:
        raise ValueError(f"Invalid MET matrix for {prefix}: {a.shape}")
    gen, reco = a[:, 2], a[:, 5]
    valid = np.all(np.isfinite(a), axis=1)
    points = []
    for tag, lo, hi in GENMET_ORDER:
        mask = valid & (gen >= lo) & (gen < hi)
        n = int(mask.sum())
        if n == 0:
            continue
        median_gen = np.median(gen[mask])
        if median_gen <= 0:
            continue
        value = np.median(reco[mask]) / median_gen
        # Ratio of medians: do not reuse the ratio-of-means uncertainty.
        error = np.nan
        points.append(((lo + hi) / 2, value, error, n))
    return tuple(np.asarray(points).T) if points else (np.array([]),) * 4


def get_response_curve(stage_hists, prefix):
    """Per-bin median(reco MET) / median(GEN MET), from saved events."""
    if not isinstance(stage_hists, PlotInputs):
        raise RuntimeError("Median MET response requires saved event-level MET data")
    return event_met_curve(stage_hists, prefix)


def get_met_resolution_curve(inputs, prefix, component):
    """Central 68% half-width (q84-q16)/2 of unbinned MET residuals, in GeV."""
    if prefix not in inputs.met_events:
        raise RuntimeError(f"Missing event-level MET data for {prefix}; rerun the updated processor and cut config")
    saved = inputs.met_events[prefix]
    a = np.asarray(saved.value if hasattr(saved, "value") else saved, dtype=float)
    if a.ndim != 2 or a.shape[1] != 6:
        raise ValueError(f"Invalid MET matrix for {prefix}: {a.shape}")
    gx, gy, gen, rx, ry, reco = a.T
    valid = np.all(np.isfinite(a), axis=1) & (gen > 0)
    safe_gen = np.where(gen > 0, gen, 1.)
    parallel = rx * (gx / safe_gen) + ry * (gy / safe_gen)
    perpendicular = -rx * (gy / safe_gen) + ry * (gx / safe_gen)
    if component == "parallel":
        residual = parallel - gen
    elif component == "perpendicular":
        residual = perpendicular
    elif component == "perpendicular_minus_gen_debug":
        residual = perpendicular - gen
    else:
        raise ValueError(f"Unknown MET component: {component}")
    valid &= np.isfinite(residual)
    points = []
    #for tag, lo, hi in GENMET_ORDER:

    # Retain the uploaded plotter's original MET resolution binning.
    resolution_edges = GENMET_EDGES[::1]
    if resolution_edges[-1] != GENMET_EDGES[-1]:
        resolution_edges = np.append(resolution_edges, GENMET_EDGES[-1])

    for lo, hi in zip(resolution_edges[:-1], resolution_edges[1:]):

        mask = valid & (gen >= lo) & (gen < hi)
        n = int(mask.sum())
        if n < 2:
            continue
        #width = np.std(residual[mask], ddof=0)
        #error = width / np.sqrt(2 * (n - 1))
        q16, q84 = np.quantile(residual[mask], [0.16, 0.84])
        width = (q84 - q16) / 2
        error = np.nan

        points.append(((lo + hi) / 2, width, error, n))
    return tuple(np.asarray(points).T) if points else (np.array([]),) * 4


# Explicit histogram suffixes: never search for a generic "resolution" key.
GENJET_ORDER = [(f"{lo}_{hi}", lo, hi) for lo, hi in
                zip([0,20,40,60,80,100,150], [20,40,60,80,100,150,200])]
METHODS_JET = [(f"{m}_jet_perf", label, marker) for m, label, marker in
               [("puppi", "Baseline PUPPI", "o"),
                ("gentarget", "Gen-target PUPPI", "s"),
                ("ml", "Neutral ML PUPPI", "^")]]


def hist_median(h):
    values, edges = values_edges(h)
    n = values.sum()
    if n <= 0:
        return np.nan
    index = min(np.searchsorted(np.cumsum(values), n / 2), len(values)-1)
    previous = values[:index].sum()
    return edges[index] + (edges[index+1]-edges[index]) * (n/2-previous)/values[index]


def component_curve(hists, prefix, order, suffix, statistic):
    points = []
    for tag, lo, hi in order:
        h = hists.get(f"{prefix}_{suffix}_{tag}")
        if h is None:
            continue
        std, n = hist_std(h)
        if n < 2:
            continue
        if statistic == "std":
            value, error = std, std / np.sqrt(2*(n-1))
        elif statistic == "mean":
            value, _ = hist_mean(h)
            error = std / np.sqrt(n)
        elif statistic == "median":
            value, error = hist_median(h), np.nan
        elif statistic == "std_over_median":
            median = hist_median(h)
            if median <= 0:
                continue
            value, error = std / median, np.nan
        else:
            raise ValueError(statistic)
        points.append(((lo+hi)/2, value, error, n))
    return tuple(np.asarray(points).T) if points else (np.array([]),)*4


def plot_curve(hists, methods, getter, title, ylabel, xlabel, path,
               xmax, ylim=None, ideal=None):
    #fig, ax = plt.subplots()
    fig, ax = plt.subplots(figsize=(10,8))
    #fig, ax = plt.subplots(figsize=(7,6))
    found = False
    for prefix, label, marker in methods:
        x, y, error, counts = getter(hists, prefix)
        if not len(x):
            print(f"Missing/empty: {title}: {prefix}")
            continue
        found = True
        # Medians and std/median have no reliable uncertainty from these
        # marginal histograms alone. Do not reuse a mean's standard error.
        yerr = error if np.all(np.isfinite(error)) else None
        color = ('magenta' if prefix.startswith('ml_') else
                 'green' if prefix.startswith('gentarget_') else 'blue')
        ax.errorbar(x, y, yerr=yerr, fmt=marker, linestyle='-',
                    color=color, capsize=2, label=label)
    if not found:
        plt.close(fig)
        return
    if ideal is not None:
        ax.axhline(ideal, color='black', linestyle='--', linewidth=1,
                   label='' if ideal == 1 else '')
    ax.set(xlabel=xlabel, ylabel=ylabel, xlim=(0, xmax))
    #ax.xaxis.label.set_size()
    #ax.yaxis.label.set_size()
    if ylim is not None:
        ax.set_ylim(*ylim)
  
    is_met = any("_met_perf" in prefix for prefix, _, _ in methods)

    text = (
        r"Central MET ($|\eta| < 2.4$)"
        if is_met else
        r"Seeded Cone Jets ($\Delta R < 0.4$)"
    )
    x_pos = 0.28 if not is_met else 0.04
    ax.text(
        x_pos, 0.1, text,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=20,
    )

    ax.legend()
    #ax.legend(fontsize=14)
    cms_label(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print('Saved', path)


def main():
    from coffea.util import load
    parser = argparse.ArgumentParser(description="MET and seeded-cone jet validation")
    parser.add_argument('--input', required=True)
    parser.add_argument('--stage', default=None)
    parser.add_argument('--outdir', default='./plots/')
    parser.add_argument('--xmax', type=float, default=200.)
    parser.add_argument('--response_ylim', nargs=2, type=float, default=[0.8,1.2])
    parser.add_argument('--resolution_ylim', nargs=2, type=float, default=[10., 50.],
                        help='MET component resolution y limits in GeV')
    parser.add_argument('--jet_xmax', type=float, default=185.)
    parser.add_argument('--jet_response_ylim', nargs=2, type=float, default=[0.62,1.01])
    parser.add_argument('--jet_resolution_ylim', nargs=2, type=float, default=[-0.03,0.62])
    parser.add_argument('--skip-distributions' , action='store_true',
                        help='Accepted for compatibility; only summary curves are produced')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    output = load(args.input)
    hists, stage = get_stage_hists(output, args.stage)
    met_events = output.get("met_events", {}).get(stage, {})
    if not met_events:
        raise RuntimeError("No event-level MET data: rerun with the updated processor.py and cut config")
    hists = PlotInputs(hists, met_events)
    print('Using stage:', stage)
    print('MET curves use saved event values (no histogram approximation). '
          'Jet statistics still use histogram bin centers and interpolated medians.')


    # Exactly five plots per comparison. Use central MET only.
    comparisons = [
        ('baseline_vs_gentarget', 'gentarget'),
        ('baseline_vs_ml', 'ml'),
    ]
    for comparison, alternative in comparisons:
        comparison_dir = os.path.join(args.outdir, comparison)
        os.makedirs(comparison_dir, exist_ok=True)
        met_methods = [method for method in METHODS_CENTRAL
                       if method[0].startswith(('puppi_', alternative + '_'))]
        jet_methods = [method for method in METHODS_JET
                       if method[0].startswith(('puppi_', alternative + '_'))]
        plot_curve(
            hists, met_methods, get_response_curve, 'Central MET',
            r'median($\mathrm{MET}^{\mathrm{reco}}$/$\mathrm{MET}^{\mathrm{gen}}$)', r'gen MET [GeV]',
            os.path.join(comparison_dir, 'met_response_ratio_of_medians_central.png'),
            args.xmax, args.response_ylim, ideal=1.,
        )
        for component, ylabel in (
            ('parallel', r'$\sigma_{\mathrm{quantile}}(\mathrm{MET}_{\parallel}^{\mathrm{reco}}-\mathrm{MET}^{\mathrm{gen}}$) [GeV]'),
            ('perpendicular', r'$\sigma_{\mathrm{quantile}}(\mathrm{MET}_{\perp}^{\mathrm{reco}}$) [GeV]'),
        ):
            getter = lambda h, p, c=component: get_met_resolution_curve(h, p, c)
            plot_curve(
                hists, met_methods, getter, f'Central MET {component}',
                ylabel, r'gen MET [GeV]',
                os.path.join(comparison_dir, f'met_resolution_{component}_central.png'),
                args.xmax, ylim=([15., 40.] if component == 'perpendicular'
                                 else args.resolution_ylim),
            )
        jet_specs = [
            ('jet_response_median', 'response_ratio', 'median',
             r'median($\mathrm{p_{T}^{reco}}$ / $\mathrm{p_{T}^{gen}}$)', 1., args.jet_response_ylim),
            ('jet_reco_pt_std_over_median', 'recojetpt_dist', 'std_over_median',
             r'$\sigma$($\mathrm{p_{T}^{reco}}$) / median($\mathrm{p_{T}^{reco}}$)', 0., args.jet_resolution_ylim),
        ]
        for name, suffix, stat, ylabel, ideal, ylim in jet_specs:
            getter = lambda h, p, s=suffix, t=stat: component_curve(h, p, GENJET_ORDER, s, t)
            plot_curve(
                hists, jet_methods, getter, 'Seeded-cone jets, R=0.4', ylabel,
                r'gen jet $\mathrm{p_{T}}$ [GeV]', os.path.join(comparison_dir, name+'.png'),
                args.jet_xmax, ylim=ylim, ideal=ideal,
            )


if __name__ == '__main__':
    main()



