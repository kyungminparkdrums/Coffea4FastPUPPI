import warnings
warnings.filterwarnings("ignore")

import coffea.processor as processor
import awkward as ak
import numpy as np
from hist import Hist, axis

from histos import fill_histo


class HistAccumulator(processor.AccumulatorABC):
    """
    Wrap hist.Hist so coffea can do identity() and merging.
    Safe for pickling (cloudpickle) on workers.
    """
    def __init__(self, h: Hist):
        # use object.__setattr__ to avoid any custom attribute logic
        object.__setattr__(self, "_h", h)

    # ---- required by coffea AccumulatorABC ----
    def identity(self):
        hnew = self._h.copy()
        view = hnew.view(flow=True)
        view[...] = 0
        try:
            var = hnew.variances(flow=True)
            if var is not None:
                var[...] = 0
        except Exception:
            pass
        return HistAccumulator(hnew)

    def add(self, other):
        if other is None:
            return
        if isinstance(other, HistAccumulator):
            self._h += other._h
        elif isinstance(other, Hist):
            self._h += other
        else:
            raise TypeError(f"Cannot add {type(other)} to HistAccumulator")

    def __iadd__(self, other):
        self.add(other)
        return self

    def fill(self, **kwargs):
        self._h.fill(**kwargs)

    def __getstate__(self):
        return {"_h": self._h}

    def __setstate__(self, state):
        object.__setattr__(self, "_h", state["_h"])

    def __getattr__(self, name):
        h = object.__getattribute__(self, "_h")  # will raise AttributeError if missing
        return getattr(h, name)

    def __repr__(self):
        return f"HistAccumulator({repr(self._h)})"


class P2L1TAnalyzer(processor.ProcessorABC):
    def __init__(self, hist_config, cut_config_module):
        self.hist_config = hist_config
        self.cut_config = cut_config_module

        stages = [name for name, _ in self.cut_config.CUTFLOW]

        # -------- build HistAccumulators per stage ----------
        hists_by_stage = {}
        for stage in stages:
            stage_hists = {}
            for hname, cfg in self.hist_config.items():
                variables = cfg["variables"]
                axes_cfg = cfg["axes"]

                axes_list = []
                for var, ax_cfg in zip(variables, axes_cfg):
                    ax_type = ax_cfg["type"]
                    label = ax_cfg.get("label", var)
                    name = ax_cfg.get("name", var)

                    if ax_type == "Regular":
                        start, stop = ax_cfg["range"][0], ax_cfg["range"][1]
                        axes_list.append(
                            axis.Regular(
                                ax_cfg["bins"],
                                start,
                                stop,
                                name=name,
                                label=label,
                            )
                        )
                    elif ax_type == "Variable":
                        axes_list.append(
                            axis.Variable(
                                ax_cfg["edges"],
                                name=name,
                                label=label,
                            )
                        )
                    else:
                        raise ValueError(f"Unsupported axis type '{ax_type}' in {hname}")

                stage_hists[hname] = HistAccumulator(Hist(*axes_list))

            hists_by_stage[stage] = processor.dict_accumulator(stage_hists)

        hists = processor.dict_accumulator(hists_by_stage)

        cutflow = processor.dict_accumulator({
            stage: processor.defaultdict_accumulator(int) for stage in stages
        })

        self._accumulator = processor.dict_accumulator({
            "hists": hists,
            "cutflow": cutflow,
        })

    @property
    def accumulator(self):
        return self._accumulator

    def _count_objects(self, obj):
        counts = {}
        for k, coll in obj.items():
            if coll is None:
                counts[k] = 0
                continue

            if hasattr(coll, "pt"):
                vals = coll.pt
                try:
                    counts[k] = int(ak.sum(ak.num(vals, axis=1)))
                except Exception:
                    counts[k] = int(ak.sum(~ak.is_none(vals)))
                continue

            try:
                counts[k] = int(ak.sum(ak.num(coll, axis=1)))
            except Exception:
                counts[k] = int(ak.sum(~ak.is_none(coll)))

        return counts

    def _count_events(self, obj):
        if not obj:
            return 0
        for _, coll in obj.items():
            if coll is None:
                continue
            try:
                return int(len(coll))
            except Exception:
                continue
        return 0

    def process(self, events):
        out = self.accumulator.identity()

        cur = self.cut_config.build_objects(events)

        for stage_name, fns in self.cut_config.CUTFLOW:
            for fn in fns:
                cur = fn(events, cur)

            n_pass = self._count_events(cur)
            out["cutflow"][stage_name]["events"] += int(n_pass)

            counts = self._count_objects(cur)
            for k, v in counts.items():
                out["cutflow"][stage_name][k] += int(v)

            stage_hists = out["hists"][stage_name]
            stage_hists = fill_histo(cur, self.hist_config, stage_hists)
            out["hists"][stage_name] = stage_hists

        return out

    def postprocess(self, accumulator):
        return accumulator

