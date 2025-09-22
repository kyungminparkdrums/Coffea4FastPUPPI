import warnings
warnings.filterwarnings("ignore") # filter out some np warning
import coffea.processor as processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import hist
from hist import Hist
import awkward as ak
import numpy as np
import utils

class P2L1TAnalyzer(processor.ProcessorABC):
    def __init__(self, hist_config, cuts=None):
        self.histos = hist_config
        self.cuts = cuts

        accumulator = {}

        for hname, cfg in self.histos.items():
            axis = hist.axis.Regular(
                cfg["bins"],
                cfg["range"][0],
                cfg["range"][1],
                name=cfg["variable"],
                label=cfg["xlabel"]
            )
            accumulator[hname] = Hist(axis)

        self._accumulator = processor.dict_accumulator(accumulator)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        out = self._accumulator

        # cut
        obj = {}
        obj['pf'] = events.L1PFCands
        obj['puppi'] = events.L1PuppiCands
        obj['gen'] = events.GenCands
        
        selected_obj = {}

        pdg = 211

        if self.cuts is not None:
            cut_mask = {}
            for typ in obj.keys():
                cut_mask[typ] = self.cuts.cut_pt(obj[typ]) & self.cuts.cut_eta(obj[typ], 0, 1.4)
                
                if typ == 'gen':
                    cut_mask[typ] = cut_mask[typ] & self.cuts.cut_pdgId(obj[typ], pdg)
                #cut_mask[typ] = cut_mask[typ] & self.cuts.cut_pdgId(obj[typ], pdg)
                selected_obj[typ] = obj[typ][cut_mask[typ]]
        else:
            selected_obj[typ] = obj[typ]

        # gen-matching
        selected_obj['matched_gen'] = utils.get_genMatched(selected_obj['gen'], selected_obj['pf'], typ='Gen')
        selected_obj['matched_pf'] = utils.get_genMatched(selected_obj['gen'], selected_obj['pf'], typ='Reco')
        selected_obj['matched_puppi'] = utils.get_genMatched(selected_obj['gen'], selected_obj['puppi'], typ='Reco')
        
        # calculate some quantities here
        multiplicity = {}
        pdgId = {}

        for typ in selected_obj.keys():
            multiplicity[typ] = ak.num(selected_obj[typ].pt)
            pdgId[typ] = np.abs(ak.flatten(selected_obj[typ].pdgId))

        # Fill histogram
        for hname, hcfg in self.histos.items():
            var = hcfg['variable']

            for typ in selected_obj.keys():
                if typ not in hname:
                    continue

                if var == "multiplicity":
                    arr = multiplicity[typ]
                elif var == "pdgId":
                    arr = pdgId[typ]
                elif var == "isMatched":
                    continue
                elif var == "isReconstructed":
                    continue
                else:
                    arr = ak.flatten(getattr(selected_obj[typ], var))

                out[hname].fill(**{var: arr})
                #print(f'Filling histogram: {hname}')

        return out

    def postprocess(self, accumulator):
        return accumulator


