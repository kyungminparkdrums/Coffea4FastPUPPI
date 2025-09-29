import warnings
warnings.filterwarnings("ignore") # filter out some np warning
import coffea.processor as processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import hist
from hist import Hist
import awkward as ak
import numpy as np
from utils import utils

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

        pdg = 130
        #pdg = 211

        if self.cuts is not None:
            cut_mask = {}
            for typ in obj.keys():
                cut_mask[typ] = self.cuts.cut_pt(obj[typ], 5) & self.cuts.cut_eta(obj[typ], 3.0, 5.0)
                #cut_mask[typ] = self.cuts.cut_pt(obj[typ], 5) & self.cuts.cut_eta(obj[typ], 2.5, 3.0)
                #cut_mask[typ] = self.cuts.cut_pt(obj[typ], 5) & self.cuts.cut_eta(obj[typ], 1.5, 2.4)
                #cut_mask[typ] = self.cuts.cut_pt(obj[typ], 5) & self.cuts.cut_eta(obj[typ], 0, 1.4)
                
                if typ == 'gen':
                    cut_mask[typ] = cut_mask[typ] & self.cuts.cut_pdgId(obj[typ], pdg)
                #elif typ == 'pf' or typ == 'puppi':
                #    cut_mask[typ] = cut_mask[typ] & self.cuts.cut_pdgId(obj[typ], pdg)
                #    cut_mask[typ] = cut_mask[typ] & self.cuts.cut_notPdgId(obj[typ], pdg)
                selected_obj[typ] = obj[typ][cut_mask[typ]]
        else:
            selected_obj[typ] = obj[typ]

        # gen-matching
        selected_obj['matched_gen'], selected_obj['nonMatched_gen'] = utils.get_genMatched(selected_obj['gen'], selected_obj['pf'], typ='Gen')
        selected_obj['matched_pf'], selected_obj['nonMatched_pf'] = utils.get_genMatched(selected_obj['gen'], selected_obj['pf'], typ='Reco')
        selected_obj['matched_puppi'], selected_obj['nonMatched_puppi'] = utils.get_genMatched(selected_obj['gen'], selected_obj['puppi'], typ='Reco')
        
        # calculate some quantities here
        multiplicity = {}
        pdgId = {}
        isGenMatched = {}
        isReconstructed = {}

        for typ in selected_obj.keys():
            multiplicity[typ] = ak.num(selected_obj[typ].pt)
            pdgId[typ] = np.abs(ak.flatten(selected_obj[typ].pdgId))
       
        isGenMatched['pf'] = ak.num(selected_obj['matched_pf'].pt) / ak.num(selected_obj['pf'].pt)
        isGenMatched['puppi'] = ak.num(selected_obj['matched_puppi'].pt) / ak.num(selected_obj['puppi'].pt)
        isReconstructed['gen'] = ak.num(selected_obj['matched_gen'].pt)/ak.num(selected_obj['gen'].pt)
        
        # Fill histogram
        for hname, hcfg in self.histos.items():
            var = hcfg['variable']

            for typ in selected_obj.keys():
                cat = hname.split("_")[0] if not 'atched' in hname else f'{hname.split("_")[0]}_{hname.split("_")[1]}'
                if typ != cat:
                    continue

                if var == "multiplicity":
                    arr = multiplicity[typ]
                elif var == "pdgId":
                    arr = pdgId[typ]
                elif typ == "pf" and var == "isMatched":
                    arr == isGenMatched['pf']
                elif typ == "puppi" and var == "isMatched":
                    arr == isGenMatched['puppi'] 
                elif typ == "gen" and var == "isReconstructed":
                    arr == isReconstructed['gen'] 
                else:
                    arr = ak.flatten(getattr(selected_obj[typ], var))

                out[hname].fill(**{var: arr})

        return out

    def postprocess(self, accumulator):
        return accumulator


