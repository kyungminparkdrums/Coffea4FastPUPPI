import warnings
warnings.filterwarnings("ignore") # filter out some np warning
import coffea.processor as processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import hist
from hist import Hist
from hist import axis
import awkward as ak
import numpy as np
from utils import utils
from histos import fill_histo

class P2L1TAnalyzer(processor.ProcessorABC):
    def __init__(self, hist_config, cuts=None, pt=5.0, genPdgId=130, recoPdgId=None, etaRegion='endcap'):
        self.histos = hist_config
        self.cuts = cuts
        self.pt = pt
        self.genPdgId = genPdgId
        self.recoPdgId = recoPdgId
        self.etaRegion = etaRegion

        # define eta range by region
        self.eta_ranges = {
            'barrel': (0.0, 1.4),
            'endcap': (1.5, 2.4),
            'endcapNoTk': (2.5, 3.0),
            'hf': (3.0, 5.0),
        }

        accumulator = {}

        for hname, cfg in self.histos.items():
            variables = cfg["variables"]
            axes_cfg = cfg["axes"]

            axes = []
            for var, ax_cfg in zip(variables, axes_cfg):
                axis_type = ax_cfg["type"]
                label = ax_cfg.get("label", var)

                if axis_type == "Regular":
                    axes.append(axis.Regular(
                        bins=ax_cfg["bins"],
                        start=ax_cfg["range"][0],
                        stop=ax_cfg["range"][1],
                        name=var,
                        label=label
                    ))
                elif axis_type == "Variable":
                    axes.append(axis.Variable(
                        edges=ax_cfg["edges"],
                        name=var,
                        label=label
                    ))
                else:
                    raise ValueError(f"Unsupported axis type '{axis_type}' in histogram '{hname}'")
            accumulator[hname] = Hist(*axes)

            self._accumulator = processor.dict_accumulator(accumulator)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        out = self._accumulator

        obj = {
            'pf': events.L1PFCands,
            'puppi': events.L1PuppiCands,
            'gen': events.GenCands,
            'pfjet': events.L1PFJets,
            'puppijet': events.L1PuppiJets
        }

        selected_obj = {}

        eta_min, eta_max = self.eta_ranges[self.etaRegion]

        if self.cuts is not None:
            cut_mask = {}
            for typ in obj.keys():
                cut_mask[typ] = (
                    self.cuts.cut_pt(obj[typ], self.pt) &
                    self.cuts.cut_eta(obj[typ], eta_min, eta_max)
                )

                if typ == 'gen':
                    cut_mask[typ] = cut_mask[typ] & self.cuts.cut_pdgId(obj[typ], self.genPdgId)

                elif typ in ['pf', 'puppi']:
                    if self.recoPdgId is not None:
                        # Apply PDG ID filter for reco objects
                        cut_mask[typ] = cut_mask[typ] & self.cuts.cut_pdgId(obj[typ], self.recoPdgId)

                selected_obj[typ] = obj[typ][cut_mask[typ]]
        else:
            selected_obj = obj

        # Matching
        selected_obj['matched_gen'], selected_obj['nonMatched_gen'] = utils.get_genMatched(selected_obj['gen'], selected_obj['pf'], typ='Gen')
        selected_obj['matched_pf'], selected_obj['nonMatched_pf'], selected_obj['matched_pfTrue'] = utils.get_genMatched(selected_obj['gen'], selected_obj['pf'], typ='Reco')
        selected_obj['matched_puppi'], selected_obj['nonMatched_puppi'], selected_obj['matched_puppiTrue'] = utils.get_genMatched(selected_obj['gen'], selected_obj['puppi'], typ='Reco')

        selected_obj['matched_genPuppi'], selected_obj['nonMatched_genPuppi'] = utils.get_genMatched(selected_obj['gen'], selected_obj['puppi'], typ='Gen')
        selected_obj['nonMatched_pfPuppi'], _, _ = utils.get_genMatched(selected_obj['puppi'], selected_obj['nonMatched_pf'], typ='Reco')

        # HGC ID cut
        puWP = {'Lpt': 0.5, 'Hpt': 0.5}
        for ptType in ['Lpt', 'Hpt']:
            selected_obj[f'pfPuIdPass{ptType}'], selected_obj[f'pfPiIdPass{ptType}'], selected_obj[f'pfEmIdPass{ptType}'] = utils.cut_hgcIdPu(selected_obj['pf'], puWP[ptType])
            selected_obj[f'matched_genPfPuIdPass{ptType}'], _= utils.get_genMatched(selected_obj['gen'], selected_obj[f'pfPuIdPass{ptType}'], typ='Gen')
            selected_obj[f'matched_pfPuIdPass{ptType}'], _, _ = utils.get_genMatched(selected_obj['gen'], selected_obj[f'pfPuIdPass{ptType}'], typ='Reco')

            selected_obj[f'matched_genPfPiIdPass{ptType}'], _= utils.get_genMatched(selected_obj['gen'], selected_obj[f'pfPiIdPass{ptType}'], typ='Gen')
            selected_obj[f'matched_pfPiIdPass{ptType}'], _, _ = utils.get_genMatched(selected_obj['gen'], selected_obj[f'pfPiIdPass{ptType}'], typ='Reco')

            selected_obj[f'gen{ptType}'] = selected_obj['gen']

        # Fill histograms
        output = fill_histo(selected_obj, self.histos, out)
       
        # Add jet stuff
        # Initialize histograms
        """
        jet_multiplicity_hist = Hist(
            "jet_multiplicity", 
            Regular(50, 0, 50, label="Jet Constituents")
        )
        
        charged_hadron_fraction_hist = Hist(
            "charged_hadron_fraction", 
            Regular(50, 0, 1, label="Fraction of Charged Hadrons")
        )
        
        neutral_hadron_fraction_hist = Hist(
            "neutral_hadron_fraction", 
            Regular(50, 0, 1, label="Fraction of Neutral Hadrons")
        )

        # Loop over jets to analyze their constituents
        jet_constituents = utils.get_jetConstituents(selected_obj['pfjet'], selected_obj['pf'])
        multiplicity = ak.num(jet_constituents)
        charged_mask = (jet_constituents.pdgId == 211)  # PDG ID for charged pions
        neutral_mask = (jet_constituents.pdgId == 130)  # PDG ID for neutral kaons

        charged_hadrons = jet_constituents[charged_mask]
        charged_fraction = ak.num(charged_hadrons) / multiplicity
        neutral_hadrons = jet_constituents[neutral_mask]
        neutral_fraction = ak.num(neutral_hadrons) / multiplicity
 
        print('multiplicity', multiplicity)
        print('charged #', ak.num(charged_hadrons))
        print('neutral #', ak.num(neutral_hadrons))
        print('charged frac', charged_fraction)
        print('neutral frac', neutral_fraction)
        """
        return output

    def postprocess(self, accumulator):
        return accumulator
