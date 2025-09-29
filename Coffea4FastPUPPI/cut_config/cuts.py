import numpy as np

def cut_pt(obj, pt_thres):
    cut = obj.pt > pt_thres

    return cut

def cut_eta(obj, eta_min, eta_max):
    eta_cut = np.abs(obj.eta) >= eta_min
    eta_cut = eta_cut & (np.abs(obj.eta) < eta_max)

    cut = eta_cut

    return cut

def cut_pdgId(obj, pdgId):
    cut = np.abs(obj.pdgId) == pdgId

    return cut

def cut_notPdgId(obj, pdgId):
    cut = np.abs(obj.pdgId) != pdgId

    return cut
