import numpy as np
import awkward as ak

def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    return (dphi + np.pi) % (2 * np.pi) - np.pi

def deltaR(eta1, phi1, eta2, phi2):
    dphi = delta_phi(phi1, phi2)
    deta = eta1 - eta2
    return np.sqrt(deta**2 + dphi**2)

def cut_hgcIdPu(pf, wp = 0.5):
    passPuId = pf.hgcIdPu < wp 
   
    passPiId = passPuId & (pf.hgcIdPi > pf.hgcIdEm)
    passEmId = passPuId & (pf.hgcIdPi < pf.hgcIdEm)

    return pf[passPuId], pf[passPiId], pf[passEmId]

def get_genMatched(gen, reco, typ='Gen', dr_cut = 0.1):
    if typ == 'Gen':
        gen_reco_pairs = ak.cartesian([gen, reco], axis=1, nested=True)
        gen_cands, reco_cands = ak.unzip(gen_reco_pairs)
    elif typ == 'Reco':
        reco_gen_pairs = ak.cartesian([reco, gen], axis=1, nested=True)
        reco_cands, gen_cands = ak.unzip(reco_gen_pairs)

    delta_eta = gen_cands.eta - reco_cands.eta
    delta_phi = np.abs(gen_cands.phi - reco_cands.phi)
    delta_phi = ak.where(delta_phi > np.pi, 2*np.pi - delta_phi, delta_phi)
    deltaR = np.sqrt(delta_eta**2 + delta_phi**2)  # shape: [events][gen][reco]

    min_deltaR = ak.min(deltaR, axis=-1)  # shape: [events][gen]
    min_deltaR = ak.fill_none(min_deltaR, 999.)
    matched_mask = min_deltaR < dr_cut  # shape: [events][gen]

    if typ == 'Gen':
        matched = gen[matched_mask]
        nonMatched = gen[~matched_mask]
    
        return matched, nonMatched

    elif typ == 'Reco':
        matched = reco[matched_mask]
        nonMatched = reco[~matched_mask]
        
        # Get indices of matched gen for each reco candidate
        min_deltaR_idx = ak.argmin(deltaR, axis=-1)  # shape: [events][reco]
        # For matched reco, get corresponding gen candidate
        matchedGen = gen[min_deltaR_idx][matched_mask]
        
        return matched, nonMatched, matchedGen

def get_jetConstituents(jet, ptcl, dr_cut=0.4):
    # Create all possible pairs between jets and particles (jets x particles)
    jet_ptcl_pairs = ak.cartesian([jet, ptcl], axis=1, nested=True)
    jet_cands, ptcl_cands = ak.unzip(jet_ptcl_pairs)

    # Calculate deltaR between each jet and particle pair
    delta_eta = jet_cands.eta - ptcl_cands.eta
    delta_phi = np.abs(jet_cands.phi - ptcl_cands.phi)
    delta_phi = ak.where(delta_phi > np.pi, 2*np.pi - delta_phi, delta_phi)
    deltaR = np.sqrt(delta_eta**2 + delta_phi**2)  # shape: [events][jets][particles]

    # Create a mask for particles that are within dR < dr_cut
    matched_mask = deltaR < dr_cut  # shape: [events][jets][particles]

    # For each jet, get the corresponding matched particles (ptcls)
    matched_ptcls = ptcl_cands[matched_mask]  # shape: [events][jets][matched particles]

    return matched_ptcls

