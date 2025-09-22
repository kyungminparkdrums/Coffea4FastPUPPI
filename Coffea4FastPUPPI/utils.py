import numpy as np
import awkward as ak

def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    return (dphi + np.pi) % (2 * np.pi) - np.pi

def deltaR(eta1, phi1, eta2, phi2):
    dphi = delta_phi(phi1, phi2)
    deta = eta1 - eta2
    return np.sqrt(deta**2 + dphi**2)

def get_genMatched(gen, reco, typ='Gen'):
    if typ == 'Gen':
        pair_genCands, pair_recoCands= ak.unzip(ak.cartesian([gen, reco], nested=True))
    elif typ == 'Reco':
        pair_recoCands, pair_genCands= ak.unzip(ak.cartesian([reco, gen], nested=True))

    deltaR = np.sqrt((pair_genCands.eta - pair_recoCands.eta)**2 + (pair_genCands.phi - pair_recoCands.phi)**2)
    matched = deltaR < 0.1

    mask_matched = ak.any(matched, axis=-1)

    if typ == 'Gen':
        matched =  gen[mask_matched]
    elif typ == 'Reco':
        matched = reco[mask_matched]

    return matched

