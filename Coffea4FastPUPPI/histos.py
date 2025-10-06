import awkward as ak
import numpy as np

def fill_histo(selected_obj, histos, out):
    multiplicity = {}
    pdgId = {}
    isGenMatched = {}
    isReconstructed = {}

    resolution = {}

    # Calculate derived quantities
    for typ in selected_obj.keys():
        if 'jet' in typ:
            continue
        multiplicity[typ] = ak.num(selected_obj[typ].pt)
        pdgId[typ] = np.abs(ak.flatten(selected_obj[typ].pdgId))

    isGenMatched['pf'] = ak.num(selected_obj['matched_pf'].pt) / ak.num(selected_obj['pf'].pt)
    isGenMatched['puppi'] = ak.num(selected_obj['matched_puppi'].pt) / ak.num(selected_obj['puppi'].pt)
    isReconstructed['gen'] = ak.num(selected_obj['matched_gen'].pt) / ak.num(selected_obj['gen'].pt)

    resolution['matched_pf'] = selected_obj['matched_pf'].pt / selected_obj['matched_pfTrue'].pt 
    resolution['matched_puppi'] = selected_obj['matched_puppi'].pt / selected_obj['matched_puppiTrue'].pt 

    for hname, hcfg in histos.items():
        variables = hcfg['variables']
        histo = histos[hname]
        for typ in selected_obj.keys():
            cat = hname.split("_")[0] if 'atched' not in hname else f'{hname.split("_")[0]}_{hname.split("_")[1]}'
            if typ != cat:
                continue
            if 'jet' in typ:
                continue

            # Build a dict of {var_name: value_array}
            data = {}

            for var in variables:
                # Post-processing plots
                if var == "multiplicity":
                    arr = ak.num(selected_obj[typ].pt)
                elif var == "pdgId":
                    arr = np.abs(ak.flatten(selected_obj[typ].pdgId))
                elif var == "isMatched":
                    arr = isGenMatched.get(typ, np.array([]))
                elif var == "isReconstructed" and typ == "gen":
                    arr = isReconstructed['gen']
                elif var == "resolution":
                    if typ == "matched_pf":
                        arr = ak.flatten(resolution["matched_pf"])
                    elif typ == "matched_puppi":
                        arr = ak.flatten(resolution["matched_puppi"])
                # Pre-calcualted plots
                else:
                    arr = ak.flatten(getattr(selected_obj[typ], var))

                data[var] = arr

            # Make sure all arrays are same length (basic check)
            if len(set(len(v) for v in data.values())) > 1:
                print(f"[Warning] Skipping histogram '{hname}' due to mismatched data lengths")
                continue
            
            # skip 2D for now, memory issue
            if len(data.keys()) > 1:
                continue

            out[hname].fill(**data)

    return out
