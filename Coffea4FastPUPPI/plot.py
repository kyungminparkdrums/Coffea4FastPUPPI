from coffea import util
import os
import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")

# import coffea
pdgId = 130
etaRegion = 'barrel'

#prefix = f'AR2025_pt5_{etaRegion}_genPdgId{pdgId}'
#prefix = f'AR2025_pt5_{etaRegion}_genPdgId{pdgId}_recoPdgId{pdgId}'
prefix = f'AR2025_pt5_{etaRegion}_genPdgId{pdgId}_recoPdgIdNot{pdgId}'

histo = util.load(f"coffea/{prefix}.coffea")

outdir = f'/eos/user/k/kypark/www/L1PF/{prefix}/'

print(f'* Plotting from: coffea/{prefix}.coffea')
print(f'* Output: {outdir}')

# Plot
list_plot = list(histo.keys())
cat_plot = [ pname.split("_")[0] if not 'atched' in pname else f'{pname.split("_")[0]}_{pname.split("_")[1]}' for pname in list_plot ]
cat_plot = list(set(cat_plot))

print("* Type of plots: ", cat_plot)
print("======== Plotting ... ========")
for cat in cat_plot:
    os.system(f'mkdir -p {outdir}/{cat}')
    os.system(f'cp utils/index.php {outdir}')
    os.system(f'cp utils/index.php {outdir}/{cat}')
    os.system(f'cp -r utils/plot-viewer {outdir}')
    os.system(f'cp -r utils/plot-viewer {outdir}/{cat}')

# Purity
for typ in ['pf', 'puppi']:
    for kinematic in ['pt', 'eta']:
        fig, ax = plt.subplots(figsize=(10,10))
        hep.histplot(histo[f'matched_{typ}_{kinematic}']/histo[f'{typ}_{kinematic}'])

        pname = f'{typ}_purity_{kinematic}'
        print(pname)
        ax.set_title(pname)

        plt.savefig(f'{outdir}/{typ}/{pname}.png')
        plt.close()

# Efficiency
for kinematic in ['pt', 'eta']:
    fig, ax = plt.subplots(figsize=(10,10))
    hep.histplot(histo[f'matched_gen_{kinematic}']/histo[f'gen_{kinematic}'])

    pname = f'gen_efficiency_{kinematic}'
    print(pname)
    ax.set_title(pname)
 
    plt.savefig(f'{outdir}/gen/{pname}.png')
    plt.close()

# Other plots
for idx, plot in enumerate(list_plot):
    print(plot)

    fig, ax = plt.subplots(figsize=(10,10))
    hep.histplot(histo[plot])
    
    ax.set_title(plot)
    
    cat = plot.split("_")[0] if not 'atched' in plot else f'{plot.split("_")[0]}_{plot.split("_")[1]}'
    plt.savefig(f'{outdir}/{cat}/{plot}.png')
    plt.close()

