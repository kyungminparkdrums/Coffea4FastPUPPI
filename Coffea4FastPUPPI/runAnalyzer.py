import os

histo_config = 'histo_config/histo_idCut.yaml'
#histo_config = 'histo_config/histo_hgcId.yaml'
#histo_config = 'histo_config/histo_baseline.yaml'

rootDir = '/eos/cms/store/group/cmst3/group/l1tr/kypark/BESTPUPPI/AR2025_PFstudies_noCutHGCid/TT_PU200/FP/140Xv0C1/'
#rootDir = '/eos/cms/store/group/cmst3/group/l1tr/kypark/BESTPUPPI/AR2025_PFstudies_NNvtx/TT_PU200/FP/140Xv0C1/'
#rootDir = '/eos/cms/store/group/cmst3/group/l1tr/kypark/BESTPUPPI/AR2025_PFstudies_HGCid/TT_PU200/FP/140Xv0C1/'
#rootDir = '/eos/cms/store/group/cmst3/group/l1tr/kypark/BESTPUPPI/Baseline_PFstudies/TT_PU200/FP/140Xv0C1/'

outPrefix = 'AR2025_hgcPuId_wp0p5'

#outPrefix = 'AR2025_noHGCid'
#outPrefix = 'AR2025_NNVtx'
#outPrefix = 'AR2025'
#outPrefix = 'AR2025_matchPUPPI'
#outPrefix = 'Baseline'
#outPrefix = 'Baseline_matchPUPPI'

pt = 5

for etaRegion in ['endcap', 'endcapNoTk', 'hf']:
#for etaRegion in ['barrel', 'endcap', 'endcapNoTk', 'hf']:
    for genPdgId in [130, 211]:
        # Run
        outName = f'{outPrefix}_pt{pt}_{etaRegion}_genPdgId{genPdgId}' 
        cmd = f'python analyzer.py --histoConfig {histo_config} --etaRegion {etaRegion} --pt {pt} --genPdgId {genPdgId} --rootDir {rootDir} --outName {outName}'
        
        #outName = f'{outPrefix}_pt{pt}_{etaRegion}_genPdgId{genPdgId}_recoPdgId{genPdgId}' 
        #cmd = f'python analyzer.py --histoConfig {histo_config} --etaRegion {etaRegion} --pt {pt} --genPdgId {genPdgId} --rootDir {rootDir} --outName {outName} --recoPdgId {genPdgId}'

        print(cmd)
        os.system(cmd)
