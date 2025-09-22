import glob
import os
import yaml
from coffea import processor
import coffea.util as util

from processor import P2L1TAnalyzer
from cut_config import cuts

"""
Customize
"""
# output setting
outdir = 'coffea'
outname = 'baseline'

# histo 
with open("histo_config/histo.yaml") as f:
    hist_config = yaml.safe_load(f)

# files
root_dir = "/eos/cms/store/group/cmst3/group/l1tr/kypark/BESTPUPPI/Baseline/TT_PU200/FP/140Xv0C1/"
files = {
    "Baseline": {
        "files": glob.glob(os.path.join(root_dir, "perfNano*.root"))
    }
}

# coffea setting
nWorkers = 4

print(f"Input file directory: {root_dir} ({len(files['Baseline']['files'])} files)")

"""
Analyzer, from processor
"""
# analyzer
processor_instance = P2L1TAnalyzer(hist_config, cuts)

runner = processor.Runner(
    executor=processor.FuturesExecutor(workers=nWorkers, merging=True),
    schema=processor.NanoAODSchema,
    # chunksize=50000,
)

result = runner(
    fileset=files,              # dict of sample: [filelist]
    treename="Events",          # name of TTree in ROOT file
    processor_instance=processor_instance,
)

# save output as binary
util.save(result, f"{outdir}/{outname}.coffea")
print(f"Saved output: {outdir}/{outname}.coffea")
