import glob
import os
import yaml
import argparse
from coffea import processor
import coffea.util as util

from processor import P2L1TAnalyzer
from cut_config import cuts

# Set up argument parser
parser = argparse.ArgumentParser(description="Run the Coffea analysis script with custom settings.")
parser.add_argument("--etaRegion", type=str, default="endcap", choices=["barrel", "endcap", "endcapNoTk", "hf"],
                    help="The eta region to analyze (default: 'endcap')")
parser.add_argument("--outName", type=str, required=True,
                    help="The output file name (e.g., 'AR2025_pt5_hf_genPdgId130')")
parser.add_argument("--rootDir", type=str, required=True,
                    help="Directory where ROOT files are stored")
parser.add_argument("--pt", type=float, default=5.0,
                    help="Pt threshold (default: 5.0)")
parser.add_argument("--genPdgId", type=int, default=130,
                    help="Generated particle PDG ID (default: 130)")
parser.add_argument("--recoPdgId", type=int, default=None,
                    help="Reconstructed particle PDG ID (default: None, meaning no cut)")
parser.add_argument("--histoConfig", type=str, default="histo_config/histo.yaml",  # Default path
                   help="Path to the histogram configuration YAML file.")

# Parse arguments
args = parser.parse_args()

# Use arguments to set values
etaRegion = args.etaRegion
outname = args.outName
rootDir = args.rootDir
pt = args.pt
genPdgId = args.genPdgId
recoPdgId = args.recoPdgId
histo_config = args.histoConfig

# Histograms
with open(histo_config) as f:
    hist_config = yaml.safe_load(f)

# Files
files = {
    "files": {
        "files": glob.glob(os.path.join(rootDir, "perfNano*.root"))
    }
}

# Coffea settings
nWorkers = 4

print(f"Input file directory: {rootDir} ({len(files['files']['files'])} files)")
print(f"Analysis settings: eta = {etaRegion}, pt threshold = {pt}, genPdgId = {genPdgId}, recoPdgId = {recoPdgId if recoPdgId else 'None'}")

# Analyzer
processor_instance = P2L1TAnalyzer(
    hist_config,
    cuts,
    pt=pt,
    genPdgId=genPdgId,
    recoPdgId=recoPdgId,
    etaRegion=etaRegion
)

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

# Save output as binary
outdir = 'coffea'
util.save(result, f"{outdir}/{outname}.coffea")
print(f"Saved output: {outdir}/{outname}.coffea")

