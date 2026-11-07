import os

cutConfig = "cut_config.cut_config_met_jet"
histo_config = "histo_config/histo_met_jet.yaml"

rootDir = "../ntuples/with_hgcID/"

outName = "torch_0p3_hgcIDscores_hdim8"

cmd = (
    f"python analyzer.py "
    f"--rootDir {rootDir} "
    f"--histoConfig {histo_config} "
    f"--cutConfig {cutConfig} "
    f"--outName {outName} "
)

print(cmd)
os.system(cmd)
