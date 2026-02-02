import os

histo_config = "histo_config/histo_pfpuppi.yaml"
cutConfig = "cut_config.cut_config_pf"

rootDir = "../TT_PU200_noIDcut.root"

outName = "test"

cmd = (
    f"python analyzer.py "
    f"--rootDir {rootDir} "
    f"--histoConfig {histo_config} "
    f"--cutConfig {cutConfig} "
    f"--outName {outName} "
)

print(cmd)
os.system(cmd)
