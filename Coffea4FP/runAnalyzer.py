import os

#histo_config = "histo_config/histo_pfpuppi.yaml"
histo_config = "histo_config/histo_tkele.yaml"

cutConfig = "cut_config.cut_config"

#rootDir = "../TT_PU200_noIDcut.root"
rootDir = "../perfNano_10044363_23.root"

#outName = "test"
outName = "test_tkele"

cmd = (
    f"python analyzer.py "
    f"--rootDir {rootDir} "
    f"--histoConfig {histo_config} "
    f"--cutConfig {cutConfig} "
    f"--outName {outName} "
)

print(cmd)
os.system(cmd)
