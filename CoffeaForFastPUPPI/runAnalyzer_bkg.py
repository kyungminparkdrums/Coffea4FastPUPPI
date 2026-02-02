import os

histo_config = "histo_config/histo_bkg_tkele.yaml"
cutConfig = "cut_config.cut_config_bkg_tkele"

rootDir = "../nugun_part.root" 
outName = "test_bkg_tkele"

cmd = (
    f"python analyzer.py "
    f"--rootDir {rootDir} "
    f"--histoConfig {histo_config} "
    f"--cutConfig {cutConfig} "
    f"--outName {outName} "
)

print(cmd)
os.system(cmd)
