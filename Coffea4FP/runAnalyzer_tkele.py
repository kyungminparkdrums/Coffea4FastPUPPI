import os

histo_config = "histo_config/histo_tkeleLeadSub.yaml"
#histo_config = "histo_config/histo_tkeleLead.yaml"
#histo_config = "histo_config/histo_tkelePair.yaml"
#histo_config = "histo_config/histo_tkele.yaml"

cutConfig = "cut_config.cut_config_tkele"

mass_points = [2, 5, 10, 15, 20, 30]

for mass in mass_points:
    rootDir = f"../m{mass}_part.root" 
    outName = f"tkele_{mass}"

    cmd = (
        f"python analyzer.py "
        f"--rootDir {rootDir} "
        f"--histoConfig {histo_config} "
        f"--cutConfig {cutConfig} "
        f"--outName {outName} "
    )

    print(cmd)
    os.system(cmd)

