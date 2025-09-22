from coffea import util

# import coffea
output = util.load("coffea/baseline.coffea")
#print(output.keys())

#h = output['matched_gen_multiplicity']
h = output['matched_puppi_pdgId']
#h = output['pf_multiplicity']

print(h)
