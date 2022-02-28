from AmplificationTimerObjects import Sample
from utility_functions import load_config
import json

"""
This script preprocesses the data
"""

path_config = "../config.json"

config = load_config(path_config)

amplifications = []
N_samples = config['summary_table'].shape[0]

for i in range(N_samples):
	samplename = config['summary_table']['samplename'][i]
	print(i, samplename)
	s = Sample(config,samplename,save=True)
	for amplification in s.amplifications:
		file_name = samplename + "_"+str(amplification.chromosome) + amplification.arm + ".json"
		path_amplification = config["path_amplifications_folder"]/file_name
		with open(path_amplification, 'w') as buff:
			json.dump(amplification.to_dict(), buff, indent=4)
