from AmplificationTimerObjects import Sample
from utility_functions import load_config
import json
# from inference import inference_t

"""
This script preprocesses the data
"""

path_config = "../config.json"

config = load_config(path_config)

amplifications = []
N_samples = config['summary_table'].shape[0]

for i in range(N_samples):
	print(i, config['summary_table']['samplename'][i])
	s = Sample(config,config['summary_table']['samplename'][i],save=True)
	amplifications.extend(a.to_dict() for a in s.amplifications)
	with open(config["path_amplifications_file"], 'w') as buff:
		json.dump(amplifications, buff, indent=4)