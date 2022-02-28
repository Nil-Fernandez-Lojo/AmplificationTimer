import json
from utility_functions import load_config, icgc
from AmplificationTimerObjects import Amplification
from inference import inference_t

path_config = "../config.json"

config = load_config(path_config, load_genome_into_config=False)

amplifications_dict = []
for file in config['path_amplifications_folder'].glob('*.json'):
	with open(file) as json_file:
		amplifications_dict.append(json.load(json_file))

for amp_dict in amplifications_dict:
	if icgc(amp_dict['clinical_data']['samplename'], config):
		amplification = Amplification(amplification_dict=amp_dict, config=config)
		for subclonality_modelled in [False, True]:
			for filter_APOBEC in [False, True]:
				for only_clock_like_SNVs in [False, True]:
					for model_index in range(1,5):
						inference_t(amplification,
									model_index,
									subclonality_modelled,
									filter_APOBEC,
									only_clock_like_SNVs,
									save=True)