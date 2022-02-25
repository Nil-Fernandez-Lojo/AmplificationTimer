import json
from utility_functions import load_config, icgc
from .AmplificationTimerObjects import Amplification
from inference import inference_t

path_config = "../config.json"

config = load_config(path_config)

N_samples = config['summary_table'].shape[0]
with open(config['path_amplifications_file']) as json_file:
	amplifications_dict = json.load(json_file)

for amp_dict in amplifications_dict:
	if icgc(amp_dict['clinical_data']['samplename'], config):
		amplification = Amplification(amplification_dict=amp_dict, config=config)
		for subclonality_modelled in [False, True]:
			for filter_APOBEC in [False, True]:
				for model_index in range(1,5):
					inference_t(amplification,model_index,subclonality_modelled,filter_APOBEC)