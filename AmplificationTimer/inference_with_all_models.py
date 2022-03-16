import os

os.environ["OMP_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"
os.environ["VECLIB_MAXIMUM_THREADS"] = "2"
os.environ["NUMEXPR_NUM_THREADS"] = "2"

import json
from utility_functions import load_config, parse_arguments_inference
from AmplificationTimerObjects import Amplification
from inference import inference_t

path_config = "../config.json"
config = load_config(path_config, load_genome_into_config=False)
args = parse_arguments_inference()

amplifications_dict = []
for file in config['path_filtered_amplifications_folder'].glob('*.json'):
    with open(file) as json_file:
        amplifications_dict.append(json.load(json_file))

for amp_dict in amplifications_dict:
    amplification = Amplification(amplification_dict=amp_dict, config=config)
    inference_t(amplification,
                args.model_index,
                args.subclonality_modelled,
                args.filter_APOBEC,
                args.only_clock_like_SNVs,
                save=True)
