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

args = parse_arguments_inference()
path_config = "../config.json"
config = load_config(path_config, load_genome_into_config=False)

amplifications_dict = []
for file in config['path_filtered_amplifications_folder'].glob('*.json'):
    with open(file) as json_file:
        amplification = Amplification(amplification_dict=json.load(json_file), config=config)

    inference_t(amplification,
                args.model_index,
                subclonality_modelled=args.subclonality_modelled,
                filter_APOBEC=args.filter_APOBEC,
                only_clock_like_SNVs = args.only_clock_like_SNVs,
                mu_ml=args.mu_ml,
                save=True)
