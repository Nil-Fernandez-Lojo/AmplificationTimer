import json
import pandas as pd
import numpy as np
from pathlib import Path
from pyensembl import EnsemblRelease

def load_config(config_file_path):
    with open(config_file_path) as json_file:
        config = json.load(json_file)
    for key, value in config.items():
        if key.startswith('path'):
            config[key] = Path(value)

    config['summary_table'] = pd.read_csv(config['path_summary_table'], sep='\t')
    config['ensembl_to_entrez'] = pd.read_csv(config['path_ensembl_to_entrez'], sep='\t') #TODO what to do with duplicates
    config['ensembl'] = EnsemblRelease(75) # TODO: I need to justify why version 75
    config['chromosome_arm_length'] = pd.read_csv(config['path_chromosome_arm_length'], sep='\t')
    config['chromosome_arm_length'].rename(columns={'chrom': 'chromosome'},inplace=True)
    config['oncogenes'] = pd.read_csv(config['path_oncognes'])
    config['cum_length_chr'] = np.zeros(25)

    cum_sum = 0
    for i,chromosome in enumerate(list(range(1,23)) + ['X','Y']):
        lengths = config['chromosome_arm_length']['length'][config['chromosome_arm_length']['chromosome'] == str(chromosome)]
        config['cum_length_chr'][i] = cum_sum
        cum_sum += lengths.sum()
    config['cum_length_chr'][-1] = cum_sum

    return config

