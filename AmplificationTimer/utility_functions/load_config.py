import json
import pandas as pd
import copy
import numpy as np
from pathlib import Path
from pyensembl import EnsemblRelease

def load_genome(fasta_file):
    with open(fasta_file, "r") as f:
        fl = f.readlines()
    genome = dict()
    chr_number = 0
    chr_name = None
    chr_seq = ""
    for line in fl:
        if line[0] == ">":
            if (chr_number!=0):
                genome[chr_name] = chr_seq
            chr_name = line[1:-1]
            chr_seq = ''
            chr_number+=1
        else:
            chr_seq+=line[:-1]
    genome[chr_name] = chr_seq
    return genome

def load_drivers(config):
    known_drivers = pd.read_csv(config['path_known_drivers_PCAWG'], sep='\t')
    cancer_types = known_drivers['ttype'].unique()
    config['drivers'] = dict()
    for c in cancer_types:
        config['drivers'][c] = {'gain':dict(), "rest":dict()}

    for index, row in known_drivers.iterrows():
        type_mutation = 'gain' if row['category'] == 'coding_amplification' else 'rest'
        dict_to_add = config['drivers'][row['ttype']][type_mutation]
        dict_to_add[row['gene']] = dict_to_add.get(row['gene'],0) + 1


def load_config(config_file_path, load_genome_into_config=True):
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

    tmp = copy.copy(config['vcf_column_types'])
    config['vcf_column_types'] = dict()
    for i in range(len(config['vcf_column_names'])):
        config['vcf_column_types'][config['vcf_column_names'][i]] = tmp[i]

    cum_sum = 0
    for i,chromosome in enumerate(list(range(1,23)) + ['X','Y']):
        lengths = config['chromosome_arm_length']['length'][config['chromosome_arm_length']['chromosome'] == str(chromosome)]
        config['cum_length_chr'][i] = cum_sum
        cum_sum += lengths.sum()
    config['cum_length_chr'][-1] = cum_sum
    if load_genome_into_config:
        config['genome'] = load_genome(config['path_fasta_genome'])

    load_drivers(config)

    return config

