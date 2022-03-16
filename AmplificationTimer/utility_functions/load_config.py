import json
import pandas as pd
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
    # TODO: drivers per cancer type
    known_drivers = pd.read_csv(config['path_known_drivers_PCAWG'], sep='\t')
    cancer_types = known_drivers['ttype'].unique()
    config['drivers'] = dict()
    for c in cancer_types:
        config['drivers'][c] = {'gain':set(), "rest":set()}

    for index, row in known_drivers.iterrows():
        if row['category'] == 'coding_amplification':
            config['drivers'][row['ttype']]['gain'].add(row['gene'])
        else:
            config['drivers'][row['ttype']]['rest'].add(row['gene'])



    # a = set(config['summary_table']['histology_pcawg'].unique())
    # b = set(config['drivers'].keys())
    # print(a.difference(b))
    # print(b.difference(a))

    # sample_ids_in_driver_file = known_drivers['sample_id'].unique()
    # for index, row in config['summary_table'].iterrows():
    #     if row['samplename'] in sample_ids_in_driver_file:
    #         cancer_type_driver = known_drivers['ttype'][known_drivers['sample_id'] == row['samplename']].unique()
    #         if len(cancer_type_driver)!=1:
    #             print('for same sample different cancer types', cancer_type_driver)
    #         else:
    #             cancer_type_summary = map_summary_to_driver_table.get(row['histology_pcawg'],row['histology_pcawg'])
    #             if cancer_type_driver[0] != cancer_type_summary:
    #                 print(row['samplename'], cancer_type_driver, cancer_type_summary)
    #


        # for summary_not, driver_not in map_summary_to_driver_table.items():
        #     samples_summary = set(
        #         config['summary_table']['samplename'][config['summary_table']['histology_pcawg'] == summary_not])
        #     samples_driver = set(known_drivers['sample_id'][known_drivers['ttype'] == driver_not].unique())
        #     print(summary_not, driver_not, 'n_samples', len(samples_summary))
        #     print(samples_summary.difference(samples_driver), len(samples_summary.difference(samples_driver)))
        #     print(samples_driver.difference(samples_summary), len(samples_driver.difference(samples_summary)))
        #     print()


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

