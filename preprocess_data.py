# author: Nil Fernandez Lojo
from decorator_equality import class_equality_attributes #TODO remove
import pandas as pd
from pathlib import Path
from sample import Sample
from chromosome import Chromosome
#from inference import inference_t
from pyensembl import EnsemblRelease
import numpy as np 
import json


#Paths input data
path_data = Path("../memoire/data")
path_chromosome_arm_length = "hg19_chr_arm_length.tsv"
path_summary_table = path_data / "summary_table_combined_annotations_v4.txt" 
path_ensembl_to_entrez= Path("ensembl_to_entrez.tsv")
path_oncognes = path_data / "Census_allTue Nov 24 18_02_37 2020.csv"
folder_cna = path_data / "consensus.20170119.somatic.cna.annotated"
folder_snv_white = path_data / "final_consensus_snv_indel_passonly_icgc.public" / "snv_mnv"
folder_snv_gray = path_data / "final_consensus_snv_indel_passonly_icgc.public" / "graylist"/ "snv_mnv"
folder_subclonal_structure = Path("data/subclonal.structures")


suffix_cnas = ".consensus.20170119.somatic.cna.annotated.txt"
suffix_snvs = ".consensus.20160830.somatic.snv_mnv.vcf"
suffix_subclonal_structure = "_subclonal_structure.txt"

vcf_column_names = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
vcf_prefix_counts = {"ref_count": "t_ref_count=", "alt_count" : "t_alt_count="}
#Path out
folder_preprocessed_data = Path('preprocessed_data')
amplifications_file = folder_preprocessed_data/"amplifications.json"

cores = 10

config = dict()
config['vcf_column_names'] = vcf_column_names
config['vcf_prefix_counts'] = vcf_prefix_counts
config['folder_cna'] = folder_cna
config['folder_snv_white'] = folder_snv_white
config['folder_snv_gray'] = folder_snv_gray
config['suffix_cnas'] = suffix_cnas
config['suffix_snvs'] = suffix_snvs
config['suffix_subclonal_structure'] = suffix_subclonal_structure
config['folder_preprocessed_data'] = folder_preprocessed_data
config['folder_subclonal_structure'] = folder_subclonal_structure
config['summary_table'] = pd.read_csv(path_summary_table, sep='\t')
config['ensembl_to_entrez'] = pd.read_csv(path_ensembl_to_entrez, sep='\t') #TODO what to do with duplicates
config['ensembl'] = EnsemblRelease(75) # TODO: I need to justify why version 75
config['chromosome_arm_length'] = pd.read_csv(path_chromosome_arm_length, sep='\t')
config['chromosome_arm_length'].rename(columns={'chrom': 'chromosome'},inplace=True)
config['oncogenes'] = pd.read_csv(path_oncognes)
N_samples = config['summary_table'].shape[0]

#TODO: should change place
config['cum_length_chr'] = np.zeros(25)
cum_sum = 0
for i,chromosome in enumerate(list(range(1,23)) + ['X','Y']):
	lengths = config['chromosome_arm_length']['length'][ config['chromosome_arm_length']['chromosome'] == str(chromosome)]
	config['cum_length_chr'][i] = cum_sum
	cum_sum += lengths.sum()
config['cum_length_chr'][-1] = cum_sum


n_MCMC_iterations = 10000
folder_models = Path("inferred_times")
amplifications = []
for i in range(N_samples):
	print(i, config['summary_table']['samplename'][i])
	s = Sample(config,config['summary_table']['samplename'][i],save=False)
	amplifications.extend(a.to_dict() for a in s.amplifications)
	with open(amplifications_file, 'w') as buff:
		json.dump(amplifications, buff, indent=4)
	continue
	if len(s.amplifications)>0 and (s.snv_type in ["white", "gray"]):
	 	for j in range(len(s.amplifications)):
	 		for subclonality_modelled in [False, True]:
	 			for model_index in range(1,5):
	 				inference_t(s,j,model_index,subclonality_modelled,n_MCMC_iterations,folder_models,cores)