# author: Nil Fernandez Lojo
from decorator_equality import class_equality_attributes #TODO remove
import pandas as pd
from pathlib import Path
from sample import Sample
from chromosome import Chromosome
from inference import inference_t
from pyensembl import EnsemblRelease
import numpy as np 
import json
import arviz as avz
import pickle

list_samplenames = ['a3914a6c-c622-11e3-bf01-24c6515278c0', '7c405ca0-c622-11e3-bf01-24c6515278c0', '98d27916-c622-11e3-bf01-24c6515278c0', '3b41cb48-c623-11e3-bf01-24c6515278c0', '7fba5aac-c622-11e3-bf01-24c6515278c0', '4a703d3e-c623-11e3-bf01-24c6515278c0', '819b4304-c622-11e3-bf01-24c6515278c0', '2572b0bc-c622-11e3-bf01-24c6515278c0', '850389d4-c622-11e3-bf01-24c6515278c0', 'd60f880a-c622-11e3-bf01-24c6515278c0', '4fdc8980-c623-11e3-bf01-24c6515278c0', 'ec5e2990-c622-11e3-bf01-24c6515278c0', '4b8943be-c623-11e3-bf01-24c6515278c0', '6622f932-c622-11e3-bf01-24c6515278c0', '030695f6-c623-11e3-bf01-24c6515278c0', '7260f57c-c623-11e3-bf01-24c6515278c0']
subclonality_modelled = True
model_index_list = [1,4]
filter_APOBEC = False
chr_amp = 8
chr_arm_amp = 'q'

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
summarry_inference_1 = dict()
summarry_inference_4 = dict()

for samplename in list_samplenames:
	print(samplename)
	s = Sample(config,samplename,save=False)
	for j,a in enumerate(s.amplifications):
		if (a.chromosome.c == chr_amp) and (a.arm == chr_arm_amp):
			for model_index in model_index_list:
				summary = inference_t(s,j,model_index,subclonality_modelled,n_MCMC_iterations,folder_models,cores,filter_APOBEC,save = False)
				if model_index == 1:
					if samplename not in summarry_inference_1.keys():
						summarry_inference_1[samplename] = dict()
					summarry_inference_1[samplename][j] = summary
					pickle.dump(summarry_inference_1,open("summary_no_filter_APOBEC_model_1.pkl", "wb" ))
				else:
					if samplename not in summarry_inference_4.keys():
						summarry_inference_4[samplename] = dict()
					summarry_inference_4[samplename][j] = summary
					pickle.dump(summarry_inference_4,open("summary_no_filter_APOBEC_model_4.pkl", "wb" ))
					