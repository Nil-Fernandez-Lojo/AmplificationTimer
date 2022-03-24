# author: Nil Fernandez Lojo
from decorator_equality import class_equality_attributes #TODO remove
import pandas as pd
from pathlib import Path
from sample import Sample
from chromosome import Chromosome
from inderence_dev import inference_t
from pyensembl import EnsemblRelease
import numpy as np 
import json

# TODO: refactor
add_snvs = True
plot_threshold_amp = False
total_cn = False
filter_APOBEC = False
list_samplenames = ['4cbe411b-b05e-46bd-bea8-126289a0866c']
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
for samplename in list_samplenames:
	print(samplename)
	s = Sample(config,samplename,save=False)
	for i,a in enumerate(s.amplifications):
		if (a.chromosome.c in [8,17]):
			margin = 0.1 * (a.segments[-1].end.position - a.segments[0].start.position)
			start_pos = a.segments[0].start.position - margin
			max_bases = a.segments[-1].end.position - a.segments[0].start.position + 2*margin

			s.plot_cn(add_snvs=add_snvs,
                      total_cn=total_cn,
                      chromosome=a.chromosome,
                      start_pos=start_pos,
                      max_bases = max_bases,
                      plot_threshold_amp = plot_threshold_amp,
                      title = samplename + ' chromosome '+str(a.chromosome.c)+' not filtered',
                      path_save = 'figures_CNP_LIRI_MYC_filtering/'+samplename + ' chromosome '+str(a.chromosome.c)+'_not_filtered.png')
			subclonality_modelled = True
			model_index = 4
			inference_t(s,i,model_index,subclonality_modelled,n_MCMC_iterations,folder_models,cores,filter_APOBEC = filter_APOBEC)

			for segment_a in a.segments:
				for segment in s.segments:
					if (segment_a.chromosome.c == segment.chromosome.c) and (segment_a.start.position == segment.start.position):
						segment.snvs = segment_a.snvs
						break

			s.plot_cn(add_snvs=add_snvs,
                      total_cn=total_cn,
                      chromosome=a.chromosome,
                      start_pos=start_pos,
                      max_bases = max_bases,
                      plot_threshold_amp = plot_threshold_amp,
                      title = samplename + ' chromosome '+str(a.chromosome.c)+' filtered',
                      path_save = 'figures_CNP_LIRI_MYC_filtering/'+samplename + ' chromosome '+str(a.chromosome.c)+'_filtered.png')
			break
	 	
