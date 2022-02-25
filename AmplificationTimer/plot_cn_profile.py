# author: Nil Fernandez Lojo
from decorator_equality import class_equality_attributes
import pandas as pd
from pathlib import Path
from sample import Sample
from chromosome import Chromosome
from pyensembl import EnsemblRelease
import numpy as np 
import json

name = 'a1e3dc5b-b81f-4890-870c-ed3b8ac36dec'
add_snvs = False
total_cn=True
chromosome=12
start_pos= 60000000
max_bases = 20000000
plot_threshold_amp = True
path_save = None
oncogene = [(69201952,69244466,'MDM2')] #position in the chr and name
#oncogene = []
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
path_amplifications_all = "preprocessed_data/amplifications_all.json"

suffix_cnas = ".consensus.20170119.somatic.cna.annotated.txt"
suffix_snvs = ".consensus.20160830.somatic.snv_mnv.vcf"
suffix_subclonal_structure = "_subclonal_structure.txt"

vcf_column_names = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
vcf_prefix_counts = {"ref_count": "t_ref_count=", "alt_count" : "t_alt_count="}
#Path out
folder_preprocessed_data = Path('preprocessed_data')
amplifications_file = folder_preprocessed_data/"amplifications.json"
suffix_snvs = ".consensus.20160830.somatic.snv_mnv.vcf"
folder_snv_white = path_data/"final_consensus_snv_indel_passonly_icgc.public" / "snv_mnv"
folder_snv_gray = path_data/"final_consensus_snv_indel_passonly_icgc.public" / "graylist" / "snv_mnv"

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

#TODO: should change place
config['cum_length_chr'] = np.zeros(25)
cum_sum = 0
for i,chrom in enumerate(list(range(1,23)) + ['X','Y']):
	lengths = config['chromosome_arm_length']['length'][ config['chromosome_arm_length']['chromosome'] == str(chrom)]
	config['cum_length_chr'][i] = cum_sum
	cum_sum += lengths.sum()
config['cum_length_chr'][-1] = cum_sum
N_samples = config['summary_table'].shape[0]


def icgc(samplename):
	file_white = folder_snv_white/(samplename+suffix_snvs)
	file_gray = folder_snv_gray/(samplename+suffix_snvs)
	return file_gray.exists() or file_white.exists()

def filer_amplifications(amplifications):
	amplifications_filtered = []
	for a in amplifications:
		samplename = a['clinical_data']['samplename']
		if icgc(samplename):
			amplifications_filtered.append(a)
	return amplifications_filtered

with open(path_amplifications_all, 'r') as fp:
	amplifications = json.load(fp)

s = Sample(config,name)
# for i,seg in enumerate(s.segments): 
# 	if seg.minor_cn == 0:
# 		seg.minor_cn = 1

# for i,seg in enumerate(s.segments):
# 	print(i,seg.chromosome.c, seg.start, seg.minor_cn)
s.plot_CN_profile(add_snvs=add_snvs,
		total_cn=total_cn,
		chromosome=chromosome,
		start_pos=start_pos,
		max_bases = max_bases,
		plot_threshold_amp = plot_threshold_amp,
		title = name,
		path_save = path_save,
		add_gene = oncogene)

# amplifications_filtered = filer_amplifications(amplifications)

# for i,a in enumerate(amplifications_filtered):
# 	name = a['clinical_data']['samplename']
# 	s = Sample(config,name)
# 	chromosome = a['segments'][0]['chromosome']
# 	margin = 0.1 * (a['segments'][-1]['end'] - a['segments'][0]['start'])
# 	start_pos = a['segments'][0]['start'] - margin
# 	max_bases = (a['segments'][-1]['end'] - a['segments'][0]['start']) + 2*margin
# 	print(i,'/',len(amplifications_filtered), name, chromosome)
# 	s.plot_CN_profile(add_snvs=add_snvs,
# 		total_cn=total_cn,
# 		chromosome=chromosome,
# 		start_pos=start_pos,
# 		max_bases = max_bases,
# 		plot_threshold_amp = plot_threshold_amp,
# 		title = name + ' chromosome '+chromosome,
# 		path_save = 'figures/plots_CNP_amps/'+name+'_chr_'+chromosome+'.png')