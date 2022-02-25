import pandas as pd
import numpy as np 
import json
from pathlib import Path


#Paths input data
path_data = Path("../memoire/data")
path_summary_table = path_data / "summary_table_combined_annotations_v4.txt" 
folder_cna = path_data / "consensus.20170119.somatic.cna.annotated"
folder_snv_white = path_data / "final_consensus_snv_indel_passonly_icgc.public" / "snv_mnv"
folder_snv_gray = path_data / "final_consensus_snv_indel_passonly_icgc.public" / "graylist"/ "snv_mnv"
suffix_snvs = ".consensus.20160830.somatic.snv_mnv.vcf"
folder_preprocessed_data = Path('preprocessed_data')
summary_table= pd.read_csv(path_summary_table, sep='\t')

def icgc(samplename):
	file_white = folder_snv_white/(samplename+suffix_snvs)
	file_gray = folder_snv_gray/(samplename+suffix_snvs)
	return file_gray.exists() or file_white.exists()

count_SNV_1_1 = 0
count_SNV_1_0 = 0
n_icgc = 0
for samplename in summary_table['samplename']:
	if icgc(samplename):
		n_icgc += 1
		path_file = folder_preprocessed_data/(samplename+'.json')
		with open(path_file, 'r') as fp:
			x = json.load(fp)
		count_SNV_1_1 += x['mutation_rate']['n_1_1']
		count_SNV_1_0 += x['mutation_rate']['n_1_0']
print(count_SNV_1_0,'count_SNV_1_0')
print(count_SNV_1_1,'count_SNV_1_1')
print(n_icgc,'n_icgc')
print('average 1 1', count_SNV_1_1/n_icgc)
print('average 1 1', count_SNV_1_0/n_icgc)
print('average both', (count_SNV_1_0+count_SNV_1_1)/n_icgc)

