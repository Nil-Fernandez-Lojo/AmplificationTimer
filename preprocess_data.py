# author: Nil Fernandez Lojo
import pandas as pd
from pathlib import Path
from Sample import Sample
from chromosome import Chromosome
from pyensembl import EnsemblRelease
from models import Model1, Model2,Model3

import numpy as np 
import arviz as az
import matplotlib.pyplot as plt


#Paths input data
path_data = Path("../memoire/data")
path_chromosome_arm_length = "hg19_chr_arm_length.tsv"
path_summary_table = path_data / "summary_table_combined_annotations_v4.txt" 
path_ensembl_to_entrez= Path("ensembl_to_entrez.tsv")
folder_cna = path_data / "consensus.20170119.somatic.cna.annotated"
folder_snv_green = path_data / "final_consensus_snv_indel_passonly_icgc.public" / "snv_mnv"
folder_snv_gray = path_data / "final_consensus_snv_indel_passonly_icgc.public" / "graylist"/ "snv_mnv"

suffix_cnas = ".consensus.20170119.somatic.cna.annotated.txt"
suffix_snvs = ".consensus.20160830.somatic.snv_mnv.vcf"

#Path out
folder_preprocessed_data = Path('preprocessed_data')


config = dict()
config['folder_cna'] = folder_cna
config['folder_snv_green'] = folder_snv_green
config['folder_snv_gray'] = folder_snv_gray
config['suffix_cnas'] = suffix_cnas
config['suffix_snvs'] = suffix_snvs
config['suffix_snvs'] = suffix_snvs
config['folder_preprocessed_data'] = folder_preprocessed_data
config['summary_table'] = pd.read_csv(path_summary_table, sep='\t')
config['ensembl_to_entrez'] = pd.read_csv(path_ensembl_to_entrez, sep='\t') #TODO what to do with duplicates
config['ensembl'] = EnsemblRelease(75) # TODO: I need to justify why version 75
config['chromosome_arm_length'] = pd.read_csv(path_chromosome_arm_length, sep='\t')
config['chromosome_arm_length'].rename(columns={'chrom': 'chromosome'},inplace=True)

N_samples = config['summary_table'].shape[0]
for i in range(13,N_samples):	
	print(i, config['summary_table']['samplename'][i])
	s = Sample(config,config['summary_table']['samplename'][i],save=False)
	x = Model1(s.amplifications[0],s.mutation_rate,s.clinical_data)
	t = np.linspace(0,1,100)
	plt.plot(t,x.get_analytical_posterior(t))
	print(x.get_approximate_MAP())
	trace = x.get_MCMC_samples_posterior(20000)
	az.plot_trace(trace)
	print(az.summary(trace, round_to=2))
	print('Model 2')
	x2 = Model2(s.amplifications[0],s.mutation_rate,s.clinical_data)
	print(x2.get_approximate_MAP())
	trace2 = x2.get_MCMC_samples_posterior(20000)
	az.plot_trace(trace2)
	print(az.summary(trace2, round_to=2))
	x3 = Model3(s.amplifications[0],s.mutation_rate,s.clinical_data)
	print(x3.get_approximate_MAP())
	trace3 = x3.get_MCMC_samples_posterior(20000)
	az.plot_trace(trace3)
	print(az.summary(trace3, round_to=2))
	plt.show()
	exit()