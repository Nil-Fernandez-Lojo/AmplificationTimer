import sys
import copy
import numpy as np
from pathlib import Path
import pandas as pd 
from pyensembl import EnsemblRelease
import pickle

from sample import Sample
from inference import inference_t
from position import Position
from chromosome import Chromosome
from snv import SNV
import arviz as avz

sample_idx = int(sys.argv[1])
model_idx = int(sys.argv[2])
file_out = "fig_2_data/sample_idx_"+str(sample_idx)+"_model_"+str(model_idx)+".tsv"
expected_mut_list = [10,100,1000,10000]
amplification_idx = 0
list_t = np.linspace(0.0, 1.0, num=21)
subclonality_modelled = True
save = False
folder_models = Path("inferred_times_simulated")

def coverage(segment,rho,nrpcc):
	cn = segment.major_cn + segment.minor_cn
	n_reads = round(nrpcc*(cn + ((1-rho)/rho)*segment.get_ploidy_healthy()))
	return n_reads

def simulation_inference_using_prior(sample_,
	amplification_idx,
	model_idx,
	folder_models,
	t = None,
	mu = None,
	min_reads_detect_SNV = 0,
	subclonality_modelled_simulation=True,
	subclonality_modelled_inference=True,
	nrpcc = 'same',
	n_MCMC_iterations = 20000,
	cores = 2,
	suffix_name = "",
	save = True):

	# TODO: usage of constant nrpcc while actually should sample from distribution with that mean
	# TODO: put rest of attributes SNV
	# TODO: keep APOBEC like mutations?
	# TODO: actually use another prior for mutation rate

	#Some preprocessing
	sample = copy.deepcopy(sample_)
	rho = sample.clinical_data['purity']
	amplification = sample.amplifications[amplification_idx]
	for segment in sample.segments:
		segment.SNVs = []
	for segment in amplification.segments:
		segment.SNVs = []

	if nrpcc == "same":
		nrpcc = sample.clinical_data["nrpcc"]
		print("nrpcc", nrpcc)

	fSNV = sample.subclonal_structure["n_snvs"].values/sample.subclonal_structure["n_snvs"].sum()
	print("fSNV", fSNV)
	subclone_ccf = sample.subclonal_structure["fraction_cancer_cells"].values
	print("subclone_ccf",subclone_ccf)

	#Sample latent variables
	if mu is None:
		mu = np.random.beta(1,10**6)
	if t is None:
		t = np.random.uniform()
	print("mu",mu)
	print("t",t)

	if model_idx ==3:
		u = np.random.uniform()
	elif model_idx ==4:
		u = np.random.uniform(size = len(sample.amplifications[amplification_idx].segments))
	else:
		u = None
	print("u", u)
	# generate SNVs on 1+0 and 1+1 segments with NRPCC
	for segment in sample.segments:
		if segment.major_cn == 1 and segment.minor_cn <=1:
			if segment.minor_cn == 0:
				p = mu
				q = rho/(rho + segment.get_ploidy_healthy()*(1-rho))
			else:
				p = 2*mu
				q = rho/(2*rho + segment.get_ploidy_healthy()*(1-rho))
			D = coverage(segment,rho,nrpcc)
			n_SNVs = np.random.binomial(n=segment.get_length(), p = p)

			for i in range(n_SNVs):
				d = np.random.binomial(n=D, p = q)
				#SNV detected
				if d>=min_reads_detect_SNV:
					segment.add_snv(SNV(segment.chromosome,
                                        Position(segment.chromosome,segment.start.position+1,sample.config),  #TODO
                                        D,
                                        D - d,
						"",  #TODO
						"")) #TODO



	# generate SNVs on Amplified segments with NRPCC
	for i,segment in enumerate(amplification.segments):
		include_u = segment.minor_cn>=2
		p_all = mu*t
		if model_idx in [1,2]:
			p_1 = segment.major_cn*mu*(1-t)
		elif model_idx == 3:
			p_1 = segment.major_cn*mu*(1-t)+segment.minor_cn*mu*(u**include_u)
		else:
			p_1 = segment.major_cn*mu*(1-t)+segment.minor_cn*mu*(u[i]**include_u)

		q_clonal_1 = rho/((segment.major_cn + segment.minor_cn)*rho+segment.get_ploidy_healthy()*(1-rho)) # 1 copy clonal
		q_clonal_all = q_clonal_1 * segment.major_cn

		if subclonality_modelled_simulation and len(fSNV) >1:
			q_subclonal_1 = q_clonal_1*subclone_ccf[1:]

			p_all = p_all*fSNV[0]
			p_1 = p_1*fSNV[0]
			p_1_subclonal = mu*fSNV[1:]
			p_0 = 1-(p_all+p_1+sum(p_1_subclonal))
			p = np.concatenate(([p_all,p_1],p_1_subclonal,[p_0]))
		else:
			q_subclonal_1 = []

			p_0 = 1-p_all-p_1
			p = np.array([p_all,p_1,p_0])

		n = np.random.multinomial(segment.get_length(),pvals = p)
		print("n",n)
		D = coverage(segment,rho,nrpcc)
		print(q_clonal_all,q_clonal_1,q_subclonal_1,[q_clonal_all, q_clonal_1],np.array([q_clonal_all, q_clonal_1]), np.array(q_subclonal_1))
		q = np.concatenate(([q_clonal_all, q_clonal_1],q_subclonal_1))
		print("q",q)

		for snv_type_idx in range(len(n)-1):
			for snv_idx in range(n[snv_type_idx]):
				d = np.random.binomial(n=D, p = q[snv_type_idx])
				if d>=min_reads_detect_SNV:
					segment.add_snv(SNV(segment.chromosome,
                                        Position(segment.chromosome,segment.start.position+1,sample.config),  #TODO
                                        D - d,
                                        d,
						"",  #TODO
						"")) #TODO

	sample.mutation_rate= sample.get_mutation_rate(sample.segments)
	return inference_t(sample,
		amplification_idx,
		model_idx,
		subclonality_modelled_inference,
		n_MCMC_iterations,
		folder_models,
		cores,
		true_t = t,
		true_u = u,
		suffix_name = suffix_name,
		save = save)

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

name = config['summary_table']['samplename'][sample_idx]
s = Sample(config,name)
length_amplification = 0
for segment in s.amplifications[amplification_idx].segments:
	length_amplification += segment.get_length()
ploidy_amp = s.amplifications[amplification_idx].get_mean_ploidy()
print("ploidy_amp",ploidy_amp)
results = pd.DataFrame(np.zeros((len(expected_mut_list)*len(list_t), 6)), 
	columns = ['sample_idx',
	'amplification_idx',
	'model_idx',
	'expected_mutations',
	't',
	'inferred_t'])
results['sample_idx'] = sample_idx
results['amplification_idx'] = amplification_idx
results['model_idx'] = model_idx

i = 0
print('length_amplification', length_amplification)
print('ploidy_amp', ploidy_amp)
print('expected_mut', expected_mut_list[0])
print(expected_mut_list[0]/(length_amplification*ploidy_amp))
for expected_mut in expected_mut_list:
	for t in list_t:
		mu = expected_mut/(length_amplification*ploidy_amp)
		trace = simulation_inference_using_prior(s,
			amplification_idx,
			model_idx,
			folder_models,
			t = t,
			mu = mu,
			subclonality_modelled_simulation = subclonality_modelled,
			subclonality_modelled_inference = subclonality_modelled,
			suffix_name = "_"+str(i),
			save = save)
		summary = avz.summary(trace)
		inferred_t = summary[['mean']].T.reset_index(drop=True)['t'].loc[0]
		results['expected_mutations'].loc[i] = expected_mut
		results['t'].loc[i] = t
		results['inferred_t'].loc[i] = inferred_t
		results.to_csv(file_out,index=False,sep='\t')

		i+=1


