import json
from pathlib import Path
import pandas as pd
import pickle
import arviz as avz
import matplotlib.pyplot as plt
import numpy as np

path_amplifications_all = "preprocessed_data/amplifications_all.json"

model_idx = 1
samplename = '1b06afe2-c623-11e3-bf01-24c6515278c0'
cancer_type = "LIRI"
oncogene = 'MYC'
path_data = Path("../memoire/data")
suffix_snvs = ".consensus.20160830.somatic.snv_mnv.vcf"
folder_snv_white = path_data/"final_consensus_snv_indel_passonly_icgc.public" / "snv_mnv"
folder_snv_gray = path_data/"final_consensus_snv_indel_passonly_icgc.public" / "graylist" / "snv_mnv"
path_summary_table = path_data / "summary_table_combined_annotations_v4.txt" 
file_statistics_amp = 'statistics_amplifications.tsv'

def icgc(samplename):
	file_white = folder_snv_white/(samplename+suffix_snvs)
	file_gray = folder_snv_gray/(samplename+suffix_snvs)
	return file_gray.exists() or file_white.exists()

def filer_amplifications(amplifications,cancer_type,oncogene=None):
	amplifications_filtered = []
	for a in amplifications:
		samplename = a['clinical_data']['samplename']
		if not a['clinical_data']['is_preferred']: 
			continue
		if cancer_type != 'all':
			if summary_table['cancer_type'].loc[summary_table['samplename']==samplename].values != cancer_type:
				continue
		if oncogene is not None:
			gene_present = False
			for gene in a['oncogenes']:
				if gene['name'] == oncogene:
					gene_present = True
					break
			if not gene_present:
				continue
		if icgc(samplename):
			amplifications_filtered.append(a)
	return amplifications_filtered

def count_patients_cancer_type(summary_table,cancer_type):
	count = 0
	for j in range(summary_table.shape[0]):
		if summary_table['cancer_type'].loc[j]!= cancer_type:
			continue
		if not summary_table['is_preferred'].loc[j]:
			continue
		if icgc(summary_table['samplename'].loc[j]):
			count +=1
	return count

def get_count_oncogene(amplifications_filtered):
	count_oncogene = {}
	for a in amplifications_filtered:
		for gene in a['oncogenes']:
			count_oncogene[gene['name']] = count_oncogene.get(gene['name'],0)+1
	return count_oncogene

def get_summary_amplifications(summary_table,amplifications, file_statistics_amp,save = True):
	summary = pd.DataFrame(index= range(len(summary_table['cancer_type'].unique())),
		columns=["cancer_type",
		"number_patients", 
		"number_patients_with_amp", 
		"most_common_oncogene", 
		"number_amp_oncogene"])

	for i,cancer_type in enumerate(summary_table['cancer_type'].unique()):
		summary["cancer_type"].loc[i] = cancer_type
		summary["number_patients"].loc[i] = count_patients_cancer_type(summary_table,cancer_type)
		amplifications_filtered = filer_amplifications(amplifications,cancer_type)
		patients_with_amp = set([a['clinical_data']['samplename'] for a in amplifications_filtered])
		summary["number_patients_with_amp"].loc[i] = len(patients_with_amp)
		count_oncogene = get_count_oncogene(amplifications_filtered)
		if len(count_oncogene) == 0:
			summary["most_common_oncogene"].loc[i] = None
			summary["number_amp_oncogene"].loc[i] = 0
		else:
			summary["most_common_oncogene"].loc[i] = max(count_oncogene, key=count_oncogene.get)
			summary["number_amp_oncogene"].loc[i] = count_oncogene[summary["most_common_oncogene"].loc[i]]

	summary = summary.astype({'number_patients_with_amp': 'int', "number_amp_oncogene": "int", "number_patients": "int"})

	summary["f_patient_amplification"] = summary["number_patients_with_amp"]/summary["number_patients"]
	summary["f_patient_oncogene"] = summary["number_amp_oncogene"]/summary["number_patients"]
	summary["f_amp_oncogene"] = summary["number_amp_oncogene"]/summary["number_patients_with_amp"]
	
	if save:
		summary.to_csv(file_statistics_amp,index=False, sep = '\t')

def get_timing_oncogene(oncogene, cancer_type,amplifications,summary_traces,summary_table,list_samples):
	timing = []
	amplifications_filtered = filer_amplifications(amplifications,cancer_type)
	amplification_id = 0
	for i,a in enumerate(amplifications_filtered):
		if i == 0 or (a['clinical_data']['samplename'] != amplifications_filtered[i-1]['clinical_data']['samplename']):
			amplification_id = 0
		else:
			amplification_id +=1
		if list_samples is not None:
			if a['clinical_data']['samplename'] not in list_samples:
				continue
		for gene in a['oncogenes']:
			if gene['name'] == oncogene:
				posterior_t = summary_traces[a['clinical_data']['samplename']][amplification_id]
				timing.append([a['clinical_data']['samplename'],
					posterior_t['mean']['t'],
					posterior_t['hdi_3%']['t'],
					posterior_t['hdi_97%']['t']])
				break
	timing = pd.DataFrame(timing, columns = ['samplename', 'mean', 'hdi_3', 'hdi_97'])
	timing = timing.sort_values(by=['mean'])
	return timing

def get_figure_3(oncogene,cancer_type,amplifications,summary_traces_4,summary_traces_1,summary_table,include_samplename=False):
	timing_4 = get_timing_oncogene(oncogene, cancer_type,amplifications,summary_traces_4,summary_table)
	timing_1 = get_timing_oncogene(oncogene, cancer_type,amplifications,summary_traces_1,summary_table)
	print(timing_1['samplename'].tolist())
	exit()
	timing_1 = timing_1.set_index('samplename')
	timing_1 = timing_1.reindex(index=timing_4['samplename'])
	timing_1 = timing_1.reset_index()
	title = cancer_type + " " +oncogene
	if include_samplename:
		y_ticks_label = timing_1['samplename']
	else:
		y_ticks_label = []
	ylabel = 'sample'
	plot_figure_3(timing_1,timing_4,title,y_ticks_label,ylabel)

def plot_figure_3(timing_1,timing_4,title,y_ticks_label,ylabel):

	plt.figure()
	plt.plot(timing_4['mean'],
		np.arange(timing_4.shape[0])+0.1,
		marker = '.', 
		linewidth=0,
		color='b', 
		markersize=5,
		label='Model 4')
	plt.plot(timing_1['mean'],
		np.arange(timing_1.shape[0])-0.1,
		marker = '.', 
		linewidth=0,
		color='r', 
		markersize=5,
		label='Model 1')

	for i in range(timing_1.shape[0]):
		plt.plot([timing_4['hdi_3'].iloc[i],timing_4['hdi_97'].iloc[i]],
			[i+0.1,i+0.1],
			color='b')
		plt.plot([timing_1['hdi_3'].iloc[i],timing_1['hdi_97'].iloc[i]],
			[i-0.1,i-0.1],
			color='r')
	plt.axvline(x=0.5,color='k',  linewidth=2)
	plt.yticks(np.arange(timing_4.shape[0]),y_ticks_label)
	plt.xticks(np.linspace(0,1,11))

	plt.title(title,fontsize=18)
	plt.xlim(0, 1)
	plt.xlabel('Inferred molecular time')
	plt.ylabel(ylabel)
	plt.text(0.2,timing_4.shape[0]-1.5,"Early",fontsize=16)
	plt.text(0.7,timing_4.shape[0]-1.5,"Late",fontsize=16)
	plt.legend(loc = 'lower right')
	plt.grid(axis='x')

def get_amplifications_sample(amplifications, samplename):
	amplifications_sample = []
	for a in amplifications:
		if a['clinical_data']['samplename'] == samplename:
			amplifications_sample.append(a)
	return amplifications_sample

def plot_timing_one_sample(samplename,summary_traces_1,summary_traces_4,amplifications_filtered):
	t_1 = summary_traces_1[samplename]
	t_4 = summary_traces_4[samplename]
	amplifications_sample = get_amplifications_sample(amplifications_filtered, samplename)

	timing_1 = []
	timing_4 = []
	for i,a in enumerate(amplifications_sample):
		chr_arm = a['chromosome']+a['arm']
		posterior_t_1 = t_1[i]
		posterior_t_4 = t_4[i]
		timing_1.append([chr_arm, 
			posterior_t_1['mean']['t'],
			posterior_t_1['hdi_3%']['t'],
			posterior_t_1['hdi_97%']['t']])
		timing_4.append([chr_arm, 
			posterior_t_4['mean']['t'],
			posterior_t_4['hdi_3%']['t'],
			posterior_t_4['hdi_97%']['t']])

	timing_1 = pd.DataFrame(timing_1, columns = ['chr_arm', 'mean', 'hdi_3', 'hdi_97'])
	timing_4 = pd.DataFrame(timing_4, columns = ['chr_arm', 'mean', 'hdi_3', 'hdi_97'])

	timing_1 = timing_1.set_index('chr_arm')
	timing_1 = timing_1.reindex(index=timing_4['chr_arm'])
	timing_1 = timing_1.reset_index()
	y_ticks_label = timing_4['chr_arm']
	ylabel = "chromosome arm"
	plot_figure_3(timing_1,timing_4,samplename,y_ticks_label,ylabel)

def get_count_amp_per_arm(amplifications_filtered):
	x = []
	for chromosome in list(range(1,23))+['X','Y']:
		for arm in ['p','q']:
			x.append([str(chromosome)+arm,0])

	count = pd.DataFrame(x,columns = ["arm", 'count'])
	for a in amplifications_filtered:
		arm = a['chromosome']+a['arm']
		count['count'].loc[count['arm'] ==arm] = count['count'].loc[count['arm'] ==arm]+1
	return count

def get_figure_3_filtering(oncogene,
	cancer_type,
	amplifications,
	summary_traces_4,
	summary_traces_4_no_filtering,
	summary_traces_1,
	summary_traces_1_no_filtering,
	summary_table,
	list_samples,
	include_samplename = True,
	include_credible = True):

	timing_4 = get_timing_oncogene(oncogene, cancer_type,amplifications,summary_traces_4,summary_table,list_samples)
	timing_1 = get_timing_oncogene(oncogene, cancer_type,amplifications,summary_traces_1,summary_table,list_samples)
	timing_4_no_filtering = get_timing_oncogene(oncogene, cancer_type,amplifications,summary_traces_4_no_filtering,summary_table,list_samples)
	timing_1_no_filtering = get_timing_oncogene(oncogene, cancer_type,amplifications,summary_traces_1_no_filtering,summary_table,list_samples)
	
	timing_1 = timing_1.set_index('samplename')
	timing_1 = timing_1.reindex(index=timing_4['samplename'])
	timing_1 = timing_1.reset_index()

	timing_4_no_filtering = timing_4_no_filtering.set_index('samplename')
	timing_4_no_filtering = timing_4_no_filtering.reindex(index=timing_4['samplename'])
	timing_4_no_filtering = timing_4_no_filtering.reset_index()

	timing_1_no_filtering = timing_1_no_filtering.set_index('samplename')
	timing_1_no_filtering = timing_1_no_filtering.reindex(index=timing_4['samplename'])
	timing_1_no_filtering = timing_1_no_filtering.reset_index()

	title = cancer_type + " " +oncogene
	if include_samplename:
		y_ticks_label = timing_1['samplename'].tolist()
		for i in range(len(y_ticks_label)):
			y_ticks_label[i] = y_ticks_label[i][:8]
	else:
		y_ticks_label = []
	ylabel = 'sample'
	plot_figure_3_filtering(timing_1,timing_4,timing_1_no_filtering,timing_4_no_filtering,title,y_ticks_label,ylabel,include_credible)

def plot_figure_3_filtering(timing_1,timing_4,timing_1_no_filtering,timing_4_no_filtering,title,y_ticks_label,ylabel,include_credible = True):
	markersize = 8
	linewidth = 1
	plt.figure()
	plt.plot(timing_4['mean'][-2:],
		np.arange(timing_1_no_filtering.shape[0]-2,timing_1_no_filtering.shape[0])+0.1,
		marker = '.', 
		linewidth=0,
		color='b', 
		markersize=markersize,
		label='Model 4')
	plt.plot(timing_1['mean'][-2:],
		np.arange(timing_1_no_filtering.shape[0]-2,timing_1_no_filtering.shape[0])-0.1,
		marker = '.', 
		linewidth=0,
		color='r', 
		markersize=markersize,
		label='Model 1')

	plt.plot(timing_4_no_filtering['mean'][-2:],
		np.arange(timing_1_no_filtering.shape[0]-2,timing_1_no_filtering.shape[0])+0.2,
		marker = 'x', 
		linewidth=0,
		color='b', 
		markersize=markersize,
		label='Model 4 (no filtering APOBEC)')
	plt.plot(timing_1_no_filtering['mean'][-2:],
		np.arange(timing_1_no_filtering.shape[0]-2,timing_1_no_filtering.shape[0])-0.2,
		marker = 'x', 
		linewidth=0,
		color='r', 
		markersize=markersize,
		label='Model 1 (no filtering APOBEC)')

	for i in range(timing_1.shape[0]-2,timing_1.shape[0]):
		if i!=timing_1.shape[0]-1:
			plt.axhline(y=i+0.5,color='k',  linewidth=0.5)
		if include_credible:
			plt.plot([timing_4['hdi_3'].iloc[i],timing_4['hdi_97'].iloc[i]],
				[i+0.1,i+0.1],
				color='b',
				linewidth=linewidth)
			plt.plot([timing_1['hdi_3'].iloc[i],timing_1['hdi_97'].iloc[i]],
				[i-0.1,i-0.1],
				color='r',
				linewidth=linewidth)
			plt.plot([timing_4_no_filtering['hdi_3'].iloc[i],timing_4_no_filtering['hdi_97'].iloc[i]],
				[i+0.2,i+0.2],
				color='b',
				linestyle = 'dashed',
				linewidth=linewidth)
			plt.plot([timing_1_no_filtering['hdi_3'].iloc[i],timing_1_no_filtering['hdi_97'].iloc[i]],
				[i-0.2,i-0.2],
				color='r',
				linestyle = 'dashed',
				linewidth=linewidth)
	plt.axvline(x=0.5,color='k',  linewidth=2)
	plt.yticks(np.arange(timing_4.shape[0]),y_ticks_label)
	plt.xticks(np.linspace(0,1,11))

	plt.title(title,fontsize=18)
	plt.xlim(0, 1)
	plt.xlabel('Inferred molecular time')
	plt.ylabel(ylabel)
	plt.text(0.2,timing_4.shape[0]-1.5,"Early",fontsize=16)
	plt.text(0.7,timing_4.shape[0]-1.5,"Late",fontsize=16)
	plt.legend(loc = 'lower left')
	plt.grid(axis='x')


summary_table = pd.read_csv(path_summary_table, sep='\t')

with open(path_amplifications_all, 'r') as fp:
	amplifications = json.load(fp)

if Path(file_statistics_amp).is_file(): 
	summary = pd.read_csv(file_statistics_amp, sep = '\t')
else:
	summary = get_summary_amplifications(summary_table,amplifications, file_statistics_amp,save = True)


path_summary_traces_1 = Path("inferred_times/traces")/"subclonality_modelled"/("Model"+str(1))/"summary.pkl"
path_summary_traces_1_no_filtering = "summary_no_filter_APOBEC_model_1.pkl"
path_summary_traces_4 = Path("inferred_times/traces")/"subclonality_modelled"/("Model"+str(4))/"summary.pkl"
path_summary_traces_4_no_filtering = "summary_no_filter_APOBEC_model_4.pkl"


summary_traces_1 = pickle.load(open(path_summary_traces_1, "rb" ))
summary_traces_1_no_filtering = pickle.load(open(path_summary_traces_1_no_filtering, "rb" ))
summary_traces_4 = pickle.load(open(path_summary_traces_4, "rb" ))
summary_traces_4_no_filtering =  pickle.load(open(path_summary_traces_4_no_filtering, "rb" ))


amplifications_filtered = filer_amplifications(amplifications,cancer_type)

print(len(summary_traces_1_no_filtering))
print(len(summary_traces_4_no_filtering))

list_samples = ['a3914a6c-c622-11e3-bf01-24c6515278c0', 
'7c405ca0-c622-11e3-bf01-24c6515278c0', 
'98d27916-c622-11e3-bf01-24c6515278c0', 
'3b41cb48-c623-11e3-bf01-24c6515278c0', 
'7fba5aac-c622-11e3-bf01-24c6515278c0', 
'4a703d3e-c623-11e3-bf01-24c6515278c0', 
'819b4304-c622-11e3-bf01-24c6515278c0', 
'2572b0bc-c622-11e3-bf01-24c6515278c0', 
'850389d4-c622-11e3-bf01-24c6515278c0', 
'd60f880a-c622-11e3-bf01-24c6515278c0', 
'4fdc8980-c623-11e3-bf01-24c6515278c0', 
'ec5e2990-c622-11e3-bf01-24c6515278c0', 
'4b8943be-c623-11e3-bf01-24c6515278c0', 
'6622f932-c622-11e3-bf01-24c6515278c0', 
'030695f6-c623-11e3-bf01-24c6515278c0', 
'7260f57c-c623-11e3-bf01-24c6515278c0']
get_figure_3_filtering(oncogene,
	cancer_type,
	amplifications,
	summary_traces_4,
	summary_traces_4_no_filtering,
	summary_traces_1,
	summary_traces_1_no_filtering,
	summary_table,
	list_samples,
	include_samplename = True,
	include_credible = True)

plt.show()


