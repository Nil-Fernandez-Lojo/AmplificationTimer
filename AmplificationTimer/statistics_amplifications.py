import json
from pathlib import Path
import pandas as pd
import pickle
import arviz as avz
import matplotlib.pyplot as plt
import numpy as np

path_amplifications_all = "preprocessed_data/amplifications_all.json"

model_idx = 1
subclonality_modelled = True
samplename = '2e76891c-b620-4cc0-9315-6f1217b09b1e'
cancer_type = "MALY"
oncogene = 'MALT1'
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

def get_timing_oncogene(oncogene, cancer_type,amplifications,summary_traces,summary_table):
	timing = []
	amplifications_filtered = filer_amplifications(amplifications,cancer_type)
	amplification_id = 0
	for i,a in enumerate(amplifications_filtered):
		if i == 0 or (a['clinical_data']['samplename'] != amplifications_filtered[i-1]['clinical_data']['samplename']):
			amplification_id = 0
		else:
			amplification_id +=1
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
	markersize = 8
	plt.figure()
	plt.plot(timing_4['mean'],
		np.arange(timing_4.shape[0])+0.1,
		marker = '.', 
		linewidth=0,
		color='b', 
		markersize=markersize,
		label='Model 4')
	plt.plot(timing_1['mean'],
		np.arange(timing_1.shape[0])-0.1,
		marker = '.', 
		linewidth=0,
		color='r', 
		markersize=markersize,
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
	title = 'MALT1 malignant lymphoma'
	plt.title(title,fontsize=18)
	plt.xlim(0, 1)
	plt.xlabel('Inferred molecular time',fontsize=14)
	plt.ylabel(ylabel,fontsize=14)
	plt.text(0.2,timing_4.shape[0]-1.5,"Early",fontsize=16)
	plt.text(0.7,timing_4.shape[0]-1.5,"Late",fontsize=16)
	plt.legend(loc = 'upper left')
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



summary_table = pd.read_csv(path_summary_table, sep='\t')

count = 0
for i in range(summary_table.shape[0]):
	if icgc(summary_table['samplename'].iloc[i]) and summary_table['is_preferred'].iloc[i]:
		count+=1
print("number samples icgc.py and is prefered",count)


with open(path_amplifications_all, 'r') as fp:
	amplifications = json.load(fp)

if Path(file_statistics_amp).is_file(): 
	summary = pd.read_csv(file_statistics_amp, sep = '\t')
else:
	summary = get_summary_amplifications(summary_table,amplifications, file_statistics_amp,save = True)

print("number of amplifications with oncogenes", )

if subclonality_modelled:
	path_summary_traces_1 = Path("inferred_times/traces")/"subclonality_modelled"/("Model"+str(1))/"summary.pkl"
	path_summary_traces_4 = Path("inferred_times/traces")/"subclonality_modelled"/("Model"+str(4))/"summary.pkl"
else:
	path_summary_traces_1 = Path("inferred_times/traces")/"subclonality_not_modelled"/("Model"+str(1))/"summary.pkl"
	path_summary_traces_4 = Path("inferred_times/traces")/"subclonality_not_modelled"/("Model"+str(4))/"summary.pkl"

summary_traces_1 = pickle.load(open(path_summary_traces_1, "rb" ))
summary_traces_4 = pickle.load(open(path_summary_traces_4, "rb" ))


print("Number samples with amplifications", len(summary_traces_1))
count_amplifications_per_sample = np.zeros(46)
for x in summary_traces_1.values():
	count_amplifications_per_sample[len(x)] +=1 

plt.figure(1)
plt.bar(np.arange(1,12), count_amplifications_per_sample[1:12], align='center')
plt.xticks(np.arange(1,12))
plt.ylabel('Number of samples')
plt.xlabel('Number of amplifications per sample')
plt.grid()

print("count_amplifications_per_sample",count_amplifications_per_sample)
print(summary)
amplifications_filtered = filer_amplifications(amplifications,cancer_type)


print("Number of amplifications", len(amplifications_filtered))

count_SNVs = np.zeros(len(amplifications_filtered))
for i,a in enumerate(amplifications_filtered):
	count = 0
	for segment in a['segments']:
		count += len(segment['snvs'])
	count_SNVs[i] = count
print("<5 : ",np.sum((count_SNVs<=5)))
print("6-10 : ",np.sum((count_SNVs>5)&(count_SNVs<=10)))
print("11-50 : ",np.sum((count_SNVs>10)&(count_SNVs<=50)))
print("51-100 : ",np.sum((count_SNVs>50)&(count_SNVs<=100)))
print("101-500 : ",np.sum((count_SNVs>100)&(count_SNVs<=500)))
print("501-1000 : ",np.sum((count_SNVs>500)&(count_SNVs<=1000)))
print(">1000 : ",np.sum((count_SNVs>1000)))

list_samples_spanning_centrome_not_filtered = ["0bfd1068-3fd3-a95b-e050-11ac0c4860c3",
"31f3ff14-7d74-447c-a5da-9ad8336c3f3f",
"34a445c2-1eb4-4a9f-8838-cddc2f82aae4",
"37522f18-77b2-4414-8df8-3c2c8048adba",
"4a703d3e-c623-11e3-bf01-24c6515278c0",
" 6ce66be0-c623-11e3-bf01-24c6515278c0",
"6e839eaf-1dbb-43f5-8846-c980e05540c7",
"9988eb07-01f6-4f83-8699-bb63e0525f08",
"a2a67c8a-c622-11e3-bf01-24c6515278c0",
"d3aff5d3-23c0-43ae-9c01-8ddd776b530b",
"da5b9926-c622-11e3-bf01-24c6515278c0"]
amplification_spanning_centromere = set()
count_amplifications_both_arms = 0
samples_amplifications_both_arms = set()
for i,a in enumerate(amplifications_filtered):
	if i>0:
		if a['clinical_data']['samplename'] == amplifications_filtered[i-1]['clinical_data']['samplename']:
			if a['chromosome'] == amplifications_filtered[i-1]['chromosome']:
				count_amplifications_both_arms+=1
				samples_amplifications_both_arms.add(a['clinical_data']['samplename'])
print('number samples with amplifications_both_arms', len(samples_amplifications_both_arms))
print("amplifications_both_arms",count_amplifications_both_arms)

for i,a in enumerate(amplifications_filtered):
	if a['clinical_data']['samplename'] in list_samples_spanning_centrome_not_filtered:
		amplification_spanning_centromere.add(a['clinical_data']['samplename'])
		if a['chromosome'] == 14:
			print(a['clinical_data']['samplename'],a['oncogenes'])
print("amplification_spanning_centromere", amplification_spanning_centromere)

amplifications_wgd_count = 0
amplifications_minor_allele_2_WGD = 0
amplifications_high_minor_cn_wgd = 0
amplifications_high_minor_cn_no_wgd = 0
amplifications_with_amplified_minor_allele = 0
amplifications_with_amplified_minor_allele_M_and_m_cp = []
amplifications_0_oncogenes = 0
amplifications_1_oncogene = 0
amplifications_multiple_oncogenes = 0
number_amplifications = 0

for a in amplifications_filtered:
	if not a['clinical_data']['is_preferred']: continue
	number_amplifications +=1
	if a["threshold_amplification"] == 10:
		amplifications_wgd_count +=1
	minor_amplified = False
	for segment in a['segments']:
		if segment['minor_cn'] > 1:
			if a["threshold_amplification"] == 10:
				amplifications_high_minor_cn_wgd +=1
			else:
				amplifications_high_minor_cn_no_wgd +=1
			if segment['minor_cn'] >= a["threshold_amplification"]:
				amplifications_with_amplified_minor_allele +=1
				minor_amplified = True
			break

	if a["threshold_amplification"] == 10:
		add_to_minor_allele_2_WGD = False
		for segment in a['segments']:
			if segment['minor_cn'] == 2:
				add_to_minor_allele_2_WGD = True
			elif segment['minor_cn'] > 2:
				add_to_minor_allele_2_WGD = False
				break
		if add_to_minor_allele_2_WGD:
			amplifications_minor_allele_2_WGD +=1


	if minor_amplified:
		major = 0
		minor = 0
		for segment in a['segments']:
			if segment['minor_cn'] > minor:
				minor = segment['minor_cn']
			if segment['major_cn'] > major:
				major = segment['major_cn'] 
		amplifications_with_amplified_minor_allele_M_and_m_cp.append((major,minor))
	if len(a["oncogenes"]) == 0:
		amplifications_0_oncogenes +=1
	elif len(a["oncogenes"]) == 1:
		amplifications_1_oncogene +=1
	else:
		amplifications_multiple_oncogenes +=1
print("number amplifications:", number_amplifications)
print("amplifications_wgd_count:", amplifications_wgd_count)
print("amplifications_high_minor_cn_wgd:", amplifications_high_minor_cn_wgd)
print("amplifications_high_minor_cn_no_wgd:", amplifications_high_minor_cn_no_wgd)
print("amplifications_with_amplified_minor_allele:", amplifications_with_amplified_minor_allele)
print('amplifications_with_amplified_minor_allele_M_and_m_cp', amplifications_with_amplified_minor_allele_M_and_m_cp)
print("amplifications_0_oncognes:", amplifications_0_oncogenes)
print("amplifications_1_oncognes:", amplifications_1_oncogene)
print("amplifications_multiple_oncognes:", amplifications_multiple_oncogenes)
print('amplifications_minor_allele_2_WGD i.e. model 3', amplifications_minor_allele_2_WGD)
n_amp_per_chr_arm = get_count_amp_per_arm(amplifications_filtered)
print("n_amp_per_chr_arm",n_amp_per_chr_arm)
count_oncogene = get_count_oncogene(amplifications_filtered)
print(count_oncogene)
print('a3914a6c-c622-11e3-bf01-24c6515278c0') 
print(summary_traces_1['a3914a6c-c622-11e3-bf01-24c6515278c0'])
print(summary_traces_4['a3914a6c-c622-11e3-bf01-24c6515278c0'])
print('7fba5aac-c622-11e3-bf01-24c6515278c0') 
print(summary_traces_1['7fba5aac-c622-11e3-bf01-24c6515278c0'])
print(summary_traces_4['7fba5aac-c622-11e3-bf01-24c6515278c0'])



#plot_timing_one_sample(samplename,summary_traces_1,summary_traces_4,amplifications)


get_figure_3(oncogene,
	cancer_type,
	amplifications,
	summary_traces_4,
	summary_traces_1,
	summary_table,
	include_samplename = False)

# if oncogene == 'all':
# 	for oncogene in count_oncogene.keys():
# 		if count_oncogene[oncogene] >=3:
# 			get_figure_3(oncogene,
# 				cancer_type,
# 				amplifications,
# 				summary_traces_4,
# 				summary_traces_1,
# 				summary_table,
# 				include_samplename = True)
# else:
# 	for gene in oncogene:
# 		if count_oncogene[gene] >=3:
# 			get_figure_3(gene,
# 				cancer_type,
# 				amplifications,
# 				summary_traces_4,
# 				summary_traces_1,
# 				summary_table,
# 				include_samplename = True)
plt.show()


