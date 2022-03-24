import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import pickle
from pathlib import Path

from utility_functions import load_config
from statistics_amplifications import filter_amplifications_dict_by_type
from AmplificationTimerObjects import Amplification

must_be_preferred = False
cancer_type = 'Liver-HCC' # put 'all' if no filter
oncogene = 'YWHAZ' # put 'all' if no filter
path_model1 = Path("../inference_results/inferred_times/traces/filter_APOBEC/subclonality_modelled/all_SNVs_types/Model1/summary")
path_model4 = Path("../inference_results/inferred_times/traces/filter_APOBEC/subclonality_modelled/all_SNVs_types/Model4/summary")

def format_timing(summary_traces):
    timing = pd.DataFrame({'samplename':list(summary_traces.keys())}, columns=['samplename', 'mean', 'hdi_3', 'hdi_97'])
    i = 0
    for amp_name, inferred_results in summary_traces.items():
        timing.loc[i] = [amp_name,
                         inferred_results['mean']['t'],
                         inferred_results['hdi_3%']['t'],
                         inferred_results['hdi_97%']['t']]
        i+=1

    timing = timing.sort_values(by=['mean'])
    return timing


def get_figure_3(summary_traces_4, summary_traces_1,cancer_type,oncogene,include_samplename=False):
    timing_4 = format_timing(summary_traces_4)
    timing_1 = format_timing(summary_traces_1)

    timing_1 = timing_1.set_index('samplename')
    timing_1 = timing_1.reindex(index=timing_4['samplename'])
    timing_1 = timing_1.reset_index()
    title = cancer_type + " " + oncogene
    if include_samplename:
        y_ticks_label = timing_1['samplename']
    else:
        y_ticks_label = []
    ylabel = 'sample'
    plot_figure_3(timing_1, timing_4, title, y_ticks_label, ylabel)


def plot_figure_3(timing_1, timing_4, title, y_ticks_label, ylabel):
    markersize = 8
    plt.figure()
    plt.plot(timing_4['mean'],
             np.arange(timing_4.shape[0]) + 0.1,
             marker='.',
             linewidth=0,
             color='b',
             markersize=markersize,
             label='Model 4')
    plt.plot(timing_1['mean'],
             np.arange(timing_1.shape[0]) - 0.1,
             marker='.',
             linewidth=0,
             color='r',
             markersize=markersize,
             label='Model 1')

    for i in range(timing_1.shape[0]):
        plt.plot([timing_4['hdi_3'].iloc[i], timing_4['hdi_97'].iloc[i]],
                 [i + 0.1, i + 0.1],
                 color='b')
        plt.plot([timing_1['hdi_3'].iloc[i], timing_1['hdi_97'].iloc[i]],
                 [i - 0.1, i - 0.1],
                 color='r')
    plt.axvline(x=0.5, color='k', linewidth=2)
    plt.yticks(np.arange(timing_4.shape[0]), y_ticks_label)
    plt.xticks(np.linspace(0, 1, 11))
    plt.title(title, fontsize=18)
    plt.xlim(0, 1)
    plt.xlabel('Inferred molecular time', fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.text(0.2, timing_4.shape[0] - 1.5, "Early", fontsize=16)
    plt.text(0.7, timing_4.shape[0] - 1.5, "Late", fontsize=16)
    plt.legend(loc='upper left')
    plt.grid(axis='x')


def get_amplifications_sample(amplifications, samplename):
    amplifications_sample = []
    for a in amplifications:
        if a['clinical_data']['samplename'] == samplename:
            amplifications_sample.append(a)
    return amplifications_sample


def plot_timing_one_sample(samplename, summary_traces_1, summary_traces_4, amplifications_filtered):
    t_1 = summary_traces_1[samplename]
    t_4 = summary_traces_4[samplename]
    amplifications_sample = get_amplifications_sample(amplifications_filtered, samplename)

    timing_1 = []
    timing_4 = []
    for i, a in enumerate(amplifications_sample):
        chr_arm = a['chromosome'] + a['arm']
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

    timing_1 = pd.DataFrame(timing_1, columns=['chr_arm', 'mean', 'hdi_3', 'hdi_97'])
    timing_4 = pd.DataFrame(timing_4, columns=['chr_arm', 'mean', 'hdi_3', 'hdi_97'])

    timing_1 = timing_1.set_index('chr_arm')
    timing_1 = timing_1.reindex(index=timing_4['chr_arm'])
    timing_1 = timing_1.reset_index()
    y_ticks_label = timing_4['chr_arm']
    ylabel = "chromosome arm"
    plot_figure_3(timing_1, timing_4, samplename, y_ticks_label, ylabel)

path_config = "../config.json"
config = load_config(path_config, load_genome_into_config=False)

amplifications_dict = []
for file in config['path_filtered_amplifications_folder'].glob('*.json'):
    with open(file) as json_file:
        amplifications_dict.append(json.load(json_file))

amplifications_dict_filtered = filter_amplifications_dict_by_type(amplifications_dict, config, cancer_type, oncogene,must_be_preferred=must_be_preferred)
print(len(amplifications_dict_filtered))
amplifications_filtered = dict()
for a in amplifications_dict_filtered:
    amp = Amplification(amplification_dict=a, config=config)
    amplifications_filtered[amp.get_name()] = amp
amp_names_filtered = set(amplifications_filtered.keys()) #TODO add _
print(amp_names_filtered)

summary_traces_1 = dict()
for file in path_model1.iterdir():
    amp_name = str(file).split('/')[-1].replace('_trace_summary.pkl','')
    if amp_name in amp_names_filtered:
        summary_traces_1[amp_name] = pickle.load(open(file, "rb"))

summary_traces_4 = dict()
for file in path_model4.iterdir():
    amp_name = str(file).split('/')[-1].replace('_trace_summary.pkl','')
    if amp_name in amp_names_filtered:
        summary_traces_4[amp_name] = pickle.load(open(file, "rb"))

get_figure_3(summary_traces_4, summary_traces_1,cancer_type,oncogene,include_samplename=True)
plt.show()