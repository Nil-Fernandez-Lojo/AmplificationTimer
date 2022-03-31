import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import pickle
from pathlib import Path

from utility_functions import load_config
from statistics_amplifications import filter_amplifications_dict_by_type


max_f_problematic_major = 1/2
must_be_preferred = True
subclonality_modelled_a = True
subclonality_modelled_b = True
mu_ml_a = False
mu_ml_b = False
model_a = 4
model_b = 1
cancer_type = 'Liver-HCC' # put 'all' if no filter
oncogene = 'YWHAZ' # put 'all' if no filter
include_samplename = True
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


def get_figure_3(summary_traces_a, summary_traces_b,cancer_type,oncogene,label_a,label_b,include_samplename=False):
    timing_b = format_timing(summary_traces_b)
    timing_a = format_timing(summary_traces_a)
    if len(timing_b) != 0:
        timing_a = timing_a.set_index('samplename')
        timing_a = timing_a.reindex(index=timing_b['samplename'])
        timing_a = timing_a.reset_index()
    title = cancer_type + " " + oncogene
    if include_samplename:
        y_ticks_label = timing_a['samplename']
    else:
        y_ticks_label = []
    ylabel = 'sample'
    plot_figure_3(timing_a, timing_b, title, y_ticks_label, ylabel,label_a,label_b)


def plot_figure_3(timing_a, timing_b, title, y_ticks_label, ylabel,label_a,label_b):
    markersize = 8
    plt.figure()
    plt.plot(timing_a['mean'],
             np.arange(timing_a.shape[0]) - 0.1,
             marker='.',
             linewidth=0,
             color='r',
             markersize=markersize,
             label=label_a)
    if len(timing_b)!= 0:
        plt.plot(timing_b['mean'],
                 np.arange(timing_b.shape[0]) + 0.1,
                 marker='.',
                 linewidth=0,
                 color='b',
                 markersize=markersize,
                 label=label_b)

    for i in range(timing_a.shape[0]):
        if len(timing_b) != 0:
            plt.plot([timing_b['hdi_3'].iloc[i], timing_b['hdi_97'].iloc[i]],
                     [i + 0.1, i + 0.1],
                     color='b')
        plt.plot([timing_a['hdi_3'].iloc[i], timing_a['hdi_97'].iloc[i]],
                 [i - 0.1, i - 0.1],
                 color='r')
    plt.axvline(x=0.5, color='k', linewidth=2)
    plt.yticks(np.arange(timing_a.shape[0]), y_ticks_label)
    plt.xticks(np.linspace(0, 1, 11))
    plt.title(title, fontsize=18)
    plt.xlim(0, 1)
    plt.xlabel('Inferred molecular time', fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.text(0.2, timing_a.shape[0] - 1.5, "Early", fontsize=16)
    plt.text(0.7, timing_a.shape[0] - 1.5, "Late", fontsize=16)
    plt.legend(loc='upper left')
    plt.grid(axis='x')


def get_amplifications_sample(amplifications, samplename):
    amplifications_sample = []
    for a in amplifications:
        if a['clinical_data']['samplename'] == samplename:
            amplifications_sample.append(a)
    return amplifications_sample


def plot_timing_one_sample(samplename, summary_traces_a, summary_traces_b, amplifications_filtered,label_a,label_b):
    t_a = summary_traces_a[samplename]
    t_b = summary_traces_b[samplename]
    amplifications_sample = get_amplifications_sample(amplifications_filtered, samplename)

    timing_a = []
    timing_b = []
    for i, a in enumerate(amplifications_sample):
        chr_arm = a['chromosome'] + a['arm']
        posterior_t_1 = t_a[i]
        posterior_t_4 = t_b[i]
        timing_a.append([chr_arm,
                         posterior_t_1['mean']['t'],
                         posterior_t_1['hdi_3%']['t'],
                         posterior_t_1['hdi_97%']['t']])
        timing_b.append([chr_arm,
                         posterior_t_4['mean']['t'],
                         posterior_t_4['hdi_3%']['t'],
                         posterior_t_4['hdi_97%']['t']])

    timing_1 = pd.DataFrame(timing_a, columns=['chr_arm', 'mean', 'hdi_3', 'hdi_97'])
    timing_4 = pd.DataFrame(timing_b, columns=['chr_arm', 'mean', 'hdi_3', 'hdi_97'])

    timing_1 = timing_1.set_index('chr_arm')
    timing_1 = timing_1.reindex(index=timing_4['chr_arm'])
    timing_1 = timing_1.reset_index()
    y_ticks_label = timing_4['chr_arm']
    ylabel = "chromosome arm"
    plot_figure_3(timing_1, timing_4, samplename, y_ticks_label, ylabel,label_a,label_b)


def get_label_model(model,subclonality_modelled,mu_ml):
    label = 'Model ' + str(model)
    if model == 1:
        if mu_ml:
            label += ' mu ml'
        else:
            label += ' mu beta prior'
    if not subclonality_modelled:
        label += ' subclonality not modelled'

    return label


def get_path(model,subclonality_modelled,mu_ml):
    path = Path("../inference_results/inferred_times/summary/filter_APOBEC")
    if subclonality_modelled:
        path = path/"subclonality_modelled"
    else:
        path = path / "subclonality_not_modelled"
    path = path/"all_SNVs_types"
    path = path/("Model"+str(model))
    if model == 1:
        if mu_ml:
            path = path / "mu_ml"
        else:
            path = path / "mu_beta_prior"
    return path


def plot(cancer_type, amplifications_dict,path_model_a,path_model_b=None):
    amplifications_dict_filtered = filter_amplifications_dict_by_type(amplifications_dict, config, cancer_type,
                                                                      oncogene, must_be_preferred=must_be_preferred)
    print(len(amplifications_dict_filtered))
    amp_names_filtered = set()
    for a in amplifications_dict_filtered:
        amp_name = a['clinical_data']['samplename'] + '_' + str(a['chromosome']) + str(a['arm'])
        amp_names_filtered.add(amp_name)
    print(amp_names_filtered)

    summary_traces_a = dict()
    for file in path_model_a.iterdir():
        amp_name = str(file).split('/')[-1].replace('_trace_summary.pkl', '')
        if amp_name in amp_names_filtered:
            summary_traces_a[amp_name] = pickle.load(open(file, "rb"))

    summary_traces_b = dict()
    if path_model_b is not None:
        for file in path_model_b.iterdir():
            amp_name = str(file).split('/')[-1].replace('_trace_summary.pkl', '')
            if amp_name in amp_names_filtered:
                summary_traces_b[amp_name] = pickle.load(open(file, "rb"))

    label_a = get_label_model(model_a, subclonality_modelled_a, mu_ml_a)
    if path_model_b is None:
        label_b = None
    else:
        label_b = get_label_model(model_b, subclonality_modelled_b, mu_ml_b)

    get_figure_3(summary_traces_a, summary_traces_b, cancer_type, oncogene, label_a, label_b,
                 include_samplename=include_samplename)

path_config = "../config.json"
config = load_config(path_config, load_genome_into_config=False)

path_model_a = get_path(model_a,subclonality_modelled_a,mu_ml_a)
path_model_b = get_path(model_b,subclonality_modelled_b,mu_ml_b)
amplifications_dict = []
for file in config['path_filtered_amplifications_folder'].glob('*.json'):
    with open(file) as json_file:
        data = json.load(json_file)
        n_bases_in_problematic_major = 0
        len_amp = 0
        for seg in data['segments']:
            len_amp += seg['end'] - seg['start']
            if seg["problematic_major"]:
                n_bases_in_problematic_major += seg['end'] - seg['start']

        if n_bases_in_problematic_major <= max_f_problematic_major * len_amp:
            amplifications_dict.append(data)

plot(cancer_type, amplifications_dict,path_model_a,path_model_b)
# list_cancer_types = config['summary_table']['histology_pcawg'].unique()
# for cancer_type in list_cancer_types:
#     plot(cancer_type, amplifications_dict,path_model_a,None)

plt.show()