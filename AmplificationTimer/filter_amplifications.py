import json
from pathlib import Path
import pandas as pd
import numpy as np
from AmplificationTimerObjects import Sample
from scipy.stats import binomtest
from utility_functions import load_config, icgc

p_value_threshold = 0.05
min_num_snvs_valid = 5
min_fraction_snvs_valid = 0.5
diff_cn_filter = 2
fraction_snvs_valid = 0.8

plot_threshold_amp = True
add_snvs= True
total_cn = False

path_config = "../config.json"
config = load_config(path_config)

def filter_icgc_preferred(amplifications):
    amplifications_icgc_preferred = []
    for a in amplifications:
        samplename = a['clinical_data']['samplename']
        if a['clinical_data']['is_preferred'] and icgc(samplename):
            amplifications_icgc_preferred.append(a)
    return amplifications_icgc_preferred

def get_number_snvs_non_intermediate(amplification,p_value_threshold):
    n = 0
    if ((amplification['clinical_data']['inferred_sex'] == 'female') or
        amplification['chromosome'] not in ['X','Y']):
        ploidy_nc = 2
    else:
        ploidy_nc = 1
    rho = amplification['clinical_data']['purity']
    for segment in amplification["segments"]:
        M = segment['major_cn']
        m = segment['minor_cn']
        for snv in segment["snvs"]:
            alt_count = snv['alt_count']
            tot_count = snv['alt_count'] + snv['ref_count']
            if tot_count ==0:
                continue
            q_clonal_one = rho / ((M+m) * rho + ploidy_nc * (1 - rho))  # 1 copy clonal
            q_clonal_all = q_clonal_one * M  # all copies clonal
            p_binom_all = binomtest(alt_count, tot_count, q_clonal_all,alternative = 'less').pvalue
            p_binom_1_clonal = binomtest(alt_count, tot_count, q_clonal_one,alternative='greater').pvalue

            if (p_binom_all> p_value_threshold ) or (p_binom_1_clonal > p_value_threshold):
                n += 1
    return n

def filter_number_SNVs(amplifications,number_snvs,min_fraction_snvs_valid,p_value_threshold):
    filtered_amps = []
    for amp in amplifications:
        tot_snvs = 0
        for segment in amp["segments"]:
            tot_snvs+= len(segment['snvs'])
        non_intermediate_snvs = get_number_snvs_non_intermediate(amp,p_value_threshold)
        if ((non_intermediate_snvs>=number_snvs) and
                (non_intermediate_snvs>= min_fraction_snvs_valid*tot_snvs)):
            filtered_amps.append(amp)
    return filtered_amps

def get_hard_assignation_multiplicity_SNV(alt_count,ref_count,tot_cn,rho,ploidy_nc):
    vaf = alt_count/(alt_count+ref_count)
    return vaf*(tot_cn + ((1-rho)/rho)*ploidy_nc)

def filter_problematic_major(amplifications,diff_cn,fraction_snvs_valid,p_value_threshold):
    filtered_amplifications = []
    for amp in amplifications:
        n_snvs_close_to_major_cn = 0
        n_snvs_pass_test_major = 0
        rho = amp['clinical_data']['purity']
        if ((amp['clinical_data']['inferred_sex'] == 'female') or amp['chromosome'] not in ['X', 'Y']):
            ploidy_nc = 2
        else:
            ploidy_nc = 1
        for segment in amp["segments"]:
            M = segment['major_cn']
            m = segment['minor_cn']
            for snv in segment["snvs"]:
                alt_count = snv['alt_count']
                ref_count = snv['ref_count']
                tot_count = alt_count+ref_count
                if tot_count == 0:
                    continue
                q_clonal_all = M*rho / ((M+m) * rho + ploidy_nc * (1 - rho))
                snv_cn = get_hard_assignation_multiplicity_SNV(alt_count, ref_count, M+m, rho, ploidy_nc)
                if M-diff_cn <= snv_cn <= M+diff_cn:
                    n_snvs_close_to_major_cn += 1
                    if (binomtest(alt_count, tot_count, q_clonal_all,alternative='less').pvalue
                            > p_value_threshold):
                        n_snvs_pass_test_major+=1
        if n_snvs_pass_test_major >= fraction_snvs_valid*n_snvs_close_to_major_cn:
            filtered_amplifications.append(amp)
    return filtered_amplifications

def plot_amplification(a):
    name = a['clinical_data']['samplename']
    s = Sample(config,name)
    chromosome = a['segments'][0]['chromosome']
    margin = 0.1 * (a['segments'][-1]['end'] - a['segments'][0]['start'])
    start_pos = a['segments'][0]['start'] - margin
    max_bases = (a['segments'][-1]['end'] - a['segments'][0]['start']) + 2*margin
    s.plot_CN_profile(add_snvs=add_snvs,
        total_cn=total_cn,
        chromosome=chromosome,
        start_pos=start_pos,
        max_bases = max_bases,
        plot_threshold_amp = plot_threshold_amp,
        title = name + ' chromosome '+chromosome)
    #path_save = 'figures/plots_CNP_amps/'+name+'_chr_'+chromosome+'.png')


amplifications_dict = []
for file in config['path_amplifications_folder'].glob('*.json'):
    with open(file) as json_file:
        amplifications_dict.append(json.load(json_file))

amplifications_icgc_preferred = filter_icgc_preferred(amplifications_dict)
amplifications_5_snvs_or_more = filter_number_SNVs(amplifications_icgc_preferred,min_num_snvs_valid,
                                                   min_fraction_snvs_valid,
                                                   p_value_threshold)
filtered_amplifications = filter_problematic_major(amplifications_5_snvs_or_more,
                                                   diff_cn_filter,
                                                   fraction_snvs_valid,
                                                   p_value_threshold)
print('Number amps 5 or more valid snvs:',len(amplifications_5_snvs_or_more))
print('Number amps pass binomial test',len(filtered_amplifications))

for i in range(20):
    plot_amplification(filtered_amplifications[i])


# print(len(amplifications))
# print(len(amplifications_icgc_preferred))