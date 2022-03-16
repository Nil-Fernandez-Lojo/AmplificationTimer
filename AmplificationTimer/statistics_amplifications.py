import pandas as pd
import json

from AmplificationTimerObjects import Amplification
from utility_functions import load_config


def get_count_oncogene(amplifications_filtered):
    count_oncogene = {}
    for a in amplifications_filtered:
        for gene in a.potential_driver_oncogenes:
            count_oncogene[gene.name] = count_oncogene.get(gene.name, 0) + 1
    return count_oncogene


def filter_amplifications_by_type(amplifications, config, cancer_type='all', oncogene='all'):
    amplifications_filtered = []
    for a in amplifications:
        samplename = a.clinical_data['samplename']
        if not a.clinical_data['is_preferred']:
            continue
        if cancer_type != 'all':
            if config['summary_table']['histology_pcawg'].loc[
                config['summary_table']['samplename'] == samplename].values != cancer_type:
                continue
        if oncogene != 'all':
            gene_present = False
            for gene in a.potential_driver_oncogenes:
                if gene.name == oncogene:
                    gene_present = True
                    break
            if not gene_present:
                continue
        amplifications_filtered.append(a)
    return amplifications_filtered

def filter_amplifications_dict_by_type(amplifications_dict, config, cancer_type='all', oncogene='all',must_be_preferred=True):
    amplifications_filtered = []
    for a_dict in amplifications_dict:
        samplename = a_dict['clinical_data']['samplename']
        if must_be_preferred and (not a_dict['clinical_data']['is_preferred']):
            continue
        if cancer_type != 'all':
            if config['summary_table']['histology_pcawg'].loc[
                config['summary_table']['samplename'] == samplename].values != cancer_type:
                continue
        if oncogene != 'all':
            gene_present = False
            for gene in a_dict['potential_driver_oncogenes']:
                if gene['name'] == oncogene:
                    gene_present = True
                    break
            if not gene_present:
                continue
        amplifications_filtered.append(a_dict)
    return amplifications_filtered


def count_patients_cancer_type(summary_table, cancer_type):
    count = 0
    for j in range(summary_table.shape[0]):
        if summary_table['histology_pcawg'].loc[j] == cancer_type and (summary_table['is_preferred'].loc[j]):
            count += 1
    return count


def get_summary_amplifications(amplifications, config):
    summary = pd.DataFrame(index=range(len(config['summary_table']['histology_pcawg'].unique())),
                           columns=["cancer_type",
                                    "number_patients",
                                    "number_patients_with_amp",
                                    "most_common_oncogene",
                                    "number_amp_oncogene"])

    for i, cancer_type in enumerate(config['summary_table']['histology_pcawg'].unique()):
        summary["cancer_type"].loc[i] = cancer_type
        summary["number_patients"].loc[i] = count_patients_cancer_type(config['summary_table'], cancer_type)
        amplifications_filtered = filter_amplifications_by_type(amplifications, config, cancer_type)
        patients_with_amp = set([a.clinical_data['samplename'] for a in amplifications_filtered])
        summary["number_patients_with_amp"].loc[i] = len(patients_with_amp)
        count_oncogene = get_count_oncogene(amplifications_filtered)
        if len(count_oncogene) == 0:
            summary["most_common_oncogene"].loc[i] = None
            summary["number_amp_oncogene"].loc[i] = 0
        else:
            summary["most_common_oncogene"].loc[i] = max(count_oncogene, key=count_oncogene.get)
            summary["number_amp_oncogene"].loc[i] = count_oncogene[summary["most_common_oncogene"].loc[i]]

    summary = summary.astype(
        {'number_patients_with_amp': 'int', "number_amp_oncogene": "int", "number_patients": "int"})

    summary["f_patient_amplification"] = summary["number_patients_with_amp"] / summary["number_patients"]
    summary["f_patient_oncogene"] = summary["number_amp_oncogene"] / summary["number_patients"]
    summary["f_amp_oncogene"] = summary["number_amp_oncogene"] / summary["number_patients_with_amp"]
    return summary


def get_count_amp_per_arm(amplifications):
    x = []
    for chromosome in list(range(1, 23)) + ['X', 'Y']:
        for arm in ['p', 'q']:
            x.append([str(chromosome) + arm, 0])

    count = pd.DataFrame(x, columns=["arm", 'count'])
    for a in amplifications:
        arm = str(a.chromosome) + a.arm
        count['count'].loc[count['arm'] == arm] = count['count'].loc[count['arm'] == arm] + 1
    return count


if __name__ == "__main__":
    path_config = "../config.json"
    config = load_config(path_config, load_genome_into_config=False)

    count = 0
    for i in range(config['summary_table'].shape[0]):
        if config['summary_table']['is_preferred'].iloc[i]:
            count += 1
    print("number samples is prefered", count)

    amplifications = []
    for file in config['path_filtered_amplifications_folder'].glob('*.json'):
        with open(file) as json_file:
            amplifications.append(Amplification(amplification_dict=json.load(json_file), config=config))

    summary = get_summary_amplifications(amplifications, config)
    summary.to_csv("../statistics_amplifications_filtered.tsv", index=False, sep='\t')
    print("Number patients with amplifications", summary['number_patients_with_amp'].sum())


    amplifications_wgd_count = 0
    amplifications_minor_allele_2_WGD = 0
    amplifications_high_minor_cn_wgd = 0
    amplifications_high_minor_cn_no_wgd = 0
    amplifications_with_amplified_minor_allele = 0
    amplifications_with_amplified_minor_allele_M_and_m_cp = []
    amplifications_0_oncogenes = 0
    amplifications_1_oncogene = 0
    amplifications_multiple_oncogenes = 0

    for a in amplifications:
        minor_amplified = False
        for segment in a.segments:
            if segment.minor_cn > 1:
                if a.threshold_amplification == 10:
                    amplifications_high_minor_cn_wgd += 1
                else:
                    amplifications_high_minor_cn_no_wgd += 1
                if segment.minor_cn >= a.threshold_amplification:
                    amplifications_with_amplified_minor_allele += 1
                    minor_amplified = True
                break

        if a.threshold_amplification == 10:
            amplifications_wgd_count += 1

            add_to_minor_allele_2_WGD = False
            for segment in a.segments:
                if segment.minor_cn == 2:
                    add_to_minor_allele_2_WGD = True
                elif segment.minor_cn > 2:
                    add_to_minor_allele_2_WGD = False
                    break
            if add_to_minor_allele_2_WGD:
                amplifications_minor_allele_2_WGD += 1

        if minor_amplified:
            major = 0
            minor = 0
            for segment in a.segments:
                if segment.minor_cn > minor:
                    minor = segment.minor_cn
                if segment.major_cn > major:
                    major = segment.major_cn
            amplifications_with_amplified_minor_allele_M_and_m_cp.append((major, minor))
        if len(a.potential_driver_oncogenes) == 0:
            amplifications_0_oncogenes += 1
        elif len(a.potential_driver_oncogenes) == 1:
            amplifications_1_oncogene += 1
        else:
            amplifications_multiple_oncogenes += 1
            print(a.clinical_data['histology_pcawg'], [g.name for g in a.potential_driver_oncogenes])
    print("amplifications_wgd_count:", amplifications_wgd_count)
    print("amplifications_high_minor_cn_wgd:", amplifications_high_minor_cn_wgd)
    print("amplifications_high_minor_cn_no_wgd:", amplifications_high_minor_cn_no_wgd)
    print("amplifications_with_amplified_minor_allele:", amplifications_with_amplified_minor_allele)
    print('amplifications_with_amplified_minor_allele_M_and_m_cp', amplifications_with_amplified_minor_allele_M_and_m_cp)
    print("amplifications_0_oncognes:", amplifications_0_oncogenes)
    print("amplifications_1_oncognes:", amplifications_1_oncogene)
    print("amplifications_multiple_oncognes:", amplifications_multiple_oncogenes)
    print('amplifications_minor_allele_2_WGD i.e. model 3', amplifications_minor_allele_2_WGD)
    n_amp_per_chr_arm = get_count_amp_per_arm(amplifications)
    print('n_amp_per_chr_arm',n_amp_per_chr_arm)




