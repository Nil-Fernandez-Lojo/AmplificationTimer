import json
from pathlib import Path
import pandas as pd
from pyensembl import EnsemblRelease
import numpy as np
from sample import Sample
from scipy.stats import binomtest

p_value_threshold = 0.05
min_num_snvs_valid = 5
min_fraction_snvs_valid = 0.5
diff_cn_filter = 2
fraction_snvs_valid = 0.8

plot_threshold_amp = True
add_snvs= True
total_cn = False

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

path_amplifications_all = "preprocessed_data/amplifications_all.json"

def icgc(samplename):
    file_white = folder_snv_white/(samplename+suffix_snvs)
    file_gray = folder_snv_gray/(samplename+suffix_snvs)
    return file_gray.exists() or file_white.exists()

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


with open(path_amplifications_all, 'r') as fp:
    amplifications = json.load(fp)

amplifications_icgc_preferred = filter_icgc_preferred(amplifications)
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