import json
from pathlib import Path

from AmplificationTimerObjects import Sample, Amplification
from scipy.stats import binomtest
from utility_functions import load_config

# Filtering parameters
p_value_threshold = 0.05
min_num_non_intermediate_snvs = 50
min_f_non_intermediate_snvs = 0.5
diff_cn_filter = 2
min_f_major_snvs_valid = 0.8

# plot parameters
folder_plot_filtered_amps = Path("../figures/CNP_filtered_amps")
folder_plot_filtered_amps_clocklike = Path("../figures/CNP_filtered_amps_clocklike")

plot_threshold_amp = True
add_snvs = True
total_cn = False
rel_margin = 0.1

path_config = "../config.json"
config = load_config(path_config)


def filter_preferred(amplifications):
    amplifications_preferred = []
    for a in amplifications:
        if a.clinical_data['is_preferred']:
            amplifications_preferred.append(a)
    return amplifications_preferred

def get_number_snvs_non_intermediate(amplification, p_value_threshold, clock_like_filter=False):
    n = 0
    if ((amplification.clinical_data['inferred_sex'] == 'female') or
            amplification.chromosome.chromosome not in ['X', 'Y']):
        ploidy_nc = 2
    else:
        ploidy_nc = 1
    rho = amplification.clinical_data['purity']
    for segment in amplification.segments:
        M = segment.major_cn
        m = segment.minor_cn
        for snv in segment.SNVs:
            if not clock_like_filter or snv.clock_like():
                alt_count = snv.alt_count
                tot_count = snv.get_tot_count()
                if tot_count == 0:
                    continue
                q_clonal_one = rho / ((M + m) * rho + ploidy_nc * (1 - rho))  # 1 copy clonal
                q_clonal_all = q_clonal_one * M  # all copies clonal
                p_binom_all = binomtest(alt_count, tot_count, q_clonal_all, alternative='less').pvalue
                p_binom_1_clonal = binomtest(alt_count, tot_count, q_clonal_one, alternative='greater').pvalue

                if (p_binom_all > p_value_threshold) or (p_binom_1_clonal > p_value_threshold):
                    n += snv.get_multiplicity(amplification.clinical_data['purity'],M+m,ploidy_nc)
    return n

def filter_number_SNVs(amplifications, number_snvs, min_fraction_snvs_valid, p_value_threshold,
                       clock_like_filter=False):
    filtered_amps = []
    for amp in amplifications:
        tot_snvs = 0
        for segment in amp.segments:
            for snv in segment.SNVs:
                if (not clock_like_filter) or snv.clock_like():
                    if snv.get_tot_count() !=0:
                        tot_snvs += snv.get_multiplicity(amp.clinical_data['purity'],segment.get_tot_cn(),segment.get_ploidy_healthy())
        non_intermediate_snvs = get_number_snvs_non_intermediate(amp, p_value_threshold, clock_like_filter)
        if ((non_intermediate_snvs >= number_snvs) and
                (non_intermediate_snvs >= min_fraction_snvs_valid * tot_snvs)):
            filtered_amps.append(amp)
    return filtered_amps


def filter_problematic_major(amplifications, diff_cn, fraction_snvs_valid, p_value_threshold, clock_like_filter=False):
    filtered_amplifications = []
    for amp in amplifications:
        n_snvs_close_to_major_cn = 0
        n_snvs_pass_test_major = 0
        rho = amp.clinical_data['purity']
        if ((amp.clinical_data['inferred_sex'] == 'female') or amp.chromosome.chromosome not in ['X', 'Y']):
            ploidy_nc = 2
        else:
            ploidy_nc = 1
        for segment in amp.segments:
            M = segment.major_cn
            m = segment.minor_cn
            for snv in segment.SNVs:
                if not clock_like_filter or snv.clock_like():
                    alt_count = snv.alt_count
                    ref_count = snv.ref_count
                    tot_count = alt_count + ref_count
                    if tot_count == 0:
                        continue
                    q_clonal_all = M * rho / ((M + m) * rho + ploidy_nc * (1 - rho))
                    snv_cn = snv.get_multiplicity(amp.clinical_data['purity'],M + m,ploidy_nc)
                    if M - diff_cn <= snv_cn <= M + diff_cn:
                        n_snvs_close_to_major_cn += 1
                        if (binomtest(alt_count, tot_count, q_clonal_all, alternative='less').pvalue
                                > p_value_threshold):
                            n_snvs_pass_test_major += 1
        if n_snvs_pass_test_major >= fraction_snvs_valid * n_snvs_close_to_major_cn:
            filtered_amplifications.append(amp)
    return filtered_amplifications


def filter_amplifications(amplifications,
           clock_like_filter,
           min_num_non_intermediate_snvs,
           min_f_non_intermediate_snvs,
           p_value_threshold,
           diff_cn_filter,
           min_f_major_snvs_valid):

    a = filter_preferred(amplifications)
    a = filter_number_SNVs(a, min_num_non_intermediate_snvs,
                           min_f_non_intermediate_snvs,
                           p_value_threshold,
                           clock_like_filter)
    a = filter_problematic_major(a,
                                 diff_cn_filter,
                                 min_f_major_snvs_valid,
                                 p_value_threshold,
                                 clock_like_filter)
    return a