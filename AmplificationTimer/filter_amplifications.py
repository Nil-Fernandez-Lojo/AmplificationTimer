def get_number_snvs_non_intermediate(amplification, clock_like_filter=False):
    n = 0
    if ((amplification.clinical_data['inferred_sex'] == 'female') or
            amplification.chromosome.chromosome not in ['X', 'Y']):
        ploidy_nc = 2
    else:
        ploidy_nc = 1
    rho = amplification.clinical_data['purity']
    for segment in amplification.segments:
        if segment.problematic_major:
            continue
        M = segment.major_cn
        m = segment.minor_cn
        ploidy_healthy = segment.get_ploidy_healthy()
        for snv in segment.snvs:
            if snv.in_kataegis:
                continue
            if not clock_like_filter or snv.clock_like():
                if snv.get_tot_count() == 0:
                    continue
                if (snv.multiplicity_accepted(rho, m + M, ploidy_healthy, 1, alternative='greater') or
                        snv.multiplicity_accepted(rho, m + M, ploidy_healthy, M, alternative='two-sided')):
                    n += snv.get_multiplicity(amplification.clinical_data['purity'], M + m, ploidy_nc)
    return n


def filter_preferred(amplification):
    return amplification.clinical_data['is_preferred']


def filter_number_SNVs(amplification, number_snvs, min_fraction_snvs_valid, clock_like_filter=False):
    tot_snvs = 0
    for segment in amplification.segments:
        if segment.problematic_major:
            continue
        for snv in segment.snvs:
            if snv.in_kataegis:
                continue
            if (not clock_like_filter) or snv.clock_like():
                if snv.get_tot_count() != 0:
                    tot_snvs += snv.get_multiplicity(amplification.clinical_data['purity'], segment.get_tot_cn(),
                                                     segment.get_ploidy_healthy())
    non_intermediate_snvs = get_number_snvs_non_intermediate(amplification, clock_like_filter)
    return ((non_intermediate_snvs >= number_snvs) and
            (non_intermediate_snvs >= min_fraction_snvs_valid * tot_snvs))


def filter_problematic_major(amplification):
    for segment in amplification.segments:
        if not segment.problematic_major:
            return True
    return False

def filter_amplification(amplification,
                         clock_like_filter,
                         min_num_non_intermediate_snvs,
                         min_f_non_intermediate_snvs):
    selected = False
    if not filter_preferred(amplification):
        reason = 'not_preferred'
    elif not filter_problematic_major(amplification):
        reason = 'problematic_major'
    elif not filter_number_SNVs(amplification, min_num_non_intermediate_snvs,
                                min_f_non_intermediate_snvs,
                                clock_like_filter):
        reason = 'too_many_intermediate_snvs'
    else:
        selected = True
        reason = 'selected'
    return selected, reason
