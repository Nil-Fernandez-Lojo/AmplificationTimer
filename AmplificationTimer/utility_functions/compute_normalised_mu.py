from AmplificationTimerObjects import Chromosome,Position

def compute_normalised_mu_one_window(segments,rho):
    n_snv = 0
    n_bases = 0
    for seg in segments:
        tot_cn = seg.get_tot_cn()
        n_bases += tot_cn * seg.get_length()
        for snv in seg.snvs:
            tot_count = snv.get_tot_count()
            if tot_count == 0:
                # TODO I have to investigate this
                continue
            vaf = snv.alt_count/tot_count
            snv_multiplicity = vaf * (tot_cn + ((1 - rho) / rho) * seg.get_ploidy_healthy())
            n_snv+=snv_multiplicity
    if len(segments) == 0:
        return float('nan')
    else:
        return n_snv/n_bases

def filter_segments(segments,chromosome,start_pos,end_pos):
    filtered_segments = []
    for s in segments:
        seg = s.subset(chromosome,start_pos,end_pos)
        if seg is not None:
            filtered_segments.append(seg)
    return filtered_segments

def compute_normalised_mu(sample,config,sliding_window_size):
    pos = []
    normalised_mu = []
    for chr_str in config['list_chromosome']:
        l = config['chromosome_arm_length']['length'][config['chromosome_arm_length']['chromosome'] == chr_str].sum()
        chromosome = Chromosome(chr_str)
        n_windows = l // sliding_window_size
        if l % sliding_window_size != 0:
            n_windows += 1
        for i in range(n_windows):
            start_pos = i*sliding_window_size
            end_pos = min(((i+1)*sliding_window_size),l)
            pos.append(Position(chromosome, (start_pos+end_pos)/2, config))
            segments = filter_segments(sample.segments,chromosome,start_pos,end_pos)
            normalised_mu.append(compute_normalised_mu_one_window(segments, sample.clinical_data['purity']))
    return pos,normalised_mu
