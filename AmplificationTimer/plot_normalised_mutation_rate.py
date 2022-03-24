from utility_functions import load_config
from AmplificationTimerObjects import Position, Chromosome, Sample, plot_normalised_mu
import matplotlib.pyplot as plt

sliding_window_size = 10**7
#samplename = 'f82d213f-caa7-fd59-e040-11ac0d483e46'
#samplename = '25224aa0-cfdd-48ec-92e5-8f3992a3e574'
#samplename = '84e601b7-dfa5-4cd5-9fef-07f03967a0d4'
samplename = 'fc950c33-faa4-0241-e040-11ac0c486786'
add_snvs = True
total_cn = True

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

def compute_normalised_mu(sample):
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



path_config = "../config.json"
config = load_config(path_config, load_genome_into_config=False)
sample = Sample(config, samplename, save=False)
pos,normalised_mu = compute_normalised_mu(sample)

f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10,5))

plot_normalised_mu(pos,normalised_mu,sample.mutation_rate.get_ml(False),config, ax=ax1)
sample.plot_cn_profile(add_snvs=add_snvs,
                  total_cn=total_cn,
                  chromosome=None,
                  start_pos=0,
                  max_bases=None,
                  plot_threshold_amp=True,
                  title=samplename,
                  path_save=None,
                  differentiate_snv_type=True,
                  ax=ax2)
plt.show()


