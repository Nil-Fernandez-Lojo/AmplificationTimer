from utility_functions import load_config
from AmplificationTimerObjects import Position, Chromosome, Sample, plot_normalised_mu, compute_normalised_mu
import matplotlib.pyplot as plt

sliding_window_size = 10**7
samplename = 'f82d213f-caa7-fd59-e040-11ac0d483e46'
#samplename = '25224aa0-cfdd-48ec-92e5-8f3992a3e574'
#samplename = '84e601b7-dfa5-4cd5-9fef-07f03967a0d4'
#samplename = '6d9d7ffc-c622-11e3-bf01-24c6515278c0'
add_snvs = True
total_cn = True

path_config = "../config.json"
config = load_config(path_config, load_genome_into_config=False)
sample = Sample(config, samplename, save=False)
pos,normalised_mu = compute_normalised_mu(sample,config,sliding_window_size)

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


