from utility_functions import load_config
from AmplificationTimerObjects import Sample, plot_rainfall
import matplotlib.pyplot as plt

#samplename = 'f82d213f-caa7-fd59-e040-11ac0d483e46'
samplename = '4e7cdeda-6dc1-4f17-b853-72a68e5aa7e1'
#samplename = '84e601b7-dfa5-4cd5-9fef-07f03967a0d4'
#samplename = 'fc950c33-faa4-0241-e040-11ac0c486786'

path_config = "../config.json"
config = load_config(path_config, load_genome_into_config=False)
sample = Sample(config, samplename, save=False)
segments = sample.segments
f,ax = plt.subplots()
plot_rainfall(segments, config, chromosome=None, start_pos=0, max_bases=None, ax=ax)
plt.show()