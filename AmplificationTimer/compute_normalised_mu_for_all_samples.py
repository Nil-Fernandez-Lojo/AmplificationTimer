from AmplificationTimerObjects import Sample
from utility_functions import load_config, compute_normalised_mu
import json
sliding_window_size = 10**7

path_config = "../config.json"
path_out = "../normalised_mu_ws"+str(sliding_window_size)+".json"
config = load_config(path_config, load_genome_into_config=False)
N_samples = config['summary_table'].shape[0]
normalised_mu = dict()
for i in range(N_samples):
    samplename = config['summary_table']['samplename'][i]
    print(i, samplename)
    s = Sample(config, samplename, save=False)
    pos, nm = compute_normalised_mu(s, config, sliding_window_size)
    normalised_mu[samplename] = nm
    with open(path_out, 'w') as buff:
        json.dump(normalised_mu, buff, indent=4)
