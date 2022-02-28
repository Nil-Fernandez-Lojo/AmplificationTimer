import sys
import numpy as np
import pandas as pd
import arviz as avz

from inference import inference_t, generate_simulated_data_under_prior

sample_idx = int(sys.argv[1])
model_idx = int(sys.argv[2])
folder_to_save =
file_out = "fig_2_data/sample_idx_"+str(sample_idx)+"_model_"+str(model_idx)+".tsv"
expected_mut_list = [10,100,1000,10000]
amplification_idx = 0
list_t = np.linspace(0.0, 1.0, num=21)
subclonality_modelled = True
filter_APOBEC = True
save = False

sample =
amplification =
length_amplification = 0
for segment in amplification.segments:
	length_amplification += segment.get_length()
ploidy_amp = amplification.get_mean_ploidy()

results = pd.DataFrame(np.zeros((len(expected_mut_list)*len(list_t), 6)),
	columns = ['sample_idx',
	'amplification_idx',
	'model_idx',
	'expected_mutations',
	't',
	'inferred_t'])
results['sample_idx'] = sample_idx
results['amplification_idx'] = amplification_idx
results['model_idx'] = model_idx

i = 0
for expected_mut in expected_mut_list:
	for t in list_t:
		mu = expected_mut/(length_amplification*ploidy_amp)
		t,u,amplification = generate_simulated_data_under_prior(sample,
                                     amplification,
                                     model_idx,
                                     t=t,
                                     mu=mu,
                                     min_reads_detect_SNV=0,
                                     subclonality_modelled_simulation=True,
                                     nrpcc='same')
		trace = inference_t(amplification,
					model_idx,
					subclonality_modelled,
					filter_APOBEC,
					only_clock_like_SNVs=False,
					true_t=t,
					true_u=u,
					save=False)

		summary = avz.summary(trace)
		inferred_t = summary[['mean']].T.reset_index(drop=True)['t'].loc[0]
		results['expected_mutations'].loc[i] = expected_mut
		results['t'].loc[i] = t
		results['inferred_t'].loc[i] = inferred_t
		results.to_csv(file_out,index=False,sep='\t')
		i+=1

