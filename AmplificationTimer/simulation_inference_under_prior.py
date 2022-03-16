import numpy as np
import pandas as pd
import arviz as avz
from pathlib import Path

from inference import inference_t, generate_simulated_data_under_prior
from utility_functions import parse_arguments_simulation, load_config
from AmplificationTimerObjects import Sample

path_config = "../config.json"

# TODO check that it works

args = parse_arguments_simulation()
# TODO: should write additional file with settings used
if not args.folder_out.exists():
    args.folder_out.mkdir(parents=True, exist_ok=True)

file_out = args.folder_out / (args.samplename + '_'+args.chromosome + args.arm+'.tsv')

# TODO not clean
path_preprocessed_data = Path("../preprocessed_data/samples") / (args.samplename + '.json')
load_genome_into_config = not path_preprocessed_data.exists()
config = load_config(path_config, load_genome_into_config=load_genome_into_config)

sample = Sample(config, args.samplename, save=False)
for a in sample.amplifications:
    if (a.chromosome.chromosome == args.chromosome) and (a.arm == args.arm):
        amplification = a

ploidy_amp = amplification.get_mean_ploidy()

results = pd.DataFrame(np.zeros((len(args.expected_mut_list) * len(args.list_t), 6)),
                       columns=['samplename',
                                'chr_arm',
                                'model_idx',
                                'expected_mutations',
                                't',
                                'inferred_t'])
results.loc['samplename'] = args.samplename
results.loc['chr_arm'] = args.chromosome + args.arm
results.loc['model_idx'] = args.model_idx

i = 0
for expected_mut in args.expected_mut_list:
    for t in args.list_t:
        mu = expected_mut / (amplification.get_length() * ploidy_amp)
        t, u, amplification = generate_simulated_data_under_prior(sample,
                                                                  amplification,
                                                                  args.model_idx,
                                                                  t=t,
                                                                  mu=mu,
                                                                  min_reads_detect_SNV=0,
                                                                  subclonality_modelled_simulation=True,
                                                                  nrpcc='same')
        n_snv = 0
        for s in amplification.segments:
            n_snv += len(s.SNVs)
        print(expected_mut, n_snv)
        continue
        trace = inference_t(amplification,
                            args.model_idx,
                            args.subclonality_modelled,
                            args.filter_APOBEC,
                            only_clock_like_SNVs=False,
                            true_t=t,
                            true_u=u,
                            save=False)

        summary = avz.summary(trace)
        inferred_t = summary[['mean']].T.reset_index(drop=True)['t'].loc[0]
        results.at[i,'expected_mutations'] = expected_mut
        results.at[i,'t'] = t
        results.at[i,'inferred_t'] = inferred_t
        results.to_csv(file_out, index=False, sep='\t')
        i += 1
