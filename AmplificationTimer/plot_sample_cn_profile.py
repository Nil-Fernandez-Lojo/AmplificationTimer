# author: Nil Fernandez Lojo

from AmplificationTimerObjects import Sample
from utility_functions import load_config
from utility_functions import parse_arguments_plot
from pathlib import Path
import matplotlib.pyplot as plt

args = parse_arguments_plot('sample')
path_config = "../config.json"

# TODO not clean
path_preprocessed_data = Path("../preprocessed_data/samples") / (args.samplename + '.json')
load_genome_into_config = not path_preprocessed_data.exists()
config = load_config(path_config,load_genome_into_config=load_genome_into_config)

f,ax = plt.subplots()
s = Sample(config, args.samplename, save=True)
s.plot_cn_profile(add_snvs=args.add_snvs,
                  total_cn=args.total_cn,
                  chromosome=args.chromosome,
                  start_pos=args.start_pos,
                  max_bases=args.max_bases,
                  plot_threshold_amp=args.plot_threshold_amp,
                  title=args.title,
                  path_save=args.path_save,
                  differentiate_snv_type= args.mark_clock_like_snvs,
                  ax=ax)
plt.show()