# author: Nil Fernandez Lojo

from AmplificationTimerObjects import Sample
from utility_functions import load_config
from utility_functions import parse_arguments

""" If sample has not been preprocessed, load_genome_into_config should be set to true when loading file"""

args = parse_arguments()
if args.title == 'default':
    title = args.samplename
    if args.chromosome is not None:
        title += ' chromosome ' + args.chromosome
else:
    title = args.title

path_config = "../config.json"

# TODO add oncogenes
#oncogene = [(69201952, 69244466, 'MDM2')]  # position in the chr and name

config = load_config(path_config,load_genome_into_config=False)

s = Sample(config, args.samplename)
s.plot_CN_profile(add_snvs=args.add_snvs,
                  total_cn=args.total_cn,
                  chromosome=args.chromosome,
                  start_pos=args.start_pos,
                  max_bases=args.max_bases,
                  plot_threshold_amp=args.plot_threshold_amp,
                  title=title,
                  path_save=args.path_save,
                  show=args.show,
                  differentiate_snv_type= args.mark_clock_like_snvs)