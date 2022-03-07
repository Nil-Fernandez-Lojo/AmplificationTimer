# author: Nil Fernandez Lojo
import json
from AmplificationTimerObjects import Amplification,Sample
from utility_functions import parse_arguments_plot, load_config

# TODO: the amplification must have been preprocessed

path_config = "../config.json"

args = parse_arguments_plot('amplification')

config = load_config(path_config,load_genome_into_config=False)

path_amp = config['path_amplifications_folder'] / (args.samplename + '_' + args.chromosome +args.arm  + '.json')
with open(path_amp) as json_file:
    amplification_dict = json.load(json_file)

amplification = Amplification(amplification_dict=amplification_dict, config=config)

if args.add_context:
    s = Sample(config, args.samplename)
    margin = args.rel_margin * (amplification.segments[-1].end.position - amplification.segments[0].start.position)
    start_pos = amplification.segments[0].start.position - margin
    max_bases = (amplification.segments[-1].end.position - amplification.segments[0].start.position) + 2 * margin

    s.plot_CN_profile(add_snvs=args.add_snvs,
                      total_cn=args.total_cn,
                      chromosome=amplification.chromosome,
                      start_pos=start_pos,
                      max_bases=max_bases,
                      plot_threshold_amp=True,
                      differentiate_snv_type=args.mark_clock_like_snvs,
                      title=args.title)
else:
    amplification.plot(add_snvs=args.add_snvs,
                       rel_margin=args.rel_margin,
                       total_cn=args.total_cn,
                       plot_threshold_amp=True,
                       title=args.title,
                       differentiate_snv_type=args.mark_clock_like_snvs)



