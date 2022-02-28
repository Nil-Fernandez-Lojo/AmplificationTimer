# author: Nil Fernandez Lojo
import json
from AmplificationTimerObjects import Amplification,Sample
from utility_functions import load_config

""" the amplification must have been preprocessed"""

path_config = "../config.json"

samplename = '00c27940-c623-11e3-bf01-24c6515278c0'
amplification_arm = '9p'
add_context = False
add_snvs = True
total_cn = False
plot_threshold_amp = True
differentiate_snv_type = True
rel_margin = 0.1
# oncogene = [(69201952, 69244466, 'MDM2')]  # position in the chr and name

config = load_config(path_config,load_genome_into_config=False)

# title = name + ' chromosome '+chromosome,
# path_save = 'figures/plots_CNP_amps/'+name+'_chr_'+chromosome+'.png')
with open(config['path_amplifications_folder'] / (samplename + '_' + amplification_arm + '.json')) as json_file:
    amplification_dict = json.load(json_file)

amplification = Amplification(amplification_dict=amplification_dict, config=config)

if add_context:
    s = Sample(config, samplename)
    margin = rel_margin * (amplification.segments[-1].end.position - amplification.segments[0].start.position)
    start_pos = amplification.segments[0].start.position - margin
    max_bases = (amplification.segments[-1].end.position - amplification.segments[0].start.position) + 2 * margin

    s.plot_CN_profile(add_snvs=add_snvs,
                      total_cn=total_cn,
                      chromosome=amplification.chromosome,
                      start_pos=start_pos,
                      max_bases=max_bases,
                      plot_threshold_amp=True,
                      differentiate_snv_type=differentiate_snv_type)
else:
    amplification.plot(add_snvs=add_snvs,
                       rel_margin=rel_margin,
                       total_cn=total_cn,
                       plot_threshold_amp=True,
                       title=None,
                       differentiate_snv_type=differentiate_snv_type)



