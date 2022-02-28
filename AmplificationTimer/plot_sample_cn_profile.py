# author: Nil Fernandez Lojo
from AmplificationTimerObjects import Sample
from utility_functions import load_config

""" If sample has not been preprocessed, load_genome_into_config should be set to true when loading file"""

path_config = "../config.json"

name = '00b9d0e6-69dc-4345-bffd-ce32880c8eef'
add_snvs = False
total_cn = False
chromosome = None
start_pos = 0
max_bases = None
plot_threshold_amp = True
path_save = None
#oncogene = [(69201952, 69244466, 'MDM2')]  # position in the chr and name

config = load_config(path_config,load_genome_into_config=False)

s = Sample(config, name)
s.plot_CN_profile(add_snvs=add_snvs,
                  total_cn=total_cn,
                  chromosome=chromosome,
                  start_pos=start_pos,
                  max_bases=max_bases,
                  plot_threshold_amp=plot_threshold_amp,
                  title=name,
                  path_save=path_save)

# amplifications_filtered = filer_amplifications(amplifications)

# for i,a in enumerate(amplifications_filtered):
# 	name = a['clinical_data']['samplename']
# 	s = Sample(config,name)
# 	chromosome = a['segments'][0]['chromosome']
# 	margin = 0.1 * (a['segments'][-1]['end'] - a['segments'][0]['start'])
# 	start_pos = a['segments'][0]['start'] - margin
# 	max_bases = (a['segments'][-1]['end'] - a['segments'][0]['start']) + 2*margin
# 	print(i,'/',len(amplifications_filtered), name, chromosome)
# 	s.plot_cn(add_snvs=add_snvs,
# 		total_cn=total_cn,
# 		chromosome=chromosome,
# 		start_pos=start_pos,
# 		max_bases = max_bases,
# 		plot_threshold_amp = plot_threshold_amp,
# 		title = name + ' chromosome '+chromosome,
# 		path_save = 'figures/plots_CNP_amps/'+name+'_chr_'+chromosome+'.png')
