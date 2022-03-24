from AmplificationTimerObjects import Amplification
from utility_functions import load_config
import json

only_clock_like_SNVs = False
out_file = '../snvs_1_copy_vs_major_cn.json'

path_config = "../config.json"
config = load_config(path_config, load_genome_into_config=False)
out_dict = []
for i,file in enumerate(config['path_filtered_amplifications_folder'].glob('*.json')):
    print(i)
    with open(file) as json_file:
        amplification = Amplification(amplification_dict=json.load(json_file), config=config)
    for segment in amplification.segments:
        out_dict.append({'samplename': amplification.clinical_data['samplename'],
                         'chr': str(amplification.chromosome),
                         'arm': str(amplification.arm),
                         'major_cn': segment.major_cn,
                         'minor_cn': segment.minor_cn,
                         'len_seg': segment.get_length(),
                         'n_snvs_1': segment.get_n_snvs_max_1_copy(only_clock_like_SNVs=only_clock_like_SNVs)})
    with open(out_file, 'w') as buff:
        json.dump(out_dict, buff, indent=4)
