import matplotlib.pyplot as plt

from AmplificationTimerObjects import Sample
from utility_functions import load_config
from filter_amplifications import filter_amplifications
import json

"""
This script preprocesses the data
"""

path_config = "../config.json"

config = load_config(path_config)

amplifications = []
N_samples = config['summary_table'].shape[0]

for i in range(N_samples):
    samplename = config['summary_table']['samplename'][i]
    print(i, samplename)
    s = Sample(config, samplename, save=True)
    filtered_amplifications = filter_amplifications(s.amplifications,
                                                    False,
                                                    config['min_num_non_intermediate_snvs'],
                                                    config['min_f_non_intermediate_snvs'],
                                                    config['p_value_threshold'],
                                                    config['diff_cn_filter'],
                                                    config['min_f_major_snvs_valid'])
    filtered_amplifications_clock_like = filter_amplifications(s.amplifications,
                                                    True,
                                                    config['min_num_non_intermediate_snvs'],
                                                    config['min_f_non_intermediate_snvs'],
                                                    config['p_value_threshold'],
                                                    config['diff_cn_filter'],
                                                    config['min_f_major_snvs_valid'])

    amp_name_filtered = [samplename + "_" + str(a.chromosome) + a.arm for a in filtered_amplifications]
    amp_name_filtered_clocklike = [samplename + "_" + str(a.chromosome) + a.arm for a in filtered_amplifications_clock_like]

    for amplification in s.amplifications:
        amplification_dict = amplification.to_dict()
        file_name = samplename + "_" + str(amplification.chromosome) + amplification.arm
        path_amplification = config["path_amplifications_folder"] / (file_name+".json")
        with open(path_amplification, 'w') as buff:
            json.dump(amplification_dict, buff, indent=4)

        paths_save_figure = [config['path_folder_plot_amps'] / (file_name + '.png')]

        if file_name in amp_name_filtered:
            path_amplification_filtered = config["path_filtered_amplifications_folder"] / (file_name+".json")
            with open(path_amplification_filtered, 'w') as buff:
                json.dump(amplification_dict, buff, indent=4)
            paths_save_figure.append(config['path_folder_plot_filtered_amps'] / (file_name + '.png'))

        if file_name in amp_name_filtered_clocklike:
            path_amplification_filtered_clocklike = config["path_filtered_amplifications_clocklike_folder"] / (file_name+".json")
            with open(path_amplification_filtered_clocklike, 'w') as buff:
                json.dump(amplification_dict, buff, indent=4)
            paths_save_figure.append(config['path_folder_plot_filtered_amps_clocklike'] / (file_name + '.png'))

        margin = config["rel_margin_plot_amplification"] * (amplification.segments[-1].end.position - amplification.segments[0].start.position)
        start_pos = amplification.segments[0].start.position - margin
        max_bases = (amplification.segments[-1].end.position - amplification.segments[
            0].start.position) + 2 * margin
        for path_save in paths_save_figure:
            fig, ax = plt.subplots()
            s.plot_cn_profile(add_snvs=True,
                              total_cn=False,
                              chromosome=amplification.chromosome,
                              start_pos=start_pos,
                              max_bases=max_bases,
                              plot_threshold_amp=True,
                              differentiate_snv_type=True,
                              title=(samplename + ' ' + amplification.chromosome.chromosome + amplification.arm),
                              path_save=path_save,
                              ax=ax)
            plt.close(fig)
