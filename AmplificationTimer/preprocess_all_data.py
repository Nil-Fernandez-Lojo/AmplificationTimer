import matplotlib.pyplot as plt

from AmplificationTimerObjects import Sample
from utility_functions import load_config
from filter_amplifications import filter_amplification
import json

"""
This script preprocesses the data
"""
# TODO: add arguments command line e.g. start sample, end sample, plot figures,...

path_config = "../config.json"
config = load_config(path_config)
N_samples = config['summary_table'].shape[0]

reasons = ['not_preferred','too_many_intermediate_snvs','problematic_major']
for r in reasons:
    for p in [config['path_folder_plot_non_selected_amps'],config['path_folder_plot_non_selected_amps_clocklike']]:
        (p/r).mkdir(parents=True, exist_ok=True)

for i in range(N_samples):
    samplename = config['summary_table']['samplename'][i]
    print(i, samplename)
    s = Sample(config, samplename, save=True)
    for amplification in s.amplifications:
        selected, reason = filter_amplification(amplification,
                                                    False,
                                                    config['min_num_non_intermediate_snvs'],
                                                    config['min_f_non_intermediate_snvs'])
        selected_clock_like, reason_clock_like = filter_amplification(amplification,
                                                          True,
                                                          config['min_num_non_intermediate_snvs'],
                                                          config['min_f_non_intermediate_snvs'])
        amplification_dict = amplification.to_dict()
        file_name = samplename + "_" + str(amplification.chromosome) + amplification.arm
        print(file_name,selected, reason)
        paths_to_save_dict = [config["path_amplifications_folder"] / (file_name+".json")]
        paths_to_save_figure = [config["path_folder_plot_amps"] / (file_name + ".png")]
        if selected:
            paths_to_save_dict.append(config["path_filtered_amplifications_folder"] / (file_name+".json"))
            paths_to_save_figure.append(config['path_folder_plot_selected_amps'] / (file_name + '.png'))
        else:
            paths_to_save_figure.append(config['path_folder_plot_non_selected_amps'] / reason/(file_name + '.png'))

        if selected_clock_like:
            paths_to_save_dict.append(config["path_filtered_amplifications_clocklike_folder"] / (file_name+".json"))
            paths_to_save_figure.append(config['path_folder_plot_selected_amps_clocklike'] / (file_name + '.png'))
        else:
            paths_to_save_figure.append(config['path_folder_plot_non_selected_amps_clocklike'] / reason_clock_like/(file_name + '.png'))

        for p in paths_to_save_dict:
            with open(p, 'w') as buff:
                json.dump(amplification_dict, buff, indent=4)
        margin = config["rel_margin_plot_amplification"] * (
                    amplification.segments[-1].end.position - amplification.segments[0].start.position)
        start_pos = amplification.segments[0].start.position - margin
        max_bases = (amplification.segments[-1].end.position - amplification.segments[
            0].start.position) + 2 * margin
        for i,path_save in enumerate(paths_to_save_figure):
            fig, ax = plt.subplots()
            s.plot_cn_profile(add_snvs=True,
                              total_cn=False,
                              chromosome=amplification.chromosome,
                              start_pos=start_pos,
                              max_bases=max_bases,
                              plot_threshold_amp=True,
                              differentiate_snv_type=True,
                              oncogenes=amplification.oncogenes,
                              title=file_name,
                              path_save=path_save,
                              ax=ax)
            plt.close(fig)
