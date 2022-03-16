import numpy as np
import matplotlib.pyplot as plt

from .chromosome import Chromosome


def plot_cn(segments,
            config,
            purity,
            add_snvs=False,
            chromosome=None,
            start_pos=0,
            max_bases=None,
            total_cn=False,
            threshold_amplification=None,
            title=None,
            path_save=None,
            differentiate_snv_type=True,
            link_segments=False,
            ax = None):
    if ax is None:
        ax = plt.gca()

    if chromosome is not None:
        if not isinstance(chromosome, Chromosome):
            chromosome = Chromosome(chromosome)
        segments_to_plot = []
        for s in segments:
            if s.chromosome == chromosome:
                segments_to_plot.append(s)
    else:
        segments_to_plot = segments
    x_start = np.zeros(len(segments_to_plot))
    x_end = np.zeros(len(segments_to_plot))
    major = np.zeros(len(segments_to_plot))
    minor = np.zeros(len(segments_to_plot))
    for i, segment in enumerate(segments_to_plot):
        if chromosome is None:
            x_start[i] = segment.start.absolute_position
            x_end[i] = segment.end.absolute_position
        else:
            x_start[i] = segment.start.position
            x_end[i] = segment.end.position
        major[i] = segment.major_cn
        minor[i] = segment.minor_cn
    if total_cn:
        second_cn = major + minor
        label_second_cn = "Total"
    else:
        second_cn = major
        label_second_cn = "Major"
    if link_segments:
        minor_to_plot = [m for m in minor for _ in range(2)]
        second_cn_to_plot = [M for M in second_cn for _ in range(2)]
        x = []
        for i in range(len(segments_to_plot)):
            x.extend([x_start[i], x_end[i]])
        plt.plot(x, minor_to_plot, label='Minor', color='tab:orange')
        plt.plot(x, second_cn_to_plot, label=label_second_cn, color='tab:blue')
    else:
        for i in range(len(segments_to_plot)):
            x_i = [x_start[i], x_end[i]]
            if i == 0:
                label_1 = "Minor"
                label_2 = label_second_cn
            else:
                label_1 = None
                label_2 = None
            plt.plot(x_i, [minor[i]] * 2, label=label_1, color='tab:orange')
            plt.plot(x_i, [second_cn[i]] * 2, label=label_2, color='tab:blue')

    if add_snvs:
        x_snvs = []
        mcn = []
        x_snvs_clock_like = []
        mcn_clock_like = []
        for s in segments_to_plot:
            for snv in s.SNVs:
                if differentiate_snv_type and snv.clock_like():
                    x = x_snvs_clock_like
                    m = mcn_clock_like
                else:
                    x = x_snvs
                    m = mcn
                if (snv.alt_count + snv.ref_count) == 0:
                    continue
                if chromosome is None:
                    x.append(snv.pos.absolute_position)
                else:
                    x.append(snv.pos.position)
                vaf = snv.alt_count / (snv.alt_count + snv.ref_count)
                tot_cn = s.major_cn + s.minor_cn
                rho = purity
                m.append(vaf * (tot_cn + ((1 - rho) / rho) * s.get_ploidy_healthy()))
        plt.plot(x_snvs, mcn, marker='.', color='k', linestyle='None', label='SNVs')
        if differentiate_snv_type:
            plt.plot(x_snvs_clock_like,
                     mcn_clock_like,
                     marker='.',
                     color='g',
                     linestyle='None',
                     label='SNVs (C->T)pG')

    plt.axhline(y=0, color='k', linestyle='--', linewidth=0.4)
    plt.axhline(y=1, color='k', linestyle='--', linewidth=0.4)
    if total_cn:
        plt.axhline(y=2, color='k', linestyle='--', linewidth=0.4)
    if threshold_amplification is not None:
        plt.axhline(y=threshold_amplification, color='k', linestyle='--', linewidth=0.4)
    if chromosome is None:
        for i in range(len(config["cum_length_chr"])):
            plt.axvline(x=config["cum_length_chr"][i], color='k', linestyle='--', linewidth=0.4)
        ax.set_xticks(config["cum_length_chr"])
        ax.set_xticklabels('')
        ax.set_xticks([(config["cum_length_chr"][i] + config["cum_length_chr"][i + 1]) / 2 for i in
                       range(len(config["cum_length_chr"]) - 1)], minor=True)
        ax.set_xticklabels(list(range(1, 23)) + ['X', 'Y'], minor=True)
        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)
            tick.label1.set_horizontalalignment('center')
        ax.set_xlabel('Chromosome')
    else:
        if start_pos != 0 or max_bases is not None:
            if max_bases is None:
                ax.set_xlim(start_pos, 1.1 * segments_to_plot[-1].position.position)
            else:
                ax.set_xlim(start_pos, start_pos + max_bases)
        ax.set_xlabel('Position')
    ax.set_ylabel('Copy number')
    plt.legend()
    if title is not None:
        plt.title(title)
    if path_save is not None:
        plt.savefig(path_save, bbox_inches='tight', dpi=300)


def plot_normalised_mu(pos,normalised_mu,mu_hat,config,chromosome=None,ax= None):
    if ax is None:
        ax = plt.gca()
    x = [p.absolute_position for p in pos]
    ax.plot(x, normalised_mu, label='normalised mu')
    if chromosome is None:
        for i in range(len(config["cum_length_chr"])):
            ax.axvline(x=config["cum_length_chr"][i], color='k', linestyle='--', linewidth=0.4)
        ax.set_xticks(config["cum_length_chr"])
        ax.set_xticklabels('')
        ax.set_xticks([(config["cum_length_chr"][i] + config["cum_length_chr"][i + 1]) / 2 for i in
                       range(len(config["cum_length_chr"]) - 1)], minor=True)
        ax.set_xticklabels(list(range(1, 23)) + ['X', 'Y'], minor=True)
        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)
            tick.label1.set_horizontalalignment('center')
        ax.set_xlabel('Chromosome')
    ax.axhline(y=mu_hat, color='k', label='estimated mu from 1+1 and 1+0')
    ax.set_ylabel('Nomalised mu')
    ax.legend()
