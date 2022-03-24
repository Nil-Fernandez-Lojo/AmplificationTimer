import numpy as np
import matplotlib.pyplot as plt
import math
from .chromosome import Chromosome


# Internal functions
def plot_segments(ax,segments_to_plot, link_segments, total_cn, unit_bases, chromosome,flag_problematic_major):
    line_style_flagged = 'dotted'
    line_style_not_flagged = 'solid'

    n = len(segments_to_plot)
    x = np.zeros(2*n)
    flagged_out = np.zeros(n,dtype=bool)
    major = np.zeros(n)
    minor = np.zeros(n)
    for i, segment in enumerate(segments_to_plot):
        if chromosome is None:
            x[2*i] = segment.start.absolute_position
            x[2*i+1] = segment.end.absolute_position
        else:
            x[2*i] = segment.start.position / unit_bases
            x[2*i+1] = segment.end.position / unit_bases
        major[i] = segment.major_cn
        minor[i] = segment.minor_cn
        flagged_out[i] = segment.problematic_major
    if total_cn:
        second_cn = major + minor
        label_second_cn = "Total"
    else:
        second_cn = major
        label_second_cn = "Major"
    minor_to_plot = [m for m in minor for _ in range(2)]
    second_cn_to_plot = [M for M in second_cn for _ in range(2)]
    for i in range(n):
        if i == 0:
            label_1 = "Minor"
            label_2 = label_second_cn
        else:
            label_1 = None
            label_2 = None
        linestyle = line_style_flagged if (flagged_out[i] and flag_problematic_major) else line_style_not_flagged
        ax.plot(x[2*i:(2*(i+1))], minor_to_plot[2*i:(2*(i+1))], linestyle = linestyle,label=label_1, color='tab:orange')
        ax.plot(x[2*i:(2*(i+1))], second_cn_to_plot[2*i:(2*(i+1))], linestyle = linestyle,label=label_2, color='tab:blue')
        if (i<(n-1)) and link_segments:
            ax.plot([x[2*i+1]]*2,[minor_to_plot[2*i], minor_to_plot[2*(i+1)]],label=None, color='tab:orange')
            ax.plot([x[2*i+1]]*2, [second_cn_to_plot[2*i], second_cn_to_plot[2*(i+1)]], label=None, color='tab:blue')


def plot_snvs(ax,
              segments_to_plot,
              chromosome,
              purity,
              differentiate_snv_type,
              mark_intermediate_snvs,
              unit_bases,
              mark_kataegis):
    pos = []
    mcn = []
    non_intermediate_snvs = []
    clock_like = []
    kataegis = []
    in_problematic_segment = []
    for s in segments_to_plot:
        tot_cn = s.major_cn + s.minor_cn
        if math.isnan(tot_cn) or tot_cn == 0:
            continue
        ploidy_healthy = s.get_ploidy_healthy()
        if not s.flagged_intermediate_snvs:
            s.flag_intermediate_snvs()

        for snv in s.snvs:
            if (snv.alt_count + snv.ref_count) == 0:
                continue
            kataegis.append(mark_kataegis and snv.in_kataegis)
            clock_like.append(differentiate_snv_type and snv.clock_like())
            non_intermediate_snvs.append((not snv.intermediate_cn) or (not mark_intermediate_snvs))
            if chromosome is None:
                p = snv.pos.absolute_position
            else:
                p = snv.pos.position / unit_bases
            pos.append(p)
            mcn.append(snv.get_multiplicity(purity, tot_cn, ploidy_healthy))
            if s.problematic_major is None:
                in_problematic_segment.append(False)
            else:
                in_problematic_segment.append(s.problematic_major)
    if len(pos) == 0:
        return
    kataegis = np.asarray(kataegis)
    non_intermediate_snvs = np.asarray(non_intermediate_snvs)
    clock_like = np.asarray(clock_like)
    pos = np.asarray(pos)
    mcn = np.asarray(mcn)
    in_problematic_segment = np.asarray(in_problematic_segment)

    m1 = non_intermediate_snvs & clock_like & (~kataegis) & (~in_problematic_segment)
    m2 = (~non_intermediate_snvs) & clock_like & (~kataegis) & (~in_problematic_segment)
    m3 = non_intermediate_snvs & (~clock_like) & (~kataegis) & (~in_problematic_segment)
    m4 = (~non_intermediate_snvs) & (~clock_like) & (~kataegis) & (~in_problematic_segment)
    m5 = non_intermediate_snvs & clock_like & (~kataegis) & in_problematic_segment
    m6 = (~non_intermediate_snvs) & clock_like & (~kataegis) & in_problematic_segment
    m7 = non_intermediate_snvs & (~clock_like) & (~kataegis) & in_problematic_segment
    m8 = (~non_intermediate_snvs) & (~clock_like) & (~kataegis) & in_problematic_segment
    m9 = kataegis & (~in_problematic_segment)
    m10 = kataegis & in_problematic_segment


    s = 20
    alpha_problematic_major = 0.3
    ax.scatter(pos[m1], mcn[m1], marker='^', color="g", ec='g', s=s, label='(C->T)pG accepted')
    ax.scatter(pos[m2], mcn[m2], marker='^', color="w", ec='r', s=s, label='(C->T)pG filtererd out')
    ax.scatter(pos[m3], mcn[m3], marker='o', color="g", ec='g', s=s, label='SNVs accepted')
    ax.scatter(pos[m4], mcn[m4], marker='o', color="w", ec='r', s=s, label='SNVs filtererd out')
    ax.scatter(pos[m9], mcn[m9], marker='P', color="r", s=s, label='SNVs in kataegis')

    ax.scatter(pos[m5], mcn[m5], marker='^', color="g", ec='g', s=s,alpha=alpha_problematic_major)
    ax.scatter(pos[m6], mcn[m6], marker='^', color="w", ec='r', s=s,alpha=alpha_problematic_major)
    ax.scatter(pos[m7], mcn[m7], marker='o', color="g", ec='g', s=s,alpha=alpha_problematic_major)
    ax.scatter(pos[m8], mcn[m8], marker='o', color="w", ec='r', s=s,alpha=alpha_problematic_major)
    ax.scatter(pos[m10], mcn[m10], marker='P', color="r", s=s,alpha=alpha_problematic_major)


def get_segments_to_plot(chromosome, segments, start_pos, max_bases):
    if chromosome is not None:
        if not isinstance(chromosome, Chromosome):
            chromosome = Chromosome(chromosome)
        segments_to_plot = []
        for s in segments:
            if ((s.chromosome == chromosome) and
                    (s.end.position > start_pos)):
                if (max_bases is None) or (s.start.position < (start_pos + max_bases)):
                    segments_to_plot.append(s)
    else:
        segments_to_plot = segments
    return segments_to_plot


def plot_horizontal_lines(ax, total_cn, threshold_amplification):
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.4)
    ax.axhline(y=1, color='k', linestyle='--', linewidth=0.4)
    if total_cn:
        ax.axhline(y=2, color='k', linestyle='--', linewidth=0.4)
    if threshold_amplification is not None:
        ax.axhline(y=threshold_amplification, color='k', linestyle='--', linewidth=0.4)


def plot_oncogenes(ax, oncogenes, chromosome, unit_bases):
    if oncogenes is not None:
        for i, gene in enumerate(oncogenes):
            if chromosome is None:
                x = (gene.start.absolute_position + gene.end.absolute_position) / 2
            else:
                x = (gene.start.position + gene.end.position) / (2 * unit_bases)
            if i == 0:
                ax.axvline(x=x, linestyle='dotted', color='r', linewidth=2, label='oncogene')
            else:
                ax.axvline(x=x, linestyle='dotted', color='r', linewidth=2)


def plot_chr_limits(ax, config, x_axis='position', cum_sum_snvs_per_chr=None):
    if x_axis == 'position':
        x = config["cum_length_chr"]
    elif x_axis == 'snv_index':
        x = cum_sum_snvs_per_chr

    for i in range(len(x)):
        ax.axvline(x=x[i], color='k', linestyle='--', linewidth=0.4)


def set_x_ticks_chromosomes(ax, config, x_axis='position', cum_sum_snvs_per_chr=None):
    if x_axis == 'position':
        x = config["cum_length_chr"]
    elif x_axis == 'snv_index':
        x = cum_sum_snvs_per_chr

    ax.set_xticks(x)
    ax.set_xticklabels('')
    ax.set_xticks([(x[i] + x[i + 1]) / 2 for i in range(len(config['list_chromosome']))], minor=True)

    ax.set_xticklabels(config['list_chromosome'], minor=True)
    for tick in ax.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0)
        tick.tick2line.set_markersize(0)
        tick.label1.set_horizontalalignment('center')


def set_x_limits(ax, segments_to_plot, start_pos, max_bases, unit_bases):
    if start_pos != 0 or max_bases is not None:
        if max_bases is None:
            ax.set_xlim(start_pos / unit_bases, 1.1 * segments_to_plot[-1].position.position / unit_bases)
        else:
            ax.set_xlim(start_pos / unit_bases, (start_pos + max_bases) / unit_bases)


def plot_snvs_rainfall(ax, snvs):
    y = []
    kataegis = []
    for i in range(len(snvs) - 1):
        kataegis.append(snvs[i].in_kataegis)
        y.append(snvs[i + 1].pos.absolute_position - snvs[i].pos.absolute_position)
    kataegis = np.asarray(kataegis)
    y = np.asarray(y)
    x = np.arange(len(snvs)-1)
    s = 20
    ax.scatter(x[~kataegis], y[~kataegis], marker='o', color="g", ec='g', s=s, label='non kataegis')
    ax.scatter(x[kataegis], y[kataegis], marker='P', color="r", s=s, label='kataegis')


def get_cum_sum_snvs_per_chr(segments_to_plot, config):
    snvs_per_chr = dict.fromkeys(config['list_chromosome'], 0)
    for segment in segments_to_plot:
        snvs_per_chr[str(segment.chromosome)] += len(segment.snvs)

    cum_sum_snvs_per_chr = [0]
    for chr_str in config['list_chromosome']:
        cum_sum_snvs_per_chr.append(cum_sum_snvs_per_chr[-1] + snvs_per_chr[chr_str])
    return cum_sum_snvs_per_chr


# Main functions
def plot_cn(segments, config, purity, add_snvs=False, chromosome=None, start_pos=0, max_bases=None, total_cn=False,
            threshold_amplification=None, title=None, path_save=None, differentiate_snv_type=True, link_segments=False,
            ax=None, oncogenes=None, mark_intermediate_snvs=True, mark_kataegis=True,flag_problematic_major=True):
    # TODO: give option to user to change this, careful also need to change label x axis in this cas
    unit_bases = 10 ** 6

    if ax is None:
        ax = plt.gca()
    segments_to_plot = get_segments_to_plot(chromosome, segments, start_pos, max_bases)
    plot_segments(ax,segments_to_plot, link_segments, total_cn, unit_bases, chromosome,flag_problematic_major)
    if add_snvs:
        plot_snvs(ax, segments_to_plot, chromosome, purity, differentiate_snv_type, mark_intermediate_snvs, unit_bases,
                  mark_kataegis)
    plot_horizontal_lines(ax, total_cn, threshold_amplification)
    plot_oncogenes(ax, oncogenes, chromosome, unit_bases)

    if chromosome is None:
        plot_chr_limits(ax, config)
        set_x_ticks_chromosomes(ax, config)
        ax.set_xlabel('Chromosome')
    else:
        set_x_limits(ax, segments_to_plot, start_pos, max_bases, unit_bases)
        ax.ticklabel_format(axis='x', style='plain')
        ax.set_xlabel('Position (Mb)')
    ax.set_ylabel('Copy number')
    #ax.legend()
    if title is not None:
        ax.set_title(title)
    if path_save is not None:
        # TODO this is dangerous
        plt.savefig(path_save, bbox_inches='tight', dpi=300)


def plot_normalised_mu(pos, normalised_mu, mu_hat, config, chromosome=None, ax=None):
    if ax is None:
        ax = plt.gca()
    x = [p.absolute_position for p in pos]
    ax.plot(x, normalised_mu, label='normalised mu')
    if chromosome is None:
        plot_chr_limits(ax, config)
        set_x_ticks_chromosomes(ax, config)
        ax.set_xlabel('Chromosome')
    ax.axhline(y=mu_hat, color='k', label='estimated mu from 1+1 and 1+0')
    ax.set_ylabel('Nomalised mu')
    ax.legend()


def plot_rainfall(segments, config, chromosome=None, start_pos=0, max_bases=None, title=None, path_save=None, ax=None):
    if ax is None:
        ax = plt.gca()
    segments_to_plot = get_segments_to_plot(chromosome, segments, start_pos, max_bases)
    snvs = []
    for s in segments_to_plot:
        snvs.extend(s.snvs)
    plot_snvs_rainfall(ax, snvs)

    if chromosome is None:
        cum_sum_snvs_per_chr = get_cum_sum_snvs_per_chr(segments_to_plot, config)
        plot_chr_limits(ax, config, x_axis='snv_index', cum_sum_snvs_per_chr=cum_sum_snvs_per_chr)
        set_x_ticks_chromosomes(ax, config, x_axis='snv_index', cum_sum_snvs_per_chr=cum_sum_snvs_per_chr)
        ax.set_xlabel('Chromosome')
    ax.axhline(y=config['mean_distance_kataegis_threshold'], color='k', linestyle='--', linewidth=0.4)
    ax.set_yscale('log')
    ax.set_ylabel('inter snv distance')
    ax.set_xlabel('SNV position')
    ax.legend()
    if title is not None:
        ax.set_title(title)
    if path_save is not None:
        # TODO this is dangerous
        plt.savefig(path_save, bbox_inches='tight', dpi=300)
