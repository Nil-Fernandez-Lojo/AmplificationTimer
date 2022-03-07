import argparse
import numpy as np
from pathlib import Path

def parse_arguments_plot(plot_type):
    parser = argparse.ArgumentParser()
    parser.add_argument("samplename", type=str, help="sample to plot")
    if plot_type == 'sample':
        parser.add_argument('--chromosome', nargs='?', default=None, choices=list(map(str, list(range(1,23))))+['X', 'Y'],
                            help='optional argument to restrict the plot to a single argument')
        parser.add_argument('--start_pos', nargs='?', default=0, type=int,
                            help='optional argument only used if --chromosome is used. It starts the plot at the given value')
        parser.add_argument('--max_bases', nargs='?', default=None, type=int,
                            help='optional argument only used if --chromosome is used. It only plots the profile for that '
                                 'number of bases')
        parser.add_argument('--plot_threshold_amp', nargs='?', default=True, type=bool,
                            help='bool that determines if the threshold used for calling an amplification is plot on the figure'
                            )
    elif plot_type == 'amplification':
        parser.add_argument("chromosome", choices=list(map(str, list(range(1,23))))+['X', 'Y'],
                            help='Chromosome of the amplification')
        parser.add_argument("arm", choices=['p', 'q'], help='Chromosome arm of the amplification')
        parser.add_argument('--no_context', dest='add_context', action='store_const',
                            const=False, default=True,
                            help='with this flag, segments non amplified in between amplified segments are not plotted')
        parser.add_argument('--rel_margin', nargs='?', default=0.1, type=float,
                            help='fraction of nucleotides before and after the amplification added to the plot')
    parser.add_argument('--path_save', nargs='?', default=None, type=str,
                        help='path to store figure, if not specified, the figure is not stored')
    parser.add_argument('--mark_clock_like_snvs', nargs='?', default=True, type=bool,
                        help='bool that determines if the (C->T)pG SNVs are marked')
    parser.add_argument('--title', nargs='?', default='default', type=str,
                        help='Title for figure, default value is sample name and chromosome if it is specified')
    parser.add_argument('-S','--add_SNVs', dest='add_snvs', action='store_const',
                        const=True, default=False,
                        help='plot SNVs in the copy number profile figure')
    parser.add_argument('--total_cn', dest='total_cn', action='store_const',
                        const=True, default=False,
                        help='plot total copy number instead of major')
    args = parser.parse_args()
    if args.title == 'default':
        args.title = args.samplename
        if args.chromosome is not None:
            args.title += ' ' + args.chromosome
        if args.arm is not None:
            args.title += args.arm
    return args

def parse_arguments_simulation():
    parser = argparse.ArgumentParser()
    parser.add_argument("samplename", type=str, help="sample on which to perform the simulation")
    # TODO: chromosome and arm should be optional arguments
    parser.add_argument("chromosome", type=str,
                        help="chromosome of the amplification on which to perform the simulation")
    parser.add_argument("arm", type=str,
                        help="chromosome arm of the amplification on which to perform the simulation")
    parser.add_argument('--model_idx', nargs='?', default=4, type=int,
                        help='model used for simulation')
    parser.add_argument('--expected_mut_list', nargs='*', default=[10, 100, 1000, 10000], type=int,
                        help='expected number of mutation on amplified segments, it defines the mutation rate used in the simulation')
    parser.add_argument('--list_t', nargs='*', default=np.linspace(0.0, 1.0, num=21), type=float,
                        help='values of t on which to perform simulation')
    parser.add_argument('--no_filter_APOBEC', dest='filter_APOBEC', action='store_const',
                        const=False, default=True,
                        help='with this flag, no filter APOBEC SNVs')
    parser.add_argument('--no_modelling_subclonality', dest='subclonality_modelled', action='store_const',
                        const=False, default=True,
                        help='with this flag, no modelling of subclones')
    parser.add_argument('--folder_out', nargs='?', default='../inference_results/simulation/', type=str,
                        help='folder to save simulation results, default: inference_results/simulation')
    args = parser.parse_args()
    args.folder_out = Path(args.folder_out)
    return args
