import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("samplename", type=str, help="sample to plot")
    parser.add_argument('-S','--add_SNVs', dest='add_snvs', action='store_const',
                        const=True, default=False,
                        help='plot SNVs in the copy number profile figure')
    parser.add_argument('--total_cn', dest='total_cn', action='store_const',
                        const=True, default=False,
                        help='plot total copy number instead of major')
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
    parser.add_argument('--path_save', nargs='?', default=None, type=str,
                        help='path to store figure, if not specified, the figure is not stored')
    parser.add_argument('--show', nargs='?', default=True, type=bool,
                        help='bool that determines if figure is shown (if set to False --path_save should be specified'
                             'otherwise this script does not do anything useful)')
    parser.add_argument('--mark_clock_like_snvs', nargs='?', default=True, type=bool,
                        help='bool that determines if the (C->T)pG SNVs are marked')
    parser.add_argument('--title', nargs='?', default='default', type=str,
                        help='Title for figure, default value is sample name and chromosome if it is specified')
    return parser.parse_args()