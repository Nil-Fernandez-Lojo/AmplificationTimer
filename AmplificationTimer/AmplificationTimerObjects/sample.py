import json
import numpy as np
import pandas as pd

from .decorator_equality import class_equality_attributes

from .snv import SNV
from .amplification import Amplification
from .chromosome import Chromosome
from .position import Position
from .segment import Segment
from .mutation_rate import MutationRate
from .plot import plot_cn


@class_equality_attributes
class Sample:

    def __init__(self, config, name, save=False):
        self.name = name
        self.config = config
        path_preprocessed_data = config['path_preprocessed_samples'] / (name + '.json')
        if path_preprocessed_data.exists():
            self.load_preprocessed_data(path_preprocessed_data)
        else:
            self.preprocess_data()
            if save:
                self.save()

    def preprocess_data(self):
        self.clinical_data = \
            self.config['summary_table'].loc[self.config['summary_table']['samplename'] == self.name].iloc[0].to_dict()

        for key in self.clinical_data.keys():  # used since json can not be saved with np.bool_ types
            if isinstance(self.clinical_data[key], np.bool_):
                self.clinical_data[key] = bool(self.clinical_data[key])
        self.segments = self.load_cnas()
        snvs, self.sample_type = self.load_snvs()  # TODO : what to do if gray list
        self.match_snvs_to_segments(snvs)

        self.major_cn_mode = self.get_major_cn_mode(self.segments)
        if self.major_cn_mode <= 1:
            self.threshold_amplification = self.config['threshold_amplification_no_WGD']
        else:
            self.threshold_amplification = self.config['threshold_amplification_WGD']
        self.mutation_rate = self.get_mutation_rate(self.segments)
        self.subclonal_structure = pd.read_csv(
            self.config["path_folder_subclonal_structure"] / (self.name + self.config['suffix_subclonal_structure']),
            sep='\t')
        self.amplifications = self.get_amplifications(self.segments, self.threshold_amplification)

    def load_cnas(self):
        df = pd.read_csv(self.config['path_folder_cna'] / (self.name + self.config['suffix_cnas']), sep='\t')
        n_segments = df.shape[0]
        list_segments = []

        for i in range(n_segments):
            list_segments.append(Segment(df.iloc[i], self.config, self.clinical_data))
            if i > 0:
                if list_segments[i].start <= list_segments[i - 1].end:
                    raise Exception('CNAs not ordered in input file')

        return list_segments

    def load_snvs(self):
        filename = self.name + self.config['suffix_snvs']
        vcf_path = self.config['path_folder_snv_white'] / filename
        if vcf_path.exists():
            return self.parse_vcf(vcf_path), 'ICGC_white'

        vcf_path = self.config['path_folder_snv_gray'] / filename
        if vcf_path.exists():
            return self.parse_vcf(vcf_path), 'ICGC_gray'

        vcf_path = self.config['path_folder_snv_TCGA'] / filename
        if vcf_path.exists():
            return self.parse_vcf(vcf_path), 'TCGA'
        else:
            return [], 'missing'

    def parse_vcf(self, vcf_path):
        df = pd.read_csv(vcf_path, sep='\t', comment='#', header=None, names=self.config['vcf_column_names'])
        snvs = []
        for index, row in df.iterrows():

            # only allow main contigs chr 1->22 + X and Y
            if str(row['CHROM']) not in self.config['list_chromosome']:
                print('SNV not loaded since chromosome:', row['CHROM'])
                continue

            chromosome = Chromosome(row['CHROM'])
            pos = Position(chromosome, row['POS'], self.config)
            info = row['INFO'].split(';')
            t_ref_count = 0
            for s in info:
                if s.startswith(self.config["vcf_prefix_counts"]["ref_count"]):
                    t_ref_count = int(s[len(self.config["vcf_prefix_counts"]["ref_count"]):])
                    break
            t_alt_count = 0
            for s in info:
                if s.startswith(self.config["vcf_prefix_counts"]["alt_count"]):
                    t_alt_count = int(s[len(self.config["vcf_prefix_counts"]["alt_count"]):])
                    break
            snvs.append(SNV(chromosome,
                            pos,
                            t_ref_count,
                            t_alt_count,
                            row['REF'],
                            row['ALT'],
                            self.config))
        snvs = sorted(snvs, key=lambda x: x.pos)
        return snvs

    def match_snvs_to_segments(self, snvs):
        i = 0  # segment index
        for snv in snvs:
            while snv.pos > self.segments[i].end:
                i += 1
                if i == len(self.segments):
                    # last SNVs not matched
                    return

            if snv.pos >= self.segments[i].start:
                self.segments[i].add_snv(snv)

    @staticmethod
    def get_major_cn_mode(segments, max_cn=5):
        # returns mode of major copy number (where all copy numbers higher than max_cn are mapped back to max_cn)
        n_bases_cn = [0 for _ in range(max_cn + 1)]
        for segment in segments:
            cn = min(max_cn, segment.major_cn)
            n_bases_cn[cn] += segment.get_length()
        return n_bases_cn.index(max(n_bases_cn))

    def get_amplifications(self, segments, threshold_amplification):
        potential_amplifications = dict()
        for chromosome in (list(range(1, 23)) + ['X', 'Y']):
            potential_amplifications[str(chromosome)] = dict()
            for arm in ['p', 'q']:
                potential_amplifications[str(chromosome)][arm] = Amplification(Chromosome(chromosome),
                                                                               arm,
                                                                               [],
                                                                               threshold_amplification,
                                                                               self.clinical_data,
                                                                               self.config,
                                                                               self.mutation_rate,
                                                                               self.subclonal_structure)

        for segment in segments:
            if segment.major_cn >= threshold_amplification:
                if segment.start.arm != segment.end.arm:
                    print('Warning, amplified segment that spans both arms', self.name, segment.start, segment.end)
                segment.match_genes()
                potential_amplifications[segment.chromosome.chromosome][segment.start.arm].add_segment(segment)

        amplifications = []
        for chromosome in (list(range(1, 23)) + ['X', 'Y']):
            for arm in ['p', 'q']:
                if len(potential_amplifications[str(chromosome)][arm].segments) != 0:
                    potential_amplifications[str(chromosome)][arm].set_oncogenes()
                    amplifications.append(potential_amplifications[str(chromosome)][arm])
        return amplifications

    @staticmethod
    def get_mutation_rate(segments):
        # TODO: maybe compute here the different subclonal snvs?
        n_1_0 = 0
        n_1_0_ctpg = 0
        l_1_0 = 0
        n_1_1 = 0
        n_1_1_ctpg = 0
        l_1_1 = 0
        for segment in segments:
            if segment.major_cn == 1:
                if segment.minor_cn == 0:
                    n_1_0 += len(segment.SNVs)
                    l_1_0 += segment.get_length()
                    for snv in segment.SNVs:
                        if snv.clock_like():
                            n_1_0_ctpg += 1
                else:
                    n_1_1 += len(segment.SNVs)
                    l_1_1 += segment.get_length()
                    for snv in segment.SNVs:
                        if snv.clock_like():
                            n_1_1_ctpg += 1

        return MutationRate(n_1_0, l_1_0, n_1_1, l_1_1, n_1_0_ctpg, n_1_1_ctpg)

    def plot_cn_profile(self,
                        add_snvs=False,
                        chromosome=None,
                        start_pos=0,
                        max_bases=None,
                        total_cn=False,
                        plot_threshold_amp=True,
                        title=None,
                        path_save=None,
                        differentiate_snv_type=True,
                        ax=None):
        if plot_threshold_amp:
            threshold_amplification = self.threshold_amplification
        else:
            threshold_amplification = None
        plot_cn(self.segments,
                self.config,
                self.clinical_data['purity'],
                add_snvs=add_snvs,
                chromosome=chromosome,
                start_pos=start_pos,
                max_bases=max_bases,
                total_cn=total_cn,
                threshold_amplification=threshold_amplification,
                title=title,
                path_save=path_save,
                differentiate_snv_type=differentiate_snv_type,
                link_segments=True,
                ax=ax)

    def load_preprocessed_data(self, path_preprocessed_data):
        with open(self.config['path_preprocessed_samples'] / (self.name + '.json'), 'r') as fp:
            data = json.load(fp)
        self.clinical_data = data['clinical_data']
        self.sample_type = data['SNVs_calls_type']
        self.segments = [Segment(segment, self.config, self.clinical_data) for segment in data['segments']]
        self.major_cn_mode = data['major_cn_mode']
        self.threshold_amplification = data['threshold_amplification']
        self.subclonal_structure = pd.DataFrame.from_dict(data['subclonal_structure'])
        self.mutation_rate = MutationRate(data['mutation_rate']['n_1_0'],
                                          data['mutation_rate']['l_1_0'],
                                          data['mutation_rate']['n_1_1'],
                                          data['mutation_rate']['l_1_1'],
                                          data['mutation_rate']['n_1_0_ctpg'],
                                          data['mutation_rate']['n_1_1_ctpg'])

        self.amplifications = []
        for amplification in data['amplifications']:
            segments = [Segment(segment, self.config, self.clinical_data) for segment in amplification['segments']]
            self.amplifications.append((Amplification(Chromosome(amplification['chromosome']),
                                                      amplification['arm'],
                                                      segments,
                                                      self.threshold_amplification,
                                                      self.clinical_data,
                                                      self.config,
                                                      self.mutation_rate,
                                                      self.subclonal_structure)))
            self.amplifications[-1].set_oncogenes()  # TODO information is there, shouldn't need to search again

    def to_dict(self):
        dic = dict()
        dic['clinical_data'] = self.clinical_data.copy()
        dic['SNVs_calls_type'] = self.sample_type
        dic['segments'] = [segment.to_dict() for segment in self.segments]
        dic['major_cn_mode'] = self.major_cn_mode
        dic['threshold_amplification'] = self.threshold_amplification
        dic['amplifications'] = [amplification.to_dict() for amplification in self.amplifications]
        dic['mutation_rate'] = self.mutation_rate.to_dict()
        dic['subclonal_structure'] = self.subclonal_structure.to_dict('records')
        return dic

    def save(self):
        if not self.config['path_preprocessed_samples'].exists():
            self.config['path_preprocessed_samples'].mkdir(parents=True, exist_ok=True)

        with open(self.config['path_preprocessed_samples'] / (self.name + '.json'), 'w') as fp:
            json.dump(self.to_dict(), fp, indent=4)
