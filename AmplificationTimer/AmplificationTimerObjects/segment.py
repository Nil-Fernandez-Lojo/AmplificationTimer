from .decorator_equality import class_equality_attributes
from .chromosome import Chromosome
from .position import Position
from .gene import Gene
from .snv import SNV
from scipy.stats import binomtest
import math
import copy



@class_equality_attributes
class Segment:
    """
    Class representing a genomic segment with constant copy number

    ...

    Attributes
    ----------
    clinical_data: dic
    config: dic
    chromosome: Chromosome
    start: Position
    end: Position
    minor_cn: int
    major_cn: int
    snvs: list of SNV
    genes: list of Gene

    Methods
    -------
    init:
    match_genes:
    add_snv:
    get_length:
    get_ploidy_healthy:
    to_dict: returns a copy of the object encoded as dictionary
    """

    def __init__(self, data, config, clinical_data, match_genes=False):
        self.clinical_data = clinical_data
        self.config = config
        self.chromosome = Chromosome(data['chromosome'])
        self.start = Position(self.chromosome, int(data['start']), config)
        self.end = Position(self.chromosome, int(data['end']), config)

        if math.isnan(data['minor_cn']):
            self.minor_cn = data['minor_cn']
        else:
            self.minor_cn = int(data['minor_cn'])
        if math.isnan(data['major_cn']):
            self.major_cn = data['minor_cn']
        else:
            self.major_cn = int(data['major_cn'])

        self.snvs = []
        for snv in data.get('snvs', []):
            chromosome = Chromosome(snv['chromosome'])
            position = Position(chromosome, snv['pos'], self.config)
            self.add_snv(SNV(chromosome,
                             position,
                             snv['ref_count'],
                             snv['alt_count'],
                             snv['ref_base'],
                             snv['alt_base'],
                             self.config,
                             snv['previous_ref_base'],
                             snv['next_ref_base'],
                             snv['in_kataegis'],
                             snv['intermediate_cn']))

        self.genes = []
        if 'genes' in data.keys():
            self.genes = [Gene(gene, config) for gene in data['genes']]
        if match_genes:
            self.match_genes()
        self.problematic_major = data.get('problematic_major', None)
        self.flagged_intermediate_snvs = data.get('flagged_intermediate_snvs', False)

    def match_genes(self):
        ensembl = self.config['ensembl']
        self.genes = []
        for gene in ensembl.genes_at_locus(contig=self.chromosome.chromosome, position=self.start.position,
                                           end=self.end.position):
            if gene.biotype == 'protein_coding':
                self.genes.append(Gene(gene, self.config))

    def add_snv(self, snv):
        self.snvs.append(snv)

    def to_dict(self):
        dic = dict()
        dic['chromosome'] = str(self.chromosome)
        dic['start'] = self.start.position
        dic['end'] = self.end.position
        dic['minor_cn'] = self.minor_cn
        dic['major_cn'] = self.major_cn
        dic['snvs'] = [snv.to_dict() for snv in self.snvs]
        if len(self.genes) > 0:
            dic['genes'] = [gene.to_dict() for gene in self.genes]
        dic['problematic_major'] = self.problematic_major
        dic['flagged_intermediate_snvs'] = self.flagged_intermediate_snvs
        return dic

    def get_length(self):
        return self.end.position - self.start.position

    def get_tot_cn(self):
        return self.minor_cn + self.major_cn

    def get_ploidy_healthy(self):
        if self.chromosome.chromosome == 'Y':
            if self.clinical_data['inferred_sex'] == 'male':
                return 1
            else:
                return 0
        elif self.chromosome.chromosome == 'X':
            if self.clinical_data['inferred_sex'] == 'female':
                return 2
            else:
                return 1
        else:
            return 2

    def check_if_problematic_major(self,rho,clock_like_filter):
        ploidy_nc = self.get_ploidy_healthy()
        M = self.major_cn
        m = self.minor_cn
        diff_cn_filter = max(self.config['min_diff_cn_problematic_major_filter'],
                             M*self.config['f_major_diff_cn_problematic_major_filter'])
        alt_count = 0
        tot_count = 0
        for snv in self.snvs:
            if ((not snv.in_kataegis)
                    and (not clock_like_filter or snv.clock_like())
                    and (snv.get_tot_count()!=0)
                    and snv.get_multiplicity(rho, M + m, ploidy_nc) > M - diff_cn_filter):
                alt_count+=snv.alt_count
                tot_count+=snv.ref_count+snv.alt_count
        if tot_count == 0:
            self.problematic_major = False
        else:
            q = M*rho / ((M+m) * rho + self.get_ploidy_healthy() * (1 - rho))
            pvalue = binomtest(alt_count, tot_count, q, alternative="two-sided").pvalue
            self.problematic_major = bool(pvalue < self.config['p_value_threshold'])
        return self.problematic_major

    def flag_intermediate_snvs(self):
        rho = self.clinical_data['purity']
        T = self.get_tot_cn()
        phi_nc = self.get_ploidy_healthy()
        for snv in self.snvs:
            if snv.get_tot_count() == 0:
                snv.intermediate_cn = True
            else:
                snv.intermediate_cn = ((not snv.multiplicity_accepted(rho, T, phi_nc, 1, alternative="greater")) and
                (not snv.multiplicity_accepted(rho, T, phi_nc, self.major_cn, alternative='two-sided')))
        self.flagged_intermediate_snvs = True

    def get_filtered_snvs(self,only_clock_like_SNVs,filter_APOBEC):
        list_snvs = []
        for snv in self.snvs:
            if snv.get_tot_count() == 0:
                continue
            elif (only_clock_like_SNVs and (not snv.clock_like())):
                continue
            elif filter_APOBEC and snv.intermediate_cn:
                continue
            else:
                list_snvs.append(snv)
        return list_snvs

    def get_n_snvs_1_copy(self,only_clock_like_SNVs):
        snvs = self.get_filtered_snvs(only_clock_like_SNVs,False)
        n = 0
        for s in snvs:
            if s.multiplicity_accepted(self.clinical_data['purity'], self.get_tot_cn(), self.get_ploidy_healthy(), 1,
                                       alternative='two-sided'):
                n += 1
        return n

    def subset(self, chromosome,start_pos, end_pos):
        if str(self.chromosome) != str(chromosome):
            return None
        elif end_pos <= self.start.position:
            return None
        elif start_pos >= self.end.position:
            return None
        else:
            segment = copy.deepcopy(self)
            segment.start = Position(segment.chromosome,max(start_pos,self.start.position),segment.config)
            segment.end = Position(segment.chromosome,min(end_pos, self.end.position),segment.config)
            snvs = copy.copy(segment.snvs)
            segment.snvs = []
            for snv in snvs:
                if segment.start <= snv.pos <= segment.end:
                    segment.add_snv(snv)
            return segment
