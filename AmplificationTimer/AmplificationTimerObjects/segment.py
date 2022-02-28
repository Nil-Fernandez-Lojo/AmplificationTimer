from .decorator_equality import class_equality_attributes
from .chromosome import Chromosome
from .position import Position
from .gene import Gene
from .snv import SNV
import math


@class_equality_attributes
class Segment():
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
	SNVs: list of SNV
	genes: list of Gene
	
	Methods
	-------
	init:
	match_genes:
	add_SNV:
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

        self.SNVs = []
        for snv in data.get('snvs', []):
            chromosome = Chromosome(snv['chromosome'])
            position = Position(chromosome, snv['pos'], self.config)
            self.add_SNV(SNV(chromosome,
                             position,
                             snv['ref_count'],
                             snv['alt_count'],
                             snv['ref_base'],
                             snv['alt_base'],
                             self.config,
                             snv['previous_ref_base'],
                             snv['next_ref_base']))

        self.genes = []
        if 'genes' in data.keys():
            self.genes = [Gene(gene, config) for gene in data['genes']]
        if match_genes:
            self.match_genes()

    def match_genes(self):
        ensembl = self.config['ensembl']
        self.genes = []
        for gene in ensembl.genes_at_locus(contig=self.chromosome.chromosome, position=self.start.position,
                                           end=self.end.position):
            if gene.biotype == 'protein_coding':
                self.genes.append(Gene(gene, self.config))

    def add_SNV(self, SNV):
        self.SNVs.append(SNV)

    def to_dict(self):
        dic = dict()
        dic['chromosome'] = str(self.chromosome)
        dic['start'] = self.start.position
        dic['end'] = self.end.position
        dic['minor_cn'] = self.minor_cn
        dic['major_cn'] = self.major_cn
        dic['snvs'] = [snv.to_dict() for snv in self.SNVs]
        if len(self.genes) > 0:
            dic['genes'] = [gene.to_dict() for gene in self.genes]
        return dic

    def get_length(self):
        return self.end.position - self.start.position

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
