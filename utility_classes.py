import math
from functools import total_ordering

from decorator_equality import class_equality_attributes

@class_equality_attributes
class Mutation_rate():
	def __init__(self,n_1_0,l_1_0,n_1_1,l_1_1):
		self.n_1_0 = n_1_0
		self.l_1_0 = l_1_0
		self.n_1_1 = n_1_1
		self.l_1_1 = l_1_1

	def get_ML(self):
		return (self.n_1_0+self.n_1_1)/(self.l_1_0+2*self.l_1_1)
	
	def get_beta_posterior_parameters(self):
		return (self.n_1_0+self.n_1_1+1,self.l_1_0+2*self.l_1_1+1)

	def to_dict(self):
		dic = dict()
		dic['n_1_0'] = self.n_1_0
		dic['l_1_0'] = self.l_1_0
		dic['n_1_1'] = self.n_1_1
		dic['l_1_1'] = self.l_1_1
		return dic

@total_ordering
@class_equality_attributes
class Chromosome:
	def __init__(self,chromosome):
		self.chromosome = str(chromosome)

		if self.chromosome == 'X':
			self.c = 23
		elif self.chromosome == 'Y':
			self.c = 24
		else:
			self.c = int(chromosome)

	def __gt__(self, other):
		return self.c > other.c

	def __str__(self):
		return self.chromosome

@total_ordering
@class_equality_attributes
class Position:
	def __init__(self,chromosome,position,config):
		self.chromosome = chromosome
		self.position = position

		chr_arm_l = config['chromosome_arm_length']
		chr_p_arm_l = chr_arm_l.loc[chr_arm_l['arm'] == 'p']
		if position <int(chr_p_arm_l['length'][chr_arm_l['chromosome'] == str(chromosome)]):
			self.arm = 'p'
		else:
			self.arm = 'q'

	def __gt__(self, other):
		if self.chromosome == other.chromosome:
			return self.position > other.position
		else:
			return self.chromosome > other.chromosome

	def __str__(self):
		return "chr: "+str(self.chromosome)+" pos: "+str(self.position)+ " arm: "+self.arm

@class_equality_attributes
class Segment():

	def __init__(self,data,config,clinical_data,match_genes = False):
		self.clinical_data = clinical_data
		self.config = config
		self.chromosome = Chromosome(data['chromosome'])
		self.start = Position(self.chromosome,int(data['start']),config)
		self.end = Position(self.chromosome, int(data['end']),config)

		if math.isnan(data['minor_cn']):
			self.minor_cn = data['minor_cn']
		else:
			self.minor_cn = int(data['minor_cn'])
		if math.isnan(data['major_cn']):
			self.major_cn = data['minor_cn']
		else:
			self.major_cn = int(data['major_cn'])
		
		self.SNVs = []
		for snv in data.get('snvs',[]):
			self.add_SNV(SNV(snv['chromosome'],
				snv['pos'],
				snv['ref_count'],
				snv['alt_count'],
				snv['ref_base'],
				snv['alt_base'],
				config))

		self.genes = []
		if 'genes' in data.keys():
			self.genes = [Gene(gene,config) for gene in data['genes']]
		if match_genes:
			 self.match_genes()

	def match_genes(self):
		ensembl = self.config['ensembl']
		self.genes = []
		for gene in ensembl.genes_at_locus(contig=self.chromosome.chromosome, position=self.start.position,end=self.end.position):
			if gene.biotype == 'protein_coding':
				self.genes.append(Gene(gene,self.config))

	def add_SNV(self,SNV):
		self.SNVs.append(SNV)

	def to_dict(self):
		dic = dict()
		dic['chromosome'] = str(self.chromosome)
		dic['start'] = self.start.position
		dic['end'] = self.end.position
		dic['minor_cn'] = self.minor_cn
		dic['major_cn'] = self.major_cn
		dic['snvs'] = [snv.to_dict() for snv in self.SNVs]
		if len(self.genes) >0:
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


@class_equality_attributes
class SNV():
	def __init__(self,chromosome,pos,ref_count,alt_count,ref_base,alt_base,config):
		self.chromosome = Chromosome(chromosome)
		self.pos = Position(Chromosome(chromosome),pos,config)
		self.ref_count = ref_count
		self.alt_count = alt_count
		self.ref_base = ref_base
		self.alt_base = str(alt_base)

	def to_dict(self):
		dic = dict()
		dic['chromosome'] = str(self.chromosome)
		dic['pos'] = self.pos.position
		dic['ref_count'] = self.ref_count
		dic['alt_count'] = self.alt_count
		dic['ref_base'] = self.ref_base
		dic['alt_base'] = self.alt_base
		return dic

	def __str__(self):
		return str(self.pos) + " alt counts: "+str(self.alt_count)+ " ref counts: "+str(self.ref_count) + " "+self.ref_base+"/"+str(self.alt_base)

@class_equality_attributes
class Amplification():
	def __init__(self,chromosome,arm,segments):
		self.chromosome = chromosome
		self.arm = arm
		self.segments = segments
		self.infered_time = None
	
	def add_segment(self,segment,match_genes = True):
		if match_genes: segment.match_genes()
		self.segments.append(segment)

	def get_genes(self):
		genes = []
		for segment in self.segments:
			genes += segment.genes
		return genes

	def to_dict(self):
		dic = dict()
		dic['chromosome'] = str(self.chromosome)
		dic['arm'] = self.arm
		dic['segments'] = [segment.to_dict() for segment in self.segments]
		return dic

@class_equality_attributes
class Gene():
	def __init__(self,gene,config):
		self.config = config
		if isinstance(gene, dict):
			self.name = gene['name']
			self.ensembl_id = gene['ensembl_id']
			self.biotype = gene['biotype']
			self.chromosome = Chromosome(gene['chromosome'])
			self.start= Position(self.chromosome,gene['start'],config)
			self.end= Position(self.chromosome,gene['end'],config)
			self.strand = gene['strand']
			self.entrez_id = gene['entrez_id']
		else:
			self.name = gene.gene_name
			self.ensembl_id = gene.gene_id
			self.biotype = gene.biotype
			self.chromosome = Chromosome(gene.contig)
			self.start= Position(self.chromosome,gene.start,config)
			self.end= Position(self.chromosome,gene.end,config)
			self.strand = gene.strand
			entrez = config['ensembl_to_entrez']['NCBI gene (formerly Entrezgene) ID'][config['ensembl_to_entrez']['Gene stable ID'] == gene.gene_id]
			if len(entrez) == 0:
				#if Ensembl not in table
				self.entrez_id = float('nan')
			elif len(entrez) == 1:
				if not math.isnan(float(entrez)):
					self.entrez_id = [int(entrez)]
				else:
					#if no entrez entry for this Ensembl id in table
					self.entrez_id = float('nan')
			else:
				self.entrez_id = []
				for i in entrez.values.tolist():
					if not math.isnan(i):
						self.entrez_id.append(int(i)) 
			

	def to_dict(self):
		dic = dict()
		dic['name'] = self.name
		dic['ensembl_id'] = self.ensembl_id
		dic['biotype'] = self.biotype
		dic['chromosome'] = self.chromosome.chromosome
		dic['start'] = self.start.position
		dic['end'] = self.end.position
		dic['strand'] = self.strand
		dic['entrez_id'] = self.entrez_id 
		return dic
