import copy

from decorator_equality import class_equality_attributes

@class_equality_attributes
class Amplification():
	"""
	Class that encodes an amplification

	...

	Attributes
	----------
	chromosome: Chromosome
		Chromosome in which the amplification took place
	arm: str
		chromosome arm ('p' or 'q')
	segments: list of Segments
		list of the segments amplified
	infered_time:


	Methods:
	init: takes as input chromosome, arm, segments
	add_segment: adds a segment to self.segments
	get_genes: returns the genes that are encoded in teh segments
	to_dict: returns a dictionary encoding the amplification
	"""
	def __init__(self,chromosome,arm,segments,threshold_amplification,clinical_data,config):
		self.chromosome = chromosome
		self.arm = arm
		self.segments = segments
		self.config = config
		self.threshold_amplification = threshold_amplification
		self.clinical_data = clinical_data
		self.config = config
		self.oncogenes = []

	def add_segment(self,segment):
		self.segments.append(segment)

	def get_genes(self):
		genes = []
		for segment in self.segments:
			for gene in segment.genes:
				if gene not in genes:
					genes.append(gene)
			genes += segment.genes
		return genes

	def set_oncogenes(self):
		self.oncogenes = []
		genes = self.get_genes()
		for gene in genes:
			if isinstance(gene.entrez_id,list):
				list_entrez_id = gene.entrez_id
			else:
				list_entrez_id = [gene.entrez_id]
			for entrez_id in list_entrez_id:
				if entrez_id in self.config['oncogenes']["Entrez GeneId"].values and (gene not in self.oncogenes):
					self.oncogenes.append(gene)				
		return copy.deepcopy(self.oncogenes)

	def get_mean_ploidy(self):
		tot_len = 0
		number_bases = 0
		for segment in self.segments:
			tot_cn = segment.major_cn + segment.minor_cn
			len_segment = segment.get_length()
			number_bases += tot_cn*len_segment
			tot_len += len_segment
		return number_bases/tot_len

	def to_dict(self):
		dic = dict()
		dic['chromosome'] = str(self.chromosome)
		dic['arm'] = self.arm
		dic['segments'] = [segment.to_dict() for segment in self.segments]
		dic['oncogenes'] = [gene.to_dict() for gene in self.oncogenes]
		dic['threshold_amplification'] = self.threshold_amplification
		dic['clinical_data'] = self.clinical_data
		# TODO: not good practice, should change this:
		# self.config is not added since it contains objects of type PATH that can't be saved in json file later
		return dic