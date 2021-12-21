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
	def __init__(self,chromosome,arm,segments):
		self.chromosome = chromosome
		self.arm = arm
		self.segments = segments
		self.infered_time = None
	
	def add_segment(self,segment):
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