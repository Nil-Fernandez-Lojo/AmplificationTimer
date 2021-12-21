from functools import total_ordering
from decorator_equality import class_equality_attributes

@total_ordering
@class_equality_attributes
class Position:
	"""
	Class that represents a genomic position (we use hg19).
	
	...

	Attributes
	----------
	chromosome: Chromosome
	position: int
	arm: str

	Methods
	-------
	init:
		takes as input:
			chromosome: Chromosome
			position: int
			config: dic

	"""
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
