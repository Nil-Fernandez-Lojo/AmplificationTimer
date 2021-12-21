from functools import total_ordering
from decorator_equality import class_equality_attributes

@total_ordering
@class_equality_attributes
class Chromosome:
	"""
	class that represents a chromosome, as input the chromosome can be given as an interger or a float. 
	Input must be in [1-23, 'X', 'Y']
	"""
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
