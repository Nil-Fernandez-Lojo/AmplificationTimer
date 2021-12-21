from decorator_equality import class_equality_attributes

@class_equality_attributes
class SNV():
	def __init__(self,chromosome,pos,ref_count,alt_count,ref_base,alt_base,config):
		self.chromosome = chromosome
		self.pos = pos
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
