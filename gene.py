from decorator_equality import class_equality_attributes
from import Chromosome
from import Position


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