from pyensembl import EnsemblRelease
ensembl = EnsemblRelease(63)
data = dict()
data['chromosome'] = 9
data['start'] = 32551000
data['end'] = 32566900
genes = ensembl.genes_at_locus(contig=data['chromosome'], position=int(data['start']),end=int(data['end']))
for gene in genes:
	print(gene.gene_id, gene.gene_name)
print()
print(ensembl.gene_by_id('ENSG00000235453'))