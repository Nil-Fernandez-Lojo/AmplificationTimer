from .decorator_equality import class_equality_attributes
from .chromosome import Chromosome
from .position import Position
import math


@class_equality_attributes
class Gene:
    def __init__(self, gene, config):
        self.config = config
        if isinstance(gene, dict):
            self.name = gene['name']
            self.ensembl_id = gene['ensembl_id']
            self.biotype = gene['biotype']
            self.chromosome = Chromosome(gene['chromosome'])
            self.start = Position(self.chromosome, gene['start'], config)
            self.end = Position(self.chromosome, gene['end'], config)
            self.strand = gene['strand']
            self.entrez_id = gene['entrez_id']
        else:
            self.name = gene.gene_name
            self.ensembl_id = gene.gene_id
            self.biotype = gene.biotype
            self.chromosome = Chromosome(gene.contig)
            self.start = Position(self.chromosome, gene.start, config)
            self.end = Position(self.chromosome, gene.end, config)
            self.strand = gene.strand
            entrez = config['ensembl_to_entrez']['NCBI gene (formerly Entrezgene) ID'][
                config['ensembl_to_entrez']['Gene stable ID'] == gene.gene_id]
            if len(entrez) == 0:
                # if Ensembl not in table
                self.entrez_id = float('nan')
            elif len(entrez) == 1:
                if not math.isnan(float(entrez)):
                    self.entrez_id = [int(entrez)]
                else:
                    # if no entrez entry for this Ensembl id in table
                    self.entrez_id = float('nan')
            else:
                self.entrez_id = []
                for i in entrez.values.tolist():
                    if not math.isnan(i):
                        self.entrez_id.append(int(i))

    def check_if_oncogene(self):
        if isinstance(self.entrez_id, list):
            list_entrez_id = self.entrez_id
        else:
            list_entrez_id = [self.entrez_id]
        for entrez_id in list_entrez_id:
            if entrez_id in self.config['oncogenes']["Entrez GeneId"].values:
                return True
        return False

    def check_if_driver_gain(self, cancer_type):
        return self.name in self.config['drivers'][cancer_type]['gain'].keys()

    def check_if_driver_rest(self, cancer_type):
        return self.name in self.config['drivers'][cancer_type]['rest'].keys()

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

def find_most_common_oncogenes(oncogene_type,genes,config,cancer_type=None):
    # only filters oncogenes for known drivers of the cancer type, not for the plain list of oncogenes
    if oncogene_type == 'oncogene':
        return genes
    else:
        type_mutation = 'gain' if oncogene_type == 'driver_gain' else 'rest'
        dict_counts = config['drivers'][cancer_type][type_mutation]
        count = 0
        most_common_oncogenes = []
        print(oncogene_type)
        print([g.name for g in genes])
        print(dict_counts)
        for g in genes:
            if dict_counts[g.name] == count:
                most_common_oncogenes.append(g)
            if dict_counts[g.name]>count:
                count = dict_counts[g.name]
                most_common_oncogenes = [g]
    return most_common_oncogenes