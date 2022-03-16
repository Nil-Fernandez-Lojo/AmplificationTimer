import pandas as pd
from .decorator_equality import class_equality_attributes
from .chromosome import Chromosome
from .mutation_rate import MutationRate
from .segment import Segment
from .plot import plot_cn


@class_equality_attributes
class Amplification:
    """
    Class that encodes an amplification

    """

    def __init__(self,
                 chromosome=None,
                 arm=None,
                 segments=None,
                 threshold_amplification=None,
                 clinical_data=None,
                 config=None,
                 mutation_rate=None,
                 subclonal_structure=None,
                 amplification_dict=None):
        if amplification_dict is not None:
            self.clinical_data = amplification_dict['clinical_data']
            self.chromosome = Chromosome(amplification_dict['chromosome'])
            self.arm = amplification_dict['arm']
            self.config = config
            self.segments = [Segment(segment, self.config, self.clinical_data) for segment in
                             amplification_dict['segments']]
            self.threshold_amplification = amplification_dict['threshold_amplification']
            self.mutation_rate = MutationRate(amplification_dict['mutation_rate']['n_1_0'],
                                              amplification_dict['mutation_rate']['l_1_0'],
                                              amplification_dict['mutation_rate']['n_1_1'],
                                              amplification_dict['mutation_rate']['l_1_1'],
                                              amplification_dict['mutation_rate']['n_1_0_ctpg'],
                                              amplification_dict['mutation_rate']['n_1_1_ctpg'])
            self.subclonal_structure = pd.DataFrame.from_dict(amplification_dict['subclonal_structure'])
        else:
            self.chromosome = chromosome
            self.arm = arm
            self.segments = segments
            self.config = config
            self.threshold_amplification = threshold_amplification
            self.clinical_data = clinical_data
            self.mutation_rate = mutation_rate
            self.subclonal_structure = subclonal_structure
        # TODO should not be recomputed each time
        self.oncogenes = []
        self.potential_driver_oncogenes = []
        self.set_oncogenes()

    def add_segment(self, segment):
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
        self.potential_driver_oncogenes = []
        genes = self.get_genes()
        cancer_type = self.config['map_cancer_type_summary_to_driver_table'].get(self.clinical_data['histology_pcawg'],
                                                                                 self.clinical_data['histology_pcawg'])
        for gene in genes:
            if isinstance(gene.entrez_id, list):
                list_entrez_id = gene.entrez_id
            else:
                list_entrez_id = [gene.entrez_id]
            for entrez_id in list_entrez_id:
                if entrez_id in self.config['oncogenes']["Entrez GeneId"].values and (gene not in self.oncogenes):
                    self.oncogenes.append(gene)
            if ((gene.name in self.config['drivers'][cancer_type]['gain'])
                    and (gene not in self.potential_driver_oncogenes)):
                self.potential_driver_oncogenes.append(gene)

    def get_mean_ploidy(self):
        tot_len = 0
        number_bases = 0
        for segment in self.segments:
            tot_cn = segment.major_cn + segment.minor_cn
            len_segment = segment.get_length()
            number_bases += tot_cn * len_segment
            tot_len += len_segment
        return number_bases / tot_len

    def get_length(self):
        length_amplification = 0
        for segment in self.segments:
            length_amplification += segment.get_length()
        return length_amplification

    def plot(self,
             add_snvs=True,
             rel_margin=0.1,
             total_cn=False,
             plot_threshold_amp=True,
             title=None,
             path_save=None,
             differentiate_snv_type=True,
             ax=None):

        if plot_threshold_amp:
            threshold_amplification = self.threshold_amplification
        else:
            threshold_amplification = None

        margin = rel_margin * (self.segments[-1].end.position - self.segments[0].start.position)
        start_pos = self.segments[0].start.position - margin
        max_bases = (self.segments[-1].end.position - self.segments[0].start.position) + 2 * margin

        plot_cn(self.segments,
                self.config,
                self.clinical_data['purity'],
                add_snvs=add_snvs,
                chromosome=self.chromosome,
                start_pos=start_pos,
                max_bases=max_bases,
                total_cn=total_cn,
                threshold_amplification=threshold_amplification,
                title=title,
                path_save=path_save,
                differentiate_snv_type=differentiate_snv_type,
                ax=ax)

    def get_name(self):
        return self.clinical_data['samplename']+str(self.chromosome)+self.arm

    def to_dict(self):
        dic = dict()
        dic['chromosome'] = str(self.chromosome)
        dic['arm'] = self.arm
        dic['segments'] = [segment.to_dict() for segment in self.segments]
        dic['oncogenes'] = [gene.to_dict() for gene in self.oncogenes]
        dic['potential_driver_oncogenes'] = [gene.to_dict() for gene in self.potential_driver_oncogenes]
        dic['threshold_amplification'] = self.threshold_amplification
        dic['clinical_data'] = self.clinical_data
        dic['mutation_rate'] = self.mutation_rate.to_dict()
        dic['subclonal_structure'] = self.subclonal_structure.to_dict('records')
        # TODO: not good practice, should change this:
        # self.config is not added since it is loaded with load_config
        return dic

