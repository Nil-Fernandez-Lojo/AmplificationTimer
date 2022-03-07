from .decorator_equality import class_equality_attributes

@class_equality_attributes
class SNV:
    def __init__(self,
                 chromosome,
                 pos,
                 ref_count,
                 alt_count,
                 ref_base,
                 alt_base,
                 config,
                 previous_ref_base=None,
                 next_ref_base=None):
        self.chromosome = chromosome
        self.pos = pos
        self.ref_count = ref_count
        self.alt_count = alt_count
        self.ref_base = ref_base
        self.alt_base = str(alt_base)
        if previous_ref_base is None:
            self.previous_ref_base = config['genome']['chr' + str(chromosome)][pos.position - 2].capitalize()
        else:
            self.previous_ref_base = previous_ref_base
        if next_ref_base is None:
            self.next_ref_base = config['genome']['chr' + str(chromosome)][pos.position].capitalize()
        else:
            self.next_ref_base = next_ref_base

    def to_dict(self):
        dic = dict()
        dic['chromosome'] = str(self.chromosome)
        dic['pos'] = self.pos.position
        dic['ref_count'] = self.ref_count
        dic['alt_count'] = self.alt_count
        dic['ref_base'] = self.ref_base
        dic['alt_base'] = self.alt_base
        dic['previous_ref_base'] = self.previous_ref_base
        dic['next_ref_base'] = self.next_ref_base
        return dic

    def clock_like(self):
        forward_strand_clock = (self.ref_base == 'C') and (self.alt_base == 'T') and (self.next_ref_base == 'G')
        backward_strand_clock = (self.ref_base == 'G') and (self.alt_base == 'A') and (self.previous_ref_base == 'C')
        return forward_strand_clock or backward_strand_clock

    def __str__(self):
        return str(self.pos) + " alt counts: " + str(self.alt_count) + " ref counts: " + str(
            self.ref_count) + " " + self.ref_base + "/" + str(self.alt_base)
