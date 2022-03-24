from .decorator_equality import class_equality_attributes


@class_equality_attributes
class MutationRate:
    """
    class that represents the mutation rate

    ...

    Attributes
    n_1_0: int
        number of mutations on 1+0 segments
    n_1_1: int
        number of mutations on 1+1 segments
    l_1_0: int
        total length of 1+0 segments
    l_1_1: int
        total length of 1+1 segments

    Methods:
    init: takes the 4 arguments in input in order
        n_1_0,l_1_0,n_1_1,l_1_1
    get_ml: gives the ML estimation of the mutation rate
    get_beta_posterior_parameters: gives the parameters of the posterior
    beta distribution if a 1,1 beta prior is used
    to_dict: returns a dictionary with the 4 attributes

    """

    def __init__(self, n_1_0, l_1_0, n_1_1, l_1_1, n_1_0_ctpg, n_1_1_ctpg):
        self.n_1_0 = n_1_0
        self.n_1_0_ctpg = n_1_0_ctpg
        self.l_1_0 = l_1_0
        self.n_1_1 = n_1_1
        self.n_1_1_ctpg = n_1_1_ctpg
        self.l_1_1 = l_1_1
        if (self.l_1_0 == 0) and (self.l_1_1 == 0):
            print("Warning, no segments of with copy number 1+1 and 1+0")

    def get_ml(self, only_clock_like):
        # TODO: change this, REALLY IMPORTANT otherwise do not trust it
        # default value if no segment of of copy number 1+1 and 1+0, is 10^-6
        if only_clock_like:
            n_1_0 = self.n_1_0_ctpg
            n_1_1 = self.n_1_1_ctpg
        else:
            n_1_0 = self.n_1_0
            n_1_1 = self.n_1_1

        if self.l_1_0 + 2 * self.l_1_1 == 0:
            return 10 ** (-6)
        elif n_1_0 + 2 * n_1_1 == 0:
            return 10 ** (-8)
        else:
            return (n_1_0 + n_1_1) / (self.l_1_0 + 2 * self.l_1_1)

    def get_beta_posterior_parameters(self, only_clock_like=False):
        if only_clock_like:
            n_1_0 = self.n_1_0_ctpg
            n_1_1 = self.n_1_1_ctpg
        else:
            n_1_0 = self.n_1_0
            n_1_1 = self.n_1_1
        return n_1_0 + n_1_1 + 1, self.l_1_0 + 2 * self.l_1_1 - (n_1_0 + n_1_1 + 1) + 1

    def to_dict(self):
        dic = dict()
        dic['n_1_0'] = self.n_1_0
        dic['l_1_0'] = self.l_1_0
        dic['n_1_1'] = self.n_1_1
        dic['l_1_1'] = self.l_1_1
        dic['n_1_0_ctpg'] = self.n_1_0_ctpg
        dic['n_1_1_ctpg'] = self.n_1_1_ctpg
        return dic
