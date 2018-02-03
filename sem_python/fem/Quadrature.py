from math import sqrt

from ..input.Parameters import Parameters


class QuadratureParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)


class Quadrature(object):

    def __init__(self, params):
        self.n_q = 2
        self.z = [-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)]
        self.w = [1.0, 1.0]
        self.Jac_divided_by_h = 0.5
        self.JxW_divided_by_h = [self.w[q] * self.Jac_divided_by_h for q in range(self.n_q)]
