from math import sqrt

from ..input.Parameters import Parameters
from ..utilities.error_utilities import error


class QuadratureParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerIntParameter("n_q_points", "Number of quadrature points per cell", 2)


class Quadrature(object):

    def __init__(self, params):
        self.n_q = params.get("n_q_points")
        if self.n_q != 2:
            error("Only 2-point quadrature has been implemented")
        self.z = [-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)]
        self.w = [1.0, 1.0]
        self.Jac_divided_by_h = 0.5
        self.JxW_divided_by_h = [self.w[q] * self.Jac_divided_by_h for q in range(self.n_q)]
