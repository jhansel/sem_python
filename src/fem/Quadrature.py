from math import sqrt

class Quadrature(object):
  def __init__(self):
    self.n_q = 2
    self.z = [-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)]
    self.w = [1.0, 1.0]
    self.Jac_divided_by_h = 0.5
    self.JxW_divided_by_h = [self.w[q] * self.Jac_divided_by_h for q in xrange(self.n_q)]
