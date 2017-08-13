from abc import ABCMeta, abstractmethod
import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

class KernelParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerIntParameter("var_index", "Index of variable")
    self.registerParameter("dof_handler", "Degree of freedom handler")

class Kernel(object):
  __metaclass__ = ABCMeta

  def __init__(self, params):
    self.var_index = params.get("var_index")
    self.dof_handler = params.get("dof_handler")
    self.i = self.dof_handler.i
    self.n = self.dof_handler.n_dof_per_cell_per_var
    self.zero = [0]

  def apply(self, data, der, r, J):
    for i_local in xrange(self.n):
      i = self.i(i_local, self.var_index)
      r[i] += sum(self.computeResidual(data, i_local))
      for j_local in xrange(self.n):
        for var_index in self.var_indices:
          j = self.i(j_local, var_index)
          J[i,j] += sum(self.computeJacobian(data, der, var_index, i_local, j_local))

  @abstractmethod
  def computeResidual(self, data, i):
    pass

  @abstractmethod
  def computeJacobian(self, data, der, var, i, j):
    pass
