import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

class BCParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerStringSelectionParameter("boundary", ["left", "right"], "Which boundary to apply boundary condition on")

class BC(object):
  def __init__(self, params, dof_handler, eos):
    self.boundary = params.get("boundary")
    self.dof_handler = dof_handler
    self.model_type = dof_handler.model_type
    self.eos = eos
    if self.boundary == "left":
      self.nx = -1.0
      self.k = 0
    else:
      self.nx = 1.0
      self.k = self.dof_handler.n_cell

  def applyStrongBCNonlinearSystem(self, U, r, J):
    pass

  def applyStrongBCLinearSystemMatrix(self, A):
    pass

  def applyStrongBCLinearSystemRHSVector(self, U_old, b):
    pass
