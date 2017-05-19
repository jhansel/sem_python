import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType, VariableName

sys.path.append(base_dir + "src/bc")
from BC import BC, BCParameters

class OnePhaseBCParameters(BCParameters):
  def __init__(self):
    BCParameters.__init__(self)
    self.registerIntParameter("phase", "Index of phase to which BC is applied")

class OnePhaseBC(BC):
  def __init__(self, params, dof_handler, eos):
    BC.__init__(self, params, dof_handler, eos)
    self.phase = params.get("phase")
    self.eos = eos[self.phase]

    # DoF indices for the conserved variables
    if (self.model_type == ModelType.TwoPhase):
      self.i_vf1 = self.dof_handler.i(self.k, VariableName.VF1)
    self.i_arho = self.dof_handler.i(self.k, VariableName.ARho, self.phase)
    self.i_arhou = self.dof_handler.i(self.k, VariableName.ARhoU, self.phase)
    self.i_arhoE = self.dof_handler.i(self.k, VariableName.ARhoE, self.phase)
