from Junction import Junction, JunctionParameters
from enums import ModelType, VariableName

class Junction1PhaseParameters(JunctionParameters):
  def __init__(self):
    JunctionParameters.__init__(self)
    self.registerIntParameter("phase", "Index of phase to which BC is applied")

## Base class for 1-phase junctions
class Junction1Phase(Junction):
  def __init__(self, params):
    Junction.__init__(self, params)
    self.phase = params.get("phase")
    self.eos = self.eos_list[self.phase]

    # get variable indices
    if (self.model_type == ModelType.TwoPhase):
      self.vf1_index = self.dof_handler.variable_index[VariableName.VF1][0]
    self.arho_index = self.dof_handler.variable_index[VariableName.ARhoA][self.phase]
    self.arhouA_index = self.dof_handler.variable_index[VariableName.ARhoUA][self.phase]
    self.arhoEA_index = self.dof_handler.variable_index[VariableName.ARhoEA][self.phase]

  def setDoFIndices(self):
    # indices for all meshes
    self.i_arhoA = [self.dof_handler.i(k, self.arho_index) for k in self.node_indices]
    self.i_arhouA = [self.dof_handler.i(k, self.arhouA_index) for k in self.node_indices]
    self.i_arhoEA = [self.dof_handler.i(k, self.arhoEA_index) for k in self.node_indices]
