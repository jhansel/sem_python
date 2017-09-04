from .Junction import Junction, JunctionParameters
from ..base.enums import ModelType, VariableName

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
    self.arho_index = self.dof_handler.variable_index[VariableName.ARho][self.phase]
    self.arhou_index = self.dof_handler.variable_index[VariableName.ARhoU][self.phase]
    self.arhoE_index = self.dof_handler.variable_index[VariableName.ARhoE][self.phase]

  def setDoFIndices(self):
    # indices for all meshes
    self.i_arho = [self.dof_handler.i(k, self.arho_index) for k in self.node_indices]
    self.i_arhou = [self.dof_handler.i(k, self.arhou_index) for k in self.node_indices]
    self.i_arhoE = [self.dof_handler.i(k, self.arhoE_index) for k in self.node_indices]
