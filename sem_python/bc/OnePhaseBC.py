from ..base.enums import ModelType, VariableName
from .BC import BC, BCParameters

class OnePhaseBCParameters(BCParameters):
  def __init__(self):
    BCParameters.__init__(self)
    self.registerIntParameter("phase", "Index of phase to which BC is applied")

class OnePhaseBC(BC):
  def __init__(self, params):
    BC.__init__(self, params)
    self.phase = params.get("phase")
    self.eos = self.eos_list[self.phase]

    # DoF indices for the conserved variables
    if (self.model_type == ModelType.TwoPhase):
      aA1_index = self.dof_handler.variable_index[VariableName.AA1][0]
      self.i_aA1 = self.dof_handler.i(self.k, aA1_index)

    arhoA_index = self.dof_handler.variable_index[VariableName.ARhoA][self.phase]
    arhouA_index = self.dof_handler.variable_index[VariableName.ARhoUA][self.phase]
    arhoEA_index = self.dof_handler.variable_index[VariableName.ARhoEA][self.phase]

    self.i_arhoA = self.dof_handler.i(self.k, arhoA_index)
    self.i_arhouA = self.dof_handler.i(self.k, arhouA_index)
    self.i_arhoEA = self.dof_handler.i(self.k, arhoEA_index)
