from enums import VariableName
from BC import BC, BCParameters

class VolumeFractionBCParameters(BCParameters):
  def __init__(self):
    BCParameters.__init__(self)

class VolumeFractionBC(BC):
  def __init__(self, params):
    BC.__init__(self, params)

    vf1_index =self.dof_handler.variable_index[VariableName.VF1][0]
    self.i_vf1 = self.dof_handler.i(self.k, vf1_index)
