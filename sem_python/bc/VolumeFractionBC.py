from ..base.enums import VariableName
from .BC import BC, BCParameters


class VolumeFractionBCParameters(BCParameters):

    def __init__(self):
        BCParameters.__init__(self)


class VolumeFractionBC(BC):

    def __init__(self, params):
        BC.__init__(self, params)

        aA1_index = self.dof_handler.variable_index[VariableName.AA1][0]
        self.i_aA1 = self.dof_handler.i(self.k, aA1_index)
