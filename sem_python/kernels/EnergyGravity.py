from .Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters
from ..base.enums import VariableName


class EnergyGravityParameters(Kernel1PhaseParameters):

    def __init__(self):
        Kernel1PhaseParameters.__init__(self)


class EnergyGravity(Kernel1Phase):

    def __init__(self, params):
        params.set("var_enum", VariableName.ARhoEA)
        Kernel1Phase.__init__(self, params)

    def computeResidual(self, data, i):
        return -data[self.arhouA] * data["g"] * data["phi"][i] * data["JxW"]

    def computeJacobian(self, data, der, var_index, i, j):
        if var_index == self.arhouA_index:
            return -data["g"] * data["phi"][j] * data["phi"][i] * data["JxW"]
        else:
            return self.zero
