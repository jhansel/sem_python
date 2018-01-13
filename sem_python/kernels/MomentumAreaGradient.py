from .Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters
from ..base.enums import VariableName


class MomentumAreaGradientParameters(Kernel1PhaseParameters):

    def __init__(self):
        Kernel1PhaseParameters.__init__(self)


class MomentumAreaGradient(Kernel1Phase):

    def __init__(self, params):
        params.set("var_enum", VariableName.ARhoUA)
        Kernel1Phase.__init__(self, params)

    def computeResidual(self, data, i):
        return -data[self.p] * data[self.vf] * data["grad_A"] * data["phi"][i] * data["JxW"]

    def computeJacobian(self, data, der, var_index, i, j):
        if var_index == self.aA1_index:
            return -(der[self.p]["aA1"] * data[self.vf] + data[self.p] * der[self.vf]["aA1"]) \
              * data["grad_A"] * data["phi"][j] * data["phi"][i] * data["JxW"]
        else:
            aux = -data[self.vf] * data["grad_A"] * data["phi"][j] * data["phi"][i] * data["JxW"]
            if var_index == self.arhoA_index:
                return der[self.p][self.arhoA] * aux
            elif var_index == self.arhouA_index:
                return der[self.p][self.arhouA] * aux
            elif var_index == self.arhoEA_index:
                return der[self.p][self.arhoEA] * aux
            else:
                return self.zero
