from .Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters
from ..base.enums import VariableName


class EnergyAdvectionParameters(Kernel1PhaseParameters):

    def __init__(self):
        Kernel1PhaseParameters.__init__(self)


class EnergyAdvection(Kernel1Phase):

    def __init__(self, params):
        params.set("var_enum", VariableName.ARhoEA)
        Kernel1Phase.__init__(self, params)

    def computeResidual(self, data, i):
        return -data[self.u] * (data[self.arhoEA] + data[self.vf] * data[self.p] * data["A"]
                                ) * data["grad_phi"][i] * data["JxW"]

    def computeJacobian(self, data, der, var_index, i, j):
        if var_index == self.aA1_index:
            return -data[self.u] * (data[self.vf] * der[self.p]["aA1"] \
              + der[self.vf]["aA1"] * data[self.p]) * data["A"] * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
        elif var_index == self.arhoA_index:
            return -(data[self.u] * (data[self.vf] * der[self.p][self.arhoA] * data["A"]) \
              + der[self.u][self.arhoA] * (data[self.arhoEA] + data[self.vf] * data[self.p] * data["A"])) * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
        elif var_index == self.arhouA_index:
            return -(data[self.u] * (data[self.vf] * der[self.p][self.arhouA] * data["A"]) \
              + der[self.u][self.arhouA] * (data[self.arhoEA] + data[self.vf] * data[self.p] * data["A"])) * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
        elif var_index == self.arhoEA_index:
            return -data[self.u] * (1 + data[self.vf] * der[self.p][self.arhoEA] * data["A"]
                                    ) * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
        else:
            return self.zero
