from .Kernel2Phase import Kernel2Phase, Kernel2PhaseParameters
from ..base.enums import VariableName


class VolumeFractionAdvectionParameters(Kernel2PhaseParameters):

    def __init__(self):
        Kernel2PhaseParameters.__init__(self)


class VolumeFractionAdvection(Kernel2Phase):

    def __init__(self, params):
        params.set("var_enum", VariableName.AA1)
        Kernel2Phase.__init__(self, params)

    def computeResidual(self, data, i):
        return data["uI"] * data["grad_vf1"] * data["A"] * data["phi"][i] * data["JxW"]

    def computeJacobian(self, data, der, var_index, i, j):
        if var_index == self.aA1_index:
            return (
                der["uI"]["aA1"] * data["grad_vf1"] * data["phi"][j] + data["uI"] * der["grad_vf1"]
                ["grad_aA1"] * data["grad_phi"][j]) * data["A"] * data["phi"][i] * data["JxW"]
        elif var_index == self.arhoA1_index:
            return der["uI"]["arhoA1"] * data["grad_vf1"] * data["A"] * data["phi"][j] * data[
                "phi"][i] * data["JxW"]
        elif var_index == self.arhouA1_index:
            return der["uI"]["arhouA1"] * data["grad_vf1"] * data["A"] * data["phi"][j] * data[
                "phi"][i] * data["JxW"]
        elif var_index == self.arhoEA1_index:
            return der["uI"]["arhoEA1"] * data["grad_vf1"] * data["A"] * data["phi"][j] * data[
                "phi"][i] * data["JxW"]
        elif var_index == self.arhoA2_index:
            return der["uI"]["arhoA2"] * data["grad_vf1"] * data["A"] * data["phi"][j] * data[
                "phi"][i] * data["JxW"]
        elif var_index == self.arhouA2_index:
            return der["uI"]["arhouA2"] * data["grad_vf1"] * data["A"] * data["phi"][j] * data[
                "phi"][i] * data["JxW"]
        elif var_index == self.arhoEA2_index:
            return der["uI"]["arhoEA2"] * data["grad_vf1"] * data["A"] * data["phi"][j] * data[
                "phi"][i] * data["JxW"]
        else:
            return self.zero
