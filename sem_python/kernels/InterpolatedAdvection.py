from .Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters


class InterpolatedAdvectionParameters(Kernel1PhaseParameters):

    def __init__(self, factory):
        Kernel1PhaseParameters.__init__(self, factory)
        self.registerStringParameter("aux_name", "Name of the aux that this kernel operates upon")


class InterpolatedAdvection(Kernel1Phase):

    def __init__(self, params):
        Kernel1Phase.__init__(self, params)
        self.aux_name = params.get("aux_name")

    def computeResidual(self, data, i):
        return -data[self.aux_name] * data["grad_phi"][i] * data["JxW"]

    def computeJacobian(self, data, der, var_index, i, j):
        if var_index == self.aA1_index:
            return -der[self.aux_name]["aA1"] * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
        elif var_index == self.arhoA_index:
            return -der[self.aux_name][self.arhoA] * data["phi"][j] * data["grad_phi"][i] * data[
                "JxW"]
        elif var_index == self.arhouA_index:
            return -der[self.aux_name][self.arhouA] * data["phi"][j] * data["grad_phi"][i] * data[
                "JxW"]
        elif var_index == self.arhoEA_index:
            return -der[self.aux_name][self.arhoEA] * data["phi"][j] * data["grad_phi"][i] * data[
                "JxW"]
        else:
            return self.zero
