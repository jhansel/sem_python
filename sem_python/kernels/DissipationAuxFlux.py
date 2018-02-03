from .Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters


class DissipationAuxFluxParameters(Kernel1PhaseParameters):

    def __init__(self, factory):
        Kernel1PhaseParameters.__init__(self, factory)
        self.registerStringParameter("flux_name", "Name of flux aux")


## Dissipation from an aux quantity flux.
#
# This dissipation is of the following form on the LHS:
# \f[
#   -\nabla\cdot \mathbf{f} ,
# \f]
# so its weak form is
# \f[
#   (-\nabla\cdot \mathbf{f}, \phi_i)_\Omega .
# \f]
# After integration by parts, this becomes
# \f[
#   (\mathbf{f}, \nabla\phi_i)_\Omega - (\mathbf{f}\cdot\mathbf{n}\phi_i)_{\partial\Omega} .
# \f]
class DissipationAuxFlux(Kernel1Phase):

    def __init__(self, params):
        Kernel1Phase.__init__(self, params)
        self.flux = params.get("flux_name")

    def computeResidual(self, data, i):
        return data[self.flux] * data["grad_phi"][i] * data["JxW"]

    def computeJacobian(self, data, der, var_index, i, j):
        if var_index == self.aA1_index:
            return (
                der[self.flux]["aA1"] * data["phi"][j] + der[self.flux]["grad_aA1"] *
                data["grad_phi"][j]) * data["grad_phi"][i] * data["JxW"]
        elif var_index == self.arhoA_index:
            return (
                der[self.flux][self.arhoA] * data["phi"][j] + der[self.flux][self.grad_arhoA] *
                data["grad_phi"][j]) * data["grad_phi"][i] * data["JxW"]
        elif var_index == self.arhouA_index:
            return (
                der[self.flux][self.arhouA] * data["phi"][j] + der[self.flux][self.grad_arhouA] *
                data["grad_phi"][j]) * data["grad_phi"][i] * data["JxW"]
        elif var_index == self.arhoEA_index:
            return (
                der[self.flux][self.arhoEA] * data["phi"][j] + der[self.flux][self.grad_arhoEA] *
                data["grad_phi"][j]) * data["grad_phi"][i] * data["JxW"]
        else:
            return self.zero
