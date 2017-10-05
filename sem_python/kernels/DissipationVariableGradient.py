from .Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters

class DissipationVariableGradientParameters(Kernel1PhaseParameters):
  def __init__(self):
    Kernel1PhaseParameters.__init__(self)

## Dissipation from a solution variable gradient.
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
class DissipationVariableGradient(Kernel1Phase):
  def __init__(self, params):
    Kernel1Phase.__init__(self, params)
    var = self.dof_handler.variableEnumToName(self.var_enum, self.phase)
    self.grad_var = "grad_" + var
    self.visccoef = "visccoef_" + var

  def computeResidual(self, data, i):
    return data[self.visccoef] * data[self.grad_var] * data["grad_phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.var_index:
      dgrad_vari_dvarj = data["grad_phi"][j]
    else:
      dgrad_vari_dvarj = self.zero

    if var_index == self.aA1_index:
      return (der[self.visccoef]["aA1"] * data["phi"][j] * data[self.grad_var] + data[self.visccoef] * dgrad_vari_dvarj) * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arhoA_index:
      return (der[self.visccoef][self.arhoA] * data["phi"][j] * data[self.grad_var] + data[self.visccoef] * dgrad_vari_dvarj) * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arhouA_index:
      return (der[self.visccoef][self.arhouA] * data["phi"][j] * data[self.grad_var] + data[self.visccoef] * dgrad_vari_dvarj) * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arhoEA_index:
      return (der[self.visccoef][self.arhoEA] * data["phi"][j] * data[self.grad_var] + data[self.visccoef] * dgrad_vari_dvarj) * data["grad_phi"][i] * data["JxW"]
    else:
      return self.zero
