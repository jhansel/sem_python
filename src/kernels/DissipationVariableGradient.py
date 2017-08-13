import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters

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

    if var_index == self.vf1_index:
      return (der[self.visccoef]["vf1"] * data["phi"][j] * data[self.grad_var] + data[self.visccoef] * dgrad_vari_dvarj) * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arho_index:
      return (der[self.visccoef][self.arho] * data["phi"][j] * data[self.grad_var] + data[self.visccoef] * dgrad_vari_dvarj) * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arhou_index:
      return (der[self.visccoef][self.arhou] * data["phi"][j] * data[self.grad_var] + data[self.visccoef] * dgrad_vari_dvarj) * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arhoE_index:
      return (der[self.visccoef][self.arhoE] * data["phi"][j] * data[self.grad_var] + data[self.visccoef] * dgrad_vari_dvarj) * data["grad_phi"][i] * data["JxW"]
    else:
      return self.zero
