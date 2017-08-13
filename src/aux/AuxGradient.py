import numpy as np

from AuxQuantity import AuxQuantity, AuxQuantityParameters

class AuxGradientParameters(AuxQuantityParameters):
  def __init__(self):
    AuxQuantityParameters.__init__(self)
    self.registerStringParameter("aux", "Name of aux for which the gradient should be computed")

## Gradient of an aux quantity
#
# Computes the gradient of an aux quantity using chain rule.
# Suppose an aux quantity \f$y\f$ is a function of \f$n\f$ variables:
# \f$y=y(u_1,\ldots,u_n)\f$. Then the gradient is computed as
# \f[
#   \nabla y = \sum\limits_{i=1}^n \frac{\partial y}{\partial u_i} \nabla u_i .
# \f]
# This aux quantity must be added after the aux quantity itself is added; this
# computation assumes that the derivatives are available in the provided
# derivatives dictionary.
#
class AuxGradient(AuxQuantity):
  def __init__(self, params):
    AuxQuantity.__init__(self, params)
    self.aux = params.get("aux")
    self.grad_aux = "grad_" + self.aux

  def compute(self, data, der):
    der_aux = der[self.aux]

    grad_aux = 0 * data[self.aux]
    grad_aux_der = dict()
    for var in der_aux:
      grad_aux += der_aux[var] * data["grad_" + var]
      grad_aux_der["grad_" + var] = der_aux[var]

    data[self.grad_aux] = grad_aux
    der[self.grad_aux] = grad_aux_der
