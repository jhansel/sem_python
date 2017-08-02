import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class LaxFriedrichsCoefficientParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)
    self.registerStringParameter("var", "Solution variable name corresponding to the " \
      + "equation to which the coefficient applies, not including the phase index")
    self.registerFloatParameter("mult", "Multiplier factor for the viscous coefficient", 1.0)

## Lax-Friedrichs artificial viscosity coefficient
#
# For phase \f$k\f$, this coefficient is computed as
# \f[
#   b_k \equiv \frac{1}{2}\Delta x (\left|u_k\right| + c_k) ,
# \f]
# where \f$\Delta x\f$ is the element width.
#
class LaxFriedrichsCoefficient(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.var = params.get("var")
    self.visccoef = "visccoef_" + self.var + self.phase
    self.mult = params.get("mult")

  def compute(self, data, der):
    b = self.mult * 0.5 * data["dx"] * (np.abs(data[self.u]) + data[self.c])
    db_du = self.mult * 0.5 * data["dx"] * np.sign(data[self.u])
    db_dc = self.mult * 0.5 * data["dx"]

    db_dvf1 = db_dc * der[self.c]["vf1"]
    db_darho = db_du * der[self.u][self.arho] + db_dc * der[self.c][self.arho]
    db_darhou = db_du * der[self.u][self.arhou] + db_dc * der[self.c][self.arhou]
    db_darhoE = db_dc * der[self.c][self.arhoE]

    data[self.visccoef] = b
    der[self.visccoef] = {"vf1" : db_dvf1, self.arho : db_darho, self.arhou : db_darhou, self.arhoE : db_darhoE}
