import numpy as np

from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters


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
        self.name = "visccoef_" + self.var + self.phase
        self.mult = params.get("mult")

    def compute(self, data, der):
        b = self.mult * 0.5 * data["dx"] * (np.abs(data[self.u]) + data[self.c])
        db_du = self.mult * 0.5 * data["dx"] * np.sign(data[self.u])
        db_dc = self.mult * 0.5 * data["dx"]

        db_daA1 = db_dc * der[self.c]["aA1"]
        db_darhoA = db_du * der[self.u][self.arhoA] + db_dc * der[self.c][self.arhoA]
        db_darhouA = db_du * der[self.u][self.arhouA] + db_dc * der[self.c][self.arhouA]
        db_darhoEA = db_dc * der[self.c][self.arhoEA]

        data[self.name] = b
        der[self.name]["aA1"] = db_daA1
        der[self.name][self.arhoA] = db_darhoA
        der[self.name][self.arhouA] = db_darhouA
        der[self.name][self.arhoEA] = db_darhoEA
