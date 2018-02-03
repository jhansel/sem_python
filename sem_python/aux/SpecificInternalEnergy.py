from numpy import vectorize

from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters
from ..utilities.error_utilities import error


def assertNonNegativeSpecificInternalEnergySingle(e, E, u):
    if e < 0:
        error(
            "Negative specific internal energy:\n" + "  e = " + str(e) + "\n" + "  E = " + str(E) +
            "\n" + "  u = " + str(u))


assertNonNegativeSpecificInternalEnergy = vectorize(assertNonNegativeSpecificInternalEnergySingle)


class SpecificInternalEnergyParameters(AuxQuantity1PhaseParameters):

    def __init__(self, factory):
        AuxQuantity1PhaseParameters.__init__(self, factory)


class SpecificInternalEnergy(AuxQuantity1Phase):

    def __init__(self, params):
        AuxQuantity1Phase.__init__(self, params)
        self.name = self.e

    def compute(self, data, der):
        u = data[self.u]
        E = data[self.E]
        data[self.name] = E - 0.5 * u * u

        assertNonNegativeSpecificInternalEnergy(data[self.name], data[self.E], data[self.u])

        de_dE = 1.0
        de_du = -u
        de_darhoA = de_dE * der[self.E][self.arhoA] + de_du * der[self.u][self.arhoA]
        de_darhouA = de_du * der[self.u][self.arhouA]
        de_darhoEA = de_dE * der[self.E][self.arhoEA]
        der[self.name][self.arhoA] = de_darhoA
        der[self.name][self.arhouA] = de_darhouA
        der[self.name][self.arhoEA] = de_darhoEA
