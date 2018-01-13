from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters


class SpecificTotalEnergyParameters(AuxQuantity1PhaseParameters):

    def __init__(self):
        AuxQuantity1PhaseParameters.__init__(self)


class SpecificTotalEnergy(AuxQuantity1Phase):

    def __init__(self, params):
        AuxQuantity1Phase.__init__(self, params)
        self.name = self.E

    def compute(self, data, der):
        arhoA = data[self.arhoA]
        arhoEA = data[self.arhoEA]
        data[self.name] = arhoEA / arhoA

        dE_darhoA = -arhoEA / arhoA / arhoA
        dE_darhoEA = 1.0 / arhoA
        der[self.name][self.arhoA] = dE_darhoA
        der[self.name][self.arhoEA] = dE_darhoEA
