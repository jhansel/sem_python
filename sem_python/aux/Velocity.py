from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters


class VelocityParameters(AuxQuantity1PhaseParameters):

    def __init__(self):
        AuxQuantity1PhaseParameters.__init__(self)


class Velocity(AuxQuantity1Phase):

    def __init__(self, params):
        AuxQuantity1Phase.__init__(self, params)
        self.name = self.u

    def compute(self, data, der):
        arhoA = data[self.arhoA]
        arhouA = data[self.arhouA]
        data[self.name] = arhouA / arhoA

        du_darhoA = -arhouA / arhoA / arhoA
        du_darhouA = 1.0 / arhoA
        der[self.name][self.arhoA] = du_darhoA
        der[self.name][self.arhouA] = du_darhouA
