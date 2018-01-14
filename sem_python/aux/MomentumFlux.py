from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters


class MomentumFluxParameters(AuxQuantity1PhaseParameters):

    def __init__(self):
        AuxQuantity1PhaseParameters.__init__(self)


class MomentumFlux(AuxQuantity1Phase):

    def __init__(self, params):
        AuxQuantity1Phase.__init__(self, params)
        self.name = self.inviscflux_arhouA

    def compute(self, data, der):
        data[
            self.name] = data[self.arhouA] * data[self.u] + data[self.vf] * data[self.p] * data["A"]

        der[self.name]["aA1"] = (
            der[self.vf]["aA1"] * data[self.p] + data[self.vf] * der[self.p]["aA1"]) * data["A"]
        der[self.name][self.arhoA] = data[self.arhouA] * der[self.u][self.arhoA] + data[self.vf] \
            * der[self.p][self.arhoA] * data["A"]
        der[self.name][self.arhouA] = data[self.u] + data[self.arhouA] * der[self.u][self.arhouA] \
            + data[self.vf] * der[self.p][self.arhouA] * data["A"]
        der[self.name][self.arhoEA] = data[self.vf] * der[self.p][self.arhoEA] * data["A"]
