from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters


class EnergyFluxParameters(AuxQuantity1PhaseParameters):

    def __init__(self, factory):
        AuxQuantity1PhaseParameters.__init__(self, factory)


class EnergyFlux(AuxQuantity1Phase):

    def __init__(self, params):
        AuxQuantity1Phase.__init__(self, params)
        self.name = self.inviscflux_arhoEA

    def compute(self, data, der):
        data[self.name] = data[self.u] * (
            data[self.arhoEA] + data[self.vf] * data[self.p] * data["A"])

        der[self.name]["aA1"] = data[self.u] * (
            der[self.vf]["aA1"] * data[self.p] + data[self.vf] * der[self.p]["aA1"]) * data["A"]
        der[self.name][self.arhoA] = der[self.u][self.arhoA] * (data[self.arhoEA] + data[self.vf] * data[self.p] * data["A"]) \
          + data[self.u] * data[self.vf] * der[self.p][self.arhoA] * data["A"]
        der[self.name][self.arhouA] = der[self.u][self.arhouA] * (data[self.arhoEA] + data[self.vf] * data[self.p] * data["A"]) \
          + data[self.u] * data[self.vf] * der[self.p][self.arhouA] * data["A"]
        der[self.name][self.arhoEA] = data[self.u] * (
            1 + data[self.vf] * der[self.p][self.arhoEA] * data["A"])
