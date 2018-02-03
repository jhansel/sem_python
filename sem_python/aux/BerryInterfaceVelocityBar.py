from .AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters


class BerryInterfaceVelocityBarParameters(AuxQuantity2PhaseParameters):

    def __init__(self, factory):
        AuxQuantity2PhaseParameters.__init__(self, factory)


class BerryInterfaceVelocityBar(AuxQuantity2Phase):

    def __init__(self, params):
        AuxQuantity2Phase.__init__(self, params)
        self.name = "uI_bar"

    def compute(self, data, der):
        u1 = data["u1"]
        u2 = data["u2"]
        z1 = data["z1"]
        z2 = data["z2"]

        numerator = (z1 * u1 + z2 * u2)
        data[self.name] = numerator / (z1 + z2)

        ddenominator = -1.0 / (z1 + z2)**2
        der[self.name]["aA1"] = (der["z1"]["aA1"] * u1 + der["z2"]["aA1"] * u2) / (z1 + z2) \
          + numerator * ddenominator * (der["z1"]["aA1"] + der["z2"]["aA1"])
        der[self.name]["arhoA1"] = (der["z1"]["arhoA1"] * u1 + z1 * der["u1"]["arhoA1"]) / (z1 + z2) \
          + numerator * ddenominator * der["z1"]["arhoA1"]
        der[self.name]["arhouA1"] = (der["z1"]["arhouA1"] * u1 + z1 * der["u1"]["arhouA1"]) / (z1 + z2) \
          + numerator * ddenominator * der["z1"]["arhouA1"]
        der[self.name]["arhoEA1"] = der["z1"]["arhoEA1"] * u1 / (
            z1 + z2) + numerator * ddenominator * der["z1"]["arhoEA1"]
        der[self.name]["arhoA2"] = (der["z2"]["arhoA2"] * u2 + z2 * der["u2"]["arhoA2"]) / (z1 + z2) \
          + numerator * ddenominator * der["z2"]["arhoA2"]
        der[self.name]["arhouA2"] = (der["z2"]["arhouA2"] * u2 + z2 * der["u2"]["arhouA2"]) / (z1 + z2) \
          + numerator * ddenominator * der["z2"]["arhouA2"]
        der[self.name]["arhoEA2"] = der["z2"]["arhoEA2"] * u2 / (
            z1 + z2) + numerator * ddenominator * der["z2"]["arhoEA2"]
