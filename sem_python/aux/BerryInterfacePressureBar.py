from .AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters


class BerryInterfacePressureBarParameters(AuxQuantity2PhaseParameters):

    def __init__(self):
        AuxQuantity2PhaseParameters.__init__(self)


class BerryInterfacePressureBar(AuxQuantity2Phase):

    def __init__(self, params):
        AuxQuantity2Phase.__init__(self, params)
        self.name = "pI_bar"

    def compute(self, data, der):
        p1 = data["p1"]
        p2 = data["p2"]
        z1 = data["z1"]
        z2 = data["z2"]

        numerator = (z1 * p2 + z2 * p1)
        denominator = z1 + z2
        data[self.name] = numerator / denominator

        dnumerator = 1.0 / denominator
        ddenominator = -numerator / denominator**2
        der[self.name]["aA1"] = dnumerator * (der["z1"]["aA1"] * p2 + z1 * der["p2"]["aA1"] + der["z2"]["aA1"] * p1 + z2 * der["p1"]["aA1"]) \
          + ddenominator * (der["z1"]["aA1"] + der["z2"]["aA1"])
        der[self.name]["arhoA1"] = dnumerator * (der["z1"]["arhoA1"] * p2 + z2 * der["p1"]["arhoA1"]) \
          + ddenominator * der["z1"]["arhoA1"]
        der[self.name]["arhouA1"] = dnumerator * (der["z1"]["arhouA1"] * p2 + z2 * der["p1"]["arhouA1"]) \
          + ddenominator * der["z1"]["arhouA1"]
        der[self.name]["arhoEA1"] = dnumerator * (der["z1"]["arhoEA1"] * p2 + z2 * der["p1"]["arhoEA1"]) \
          + ddenominator * der["z1"]["arhoEA1"]
        der[self.name]["arhoA2"] = dnumerator * (z1 * der["p2"]["arhoA2"] + der["z2"]["arhoA2"] * p1) \
          + ddenominator * der["z2"]["arhoA2"]
        der[self.name]["arhouA2"] = dnumerator * (z1 * der["p2"]["arhouA2"] + der["z2"]["arhouA2"] * p1) \
          + ddenominator * der["z2"]["arhouA2"]
        der[self.name]["arhoEA2"] = dnumerator * (z1 * der["p2"]["arhoEA2"] + der["z2"]["arhoEA2"] * p1) \
          + ddenominator * der["z2"]["arhoEA2"]
