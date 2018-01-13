from .AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters


class BerryPressureRelaxationCoefParameters(AuxQuantity2PhaseParameters):

    def __init__(self):
        AuxQuantity2PhaseParameters.__init__(self)


class BerryPressureRelaxationCoef(AuxQuantity2Phase):

    def __init__(self, params):
        AuxQuantity2Phase.__init__(self, params)
        self.name = "p_relax"

    def compute(self, data, der):
        a_int = data["a_int"]
        z1 = data["z1"]
        z2 = data["z2"]

        data[self.name] = a_int / (z1 + z2)

        ddenominator = -a_int / (z1 + z2)**2
        der[self.name]["aA1"] = der["a_int"]["aA1"] / (z1 + z2) + ddenominator * (
            der["z1"]["aA1"] + der["z2"]["aA1"])
        der[self.name]["arhoA1"] = ddenominator * der["z1"]["arhoA1"]
        der[self.name]["arhouA1"] = ddenominator * der["z1"]["arhouA1"]
        der[self.name]["arhoEA1"] = ddenominator * der["z1"]["arhoEA1"]
        der[self.name]["arhoA2"] = ddenominator * der["z2"]["arhoA2"]
        der[self.name]["arhouA2"] = ddenominator * der["z2"]["arhouA2"]
        der[self.name]["arhoEA2"] = ddenominator * der["z2"]["arhoEA2"]
