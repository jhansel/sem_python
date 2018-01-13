from numpy import sign

from .AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters


class BerryInterfacePressureParameters(AuxQuantity2PhaseParameters):

    def __init__(self):
        AuxQuantity2PhaseParameters.__init__(self)


class BerryInterfacePressure(AuxQuantity2Phase):

    def __init__(self, params):
        AuxQuantity2Phase.__init__(self, params)
        self.name = "pI"

    def compute(self, data, der):
        pI_bar = data["pI_bar"]
        grad_aA1 = data["grad_aA1"]
        u1 = data["u1"]
        u2 = data["u2"]
        z1 = data["z1"]
        z2 = data["z2"]

        sign_grad_aA1 = sign(grad_aA1)
        numerator = sign_grad_aA1 * z1 * z2 * (u2 - u1)
        denominator = z1 + z2
        data[self.name] = pI_bar + numerator / denominator

        dnumerator = 1.0 / denominator
        ddenominator = -numerator / denominator**2
        der[self.name]["aA1"] = der["pI_bar"]["aA1"] + dnumerator * sign_grad_aA1 * (der["z1"]["aA1"] * z2 + z1 * der["z2"]["aA1"]) * (u2 - u1) \
          + ddenominator * (der["z1"]["aA1"] + der["z2"]["aA1"])
        der[self.name]["arhoA1"] = der["pI_bar"]["arhoA1"] + dnumerator * sign_grad_aA1 * z2 * (der["z1"]["arhoA1"] * (u2 - u1) - z1 * der["u1"]["arhoA1"]) \
          + ddenominator * der["z1"]["arhoA1"]
        der[self.name]["arhouA1"] = der["pI_bar"]["arhouA1"] + dnumerator * sign_grad_aA1 * z2 * (der["z1"]["arhouA1"] * (u2 - u1) - z1 * der["u1"]["arhouA1"]) \
          + ddenominator * der["z1"]["arhouA1"]
        der[self.name]["arhoEA1"] = der["pI_bar"]["arhoEA1"] + dnumerator * sign_grad_aA1 * der["z1"]["arhoEA1"] * z2 * (u2 - u1) \
          + ddenominator * der["z1"]["arhoEA1"]
        der[self.name]["arhoA2"] = der["pI_bar"]["arhoA2"] + dnumerator * sign_grad_aA1 * z1 * (der["z2"]["arhoA2"] * (u2 - u1) + z2 * der["u2"]["arhoA2"]) \
          + ddenominator * der["z2"]["arhoA2"]
        der[self.name]["arhouA2"] = der["pI_bar"]["arhouA2"] + dnumerator * sign_grad_aA1 * z1 * (der["z2"]["arhouA2"] * (u2 - u1) + z2 * der["u2"]["arhouA2"]) \
          + ddenominator * der["z2"]["arhouA2"]
        der[self.name]["arhoEA2"] = der["pI_bar"]["arhoEA2"] + dnumerator * sign_grad_aA1 * z1 * der["z2"]["arhoEA2"] * (u2 - u1) \
          + ddenominator * der["z2"]["arhoEA2"]
