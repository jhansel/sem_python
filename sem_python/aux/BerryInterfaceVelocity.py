from numpy import sign

from .AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters


class BerryInterfaceVelocityParameters(AuxQuantity2PhaseParameters):

    def __init__(self):
        AuxQuantity2PhaseParameters.__init__(self)


class BerryInterfaceVelocity(AuxQuantity2Phase):

    def __init__(self, params):
        AuxQuantity2Phase.__init__(self, params)
        self.name = "uI"

    def compute(self, data, der):
        uI_bar = data["uI_bar"]
        grad_aA1 = data["grad_aA1"]
        p1 = data["p1"]
        p2 = data["p2"]
        z1 = data["z1"]
        z2 = data["z2"]

        sign_grad_aA1 = sign(grad_aA1)
        numerator = sign_grad_aA1 * (p2 - p1)
        denominator = z1 + z2
        data[self.name] = uI_bar + numerator / denominator

        dnumerator = 1.0 / denominator
        ddenominator = -numerator / denominator**2
        der[self.name]["aA1"] = der["uI_bar"]["aA1"] + dnumerator * sign_grad_aA1 * (der["p2"]["aA1"] - der["p1"]["aA1"]) \
          + ddenominator * (der["z1"]["aA1"] + der["z2"]["aA1"])
        der[self.name][
            "arhoA1"] = der["uI_bar"]["arhoA1"] - dnumerator * sign_grad_aA1 * der["p1"]["arhoA1"] + ddenominator * der["z1"]["arhoA1"]
        der[self.name][
            "arhouA1"] = der["uI_bar"]["arhouA1"] - dnumerator * sign_grad_aA1 * der["p1"]["arhouA1"] + ddenominator * der["z1"]["arhouA1"]
        der[self.name][
            "arhoEA1"] = der["uI_bar"]["arhoEA1"] - dnumerator * sign_grad_aA1 * der["p1"]["arhoEA1"] + ddenominator * der["z1"]["arhoEA1"]
        der[self.name][
            "arhoA2"] = der["uI_bar"]["arhoA2"] + dnumerator * sign_grad_aA1 * der["p2"]["arhoA2"] + ddenominator * der["z2"]["arhoA2"]
        der[self.name][
            "arhouA2"] = der["uI_bar"]["arhouA2"] + dnumerator * sign_grad_aA1 * der["p2"]["arhouA2"] + ddenominator * der["z2"]["arhouA2"]
        der[self.name][
            "arhoEA2"] = der["uI_bar"]["arhoEA2"] + dnumerator * sign_grad_aA1 * der["p2"]["arhoEA2"] + ddenominator * der["z2"]["arhoEA2"]
