from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters


class TemperatureParameters(AuxQuantity1PhaseParameters):

    def __init__(self):
        AuxQuantity1PhaseParameters.__init__(self)
        self.registerFunctionParameter("T_function", "Temperature function from EOS")


class Temperature(AuxQuantity1Phase):

    def __init__(self, params):
        AuxQuantity1Phase.__init__(self, params)
        self.name = self.T
        self.T_function = params.get("T_function")

    def compute(self, data, der):
        T, dT_dv, dT_de = self.T_function(data[self.v], data[self.e])
        data[self.name] = T

        dT_daA1 = dT_dv * der[self.v]["aA1"]
        dT_darhoA = dT_dv * der[self.v][self.arhoA] + dT_de * der[self.e][self.arhoA]
        dT_darhouA = dT_de * der[self.e][self.arhouA]
        dT_darhoEA = dT_de * der[self.e][self.arhoEA]
        der[self.name]["aA1"] = dT_daA1
        der[self.name][self.arhoA] = dT_darhoA
        der[self.name][self.arhouA] = dT_darhouA
        der[self.name][self.arhoEA] = dT_darhoEA
