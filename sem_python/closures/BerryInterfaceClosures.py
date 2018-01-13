from .InterfaceClosures import InterfaceClosures, InterfaceClosuresParameters


class BerryInterfaceClosuresParameters(InterfaceClosuresParameters):

    def __init__(self):
        InterfaceClosuresParameters.__init__(self)
        self.registerFloatParameter("a_int_min", "Minimum interfacial area density", 0)
        self.registerFloatParameter("a_int_max", "Maximum interfacial area density")


class BerryInterfaceClosures(InterfaceClosures):

    def __init__(self, params):
        InterfaceClosures.__init__(self, params)
        self.a_int_min = params.get("a_int_min")
        self.a_int_max = params.get("a_int_max")

    def createIndependentPhaseAuxQuantities(self, phase):
        aux_list = list()
        params = {"phase": phase}
        aux_list.append(self.factory.createObject("AcousticImpedance", params))
        return aux_list

    def createAuxQuantities(self):
        aux_list = self.createIndependentPhaseAuxQuantities(
            0) + self.createIndependentPhaseAuxQuantities(1)
        aux_names = [
            "BerryInterfacialAreaDensity", "BerryPressureRelaxationCoef",
            "BerryVelocityRelaxationCoef", "BerryInterfacePressureBar", "BerryInterfaceVelocityBar",
            "BerryInterfacePressure", "BerryInterfaceVelocity"
        ]
        for aux_name in aux_names:
            params = dict()
            if aux_name == "BerryInterfacialAreaDensity":
                params["a_int_min"] = self.a_int_min
                params["a_int_max"] = self.a_int_max
            aux_list.append(self.factory.createObject(aux_name, params))

        return aux_list
