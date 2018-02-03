from .InterfaceClosures import InterfaceClosures, InterfaceClosuresParameters


class AmbrosoInterfaceClosuresParameters(InterfaceClosuresParameters):

    def __init__(self, factory):
        InterfaceClosuresParameters.__init__(self, factory)
        self.registerFloatParameter("chi", "Weight fraction for phase 1", 0.5)
        self.registerFloatParameter("pressure_relaxation_time", "Relaxation time for pressures")


class AmbrosoInterfaceClosures(InterfaceClosures):

    def __init__(self, params):
        InterfaceClosures.__init__(self, params)
        self.pressure_relaxation_time = params.get("pressure_relaxation_time")
        self.chi = params.get("chi")  # should be in (0,1)

    def createAuxQuantities(self):
        interaction_aux_names = [
            "AmbrosoBeta", "AmbrosoVelocityRelaxationCoef", "AmbrosoPressureRelaxationCoef",
            "AmbrosoInterfaceVelocity", "AmbrosoInterfacePressure"
        ]
        interaction_aux = list()
        for aux_name in interaction_aux_names:
            params = {"size": self.n_q}
            if aux_name == "AmbrosoBeta":
                params["chi"] = self.chi
            elif aux_name == "AmbrosoPressureRelaxationCoef":
                params["pressure_relaxation_time"] = self.pressure_relaxation_time
            interaction_aux.append(self.factory.createObject(aux_name, params))

        params = {"original_aux": "pI", "copy_aux": "pI_bar", "size": self.n_q}
        interaction_aux.append(self.factory.createObject("IdenticalAux", params))

        return interaction_aux
