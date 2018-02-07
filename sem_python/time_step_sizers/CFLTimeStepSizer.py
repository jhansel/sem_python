from .TimeStepSizer import TimeStepSizer, TimeStepSizerParameters


class CFLTimeStepSizerParameters(TimeStepSizerParameters):

    def __init__(self, factory):
        TimeStepSizerParameters.__init__(self, factory)
        self.registerFloatParameter("cfl", "CFL number to compute time step size")
        self.registerParameter("meshes", "List of meshes")
        self.registerParameter("model_type", "Model type")
        self.registerParameter("eos_list", "List of EoS objects")
        self.registerParameter("dof_handler", "DoF handler")


class CFLTimeStepSizer(TimeStepSizer):

    def __init__(self, params):
        TimeStepSizer.__init__(self, params)
        self.cfl = params.get("cfl")
        self.model_type = params.get("model_type")
        self.eos_list = params.get("eos_list")
        self.dof_handler = params.get("dof_handler")
        self.factory = params.get("factory")

        # create the necessary aux quantities
        self.aux_list = self.createCFLAuxQuantities()
        self.aux_names = [aux.name for aux in self.aux_list]

        # determine minimum cell size
        meshes = params.get("meshes")
        dx_min = meshes[0].getMinimumCellWidth()
        for mesh in meshes:
            dx_min = min(dx_min, mesh.getMinimumCellWidth())
        self.dx_min = dx_min

    def createCFLAuxQuantities(self):
        aux_list = list()

        if self.model_type == ModelType.OnePhase:
            vf_classes = ["VolumeFraction1Phase"]
            self.n_phases = 1
        else:
            vf_classes = ["VolumeFractionPhase1", "VolumeFractionPhase2"]
            self.n_phases = 2

        for phase in range(self.n_phases):
            names = [vf_classes[phase]] + [
                "Density", "SpecificVolume", "Velocity", "SpecificTotalEnergy",
                "SpecificInternalEnergy", "Pressure", "SoundSpeed"
            ]
            for name in names:
                if name == "Pressure":
                    params = {"phase": phase, "p_function": self.eos_list[phase].p}
                elif name == "SoundSpeed":
                    params = {"phase": phase, "c_function": self.eos_list[phase].c}
                else:
                    params = {"phase": phase}
                aux_list.append(self.factory.createObject(name, params))

        return aux_list

    def getTimeStepSizeInternal(self, U):
        # determine maximum wave speed
        data = dict()
        der = dict()

        for phase in range(self.n_phases):
            phase_str = str(phase + 1)
            vf = "vf" + phase_str
            arhoA = "arhoA" + phase_str
            arhouA = "arhouA" + phase_str
            arhoEA = "arhoEA" + phase_str
            data[vf], data[arhoA], data[arhouA], data[arhoEA] = self.dof_handler.getPhaseSolution(
                U, phase)

        for aux in self.aux_list:
            aux.compute(data, der)

        wave_speed_max = 0.0
        for phase in range(self.n_phases):
            phase_str = str(phase + 1)
            u = "u" + phase_str
            c = "c" + phase_str
            wave_speed = "wave_speed" + phase_str

            data[wave_speed] = abs(data[u]) + data[c]
            wave_speed_max = max(wave_speed_max, max(data[wave_speed]))

        # compute CFL time step size
        dt_CFL = self.dx_min / wave_speed_max

        return self.cfl * dt_CFL
