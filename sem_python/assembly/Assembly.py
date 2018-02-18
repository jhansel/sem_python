from abc import ABCMeta, abstractmethod

from ..input.Parameters import Parameters


class AssemblyParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerParameter("model", "Model")
        self.registerParameter("bcs", "Boundary conditions")
        self.registerParameter("junctions", "List of junctions")
        self.registerParameter("eos_list", "List of equations of state")
        self.registerParameter("interface_closures", "Interface closures")
        self.registerFloatListParameter("gravity", "3-D gravitational acceleration vector")
        self.registerParameter("ht_data", "List of HeatTransferData objects")
        self.registerParameter("dof_handler", "Degree of freedom handler")
        self.registerParameter("meshes", "List of meshes")
        self.registerParameter("factory", "Factory")

class Assembly(object, metaclass=ABCMeta):

    def __init__(self, params):
        self.model = params.get("model")
        self.model_type = self.model.model_type
        self.bcs = params.get("bcs")
        self.junctions = params.get("junctions")
        self.eos_list = params.get("eos_list")
        self.interface_closures = params.get("interface_closures")
        self.gravity = params.get("gravity")
        self.ht_data = params.get("ht_data")
        self.dof_handler = params.get("dof_handler")
        self.meshes = params.get("meshes")
        self.factory = params.get("factory")

        # create aux quantities
        self.aux_list = self.createAuxQuantities()

        # create list of source kernels
        self.source_kernels = self.createSourceKernels()

    ##
    # Creates aux quantities
    #
    # @returns list of aux quantities
    #
    def createAuxQuantities(self):
        aux_list = list()
        for phase in self.model.phases:
            # create list of aux quantities to create
            aux_names_phase = list()
            if phase == 0:
                if self.model.n_phases == 1:
                    aux_names_phase.append("VolumeFraction1Phase")
                else:
                    aux_names_phase.append("VolumeFractionPhase1")
            else:
                aux_names_phase.append("VolumeFractionPhase2")
            aux_names_phase += ["Velocity", "SpecificTotalEnergy", "Density", \
              "SpecificVolume", "SpecificInternalEnergy", "Pressure", "Temperature", "SoundSpeed"]

            # create the aux quantities for this phase
            for aux_name in aux_names_phase:
                params = {"phase": phase}
                if aux_name == "Pressure":
                    params["p_function"] = self.eos_list[phase].p
                elif aux_name == "Temperature":
                    params["T_function"] = self.eos_list[phase].T
                elif aux_name == "SoundSpeed":
                    params["c_function"] = self.eos_list[phase].c
                aux_list.append(self.factory.createObject(aux_name, params))

        if self.model.phase_interaction:
            params = {"aux": "vf1", "variable_names": ["aA1", "A"]}
            aux_list.append(self.factory.createObject("AuxGradient", params))

            # add auxes from interfacial closures object
            aux_list += self.interface_closures.createAuxQuantities()

        return aux_list

    ##
    # Creates source kernels
    #
    # @returns list of source kernels
    #
    def createSourceKernels(self):
        kernels = list()
        for phase in self.model.phases:
            params = {"phase": phase, "dof_handler": self.dof_handler}
            kernel_names = ["MomentumGravity", "EnergyGravity", "EnergyHeatTransfer"]
            kernels += [self.factory.createObject(kernel_name, params) for kernel_name in kernel_names]

        if self.model.phase_interaction:
            params1 = {"phase": 0, "dof_handler": self.dof_handler}
            params2 = {"phase": 1, "dof_handler": self.dof_handler}

            kernels.append(self.factory.createObject("VolumeFractionPressureRelaxation", params1))
            kernels.append(self.factory.createObject("EnergyPressureRelaxation", params1))
            kernels.append(self.factory.createObject("EnergyPressureRelaxation", params2))

        return kernels

    @abstractmethod
    def performTransientSetup(self):
        pass

    @abstractmethod
    def assembleSteadyStateSystemWithoutConstraints(self, U):
        pass

    @abstractmethod
    def assembleTransientSystem(self, U):
        pass

    @abstractmethod
    def applyConstraintsToNonlinearSystem(self, U, r, J):
        pass

    @abstractmethod
    def applyConstraintsToLinearSystemMatrix(self, A):
        pass

    @abstractmethod
    def applyConstraintsToLinearSystemRHSVector(self, U, b):
        pass
