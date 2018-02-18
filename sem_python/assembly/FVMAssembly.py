import numpy as np

from .Assembly import Assembly, AssemblyParameters
from ..base.enums import ModelType, VariableName
from ..utilities.assembly_utilities import initializeDerivativeData


class FVMAssemblyParameters(AssemblyParameters):

    def __init__(self, factory):
        AssemblyParameters.__init__(self, factory)


class FVMAssembly(Assembly):

    def __init__(self, params):
        Assembly.__init__(self, params)

        # create aux quantities
        self.aux_list = self.createIndependentPhaseAuxQuantities(0)
        if self.model_type != ModelType.OnePhase:
            self.aux_list += self.createIndependentPhaseAuxQuantities(1)
        if self.model_type == ModelType.TwoPhase:
            self.aux_list += self.createPhaseInteractionAuxQuantities(
            ) + interface_closures.createAuxQuantities()
        self.aux_list += stabilization.createAuxQuantities()

        # get list of aux quantities
        self.aux_names = [aux.name for aux in self.aux_list]

        # create nodal aux quantities if needed by group FVM
        if self.group_fem:
            self.nodal_aux_list = self.createNodalAuxQuantities(0)
            if self.model_type != ModelType.OnePhase:
                self.nodal_aux_list = self.createNodalAuxQuantities(1)
            self.nodal_aux_names = [aux.name for aux in self.nodal_aux_list]

            self.group_fem_interp_aux_list = self.createInterpolatedAux(0)
            if self.model_type != ModelType.OnePhase:
                self.group_fem_interp_aux_list += self.createInterpolatedAux(1)
            self.group_fem_aux_names = [aux.name for aux in self.group_fem_interp_aux_list]
        else:
            self.nodal_aux_list = list()
            self.nodal_aux_names = list()
            self.group_fem_interp_aux_list = list()
            self.group_fem_aux_names = list()

        # create list of source kernels
        self.source_kernels = self.createIndependentPhaseSourceKernels(0)
        if self.model_type != ModelType.OnePhase:
            self.source_kernels += self.createIndependentPhaseSourceKernels(1)
        if self.model_type == ModelType.TwoPhase:
            self.source_kernels += self.createPhaseInteractionSourceKernels()

        # create advection kernels
        self.fem_kernels = self.createIndependentPhaseAdvectionKernels(
            0) + stabilization.createIndependentPhaseKernels(0)
        if self.model_type != ModelType.OnePhase:
            self.fem_kernels += self.createIndependentPhaseAdvectionKernels(
                1) + stabilization.createIndependentPhaseKernels(1)
        if self.model_type == ModelType.TwoPhase:
            self.fem_kernels += self.createPhaseInteractionAdvectionKernels(
            ) + stabilization.createPhaseInteractionKernels()

        # add source kernels to kernel lists
        self.fem_kernels += self.source_kernels

    def createIndependentPhaseAuxQuantities(self, phase):
        # create list of aux quantities to create
        aux_names_phase = list()
        if phase == 0:
            if self.model_type == ModelType.OnePhase:
                aux_names_phase.append("VolumeFraction1Phase")
            else:
                aux_names_phase.append("VolumeFractionPhase1")
        else:
            aux_names_phase.append("VolumeFractionPhase2")
        aux_names_phase += ["Velocity", "SpecificTotalEnergy", "Density", \
          "SpecificVolume", "SpecificInternalEnergy", "Pressure", "Temperature", "SoundSpeed"]

        # create the aux quantities for this phase
        aux_list = list()
        for aux_name in aux_names_phase:
            params = {"phase": phase, "size": self.quadrature.n_q}
            if aux_name == "Pressure":
                params["p_function"] = self.eos_list[phase].p
            elif aux_name == "Temperature":
                params["T_function"] = self.eos_list[phase].T
            elif aux_name == "SoundSpeed":
                params["c_function"] = self.eos_list[phase].c
            aux_list.append(self.factory.createObject(aux_name, params))

        return aux_list

    def createPhaseInteractionAuxQuantities(self):
        aux_list = list()
        params = {"aux": "vf1", "variable_names": ["aA1", "A"], "size": self.quadrature.n_q}
        aux_list.append(self.factory.createObject("AuxGradient", params))

        return aux_list

    def createNodalAuxQuantities(self, phase):
        # create list of aux quantities to create
        aux_names = list()
        if phase == 0:
            if self.model_type == ModelType.OnePhase:
                aux_names.append("VolumeFraction1Phase")
            else:
                aux_names.append("VolumeFractionPhase1")
        else:
            aux_names.append("VolumeFractionPhase2")
        aux_names += [
            "Velocity", "Density", "SpecificVolume", "SpecificTotalEnergy",
            "SpecificInternalEnergy", "Pressure", "MassFlux", "MomentumFlux", "EnergyFlux"
        ]

        # create the aux quantities for this phase
        aux_list = list()
        for aux_name in aux_names:
            params = {"phase": phase, "size": self.dof_handler.n_dof_per_cell_per_var}
            if aux_name == "Pressure":
                params["p_function"] = self.eos_list[phase].p
            aux_list.append(self.factory.createObject(aux_name, params))

        return aux_list

    def createInterpolatedAux(self, phase):
        aux_list = list()

        phase_str = str(phase + 1)
        arhoA = "arhoA" + phase_str
        arhouA = "arhouA" + phase_str
        arhoEA = "arhoEA" + phase_str

        # mass
        var = "inviscflux_arhoA" + phase_str
        params = {
            "variable": var,
            "dependencies": [arhouA],
            "size": self.dof_handler.n_dof_per_cell_per_var
        }
        aux_list.append(self.factory.createObject("InterpolatedAux", params))

        # momentum
        var = "inviscflux_arhouA" + phase_str
        params = {
            "variable": var,
            "dependencies": ["aA1", arhoA, arhouA, arhoEA],
            "size": self.dof_handler.n_dof_per_cell_per_var
        }
        aux_list.append(self.factory.createObject("InterpolatedAux", params))

        # energy
        var = "inviscflux_arhoEA" + phase_str
        params = {
            "variable": var,
            "dependencies": ["aA1", arhoA, arhouA, arhoEA],
            "size": self.dof_handler.n_dof_per_cell_per_var
        }
        aux_list.append(self.factory.createObject("InterpolatedAux", params))

        return aux_list

    def createIndependentPhaseAdvectionKernels(self, phase):
        kernels = list()
        if self.group_fem:
            phase_str = str(phase + 1)
            var_enums = [VariableName.ARhoA, VariableName.ARhoUA, VariableName.ARhoEA]
            aux_vars = [
                "inviscflux_arhoA" + phase_str, "inviscflux_arhouA" + phase_str,
                "inviscflux_arhoEA" + phase_str
            ]
            for var_enum, aux_var in zip(var_enums, aux_vars):
                params = {
                    "var_enum": var_enum,
                    "aux_name": aux_var,
                    "phase": phase,
                    "dof_handler": self.dof_handler
                }
                kernels.append(self.factory.createObject("InterpolatedAdvection", params))
        else:
            params = {"phase": phase, "dof_handler": self.dof_handler}
            kernel_name_list = [
                "MassAdvection", "MomentumAdvection", "MomentumAreaGradient", "EnergyAdvection"
            ]
            for kernel_name in kernel_name_list:
                kernels.append(self.factory.createObject(kernel_name, params))
        return kernels

    def createIndependentPhaseSourceKernels(self, phase):
        params = {"phase": phase, "dof_handler": self.dof_handler}
        kernel_names = ["MomentumGravity", "EnergyGravity", "EnergyHeatTransfer"]
        kernels = [self.factory.createObject(kernel_name, params) for kernel_name in kernel_names]
        return kernels

    def createPhaseInteractionAdvectionKernels(self):
        params1 = {"phase": 0, "dof_handler": self.dof_handler}
        params2 = {"phase": 1, "dof_handler": self.dof_handler}

        kernels = list()
        kernels.append(self.factory.createObject("VolumeFractionAdvection", params1))
        kernels.append(self.factory.createObject("MomentumVolumeFractionGradient", params1))
        kernels.append(self.factory.createObject("MomentumVolumeFractionGradient", params2))
        kernels.append(self.factory.createObject("EnergyVolumeFractionGradient", params1))
        kernels.append(self.factory.createObject("EnergyVolumeFractionGradient", params2))

        return kernels

    def createPhaseInteractionSourceKernels(self):
        params1 = {"phase": 0, "dof_handler": self.dof_handler}
        params2 = {"phase": 1, "dof_handler": self.dof_handler}

        kernels = list()
        kernels.append(self.factory.createObject("VolumeFractionPressureRelaxation", params1))
        kernels.append(self.factory.createObject("EnergyPressureRelaxation", params1))
        kernels.append(self.factory.createObject("EnergyPressureRelaxation", params2))

        return kernels

    # computes the steady-state residual and Jacobian without applying strong BC
    def assembleSteadyStateSystemWithoutConstraints(self, U):
        r = np.zeros(self.dof_handler.n_dof)
        J = np.zeros(shape=(self.dof_handler.n_dof, self.dof_handler.n_dof))

        # volumetric terms
        self.addVolumetricTerms(U, r, J)

        # BCs
        for bc in self.bcs:
            bc.applyWeakBC(U, r, J)

        # junctions
        for junction in self.junctions:
            junction.applyWeaklyToNonlinearSystem(U, r, J)

        return (r, J)

    def performTransientSetup(self):
        self.M = self.computeMassMatrix()

    def computeMassMatrix(self):
        M = np.zeros(shape=(self.dof_handler.n_dof, self.dof_handler.n_dof))

        self.addMassMatrixPhase(M, 0)
        if (self.model_type != ModelType.OnePhase):
            self.addMassMatrixPhase(M, 1)
        if (self.model_type == ModelType.TwoPhase):
            self.addMassMatrixVolumeFraction(M)

        return M

    def addMassMatrixPhase(self, M, phase):
        arhoA_index = self.dof_handler.variable_index[VariableName.ARhoA][phase]
        arhouA_index = self.dof_handler.variable_index[VariableName.ARhoUA][phase]
        arhoEA_index = self.dof_handler.variable_index[VariableName.ARhoEA][phase]

        for e in range(self.dof_handler.n_cell):
            dx = self.dof_handler.dx(e)

            i_arhoA = self.dof_handler.i(e, arhoA_index)
            i_arhouA = self.dof_handler.i(e, arhouA_index)
            i_arhoEA = self.dof_handler.i(e, arhoEA_index)

            M[i_arhoA, i_arhoA] = dx
            M[i_arhouA, i_arhouA] = dx
            M[i_arhoEA, i_arhoEA] = dx

    def addMassMatrixVolumeFraction(self, M):
        aA1_index = self.dof_handler.variable_index[VariableName.AA1][0]

        for e in range(self.dof_handler.n_cell):
            dx = self.dof_handler.dx(e)
            i_aA1 = self.dof_handler.i(e, aA1_index)
            M[i_aA1, i_aA1] = dx

    def assembleTransientSystem(self, U, U_old):
        M_dU = np.matmul(self.M, U - U_old)
        return (M_dU, self.M)

    ## Applies constraints to a nonlinear system solved with Newton's method
    # @param[in] U  implicit solution vector
    # @param[in] r  nonlinear system residual vector
    # @param[in] J  nonlinear system Jacobian matrix
    def applyConstraintsToNonlinearSystem(self, U, r, J):
        # BCs
        for bc in self.bcs:
            bc.applyStrongBCNonlinearSystem(U, r, J)

        # junctions
        for junction in self.junctions:
            junction.applyStronglyToNonlinearSystem(U, r, J)

    ## Applies constraints to a linear system matrix.
    #
    # This is separated from the corresponding RHS vector modification function
    # because the matrix needs to be modified only once; the RHS vector might
    # depend on time or the solution vector.
    #
    # @param[in] A      linear system matrix
    def applyConstraintsToLinearSystemMatrix(self, A):
        # BCs
        for bc in self.bcs:
            bc.applyStrongBCLinearSystemMatrix(A)

        # junctions
        for junction in self.junctions:
            junction.applyStronglyToLinearSystemMatrix(A)

    ## Applies strong constraints to a linear system RHS vector.
    # @param[in] U   solution, needed if Dirichlet values are solution-dependent
    # @param[in] b   linear system RHS vector
    def applyConstraintsToLinearSystemRHSVector(self, U, b):
        # BCs
        for bc in self.bcs:
            bc.applyStrongBCLinearSystemRHSVector(U, b)

        # junctions
        for junction in self.junctions:
            junction.applyStronglyToLinearSystemRHSVector(U, b)

    ## Computes the steady-state residual and Jacobian
    def addVolumetricTerms(self, U, r, J):




        data = dict()
        der = initializeDerivativeData(
            self.aux_names + self.group_fem_aux_names, self.quadrature.n_q)
        nodal_data = dict()
        nodal_der = initializeDerivativeData(
            self.nodal_aux_names, self.dof_handler.n_dof_per_cell_per_var)

        data["phi"] = self.fe_values.get_phi()
        for elem in range(self.dof_handler.n_cell):
            r_cell = np.zeros(self.dof_handler.n_dof_per_cell)
            J_cell = np.zeros(
                shape=(self.dof_handler.n_dof_per_cell, self.dof_handler.n_dof_per_cell))

            i_mesh = self.dof_handler.elem_to_mesh_index[elem]

            data["g"] = np.dot(self.meshes[i_mesh].orientation, self.gravity)
            data["T_wall"] = self.ht_data[i_mesh].T_wall
            data["htc_wall"] = self.ht_data[i_mesh].htc_wall
            data["P_heat"] = self.ht_data[i_mesh].P_heat

            # compute elemental solution
            self.computeLocalCellSolution(U, elem, data)

            # compute elemental auxiliary quantities
            for aux in self.aux_list:
                aux.compute(data, der)

            # compute the local residual and Jacobian
            for kernel in self.kernels:
                kernel.apply(data, der, r_cell, J_cell)

            # aggregate cell residual and matrix into global residual and matrix
            self.dof_handler.aggregateLocalCellVector(r, r_cell, elem)
            self.dof_handler.aggregateLocalCellMatrix(J, J_cell, elem)

    ## Computes the local cell solution and gradients for 1-phase flow
    def computeLocalCellSolutionOnePhase(self, U, elem, data):
        data["A"] = self.fe_values.computeLocalArea(elem)
        data["aA1"] = self.fe_values.computeLocalVolumeFractionSolution(U, elem)
        data["arhoA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoA, 0, elem)
        data["arhouA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoUA, 0, elem)
        data["arhoEA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoEA, 0, elem)
        data["grad_A"] = self.fe_values.computeLocalAreaGradient(elem)
        if self.need_solution_gradients:
            data["grad_aA1"] = self.fe_values.computeLocalVolumeFractionSolutionGradient(U, elem)
            data["grad_arhoA1"] = self.fe_values.computeLocalSolutionGradient(
                U, VariableName.ARhoA, 0, elem)
            data["grad_arhouA1"] = self.fe_values.computeLocalSolutionGradient(
                U, VariableName.ARhoUA, 0, elem)
            data["grad_arhoEA1"] = self.fe_values.computeLocalSolutionGradient(
                U, VariableName.ARhoEA, 0, elem)

    ## Computes the local node solution for 1-phase flow
    def computeLocalNodeSolutionOnePhase(self, U, k, data):
        arhoA1_index = self.dof_handler.arhoA_index[0]
        arhouA1_index = self.dof_handler.arhouA_index[0]
        arhoEA1_index = self.dof_handler.arhoEA_index[0]
        data["A"] = self.dof_handler.A[k]
        data["aA1"] = self.dof_handler.aA1(U, k)
        data["arhoA1"] = U[self.dof_handler.i(k, arhoA1_index)]
        data["arhouA1"] = U[self.dof_handler.i(k, arhouA1_index)]
        data["arhoEA1"] = U[self.dof_handler.i(k, arhoEA1_index)]

    ## Computes the local cell solution and gradients for 2-phase flow
    def computeLocalCellSolutionTwoPhase(self, U, elem, data):
        data["A"] = self.fe_values.computeLocalArea(elem)
        data["aA1"] = self.fe_values.computeLocalVolumeFractionSolution(U, elem)
        data["arhoA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoA, 0, elem)
        data["arhouA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoUA, 0, elem)
        data["arhoEA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoEA, 0, elem)
        data["arhoA2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoA, 1, elem)
        data["arhouA2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoUA, 1, elem)
        data["arhoEA2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoEA, 1, elem)
        data["grad_A"] = self.fe_values.computeLocalAreaGradient(elem)
        data["grad_aA1"] = self.fe_values.computeLocalVolumeFractionSolutionGradient(U, elem)
        if self.need_solution_gradients:
            data["grad_arhoA1"] = self.fe_values.computeLocalSolutionGradient(
                U, VariableName.ARhoA, 0, elem)
            data["grad_arhouA1"] = self.fe_values.computeLocalSolutionGradient(
                U, VariableName.ARhoUA, 0, elem)
            data["grad_arhoEA1"] = self.fe_values.computeLocalSolutionGradient(
                U, VariableName.ARhoEA, 0, elem)
            data["grad_arhoA2"] = self.fe_values.computeLocalSolutionGradient(
                U, VariableName.ARhoA, 1, elem)
            data["grad_arhouA2"] = self.fe_values.computeLocalSolutionGradient(
                U, VariableName.ARhoUA, 1, elem)
            data["grad_arhoEA2"] = self.fe_values.computeLocalSolutionGradient(
                U, VariableName.ARhoEA, 1, elem)

    ## Computes the local node solution for 2-phase flow
    def computeLocalNodeSolutionTwoPhase(self, U, k, data):
        arhoA1_index = self.dof_handler.arhoA_index[0]
        arhouA1_index = self.dof_handler.arhouA_index[0]
        arhoEA1_index = self.dof_handler.arhoEA_index[0]
        arhoA2_index = self.dof_handler.arhoA_index[1]
        arhouA2_index = self.dof_handler.arhouA_index[1]
        arhoEA2_index = self.dof_handler.arhoEA_index[1]
        data["A"] = self.dof_handler.A[k]
        data["aA1"] = self.dof_handler.aA1(U, k)
        data["arhoA1"] = U[self.dof_handler.i(k, arhoA1_index)]
        data["arhouA1"] = U[self.dof_handler.i(k, arhouA1_index)]
        data["arhoEA1"] = U[self.dof_handler.i(k, arhoEA1_index)]
        data["arhoA2"] = U[self.dof_handler.i(k, arhoA2_index)]
        data["arhouA2"] = U[self.dof_handler.i(k, arhouA2_index)]
        data["arhoEA2"] = U[self.dof_handler.i(k, arhoEA2_index)]

    def getLocalNodalSolutionOnePhase(self, U, elem, data):
        data["A"] = self.fe_values.getLocalNodalArea(elem)
        data["arhoA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoA, 0, elem)
        data["arhouA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoUA, 0, elem)
        data["arhoEA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoEA, 0, elem)

    def getLocalNodalSolutionTwoPhase(self, U, elem, data):
        data["A"] = self.fe_values.getLocalNodalArea(elem)
        data["aA1"] = self.fe_values.getLocalNodalVolumeFractionSolution(U, elem)
        data["arhoA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoA, 0, elem)
        data["arhouA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoUA, 0, elem)
        data["arhoEA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoEA, 0, elem)
        data["arhoA2"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoA, 1, elem)
        data["arhouA2"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoUA, 1, elem)
        data["arhoEA2"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoEA, 1, elem)
