from copy import deepcopy
import numpy as np
from termcolor import colored

from ..base.enums import ModelType, VariableName
from .Executioner import Executioner, ExecutionerParameters
from ..utilities.error_utilities import error
from ..utilities.assembly_utilities import initializeDerivativeData


class TransientExecutionerParameters(ExecutionerParameters):

    def __init__(self, factory):
        ExecutionerParameters.__init__(self, factory)
        self.registerParameter("time_step_sizer", "Time step sizer")
        self.registerBoolParameter("lump_mass_matrix", "Lump the mass matrix?", False)
        self.registerBoolParameter("multiply_by_dt", "Multiply the nonlinear system by dt?", True)
        self.registerFloatParameter("ss_tol", "Tolerance for steady-state check")


class TransientExecutioner(Executioner):

    def __init__(self, params):
        Executioner.__init__(self, params)

        self.time_step_sizer = params.get("time_step_sizer")
        self.lump_mass_matrix = params.get("lump_mass_matrix")
        self.multiply_by_dt = params.get("multiply_by_dt")

        if params.has("ss_tol"):
            self.check_ss = True
            self.ss_tol = params.get("ss_tol")
        else:
            self.check_ss = False

        self.U_old = deepcopy(self.U)

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
        phi = self.fe_values.get_phi()

        arhoA_index = self.dof_handler.variable_index[VariableName.ARhoA][phase]
        arhouA_index = self.dof_handler.variable_index[VariableName.ARhoUA][phase]
        arhoEA_index = self.dof_handler.variable_index[VariableName.ARhoEA][phase]

        for e in range(self.dof_handler.n_cell):
            M_cell = np.zeros(
                shape=(self.dof_handler.n_dof_per_cell, self.dof_handler.n_dof_per_cell))

            JxW = self.fe_values.get_JxW(e)
            for q in range(self.quadrature.n_q):
                for k_local in range(self.dof_handler.n_dof_per_cell_per_var):
                    i_arhoA = self.dof_handler.i(k_local, arhoA_index)
                    i_arhouA = self.dof_handler.i(k_local, arhouA_index)
                    i_arhoEA = self.dof_handler.i(k_local, arhoEA_index)
                    for l_local in range(self.dof_handler.n_dof_per_cell_per_var):
                        if self.lump_mass_matrix:
                            M_cell[i_arhoA, i_arhoA] += phi[k_local, q] * phi[l_local, q] * JxW[q]
                            M_cell[i_arhouA, i_arhouA] += phi[k_local, q] * phi[l_local, q] * JxW[q]
                            M_cell[i_arhoEA, i_arhoEA] += phi[k_local, q] * phi[l_local, q] * JxW[q]
                        else:
                            j_arhoA = self.dof_handler.i(l_local, arhoA_index)
                            j_arhouA = self.dof_handler.i(l_local, arhouA_index)
                            j_arhoEA = self.dof_handler.i(l_local, arhoEA_index)

                            M_cell[i_arhoA, j_arhoA] += phi[k_local, q] * phi[l_local, q] * JxW[q]
                            M_cell[i_arhouA, j_arhouA] += phi[k_local, q] * phi[l_local, q] * JxW[q]
                            M_cell[i_arhoEA, j_arhoEA] += phi[k_local, q] * phi[l_local, q] * JxW[q]

            # aggregate cell matrix into global matrix
            self.dof_handler.aggregateLocalCellMatrix(M, M_cell, e)

    def addMassMatrixVolumeFraction(self, M):
        phi = self.fe_values.get_phi()

        aA1_index = self.dof_handler.variable_index[VariableName.AA1][0]

        for e in range(self.dof_handler.n_cell):
            M_cell = np.zeros(
                shape=(self.dof_handler.n_dof_per_cell, self.dof_handler.n_dof_per_cell))

            JxW = self.fe_values.get_JxW(e)
            for q in range(self.quadrature.n_q):
                for k_local in range(self.dof_handler.n_dof_per_cell_per_var):
                    i_aA1 = self.dof_handler.i(k_local, aA1_index)
                    for l_local in range(self.dof_handler.n_dof_per_cell_per_var):
                        if self.lump_mass_matrix:
                            M_cell[i_aA1, i_aA1] += phi[k_local, q] * phi[l_local, q] * JxW[q]
                        else:
                            j_vf1 = self.dof_handler.i(l_local, aA1_index)
                            M_cell[i_aA1, j_vf1] += phi[k_local, q] * phi[l_local, q] * JxW[q]

            # aggregate cell matrix into global matrix
            self.dof_handler.aggregateLocalCellMatrix(M, M_cell, e)

    def assembleTransientSystem(self, U):
        M_dU = np.matmul(self.M, U - self.U_old)
        return (M_dU, self.M)

    ## Integrates the source terms when source-splitting is used
    def takeSourceStep(self, U, dt):
        data = dict()
        der = initializeDerivativeData(self.aux_names, 1)

        data["phi"] = np.zeros(shape=(self.dof_handler.n_dof_per_cell_per_var, 1))
        data["grad_phi"] = np.zeros(shape=(self.dof_handler.n_dof_per_cell_per_var, 1))
        data["phi"].fill(1.0)
        data["grad_phi"].fill(float("NaN"))
        data["JxW"] = 1.0

        n_var = self.dof_handler.n_var

        # loop over nodes
        for k in range(self.dof_handler.n_node):
            r_node = np.zeros(n_var)
            J_node = np.zeros(shape=(n_var, n_var))

            i_mesh = self.dof_handler.node_to_mesh_index[k]

            data["g"] = np.dot(self.meshes[i_mesh].orientation, self.gravity)
            data["T_wall"] = self.ht_data[i_mesh].T_wall
            data["htc_wall"] = self.ht_data[i_mesh].htc_wall
            data["P_heat"] = self.ht_data[i_mesh].P_heat

            # compute solution
            self.computeLocalNodeSolution(U, k, data)

            # compute auxiliary quantities
            for aux in self.aux_list:
                aux.compute(data, der)

            # compute the local residual and Jacobian
            for kernel in self.source_kernels:
                kernel.apply(data, der, r_node, J_node)

            # add integrated source (recall kernels assume LHS)
            self.dof_handler.aggregateLocalNodeVector(U, -dt * r_node, k)

    def run(self):
        while (self.time_step_sizer.transientIncomplete()):
            # compute time step size
            self.dt = self.time_step_sizer.getTimeStepSize(self.U)
            if self.verbose:
                self.time_step_sizer.printTimeStepInfo()

            # solve the time step
            self.solve()

            # perform source term integration
            if self.split_source:
                self.takeSourceStep(self.U, self.dt)

            # check for steady-state
            if self.check_ss:
                dU_dt = (self.U - self.U_old) / self.dt
                dU_dt_norm = np.linalg.norm(dU_dt, 2)
                U_old_norm = np.linalg.norm(self.U_old, 2)
                U_change_norm = dU_dt_norm / U_old_norm
                if self.verbose:
                    print("Relative solution change: %e" % (U_change_norm))
                if U_change_norm < self.ss_tol:
                    if self.verbose:
                        print(colored("\nConverged to steady-state!\n", "green"))
                    return self.U

            # save old solution and increment time step index
            self.U_old = deepcopy(self.U)

        if self.verbose:
            print("")

        return self.U

    def solve(self):
        self.nonlinear_solver.solve(self.U)
