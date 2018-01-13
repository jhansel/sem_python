from copy import deepcopy
import numpy as np

from sem_python.base.enums import ModelType
from sem_python.base.Factory import Factory
from sem_python.utilities.numeric_utilities import computeRelativeDifference


class BCTester(object):

    def checkJacobian(
            self,
            bc_name,
            model_type=ModelType.OnePhase,
            boundary="right",
            bc_params=dict(),
            fd_eps=1e-8):
        self.model_type = model_type

        # factory
        factory = Factory()

        # mesh
        params = {"n_cell": 1}
        meshes = [factory.createObject("UniformMesh", params)]

        # test node index
        if boundary == "left":
            k_test = 0
        else:
            k_test = 1

        # DoF handler
        dof_handler_params = {"meshes": meshes}
        if self.model_type == ModelType.OnePhase:
            ic_params = {"mesh_name": meshes[0].name, "A": "0.2", "rho": "1", "u": "1", "p": "1"}
            ics = [factory.createObject("InitialConditions1Phase", ic_params)]
            dof_handler_class = "DoFHandler1Phase"
        else:
            ic_params = {
                "mesh_name": meshes[0].name,
                "A": "0.2",
                "vf1": "0.3",
                "rho1": "1",
                "u1": "1",
                "p1": "1",
                "rho2": "1",
                "u2": "1",
                "p2": "1"
            }
            ics = [factory.createObject("InitialConditions2Phase", ic_params)]
            if self.model_type == ModelType.TwoPhaseNonInteracting:
                dof_handler_class = "DoFHandler2PhaseNonInteracting"
            elif self.model_type == ModelType.TwoPhase:
                dof_handler_class = "DoFHandler2Phase"
        dof_handler_params["ics"] = ics
        dof_handler = factory.createObject(dof_handler_class, dof_handler_params)
        n_dof = dof_handler.n_dof
        n_var = dof_handler.n_var

        # equation of state
        eos_params = dict()
        eos_params["gamma"] = 1.4
        eos_params["R"] = 2.0
        eos_list = [factory.createObject("IdealGasEoS", eos_params)]

        # BC
        bc_params["mesh_name"] = meshes[0].name
        bc_params["boundary"] = boundary
        bc_params["dof_handler"] = dof_handler
        bc_params["eos_list"] = eos_list
        bc = factory.createObject(bc_name, bc_params)

        # compute base solution
        U = np.zeros(n_dof)
        for i in range(n_dof):
            U[i] = i + 1.0

        # base calculation
        r = np.zeros(n_dof)
        J_hand_coded = np.zeros(shape=(n_dof, n_dof))
        bc.applyWeakBC(U, r, J_hand_coded)

        # finite difference Jacobians
        rel_diffs = np.zeros(shape=(n_var, n_var))
        J_fd = np.zeros(shape=(n_var, n_var))
        for var_index_j in range(dof_handler.n_var):
            # solution index to perturb
            j = dof_handler.i(k_test, var_index_j)

            # perturb solution
            U_perturbed = deepcopy(U)
            U_perturbed[j] += fd_eps

            # compute finite difference Jacobian
            r_perturbed = np.zeros(n_dof)
            J_perturbed = np.zeros(shape=(n_dof, n_dof))
            bc.applyWeakBC(U_perturbed, r_perturbed, J_perturbed)
            for var_index_i in range(n_var):
                # residual index tested
                i = dof_handler.i(k_test, var_index_i)

                J_fd[var_index_i][var_index_j] = (r_perturbed[i] - r[i]) / fd_eps
                rel_diffs[var_index_i][var_index_j] = computeRelativeDifference(
                    J_hand_coded[i][j], J_fd[var_index_i][var_index_j])

        # print results
        for var_index_i in range(n_var):
            i = dof_handler.i(k_test, var_index_i)
            var_i = dof_handler.variable_names[var_index_i]
            print("\nEquation variable:", var_i)
            for var_index_j in range(n_var):
                j = dof_handler.i(k_test, var_index_j)
                var_j = dof_handler.variable_names[var_index_j]
                print("\n  Derivative variable:", var_j)
                print("    Hand-coded        =", J_hand_coded[i][j])
                print("    Finite difference =", J_fd[var_index_i][var_index_j])
                print("    Rel. difference   =", rel_diffs[var_index_i][var_index_j])

        # take the absolute value of the relative differences
        return abs(rel_diffs)
