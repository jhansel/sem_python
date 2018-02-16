from abc import ABCMeta, abstractmethod
import numpy as np

from ..base.enums import ModelType, VariableName
from ..input.Parameters import Parameters


class ExecutionerParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerParameter("assembly", "Assembly")
        self.registerParameter("model", "Model")
        self.registerParameter("ics", "Initial conditions")
        self.registerParameter("junctions", "List of junctions")
        self.registerParameter("eos_list", "List of equations of state")
        self.registerParameter("dof_handler", "Degree of freedom handler")
        self.registerParameter("meshes", "List of meshes")
        self.registerBoolParameter("verbose", "Print execution information?", True)


class Executioner(object, metaclass=ABCMeta):

    def __init__(self, params):
        self.assembly = params.get("assembly")
        self.model = params.get("model")
        ics = params.get("ics")
        self.model_type = self.model.model_type
        self.junctions = params.get("junctions")
        self.eos_list = params.get("eos_list")
        self.dof_handler = params.get("dof_handler")
        self.meshes = params.get("meshes")
        self.verbose = params.get("verbose")

        # initialize the solution
        self.U = np.zeros(self.dof_handler.n_dof)
        self.initializePhaseSolution(ics, 0)
        if self.model_type != ModelType.OnePhase:
            self.initializePhaseSolution(ics, 1)
        if self.model_type == ModelType.TwoPhase:
            self.initializeVolumeFractionSolution(ics)

        # initialize the junction constraint variables
        for junction in self.junctions:
            junction.initializeConstraintVariables(self.U)

    @abstractmethod
    def run(self):
        pass

    def initializePhaseSolution(self, ics, phase):
        eos_phase = self.eos_list[phase]
        arhoA_index = self.dof_handler.variable_index[VariableName.ARhoA][phase]
        arhouA_index = self.dof_handler.variable_index[VariableName.ARhoUA][phase]
        arhoEA_index = self.dof_handler.variable_index[VariableName.ARhoEA][phase]

        for ic in ics:
            # get corresponding mesh
            mesh_name = ic.mesh_name
            i_mesh = self.dof_handler.mesh_name_to_mesh_index[mesh_name]
            mesh = self.meshes[i_mesh]

            # get appropriate volume fraction function
            if self.model_type == ModelType.OnePhase:

                def initial_vf(x):
                    return 1
            else:
                initial_vf1 = ic.vf1
                if phase == 0:
                    initial_vf = initial_vf1
                else:

                    def initial_vf(x):
                        return 1 - initial_vf1(x)

            # get relevant IC functions
            A_function = ic.A
            initial_p = ic.p[phase]
            initial_u = ic.u[phase]
            if ic.specified_rho:
                initial_rho = ic.rho[phase]
            else:
                initial_T = ic.T[phase]

            # compute IC
            for k_mesh in range(mesh.n_node):
                k = self.dof_handler.k_from_k_mesh(k_mesh, i_mesh)

                x = self.dof_handler.x[k]

                A = A_function(x)
                vf = initial_vf(x)
                p = initial_p(x)
                u = initial_u(x)
                if ic.specified_rho:
                    rho = initial_rho(x)
                else:
                    T = initial_T(x)
                    rho, _, _ = eos_phase.rho(p, T)
                e = eos_phase.e(1.0 / rho, p)[0]
                E = e + 0.5 * u * u
                self.U[self.dof_handler.i(k, arhoA_index)] = vf * rho * A
                self.U[self.dof_handler.i(k, arhouA_index)] = vf * rho * u * A
                self.U[self.dof_handler.i(k, arhoEA_index)] = vf * rho * E * A

    def initializeVolumeFractionSolution(self, ics):
        aA1_index = self.dof_handler.variable_index[VariableName.AA1][0]

        for ic in ics:
            # get corresponding mesh
            mesh_name = ic.mesh_name
            i_mesh = self.dof_handler.mesh_name_to_mesh_index[mesh_name]
            mesh = self.meshes[i_mesh]

            for k_mesh in range(mesh.n_node):
                k = self.dof_handler.k_from_k_mesh(k_mesh, i_mesh)
                x = self.dof_handler.x[k]
                self.U[self.dof_handler.i(k, aA1_index)] = ic.vf1(x) * ic.A(x)
