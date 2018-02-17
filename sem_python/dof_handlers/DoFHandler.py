from abc import ABCMeta, abstractmethod
import numpy as np
from operator import add

from ..base.enums import ModelType, VariableName
from ..closures.thermodynamic_functions import computeVolumeFraction
from ..input.Parameters import Parameters


class DoFHandlerParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerParameter("meshes", "List of meshes")
        self.registerParameter("ics", "List of initial conditions")
        self.registerParameter("model_type", "Model type")


class DoFHandler(object, metaclass=ABCMeta):

    def __init__(self, params):
        self.meshes = params.get("meshes")
        self.ics = params.get("ics")
        self.model_type = params.get("model_type")

        # number of cells and DoF nodes
        self.n_cell = 0
        self.n_node = 0
        for mesh in self.meshes:
            self.n_cell += mesh.n_cell
            self.n_node += self.computeNumberOfNodesInMesh(mesh.n_cell)

        # create various mesh and indexing quantities
        self.n_meshes = len(self.meshes)
        self.x = np.zeros(self.n_node)
        self.h = np.zeros(self.n_cell)
        self.mesh_name_to_mesh_index = dict()
        self.elem_to_mesh_index = [0] * self.n_cell
        self.node_to_mesh_index = [0] * self.n_node
        self.n_nodes_before_mesh = [0] * self.n_meshes
        elem_begin = 0
        k_begin = 0
        for i_mesh, mesh in enumerate(self.meshes):
            self.mesh_name_to_mesh_index[mesh.name] = i_mesh
            mesh_n_node = self.computeNumberOfNodesInMesh(mesh.n_cell)
            elem_end = elem_begin + mesh.n_cell - 1
            k_end = k_begin + mesh_n_node - 1
            self.x[k_begin:k_end + 1] = self.getNodalPositionsFromMesh(mesh)
            self.h[elem_begin:elem_end + 1] = mesh.h
            self.elem_to_mesh_index[elem_begin:elem_end + 1] = [i_mesh] * mesh.n_cell
            self.node_to_mesh_index[k_begin:k_end + 1] = [i_mesh] * mesh_n_node
            if i_mesh != self.n_meshes - 1:
                for j_mesh in range(i_mesh + 1, self.n_meshes):
                    self.n_nodes_before_mesh[j_mesh] += mesh_n_node
            elem_begin += mesh.n_cell
            k_begin += mesh_n_node

        # number of DoFs per cell per variable (2 for linear FEM)
        self.n_dof_per_cell_per_var = self.getNumberOfDoFsPerCellPerVariable()

        # initialize number of constraints to zero
        self.n_constraints = [0] * self.n_meshes

        # compute area at each node
        self.A = np.zeros(self.n_node)
        for ic in self.ics:
            # get corresponding mesh
            mesh_name = ic.mesh_name
            i_mesh = self.mesh_name_to_mesh_index[mesh_name]
            mesh = self.meshes[i_mesh]

            # compute area for each node on mesh
            A = ic.A
            mesh_n_node = self.computeNumberOfNodesInMesh(mesh.n_cell)
            for k_mesh in range(mesh_n_node):
                k = self.k_from_k_mesh(k_mesh, i_mesh)
                self.A[k] = A(self.x[k])

        # set model-dependent data
        if self.model_type == ModelType.OnePhase:
            self.n_phases = 1
            self.n_vf_equations = 0

            self.aA1 = self.aA1_1phase

        elif self.model_type == ModelType.TwoPhaseNonInteracting:
            self.n_phases = 2
            self.n_vf_equations = 0

            # create array for volume fraction
            self.vf1 = np.zeros(self.n_node)
            for ic in self.ics:
                # get corresponding mesh
                mesh_name = ic.mesh_name
                i_mesh = self.mesh_name_to_mesh_index[mesh_name]
                mesh = self.meshes[i_mesh]

                # compute volume fraction for each node on mesh
                vf1 = ic.vf1
                mesh_n_node = self.computeNumberOfNodesInMesh(mesh.n_cell)
                for k_mesh in range(mesh_n_node):
                    k = self.k_from_k_mesh(k_mesh, i_mesh)
                    self.vf1[k] = vf1(self.x[k])

            self.aA1 = self.aA1_2phase_noninteracting

        elif self.model_type == ModelType.TwoPhase:
            self.n_phases = 2
            self.n_vf_equations = 1

            self.aA1 = self.aA1_2phase

        # perform model-dependent setup
        self.setup()

    ##
    # Computes the number of nodes for a mesh from its number of cells
    #
    # @param[in] n_cells   number of cells in the mesh
    # @returns number of DoF nodes in the mesh
    #
    @abstractmethod
    def computeNumberOfNodesInMesh(self, n_cells):
        pass

    ##
    # Gets the nodal DoF positions in a mesh
    #
    # @param[in] mesh   mesh
    # @returns nodal DoF positions in the mesh
    #
    @abstractmethod
    def getNodalPositionsFromMesh(self, mesh):
        pass

    ##
    # Gets the number of DoFs per cell per variable
    #
    # @returns the number of DoFs per cell per variable
    #
    @abstractmethod
    def getNumberOfDoFsPerCellPerVariable(self):
        pass

    def aA1_1phase(self, U, k):
        return self.A[k]

    def aA1_2phase_noninteracting(self, U, k):
        return self.vf1[k] * self.A[k]

    def aA1_2phase(self, U, k):
        return U[self.i(k, self.aA1_index[0])]

    def updateWithJunctionConstraints(self, junctions):
        # add the number of constraints from each junction
        for junction in junctions:
            # get corresponding mesh index
            mesh_names = junction.mesh_names
            mesh_indices = [self.mesh_name_to_mesh_index[name] for name in mesh_names]
            mesh_index_min = min(mesh_indices)

            n_constraints = junction.n_constraints
            self.n_constraints[mesh_index_min + 1] += n_constraints

        # give each junction its constraint DoF indices
        current_local_constraint_index = [1] * self.n_meshes
        for junction in junctions:
            # get corresponding mesh index
            mesh_names = junction.mesh_names
            mesh_indices = [self.mesh_name_to_mesh_index[name] for name in mesh_names]
            mesh_index_min = min(mesh_indices)

            # get the index of the DoF before this junction's constraint DoFs
            i_previous = 0
            # add the DoFs for all non-constraint variables
            k_previous = self.getNodeIndexFromRight(self.meshes[mesh_index_min].name, 0)
            i_previous += (k_previous + 1) * self.n_var - 1
            # add the constraint DoFs
            for mesh_index in range(mesh_index_min + 1):
                i_previous += self.n_constraints[mesh_index]

            # finish computation of the constraint DoF indices
            n_constraints = junction.n_constraints
            i_local_begin = current_local_constraint_index[mesh_index_min]
            local_constraint_dof_indices = list(range(i_local_begin, i_local_begin + n_constraints))
            constraint_dof_indices = list(
                map(add, [i_previous] * n_constraints, local_constraint_dof_indices))

            # set the constraint DoF indices for the junction
            junction.setConstraintDoFIndices(constraint_dof_indices)

            # update the current local constraint index
            current_local_constraint_index[mesh_index_min] += n_constraints

        # update total number of DoFs
        self.n_dof += sum(self.n_constraints)

    def setup(self):
        arhoA_index_phase = 0
        arhouA_index_phase = 2
        arhoEA_index_phase = 1
        self.n_var = self.n_vf_equations + self.n_phases * 3
        self.arhoA_index = list()
        self.arhouA_index = list()
        self.arhoEA_index = list()
        self.variable_names = [""] * self.n_var
        self.index_to_variable = dict()
        self.index_to_phase = dict()
        for phase in range(self.n_phases):
            self.arhoA_index.append(self.n_vf_equations + phase * 3 + arhoA_index_phase)
            self.arhouA_index.append(self.n_vf_equations + phase * 3 + arhouA_index_phase)
            self.arhoEA_index.append(self.n_vf_equations + phase * 3 + arhoEA_index_phase)

            self.variable_names[self.arhoA_index[phase]] = "arhoA" + str(phase + 1)
            self.variable_names[self.arhouA_index[phase]] = "arhouA" + str(phase + 1)
            self.variable_names[self.arhoEA_index[phase]] = "arhoEA" + str(phase + 1)

            self.index_to_variable[self.arhoA_index[phase]] = VariableName.ARhoA
            self.index_to_variable[self.arhouA_index[phase]] = VariableName.ARhoUA
            self.index_to_variable[self.arhoEA_index[phase]] = VariableName.ARhoEA

            self.index_to_phase[self.arhoA_index[phase]] = phase
            self.index_to_phase[self.arhouA_index[phase]] = phase
            self.index_to_phase[self.arhoEA_index[phase]] = phase

        self.variable_index = {
            VariableName.ARhoA: self.arhoA_index,
            VariableName.ARhoUA: self.arhouA_index,
            VariableName.ARhoEA: self.arhoEA_index
        }

        if self.model_type == ModelType.TwoPhase:
            self.aA1_index = [0]
            phase = 0
            self.variable_index[VariableName.AA1] = self.aA1_index
            self.variable_names[self.aA1_index[phase]] = "aA1"
            self.index_to_variable[self.aA1_index[phase]] = VariableName.AA1
            self.index_to_phase[self.aA1_index[phase]] = phase

        # total number of DoFs
        self.n_dof_per_cell = self.n_dof_per_cell_per_var * self.n_var
        self.n_dof = self.n_node * self.n_var

    ## Returns global DoF index corresponding to a node and variable
    # @param[in] k  global node index
    # @param[in] var_index  variable index
    def i(self, k, var_index):
        # get mesh index from node index
        m = self.node_to_mesh_index[k]

        return k * self.n_var + var_index + sum(self.n_constraints[0:m + 1])

    ## Returns global node index from a node index on a mesh
    # @param[in] k_mesh  node index on a mesh
    # @param[in] i_mesh  mesh index
    def k_from_k_mesh(self, k_mesh, i_mesh):
        return k_mesh + self.n_nodes_before_mesh[i_mesh]

    ## Returns a node index for a mesh, counted from the left
    # @param[in] mesh_name  name of the mesh
    # @param[in] k_i  node index of a mesh, counted from the left
    def getNodeIndexFromLeft(self, mesh_name, k_i):
        i_mesh = self.mesh_name_to_mesh_index[mesh_name]
        k = 0
        for i in range(i_mesh):
            k += self.computeNumberOfNodesInMesh(self.meshes[i].n_cell)
        k += k_i
        return k

    ## Returns a node index for a mesh, counted from the right
    # @param[in] mesh_name  name of the mesh
    # @param[in] k_i  node index of a mesh, counted from the right
    def getNodeIndexFromRight(self, mesh_name, k_i):
        i_mesh = self.mesh_name_to_mesh_index[mesh_name]
        k = 0
        for i in range(i_mesh + 1):
            k += self.computeNumberOfNodesInMesh(self.meshes[i].n_cell)
        k -= (k_i + 1)
        return k

    ## Converts variable enum to its string name with phase index
    # @param[in] var  variable enum
    # @param[in] phase  phase
    def variableEnumToName(self, var, phase):
        index = self.variable_index[var][phase]
        return self.variable_names[index]

    def getSolution(self, U, variable_name, phase):
        var_index = self.variable_index[variable_name][phase]
        return np.array([U[self.i(k, var_index)] for k in range(self.n_node)])

    def getPhaseSolution(self, U, phase):
        aA1 = np.array([self.aA1(U, k) for k in range(self.n_node)])
        vf, _ = computeVolumeFraction(aA1, self.A, phase, self.model_type)
        arhoA = self.getSolution(U, VariableName.ARhoA, phase)
        arhouA = self.getSolution(U, VariableName.ARhoUA, phase)
        arhoEA = self.getSolution(U, VariableName.ARhoEA, phase)
        return (vf, arhoA, arhouA, arhoEA)

    def separateNodalQuantityByMesh(self, y):
        y_by_mesh = [0] * self.n_meshes
        for i_mesh in range(self.n_meshes):
            k_min = self.n_nodes_before_mesh[i_mesh]
            k_max = k_min + self.computeNumberOfNodesInMesh(self.meshes[i_mesh].n_cell)
            y_by_mesh[i_mesh] = y[k_min:k_max]
        return y_by_mesh

    ## Applies scaling factors to the nonlinear residual based on variable
    # @param[in,out] r  nonlinear residual vector to modify
    # @param[in] scaling  dictionary of variable name and phase to scaling factor
    def applyScalingFactors(self, r, scaling):
        for m in range(self.n_var):
            variable_m = self.index_to_variable[m]
            phase_m = self.index_to_phase[m]
            scaling_m = scaling[variable_m][phase_m]
            for k in range(self.n_node):
                r[self.i(k, m)] *= scaling_m
