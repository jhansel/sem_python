from abc import ABCMeta, abstractmethod
import numpy as np
from operator import add

from ..base.enums import ModelType, VariableName
from ..closures.thermodynamic_functions import computeVolumeFraction
from ..input.Parameters import Parameters


class DoFHandlerParameters(Parameters):

    def __init__(self):
        Parameters.__init__(self)
        self.registerParameter("meshes", "List of meshes")
        self.registerParameter("ics", "List of initial conditions")


class DoFHandler(object, metaclass=ABCMeta):

    def __init__(self, params):
        self.meshes = params.get("meshes")
        self.ics = params.get("ics")

        # number of cells and nodes
        self.n_cell = 0
        self.n_node = 0
        for mesh in self.meshes:
            self.n_cell += mesh.n_cell
            self.n_node += mesh.n_cell + 1

        # create various mesh and indexing quantities
        self.n_meshes = len(self.meshes)
        self.mesh_name_to_mesh_index = dict()
        self.x = np.zeros(self.n_node)
        self.h = np.zeros(self.n_cell)
        self.elem_to_mesh_index = [0] * self.n_cell
        self.node_to_mesh_index = [0] * self.n_node
        self.n_nodes_before_mesh = [0] * self.n_meshes
        k_begin = 0
        elem_begin = 0
        for i_mesh, mesh in enumerate(self.meshes):
            self.mesh_name_to_mesh_index[mesh.name] = i_mesh
            mesh_n_node = mesh.n_node
            k_end = k_begin + mesh_n_node - 1
            elem_end = elem_begin + mesh.n_cell - 1
            self.x[k_begin:k_end + 1] = mesh.x
            self.h[elem_begin:elem_end + 1] = mesh.h
            self.elem_to_mesh_index[elem_begin:elem_end + 1] = [i_mesh] * mesh.n_cell
            self.node_to_mesh_index[k_begin:k_end + 1] = [i_mesh] * mesh_n_node
            if i_mesh != self.n_meshes - 1:
                for j_mesh in range(i_mesh + 1, self.n_meshes):
                    self.n_nodes_before_mesh[j_mesh] += mesh_n_node
            k_begin += mesh_n_node
            elem_begin += mesh.n_cell

        # number of DoFs per cell per variable (2 for linear FEM)
        self.n_dof_per_cell_per_var = 2

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
            for k_mesh in range(mesh.n_node):
                k = self.k_from_k_mesh(k_mesh, i_mesh)
                self.A[k] = A(mesh.x[k_mesh])

        # data to be defined by derived classes
        self.n_vf_equations = None
        self.n_phases = None
        self.model_type = None

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

    ## Returns global node index for an element index and local node index
    # @param[in] e  element index
    # @param[in] k_local  local node index
    def k(self, e, k_local):
        return e + k_local + self.elem_to_mesh_index[e]

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
            k += self.meshes[i].n_cell + 1
        k += k_i
        return k

    ## Returns a node index for a mesh, counted from the right
    # @param[in] mesh_name  name of the mesh
    # @param[in] k_i  node index of a mesh, counted from the right
    def getNodeIndexFromRight(self, mesh_name, k_i):
        i_mesh = self.mesh_name_to_mesh_index[mesh_name]
        k = 0
        for i in range(i_mesh + 1):
            k += self.meshes[i].n_cell + 1
        k -= (k_i + 1)
        return k

    ## Converts variable enum to its string name with phase index
    # @param[in] var  variable enum
    # @param[in] phase  phase
    def variableEnumToName(self, var, phase):
        index = self.variable_index[var][phase]
        return self.variable_names[index]

    @abstractmethod
    def aA1(self, U, k):
        pass

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
            k_max = k_min + self.meshes[i_mesh].n_node
            y_by_mesh[i_mesh] = y[k_min:k_max]
        return y_by_mesh

    # aggregates local cell vector into global vector
    def aggregateLocalCellVector(self, r, r_cell, e):
        i_min = self.i(self.k(e, 0), 0)
        i_max = self.i(self.k(e, 1), self.n_var - 1)
        r[i_min:i_max + 1] += r_cell

    # aggregates local node vector into global vector
    def aggregateLocalNodeVector(self, r, r_node, k):
        i_min = k * self.n_var
        i_max = (k + 1) * self.n_var - 1
        r[i_min:i_max + 1] += r_node

    # aggregates local cell matrix into global matrix
    def aggregateLocalCellMatrix(self, J, J_cell, e):
        i_min = self.i(self.k(e, 0), 0)
        i_max = self.i(self.k(e, 1), self.n_var - 1)
        J[i_min:i_max + 1, i_min:i_max + 1] += J_cell

    # aggregates local node matrix into global matrix
    def aggregateLocalNodeMatrix(self, J, J_node, k):
        i_min = k * self.n_var
        i_max = (k + 1) * self.n_var - 1
        J[i_min:i_max + 1, i_min:i_max + 1] += J_node

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
