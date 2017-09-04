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

class DoFHandler(object):
  __metaclass__ = ABCMeta
  def __init__(self, params):
    self.meshes = params.get("meshes")

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
    k_begin = 0
    elem_begin = 0
    for i_mesh, mesh in enumerate(self.meshes):
      self.mesh_name_to_mesh_index[mesh.name] = i_mesh
      mesh_n_node = mesh.n_cell + 1
      k_end = k_begin + mesh_n_node - 1
      elem_end = elem_begin + mesh.n_cell - 1
      self.x[k_begin:k_end+1] = mesh.x
      self.h[elem_begin:elem_end+1] = mesh.h
      self.elem_to_mesh_index[elem_begin:elem_end+1] = [i_mesh] * mesh.n_cell
      self.node_to_mesh_index[k_begin:k_end+1] = [i_mesh] * (mesh.n_cell + 1)
      k_begin += mesh_n_node
      elem_begin += mesh.n_cell

    # determine min and max x-positions
    self.x_min = self.meshes[0].x_min
    self.x_max = self.meshes[0].x_max
    for mesh in self.meshes:
      self.x_min = min(self.x_min, mesh.x_min)
      self.x_max = max(self.x_max, mesh.x_max)

    # number of DoFs per cell per variable (2 for linear FEM)
    self.n_dof_per_cell_per_var = 2

    # initialize number of constraints to zero
    self.n_constraints = [0] * self.n_meshes

  def updateWithJunctionConstraints(self, junctions):
    # add the number of constraints from each junction
    for junction in junctions:
      # get corresponding mesh index
      mesh_names = junction.mesh_names
      mesh_indices = [self.mesh_name_to_mesh_index[name] for name in mesh_names]
      mesh_index_min = min(mesh_indices)

      n_constraints = junction.n_constraints
      self.n_constraints[mesh_index_min+1] += n_constraints

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
      for mesh_index in xrange(mesh_index_min):
        i_previous += self.n_constraints[mesh_index_min]

      # finish computation of the constraint DoF indices
      n_constraints = junction.n_constraints
      i_local_begin = current_local_constraint_index[mesh_index_min]
      local_constraint_dof_indices = range(i_local_begin, i_local_begin + n_constraints)
      constraint_dof_indices = map(add, [i_previous] * n_constraints, local_constraint_dof_indices)

      # set the constraint DoF indices for the junction
      junction.setConstraintDoFIndices(constraint_dof_indices)

      # update the current local constraint index
      current_local_constraint_index[mesh_index_min] += n_constraints

    # update total number of DoFs
    self.n_dof += sum(self.n_constraints)

  def setup(self):
    arho_index_phase = 0
    arhou_index_phase = 2
    arhoE_index_phase = 1
    self.n_var = self.n_vf_equations + self.n_phases * 3
    self.arho_index = list()
    self.arhou_index = list()
    self.arhoE_index = list()
    self.variable_names = [""] * self.n_var
    self.index_to_variable = dict()
    self.index_to_phase = dict()
    for phase in xrange(self.n_phases):
      self.arho_index.append(self.n_vf_equations + phase * 3 + arho_index_phase)
      self.arhou_index.append(self.n_vf_equations + phase * 3 + arhou_index_phase)
      self.arhoE_index.append(self.n_vf_equations + phase * 3 + arhoE_index_phase)

      self.variable_names[self.arho_index[phase]] = "arho" + str(phase+1)
      self.variable_names[self.arhou_index[phase]] = "arhou" + str(phase+1)
      self.variable_names[self.arhoE_index[phase]] = "arhoE" + str(phase+1)

      self.index_to_variable[self.arho_index[phase]] = VariableName.ARho
      self.index_to_variable[self.arhou_index[phase]] = VariableName.ARhoU
      self.index_to_variable[self.arhoE_index[phase]] = VariableName.ARhoE

      self.index_to_phase[self.arho_index[phase]] = phase
      self.index_to_phase[self.arhou_index[phase]] = phase
      self.index_to_phase[self.arhoE_index[phase]] = phase

    self.variable_index = {
      VariableName.ARho: self.arho_index,
      VariableName.ARhoU: self.arhou_index,
      VariableName.ARhoE: self.arhoE_index}

    if self.model_type == ModelType.TwoPhase:
      self.vf1_index = [0]
      phase = 0
      self.variable_index[VariableName.VF1] = self.vf1_index
      self.variable_names[self.vf1_index[phase]] = "vf1"
      self.index_to_variable[self.vf1_index[phase]] = VariableName.VF1
      self.index_to_phase[self.vf1_index[phase]] = phase

    # total number of DoFs
    self.n_dof_per_cell = self.n_dof_per_cell_per_var * self.n_var
    self.n_dof = self.n_node * self.n_var

  ## Returns global DoF index corresponding to a node and variable
  # @param[in] k  global node index
  # @param[in] var_index  variable index
  def i(self, k, var_index):
    return k * self.n_var + var_index + self.n_constraints[self.node_to_mesh_index[k]]

  ## Returns global node index for an element index and local node index
  # @param[in] e  element index
  # @param[in] k_local  local node index
  def k(self, e, k_local):
    return e + k_local + self.elem_to_mesh_index[e]

  ## Returns a node index for a mesh, counted from the left
  # @param[in] mesh_name  name of the mesh
  # @param[in] k_i  node index of a mesh, counted from the left
  def getNodeIndexFromLeft(self, mesh_name, k_i):
    i_mesh = self.mesh_name_to_mesh_index[mesh_name]
    k = 0
    for i in xrange(i_mesh):
      k += self.meshes[i].n_cell + 1
    k += k_i
    return k

  ## Returns a node index for a mesh, counted from the right
  # @param[in] mesh_name  name of the mesh
  # @param[in] k_i  node index of a mesh, counted from the right
  def getNodeIndexFromRight(self, mesh_name, k_i):
    i_mesh = self.mesh_name_to_mesh_index[mesh_name]
    k = 0
    for i in xrange(i_mesh + 1):
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
  def getVolumeFraction(self, U, k):
    pass

  def getSolution(self, U, variable_name, phase):
    var_index = self.variable_index[variable_name][phase]
    return np.array([U[self.i(k, var_index)] for k in xrange(self.n_node)])

  def getPhaseSolution(self, U, phase):
    vf1 = np.array([self.getVolumeFraction(U, k) for k in xrange(self.n_node)])
    vf, _ = computeVolumeFraction(vf1, phase, self.model_type)
    arho = self.getSolution(U, VariableName.ARho, phase)
    arhou = self.getSolution(U, VariableName.ARhoU, phase)
    arhoE = self.getSolution(U, VariableName.ARhoE, phase)
    return (vf, arho, arhou, arhoE)

  # aggregates local cell vector into global vector
  def aggregateLocalCellVector(self, r, r_cell, e):
    i_min = self.i(self.k(e, 0), 0)
    i_max = self.i(self.k(e, 1), self.n_var - 1)
    r[i_min:i_max+1] += r_cell

  # aggregates local node vector into global vector
  def aggregateLocalNodeVector(self, r, r_node, k):
    i_min = k * self.n_var
    i_max = (k+1) * self.n_var - 1
    r[i_min:i_max+1] += r_node

  # aggregates local cell matrix into global matrix
  def aggregateLocalCellMatrix(self, J, J_cell, e):
    i_min = self.i(self.k(e, 0), 0)
    i_max = self.i(self.k(e, 1), self.n_var - 1)
    J[i_min:i_max+1, i_min:i_max+1] += J_cell

  # aggregates local node matrix into global matrix
  def aggregateLocalNodeMatrix(self, J, J_node, k):
    i_min = k * self.n_var
    i_max = (k+1) * self.n_var - 1
    J[i_min:i_max+1, i_min:i_max+1] += J_node

  ## Applies scaling factors to the nonlinear residual based on variable
  # @param[in,out] r  nonlinear residual vector to modify
  # @param[in] scaling  dictionary of variable name and phase to scaling factor
  def applyScalingFactors(self, r, scaling):
    for m in xrange(self.n_var):
      variable_m = self.index_to_variable[m]
      phase_m = self.index_to_phase[m]
      scaling_m = scaling[variable_m][phase_m]
      for k in xrange(self.n_node):
        r[self.i(k, m)] *= scaling_m

  ## Initializes derivative data to zero for each solution variable
  # @param[in] names  list of aux quantity names
  def initializeDerivativeData(self, names):
    variable_names = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]
    der = dict()
    for name in names:
      der[name] = dict()
      for var in variable_names:
        der[name][var] = 0
        der[name]["grad_" + var] = 0
    return der
