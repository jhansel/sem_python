import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType, PhaseType, VariableName

class DoFHandler(object):
  def __init__(self, mesh, model_type, ics):
    # counts
    self.n_cell = mesh.n_cell
    self.n_dof_per_var = self.n_cell + 1
    self.n_node = self.n_dof_per_var
    self.n_dof_per_cell_per_var = 2

    self.model_type = model_type

    # variable ordering
    if (self.model_type == ModelType.OnePhase):
      self.arho_index = {PhaseType.First: 0}
      self.arhou_index = {PhaseType.First: 2}
      self.arhoE_index = {PhaseType.First: 1}
      self.variable_index = {
        VariableName.ARho: self.arho_index,
        VariableName.ARhoU: self.arhou_index,
        VariableName.ARhoE: self.arhoE_index}
      self.n_var = 3
      self.variable_names = [None] * self.n_var
      self.variable_names[self.arho_index[PhaseType.First]] = "arho"
      self.variable_names[self.arhou_index[PhaseType.First]] = "arhou"
      self.variable_names[self.arhoE_index[PhaseType.First]] = "arhoE"
      self.index_to_variable = [None] * self.n_var
      self.index_to_variable[self.arho_index[PhaseType.First]] = VariableName.ARho
      self.index_to_variable[self.arhou_index[PhaseType.First]] = VariableName.ARhoU
      self.index_to_variable[self.arhoE_index[PhaseType.First]] = VariableName.ARhoE
      self.index_to_phase = [None] * self.n_var
      self.index_to_phase[self.arho_index[PhaseType.First]] = PhaseType.First
      self.index_to_phase[self.arhou_index[PhaseType.First]] = PhaseType.First
      self.index_to_phase[self.arhoE_index[PhaseType.First]] = PhaseType.First
    elif (self.model_type == ModelType.TwoPhaseNonInteracting):
      self.arho_index = {PhaseType.First: 0, PhaseType.Second: 3}
      self.arhou_index = {PhaseType.First: 2, PhaseType.Second: 5}
      self.arhoE_index = {PhaseType.First: 1, PhaseType.Second: 4}
      self.variable_index = {
        VariableName.ARho: self.arho_index,
        VariableName.ARhoU: self.arhou_index,
        VariableName.ARhoE: self.arhoE_index}
      self.n_var = 6
      self.variable_names = [None] * self.n_var
      self.variable_names[self.arho_index[PhaseType.First]] = "arho1"
      self.variable_names[self.arhou_index[PhaseType.First]] = "arhou1"
      self.variable_names[self.arhoE_index[PhaseType.First]] = "arhoE1"
      self.variable_names[self.arho_index[PhaseType.Second]] = "arho2"
      self.variable_names[self.arhou_index[PhaseType.Second]] = "arhou2"
      self.variable_names[self.arhoE_index[PhaseType.Second]] = "arhoE2"
      self.index_to_variable = [None] * self.n_var
      self.index_to_variable[self.arho_index[PhaseType.First]] = VariableName.ARho
      self.index_to_variable[self.arhou_index[PhaseType.First]] = VariableName.ARhoU
      self.index_to_variable[self.arhoE_index[PhaseType.First]] = VariableName.ARhoE
      self.index_to_variable[self.arho_index[PhaseType.Second]] = VariableName.ARho
      self.index_to_variable[self.arhou_index[PhaseType.Second]] = VariableName.ARhoU
      self.index_to_variable[self.arhoE_index[PhaseType.Second]] = VariableName.ARhoE
      self.index_to_phase = [None] * self.n_var
      self.index_to_phase[self.arho_index[PhaseType.First]] = PhaseType.First
      self.index_to_phase[self.arhou_index[PhaseType.First]] = PhaseType.First
      self.index_to_phase[self.arhoE_index[PhaseType.First]] = PhaseType.First
      self.index_to_phase[self.arho_index[PhaseType.Second]] = PhaseType.Second
      self.index_to_phase[self.arhou_index[PhaseType.Second]] = PhaseType.Second
      self.index_to_phase[self.arhoE_index[PhaseType.Second]] = PhaseType.Second
    elif (self.model_type == ModelType.TwoPhase):
      self.vf1_index = {PhaseType.First: 0}
      self.arho_index = {PhaseType.First: 1, PhaseType.Second: 4}
      self.arhou_index = {PhaseType.First: 3, PhaseType.Second: 6}
      self.arhoE_index = {PhaseType.First: 2, PhaseType.Second: 5}
      self.variable_index = {
        VariableName.VF1: self.vf1_index,
        VariableName.ARho: self.arho_index,
        VariableName.ARhoU: self.arhou_index,
        VariableName.ARhoE: self.arhoE_index}
      self.n_var = 7
      self.variable_names = [None] * self.n_var
      self.variable_names[self.vf1_index[PhaseType.First]] = "vf1"
      self.variable_names[self.arho_index[PhaseType.First]] = "arho1"
      self.variable_names[self.arhou_index[PhaseType.First]] = "arhou1"
      self.variable_names[self.arhoE_index[PhaseType.First]] = "arhoE1"
      self.variable_names[self.arho_index[PhaseType.Second]] = "arho2"
      self.variable_names[self.arhou_index[PhaseType.Second]] = "arhou2"
      self.variable_names[self.arhoE_index[PhaseType.Second]] = "arhoE2"
      self.index_to_variable = [None] * self.n_var
      self.index_to_variable[self.vf1_index[PhaseType.First]] = VariableName.VF1
      self.index_to_variable[self.arho_index[PhaseType.First]] = VariableName.ARho
      self.index_to_variable[self.arhou_index[PhaseType.First]] = VariableName.ARhoU
      self.index_to_variable[self.arhoE_index[PhaseType.First]] = VariableName.ARhoE
      self.index_to_variable[self.arho_index[PhaseType.Second]] = VariableName.ARho
      self.index_to_variable[self.arhou_index[PhaseType.Second]] = VariableName.ARhoU
      self.index_to_variable[self.arhoE_index[PhaseType.Second]] = VariableName.ARhoE
      self.index_to_phase = [None] * self.n_var
      self.index_to_phase[self.vf1_index[PhaseType.First]] = PhaseType.First
      self.index_to_phase[self.arho_index[PhaseType.First]] = PhaseType.First
      self.index_to_phase[self.arhou_index[PhaseType.First]] = PhaseType.First
      self.index_to_phase[self.arhoE_index[PhaseType.First]] = PhaseType.First
      self.index_to_phase[self.arho_index[PhaseType.Second]] = PhaseType.Second
      self.index_to_phase[self.arhou_index[PhaseType.Second]] = PhaseType.Second
      self.index_to_phase[self.arhoE_index[PhaseType.Second]] = PhaseType.Second
    else:
      raise NotImplementedError("Selected model type not implemented")

    # total number of DoFs
    self.n_dof_per_cell = self.n_dof_per_cell_per_var * self.n_var
    self.n_dof = self.n_dof_per_var * self.n_var

    # create array for volume fraction if two-phase non-interacting
    if self.model_type == ModelType.TwoPhaseNonInteracting:
      self.vf1 = np.zeros(self.n_node)
      for k in xrange(self.n_node):
        self.vf1[k] = ics.vf0(mesh.x[k])

  # DoF index
  def i(self, k, variable_name, phase=PhaseType.First):
    return k * self.n_var + self.variable_index[variable_name][phase]

  # global node index
  def k(self, e, k_local):
    return e + k_local

  def getVolumeFraction(self, U, k, phase):
    if (self.model_type == ModelType.TwoPhase):
      vf1 = U[self.i(k, VariableName.VF1)]
      if (phase == PhaseType.First):
        vf = vf1
        dvf_dvf1 = 1.0
      else:
        vf = 1.0 - vf1
        dvf_dvf1 = -1.0
    elif (self.model_type == ModelType.TwoPhaseNonInteracting):
      vf1 = self.vf1[k]
      if phase == PhaseType.First:
        vf = vf1
        dvf_dvf1 = float("NaN")
      else:
        vf = 1 - vf1
        dvf_dvf1 = float("NaN")
    elif self.model_type == ModelType.OnePhase:
      vf = 1
      dvf_dvf1 = float("NaN")
    return (vf, dvf_dvf1)

  def getSolution(self, U, variable_name, phase):
    return np.array([U[self.i(k, variable_name, phase)] for k in xrange(self.n_node)])

  def getPhaseSolution(self, U, phase):
    vf = np.array([self.getVolumeFraction(U, k, phase)[0] for k in xrange(self.n_node)])
    arho = self.getSolution(U, VariableName.ARho, phase)
    arhou = self.getSolution(U, VariableName.ARhoU, phase)
    arhoE = self.getSolution(U, VariableName.ARhoE, phase)
    return (vf, arho, arhou, arhoE)

  def getVolumeFractionSolution(self, U, phase):
    return np.array([self.getVolumeFraction(U, k, phase)[0] for k in xrange(self.n_node)])

  # aggregates local vector into global vector
  def aggregateLocalVector(self, r, r_cell, e):
    r[e * self.n_var : (e+2) * self.n_var] += r_cell

  # aggregates local matrix into global matrix
  def aggregateLocalMatrix(self, J, J_cell, e):
    J[e * self.n_var : (e+2) * self.n_var, e * self.n_var : (e+2) * self.n_var] += J_cell

  ## Applies scaling factors to the nonlinear residual based on variable
  # @param[in,out] r  nonlinear residual vector to modify
  # @param[in] scaling  dictionary of variable name and phase to scaling factor
  def applyScalingFactors(self, r, scaling):
    for m in xrange(self.n_var):
      variable_m = self.index_to_variable[m]
      phase_m = self.index_to_phase[m]
      scaling_m = scaling[variable_m][phase_m]
      for k in xrange(self.n_node):
        r[self.i(k, variable_m, phase_m)] *= scaling_m
