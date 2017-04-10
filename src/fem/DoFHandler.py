import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType, PhaseType, VariableName

class DoFHandler(object):
  def __init__(self, n_cell, model_type):
    # counts
    self.n_cell = n_cell
    self.n_dof_per_var = n_cell + 1
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
    elif (self.model_type == ModelType.TwoPhaseNonInteracting):
      self.arho_index = {PhaseType.First: 0, PhaseType.Second: 3}
      self.arhou_index = {PhaseType.First: 2, PhaseType.Second: 5}
      self.arhoE_index = {PhaseType.First: 1, PhaseType.Second: 4}
      self.variable_index = {
        VariableName.ARho: self.arho_index,
        VariableName.ARhoU: self.arhou_index,
        VariableName.ARhoE: self.arhoE_index}
      self.n_var = 6
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
    else:
      raise NotImplementedError("Selected model type not implemented")

    # total number of DoFs
    self.n_dof_per_cell = self.n_dof_per_cell_per_var * self.n_var
    self.n_dof = self.n_dof_per_var * self.n_var

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
    else:
      vf = 1.0
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
