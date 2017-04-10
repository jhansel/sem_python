import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import PhaseType

class FEValues(object):
  def __init__(self, quadrature, dof_handler, mesh):
    self.quadrature = quadrature
    self.dof_handler = dof_handler
    self.mesh = mesh
    self.phi = np.zeros(shape=(dof_handler.n_dof_per_cell_per_var, quadrature.n_q))
    self.grad_phi = np.zeros(shape=(dof_handler.n_dof_per_cell_per_var, quadrature.n_q))
    for q in xrange(quadrature.n_q):
      self.phi[0, q] = self.phi_left(quadrature.z[q])
      self.phi[1, q] = self.phi_right(quadrature.z[q])
      self.grad_phi[0, q] = -0.5
      self.grad_phi[1, q] = 0.5

  def phi_left(self, z):
    return - 0.5 * z + 0.5

  def phi_right(self, z):
    return 0.5 * z + 0.5

  def get_phi(self):
    return self.phi

  def get_grad_phi(self, e):
    Jac = self.quadrature.Jac_divided_by_h * self.mesh.h[e]
    return self.grad_phi / Jac

  def get_JxW(self, e):
    return [self.quadrature.JxW_divided_by_h[q] * self.mesh.h[e] for q in xrange(self.quadrature.n_q)]

  def computeLocalVolumeFractionSolution(self, U, phase, e):
    vf = np.zeros(self.quadrature.n_q)
    dvf_dvf1 = np.zeros(self.quadrature.n_q)
    for q in xrange(self.quadrature.n_q):
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        k = self.dof_handler.k(e, k_local)
        vf_k, dvf_dvf1_k = self.dof_handler.getVolumeFraction(U, k, phase)
        vf[q] += vf_k * self.phi[k_local, q]
        dvf_dvf1[q] += dvf_dvf1_k * self.phi[k_local, q]
    return (vf, dvf_dvf1)

  def computeLocalVolumeFractionSolutionGradient(self, U, phase, e):
    dvf_dx = np.zeros(self.quadrature.n_q)
    dvf_dx_dvf1 = np.zeros(self.quadrature.n_q)
    Jac = self.quadrature.Jac_divided_by_h * self.mesh.h[e]
    for q in xrange(self.quadrature.n_q):
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        k = self.dof_handler.k(e, k_local)
        vf_k, dvf_dvf1_k = self.dof_handler.getVolumeFraction(U, k, phase)
        dvf_dx[q] += vf_k * self.grad_phi[k_local, q] / Jac
        dvf_dx_dvf1[q] += dvf_dvf1_k * self.grad_phi[k_local, q] / Jac
    return (dvf_dx, dvf_dx_dvf1)

  def computeLocalSolution(self, U, variable_name, phase, e):
    solution = np.zeros(self.quadrature.n_q)
    for q in xrange(self.quadrature.n_q):
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        k = self.dof_handler.k(e, k_local)
        i = self.dof_handler.i(k, variable_name, phase)
        solution[q] += U[i] * self.phi[k_local, q]
    return solution

  def computeLocalSolutionGradient(self, U, variable_name, phase, e):
    solution_grad = np.zeros(self.quadrature.n_q)
    Jac = self.quadrature.Jac_divided_by_h * self.mesh.h[e]
    for q in xrange(self.quadrature.n_q):
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        k = self.dof_handler.k(e, k_local)
        i = self.dof_handler.i(k, variable_name, phase)
        # compute gradient in real space, not reference space
        solution_grad[q] += U[i] * self.grad_phi[k_local, q] / Jac
    return solution_grad
