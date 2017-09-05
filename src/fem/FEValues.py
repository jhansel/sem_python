import numpy as np

from Parameters import Parameters

class FEValuesParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerParameter("quadrature", "Quadrature")
    self.registerParameter("dof_handler", "Degree of freedom handler")
    self.registerParameter("meshes", "List of meshes")

class FEValues(object):
  def __init__(self, params):
    self.quadrature = params.get("quadrature")
    self.dof_handler = params.get("dof_handler")
    self.meshes = params.get("meshes")

    self.phi = np.zeros(shape=(self.dof_handler.n_dof_per_cell_per_var, self.quadrature.n_q))
    self.grad_phi = np.zeros(shape=(self.dof_handler.n_dof_per_cell_per_var, self.quadrature.n_q))
    for q in xrange(self.quadrature.n_q):
      self.phi[0, q] = self.phi_left(self.quadrature.z[q])
      self.phi[1, q] = self.phi_right(self.quadrature.z[q])
      self.grad_phi[0, q] = -0.5
      self.grad_phi[1, q] = 0.5

  def phi_left(self, z):
    return - 0.5 * z + 0.5

  def phi_right(self, z):
    return 0.5 * z + 0.5

  def get_phi(self):
    return self.phi

  def get_grad_phi(self, e):
    Jac = self.quadrature.Jac_divided_by_h * self.dof_handler.h[e]
    return self.grad_phi / Jac

  def get_JxW(self, e):
    return [self.quadrature.JxW_divided_by_h[q] * self.dof_handler.h[e] for q in xrange(self.quadrature.n_q)]

  def computeLocalArea(self, e):
    A = np.zeros(self.quadrature.n_q)
    for q in xrange(self.quadrature.n_q):
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        k = self.dof_handler.k(e, k_local)
        A[q] += self.dof_handler.A[k] * self.phi[k_local, q]
    return A

  def computeLocalAreaGradient(self, e):
    grad_A = np.zeros(self.quadrature.n_q)
    Jac = self.quadrature.Jac_divided_by_h * self.dof_handler.h[e]
    for q in xrange(self.quadrature.n_q):
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        k = self.dof_handler.k(e, k_local)
        grad_A[q] += self.dof_handler.A[k] * self.grad_phi[k_local, q] / Jac
    return grad_A

  def computeLocalVolumeFractionSolution(self, U, e):
    aA1 = np.zeros(self.quadrature.n_q)
    for q in xrange(self.quadrature.n_q):
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        k = self.dof_handler.k(e, k_local)
        aA1_k = self.dof_handler.aA1(U, k)
        aA1[q] += aA1_k * self.phi[k_local, q]
    return aA1

  def computeLocalVolumeFractionSolutionGradient(self, U, e):
    aA1_grad = np.zeros(self.quadrature.n_q)
    Jac = self.quadrature.Jac_divided_by_h * self.dof_handler.h[e]
    for q in xrange(self.quadrature.n_q):
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        k = self.dof_handler.k(e, k_local)
        aA1_k = self.dof_handler.aA1(U, k)
        aA1_grad[q] += aA1_k * self.grad_phi[k_local, q] / Jac
    return aA1_grad

  def computeLocalSolution(self, U, variable_name, phase, e):
    var_index = self.dof_handler.variable_index[variable_name][phase]
    solution = np.zeros(self.quadrature.n_q)
    for q in xrange(self.quadrature.n_q):
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        k = self.dof_handler.k(e, k_local)
        i = self.dof_handler.i(k, var_index)
        solution[q] += U[i] * self.phi[k_local, q]
    return solution

  def computeLocalSolutionGradient(self, U, variable_name, phase, e):
    var_index = self.dof_handler.variable_index[variable_name][phase]
    solution_grad = np.zeros(self.quadrature.n_q)
    Jac = self.quadrature.Jac_divided_by_h * self.dof_handler.h[e]
    for q in xrange(self.quadrature.n_q):
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        k = self.dof_handler.k(e, k_local)
        i = self.dof_handler.i(k, var_index)
        # compute gradient in real space, not reference space
        solution_grad[q] += U[i] * self.grad_phi[k_local, q] / Jac
    return solution_grad
