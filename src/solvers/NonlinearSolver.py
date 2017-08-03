from copy import deepcopy
import numpy as np
from numpy.linalg import matrix_rank, inv
from termcolor import colored

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import VariableName

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

sys.path.append(base_dir + "src/utilities")
from display_utilities import computeRelativeDifferenceMatrix, printMatrix, printMatrixDifference
from error_utilities import errorNoTraceback

class NonlinearSolverParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFloatParameter("absolute_tolerance", "Absolute tolerance for nonlinear solve", 1e-6)
    self.registerIntParameter("max_iterations", "Maximum number of nonlinear iterations", 10)
    self.registerBoolParameter("print_variable_residual_norms", "Option to print the individual variable residual norms", False)
    self.registerBoolParameter("print_residual", "Option to print the residual vector", False)
    self.registerBoolParameter("debug_jacobian", "Option to debug the Jacobian", False)
    self.registerFloatParameter("finite_difference_eps", "Parameter for the FD debug Jacobian", 1e-4)
    self.registerBoolParameter("print_jacobian_inverse", "Option to print the inverse of the Jacobian", False)
    self.registerBoolParameter("verbose", "Option to print out iterations", True)

    self.registerFloatParameter("scaling_vf1", "Scaling factor for vf1", 1)
    self.registerFloatParameter("scaling_arho1", "Scaling factor for arho1", 1)
    self.registerFloatParameter("scaling_arhou1", "Scaling factor for arhou1", 1)
    self.registerFloatParameter("scaling_arhoE1", "Scaling factor for arhoE1", 1)
    self.registerFloatParameter("scaling_arho2", "Scaling factor for arho2", 1)
    self.registerFloatParameter("scaling_arhou2", "Scaling factor for arhou2", 1)
    self.registerFloatParameter("scaling_arhoE2", "Scaling factor for arhoE2", 1)

class NonlinearSolver(object):
  def __init__(self, params, assembleSystem, dof_handler):
    self.assembleSystem = assembleSystem
    self.dof_handler = dof_handler

    self.max_iterations = params.get("max_iterations")
    self.absolute_tol = params.get("absolute_tolerance")
    self.print_variable_residual_norms = params.get("print_variable_residual_norms")
    self.print_residual = params.get("print_residual")
    self.debug_jacobian = params.get("debug_jacobian")
    self.fd_eps = params.get("finite_difference_eps")
    self.print_jacobian_inverse = params.get("print_jacobian_inverse")
    self.verbose = params.get("verbose")

    self.scaling = dict()
    self.scaling[VariableName.VF1] = [params.get("scaling_vf1")]
    self.scaling[VariableName.ARho] = [params.get("scaling_arho1"), params.get("scaling_arho2")]
    self.scaling[VariableName.ARhoU] = [params.get("scaling_arhou1"), params.get("scaling_arhou2")]
    self.scaling[VariableName.ARhoE] = [params.get("scaling_arhoE1"), params.get("scaling_arhoE2")]

  def solve(self, U, residual_factor=1.0):
    # begin Newton solve
    it = 1
    converged = False
    r_norm_old = 1e15
    while it <= self.max_iterations:
      # compute the residual and Jacobian
      r, J = self.assembleSystem(U)

      # print the Jacobian inverse if specified
      if (self.print_jacobian_inverse):
        print("\nJacobian:")
        printMatrix(J)
        print("\nJacobian inverse:")
        printMatrix(inv(J))

        sys.exit()

      # compare Jacobian to finite difference Jacobian if in debug mode
      if (self.debug_jacobian):
        # compute finite difference Jacobian
        n = r.size
        J_fd = np.zeros(shape=(n, n))
        for j in xrange(n):
          U_forward = deepcopy(U)
          U_eps = max(self.fd_eps, abs(U[j] * self.fd_eps))
          U_forward[j] += U_eps
          r_forward, J_unused = self.assembleSystem(U_forward)
          J_fd[:,j] = (r_forward - r) / U_eps

        # print the matrices
        print("\nHand-coded Jacobian:")
        printMatrix(J)
        print("\nFinite Difference Jacobian:")
        printMatrix(J_fd)
        print("\nAbsolute Difference:")
        printMatrixDifference(J - J_fd)
        print("\nRelative Difference:")
        J_relative_difference = computeRelativeDifferenceMatrix(J, J_fd)
        printMatrixDifference(J_relative_difference, 1e-1, 1e-3)

        # exit
        sys.exit()

      # apply residual factor
      r_scaled = r / residual_factor

      # apply scaling factors
      self.dof_handler.applyScalingFactors(r_scaled, self.scaling)

      # report nonlinear residual
      r_norm = np.linalg.norm(r_scaled, 2)
      if self.verbose:
        if (r_norm < r_norm_old):
          color = "green"
        else:
          color = "red"
        sys.stdout.write("Iter %2i: " % (it)
                         + colored("res = %.3e\n" % (r_norm), color))
        # print individual variable residual norms
        if self.print_variable_residual_norms:
          for m in xrange(self.dof_handler.n_var):
            r_m = np.zeros(self.dof_handler.n_node)
            for k in xrange(self.dof_handler.n_node):
              r_m[k] = r_scaled[k * self.dof_handler.n_var + m]
            r_m_norm = np.linalg.norm(r_m, 2)
            print "%7s: res = %.3e" % (self.dof_handler.variable_names[m], r_m_norm)
          print ""
        # print residual vector
        if self.print_residual:
          header_items = ("i",) + tuple(self.dof_handler.variable_names)
          header_format = "%4s" + " %12s" * self.dof_handler.n_var
          entry_format = "%4i" + " %12.3e" * self.dof_handler.n_var
          print header_format % header_items
          for k in xrange(self.dof_handler.n_node):
            res_k = list()
            for m in xrange(self.dof_handler.n_var):
              res_k.append(r_scaled[k * self.dof_handler.n_var + m])
            entry_items = (k,) + tuple(res_k)
            print entry_format % entry_items
          print ""

      # check for convergence
      if (r_norm <= self.absolute_tol):
        if self.verbose:
          print colored("Solution converged!\n", "green")
        converged = True
        break

      # solve the linear system J * dU = - r
      try:
        dU = np.linalg.solve(J, -r)
      except:
        n = r.size
        Jr = np.concatenate((J, -r.reshape((n, 1))), axis=1)
        if (matrix_rank(Jr) == matrix_rank(J)):
          errorNoTraceback("Infinitely many solutions!\n")
        elif (matrix_rank(Jr) == matrix_rank(J) + 1):
          errorNoTraceback("No solutions!\n")
        else:
          errorNoTraceback("Unknown number of solutions!\n")

      # update solution
      U += dU

      # increment iteration index
      it += 1

      # save old residual norm
      r_norm_old = r_norm

    if (not converged):
      errorNoTraceback("Solution did not converge in " + str(self.max_iterations) + " iterations.\n")

    return U
