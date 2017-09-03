from copy import deepcopy
import numpy as np
from numpy.linalg import matrix_rank, inv
import sys
from termcolor import colored

from enums import VariableName
from Parameters import Parameters
from display_utilities import computeRelativeDifferenceMatrix, printMatrix, \
  printRelativeMatrixDifference, printDoFVector
from error_utilities import errorNoTraceback

class NonlinearSolverParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFloatParameter("absolute_tolerance", "Absolute tolerance for nonlinear solve", 1e-6)
    self.registerFloatParameter("relative_tolerance", "Relative tolerance for nonlinear solve", 1e-6)
    self.registerIntParameter("max_iterations", "Maximum number of nonlinear iterations", 10)
    self.registerBoolParameter("print_variable_residual_norms", "Option to print the individual variable residual norms", False)
    self.registerBoolParameter("print_residual", "Option to print the residual vector", False)
    self.registerBoolParameter("debug_jacobian", "Option to debug the Jacobian", False)
    self.registerFloatParameter("finite_difference_eps", "Parameter for the FD debug Jacobian", 1e-4)
    self.registerBoolParameter("print_jacobian_inverse", "Option to print the inverse of the Jacobian", False)
    self.registerBoolParameter("verbose", "Option to print out iterations", True)

    self.registerFloatParameter("scaling_vf1", "Scaling factor for vf1", 1)
    self.registerFloatParameter("scaling_arhoA1", "Scaling factor for arhoA1", 1)
    self.registerFloatParameter("scaling_arhouA1", "Scaling factor for arhouA1", 1)
    self.registerFloatParameter("scaling_arhoEA1", "Scaling factor for arhoEA1", 1)
    self.registerFloatParameter("scaling_arhoA2", "Scaling factor for arhoA2", 1)
    self.registerFloatParameter("scaling_arhouA2", "Scaling factor for arhouA2", 1)
    self.registerFloatParameter("scaling_arhoEA2", "Scaling factor for arhoEA2", 1)

    self.registerParameter("assemble_system_function", "System assembly function")
    self.registerParameter("dof_handler", "Degree of freedom handler")

class NonlinearSolver(object):
  def __init__(self, params):
    self.assembleSystem = params.get("assemble_system_function")
    self.dof_handler = params.get("dof_handler")

    self.max_iterations = params.get("max_iterations")
    self.absolute_tol = params.get("absolute_tolerance")
    self.relative_tol = params.get("relative_tolerance")
    self.print_variable_residual_norms = params.get("print_variable_residual_norms")
    self.print_residual = params.get("print_residual")
    self.debug_jacobian = params.get("debug_jacobian")
    self.fd_eps = params.get("finite_difference_eps")
    self.print_jacobian_inverse = params.get("print_jacobian_inverse")
    self.verbose = params.get("verbose")

    self.scaling = dict()
    self.scaling[VariableName.VF1] = [params.get("scaling_vf1")]
    self.scaling[VariableName.ARhoA] = [params.get("scaling_arhoA1"), params.get("scaling_arhoA2")]
    self.scaling[VariableName.ARhoUA] = [params.get("scaling_arhouA1"), params.get("scaling_arhouA2")]
    self.scaling[VariableName.ARhoEA] = [params.get("scaling_arhoEA1"), params.get("scaling_arhoEA2")]

  def solve(self, U, residual_factor=1.0):
    # begin Newton solve
    it = 1
    converged = False
    r_norm_abs_old = 1e15
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
        print("\nRelative Difference:")
        J_relative_difference = computeRelativeDifferenceMatrix(J, J_fd)
        printRelativeMatrixDifference(J_relative_difference, J - J_fd, 1e-1, 1e-3)

        # exit
        sys.exit()

      # apply residual factor
      r_scaled = r / residual_factor

      # apply scaling factors
      self.dof_handler.applyScalingFactors(r_scaled, self.scaling)

      # compute absolute nonlinear residual norm
      r_norm_abs = np.linalg.norm(r_scaled, 2)
      if it == 1:
        r_norm_abs_initial = r_norm_abs

      # compute relative nonlinear residual norm
      if abs(r_norm_abs_initial) < 1e-15:
        r_norm_rel = r_norm_abs
      else:
        r_norm_rel = r_norm_abs / r_norm_abs_initial

      # report residual norms
      if self.verbose:
        if (r_norm_abs < r_norm_abs_old):
          color = "green"
        else:
          color = "red"
        sys.stdout.write("Iter %2i: " % (it)
                         + colored("abs = %.3e, rel = %.3e\n" % (r_norm_abs, r_norm_rel), color))
        # print individual variable residual norms
        if self.print_variable_residual_norms:
          for m in xrange(self.dof_handler.n_var):
            r_m = np.zeros(self.dof_handler.n_node)
            for k in xrange(self.dof_handler.n_node):
              r_m[k] = r_scaled[k * self.dof_handler.n_var + m]
            r_m_norm = np.linalg.norm(r_m, 2)
            print "%7s: abs = %.3e" % (self.dof_handler.variable_names[m], r_m_norm)
          print ""
        # print residual vector
        if self.print_residual:
          printDoFVector(r_scaled, self.dof_handler)

      # check for convergence
      if r_norm_abs <= self.absolute_tol or r_norm_rel <= self.relative_tol:
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
      r_norm_abs_old = r_norm_abs

    if (not converged):
      errorNoTraceback("Solution did not converge in " + str(self.max_iterations) + " iterations.\n")

    return U
