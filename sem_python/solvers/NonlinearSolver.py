from copy import deepcopy
import numpy as np
from numpy.linalg import matrix_rank, inv, eigvals, svd
import sys
from termcolor import colored

from ..base.enums import VariableName
from ..input.Parameters import Parameters
from ..utilities.display_utilities import computeRelativeDifferenceMatrix, printMatrix, \
  printRelativeMatrixDifference, printDoFVector, printMatrixSparsity
from ..utilities.error_utilities import errorNoTraceback

class NonlinearSolverParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFloatParameter("absolute_tolerance", "Absolute tolerance for nonlinear solve", 1e-6)
    self.registerFloatParameter("relative_tolerance", "Relative tolerance for nonlinear solve", 1e-6)
    self.registerIntParameter("max_iterations", "Maximum number of nonlinear iterations", 10)

    self.registerBoolParameter("print_variable_residual_norms", "Option to print the individual variable residual norms", False)
    self.registerBoolParameter("print_residual", "Option to print the residual vector", False)
    self.registerBoolParameter("debug_jacobian", "Option to debug the Jacobian", False)
    self.registerBoolParameter("visualize_sparsity", "Option to visualize the Jacobian sparsity pattern", False)
    self.registerFloatParameter("finite_difference_eps", "Parameter for the FD debug Jacobian", 1e-4)
    self.registerBoolParameter("print_jacobian_inverse", "Option to print the inverse of the Jacobian", False)
    self.registerBoolParameter("verbose", "Option to print out iterations", True)
    self.registerBoolParameter("print_svd", "Option to print singular value decomposition", False)
    self.registerBoolParameter("print_eigenvalues", "Option to print eigenvalues", False)
    self.registerBoolParameter("save_linear_system", "Option to save nonlinear residual and jacobian to CSV files", False)

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
    self.visualize_sparsity = params.get("visualize_sparsity")
    self.fd_eps = params.get("finite_difference_eps")
    self.print_jacobian_inverse = params.get("print_jacobian_inverse")
    self.verbose = params.get("verbose")
    self.print_eigenvalues = params.get("print_eigenvalues")
    self.print_svd = params.get("print_svd")
    self.save_linear_system = params.get("save_linear_system")

    self.scaling = dict()
    self.scaling[VariableName.AA1] = [params.get("scaling_vf1")]
    self.scaling[VariableName.ARhoA] = [params.get("scaling_arhoA1"), params.get("scaling_arhoA2")]
    self.scaling[VariableName.ARhoUA] = [params.get("scaling_arhouA1"), params.get("scaling_arhouA2")]
    self.scaling[VariableName.ARhoEA] = [params.get("scaling_arhoEA1"), params.get("scaling_arhoEA2")]

  def solve(self, U):
    # begin Newton solve
    it = 1
    converged = False
    r_norm_abs_old = 1e15
    while it <= self.max_iterations:
      # compute the residual and Jacobian
      r, J = self.assembleSystem(U)

      # print eigenvalues
      if self.print_eigenvalues:
        print eigvals(J)

      # print SVD
      if self.print_svd:
        _, svds, _ = svd(J, full_matrices=True)
        print np.sort(svds)

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

      # visualize sparsity pattern of Jacobian
      if (self.visualize_sparsity):
        printMatrixSparsity(J)
        sys.exit()

      # save the linear system to CSV files
      if self.save_linear_system:
        np.savetxt("residual.csv", r, delimiter=",")
        np.savetxt("jacobian.csv", J, delimiter=",")
        sys.exit()

      # apply scaling factors
      r_scaled = deepcopy(r)
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
          print colored("Solution converged!", "green")
        converged = True
        break

      # solve the linear system J * dU = - r
      try:
        dU = np.linalg.solve(J, -r)
      except:
        # If the linear solve fails, it is assumed that the matrix is singular
        # and thus is rank-deficient. If the concatenated matrix J|r resolves
        # a deficiency, i.e., rank(J|r) > rank(J), then this implies that there
        # are contradictions in the system and thus there are no solutions.
        # Else the ranks are equal, which implies that there is a redundancy
        # in the system. For example, suppose the following equations are in
        # the system:
        #   x1 = A
        #   x1 = B
        # Assuming B does not equal A, this implies a contradiction and thus no
        # solution exists. Assuming the other equations are linearly independent,
        # the Jacobian matrix J will have rank N-1. However,tThe matrix J|r has
        # rank N. If B were equal to A, then there would be no contradiction,
        # just a redundancy, so there would be infinite solutions. In this case,
        # the matrix J|r would have the same rank as J: N-1.
        n = r.size
        Jr = np.concatenate((J, -r.reshape((n, 1))), axis=1)
        if matrix_rank(Jr) == matrix_rank(J):
          errorNoTraceback("Infinitely many solutions!\n")
        elif matrix_rank(Jr) > matrix_rank(J):
          errorNoTraceback("No solutions!\n")
        else:
          # This should not be possible
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
