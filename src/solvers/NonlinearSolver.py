from copy import deepcopy
import numpy as np
from numpy.linalg import matrix_rank, inv
from termcolor import colored

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

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
    self.registerBoolParameter("debug_jacobian", "Option to debug the Jacobian", False)
    self.registerFloatParameter("finite_difference_eps", "Parameter for the FD debug Jacobian", 1e-4)
    self.registerBoolParameter("print_jacobian_inverse", "Option to print the inverse of the Jacobian", False)
    self.registerBoolParameter("verbose", "Option to print out iterations", True)

class NonlinearSolver(object):
  def __init__(self, params, assembleSystem):
    self.assembleSystem = assembleSystem
    self.max_iterations = params.get("max_iterations")
    self.absolute_tol = params.get("absolute_tolerance")
    self.debug_jacobian = params.get("debug_jacobian")
    self.fd_eps = params.get("finite_difference_eps")
    self.print_jacobian_inverse = params.get("print_jacobian_inverse")
    self.verbose = params.get("verbose")

  def solve(self, U):
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

      # report nonlinear residual
      r_norm = np.linalg.norm(r, 2)
      if self.verbose:
        if (r_norm < r_norm_old):
          color = "green"
        else:
          color = "red"
        sys.stdout.write("Iter %2i: " % (it)
                         + colored("res = %.3e\n" % (r_norm), color))

      # check for convergence
      if (r_norm <= self.absolute_tol):
        if self.verbose:
          print colored("Solution converged!", "green")
        converged = True
        break

      # solve the linear system J * dU = - r
      try:
        dU = np.linalg.solve(J, -r)
      except:
        n = r.size
        Jr = np.concatenate((J, -r.reshape((n, 1))), axis=1)
        if (matrix_rank(Jr) == matrix_rank(J)):
          errorNoTraceback("Infinitely many solutions!")
        elif (matrix_rank(Jr) == matrix_rank(J) + 1):
          errorNoTraceback("No solutions!")
        else:
          errorNoTraceback("Unknown number of solutions!")

      # update solution
      U += dU

      # increment iteration index
      it += 1

      # save old residual norm
      r_norm_old = r_norm

    if (not converged):
      errorNoTraceback("Solution did not converge in " + str(self.max_iterations) + " iterations.")

    return U
