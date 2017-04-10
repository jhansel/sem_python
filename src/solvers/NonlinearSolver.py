from copy import deepcopy
import numpy as np
from numpy.linalg import matrix_rank

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

sys.path.append(base_dir + "src/utilities")
from display_utilities import computeRelativeDifferenceMatrix, printMatrix

class NonlinearSolverParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFloatParameter("absolute_tolerance", "Absolute tolerance for nonlinear solve", 1e-6)
    self.registerIntParameter("max_iterations", "Maximum number of nonlinear iterations", 10)
    self.registerBoolParameter("debug_jacobian", "Option to debug the Jacobian", False)
    self.registerFloatParameter("finite_difference_eps", "Parameter for the FD debug Jacobian", 1e-8)
    self.registerBoolParameter("verbose", "Option to print out iterations", True)

class NonlinearSolver(object):
  def __init__(self, params, assembleSystem):
    self.assembleSystem = assembleSystem
    self.max_iterations = params.get("max_iterations")
    self.absolute_tol = params.get("absolute_tolerance")
    self.debug_jacobian = params.get("debug_jacobian")
    self.fd_eps = params.get("finite_difference_eps")
    self.verbose = params.get("verbose")

  def solve(self, U):
    # begin Newton solve
    it = 1
    converged = False
    while it <= self.max_iterations:
      # compute the residual and Jacobian
      r, J = self.assembleSystem(U)

      # compare Jacobian to finite difference Jacobian if in debug mode
      if (self.debug_jacobian):
        # compute finite difference Jacobian
        n = r.size
        J_fd = np.zeros(shape=(n, n))
        for i in xrange(n):
          U_forward = deepcopy(U)
          U_forward[i] *= (1 + self.fd_eps)
          r_forward, J_unused = self.assembleSystem(U_forward)
          J_fd[:,i] = (r_forward - r) / (self.fd_eps * U[i])

        # print the matrices
        print("\nHand-coded Jacobian:")
        printMatrix(J)
        print("\nFinite Difference Jacobian:")
        printMatrix(J_fd)
        print("\nAbsolute Difference:")
        printMatrix(J - J_fd)
        print("\nRelative Difference:")
        J_relative_difference = computeRelativeDifferenceMatrix(J, J_fd)
        printMatrix(J_relative_difference, 1e-3, 1e-7)

        # exit
        sys.exit()

      # report nonlinear residual
      r_norm = np.linalg.norm(r, 2)
      if self.verbose:
        print "Iter %i: res = %.3e" % (it, r_norm)

      # check for convergence
      if (r_norm <= self.absolute_tol):
        if self.verbose:
          print "Solution converged!"
        converged = True
        break

      # solve the linear system J * dU = - r
      try:
        dU = np.linalg.solve(J, -r)
      except:
        n = r.size
        Jr = np.concatenate((J, -r.reshape((n, 1))), axis=1)
        if (matrix_rank(Jr) == matrix_rank(J)):
          print('Infinitely many solutions!')
        elif (matrix_rank(Jr) == matrix_rank(J) + 1):
          print('No solutions!')
        else:
          print('Unknown number of solutions!')
        sys.exit()

      # update solution
      U += dU

      # increment iteration index
      it += 1

    if (not converged):
      print("\nSolution did not converge in " + str(self.max_iterations) + " iterations.")
      sys.exit()

    return U
