from copy import deepcopy

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "src/utilities")
from numeric_utilities import computeRelativeDifference

class AuxDerivativesTester(object):
  def __init__(self, verbose=False):
    self.verbose = verbose

  def checkDerivatives(self, test_aux, test_var, other_aux, other_vars, root_vars, constant_data=dict(), fd_eps=1e-8):
    # setup input data
    data = constant_data
    for i,x in enumerate(root_vars):
      data[x] = i + 2.0
    der = dict()

    # base computation
    for var in other_vars:
      other_aux[var].compute(data, der)
    test_aux.compute(data, der)
    base = data[test_var]
    hand_der = deepcopy(der)

    # compute finite difference derivatives and compute relative differences
    fd_der = dict()
    rel_diffs = dict()
    for x in root_vars:
      data_perturbed = deepcopy(data)
      data_perturbed[x] += fd_eps
      for var in other_vars:
        other_aux[var].compute(data_perturbed, der)
      test_aux.compute(data_perturbed, der)
      fd_der[x] = (data_perturbed[test_var] - base) / fd_eps
      rel_diffs[x] = computeRelativeDifference(hand_der[test_var][x], fd_der[x])

    # print results
    if self.verbose:
      print "\nTest quantity:", test_var
      for x in root_vars:
        print "\nDerivative:", x
        print "  Hand-coded        =", hand_der[test_var][x]
        print "  Finite difference =", fd_der[x]
        print "  Rel. difference   =", rel_diffs[x]

    # take the absolute value of the relative differences
    for x in rel_diffs:
      rel_diffs[x] = abs(rel_diffs[x])
    return rel_diffs
