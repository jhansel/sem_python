from copy import deepcopy

from TestAux import TestAux, TestAuxParameters
from numeric_utilities import computeRelativeDifference

class AuxDerivativesTester(object):
  def __init__(self, verbose=False):
    self.verbose = verbose

  def checkDerivatives(self, test_aux, other_aux, root_vars, constant_data=dict(), fd_eps=1e-8):
    # name of test aux
    test_var = test_aux.name

    # setup input data
    data = constant_data
    for i,x in enumerate(root_vars):
      data[x] = i + 2.0

    # initialize derivatives to zero
    derivative_list = deepcopy(root_vars)
    if "vf1" not in derivative_list:
      derivative_list.append("vf1")
    der = dict()
    for aux in other_aux + [test_aux]:
      der[aux.name] = dict()
      for var in derivative_list:
        der[aux.name][var] = 0

    # base computation
    for aux in other_aux:
      aux.compute(data, der)
    test_aux.compute(data, der)
    base = data[test_var]
    hand_der = deepcopy(der)

    # compute finite difference derivatives and compute relative differences
    fd_der = dict()
    rel_diffs = dict()
    for x in root_vars:
      data_perturbed = deepcopy(data)
      data_perturbed[x] += fd_eps
      for aux in other_aux:
        aux.compute(data_perturbed, der)
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
