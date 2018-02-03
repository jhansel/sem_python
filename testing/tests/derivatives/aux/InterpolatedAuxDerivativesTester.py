from copy import deepcopy
import numpy as np
import unittest

from sem_python.base.Factory import Factory
from sem_python.utilities.numeric_utilities import computeRelativeDifference

factory = Factory()

# test aux
params = dict()
test_var = "testvar"
params["variable"] = test_var
params["dependencies"] = ["dep1", "dep2"]
n_dof_per_cell_per_var = 2
params["size"] = n_dof_per_cell_per_var
test_aux = factory.createObject("InterpolatedAux", params)

# nodal test aux
params = dict()
params["var"] = test_var
params["other_vars"] = ["dep1", "dep2"]
params["coefs"] = [2.0, 3.0]
params["size"] = 2
nodal_test_aux = factory.createObject("TestAux", params)


class InterpolatedAuxDerivativesTester(unittest.TestCase):

    def test(self):
        root_vars = ["dep1", "dep2"]

        # setup nodal input data
        nodal_data = dict()
        for i, x in enumerate(root_vars):
            nodal_data[x] = np.array([i + 2.0, i + 1.3])

        # setup elemental input data
        elem_data = {"phi": np.array([[2.4, 1.2], [1.1, 2.2]])}

        # initialize nodal derivatives to zero
        nodal_der = {test_var: dict()}
        for var in root_vars:
            nodal_der[test_var][var] = np.zeros(n_dof_per_cell_per_var)

        # initialize elemental derivatives to zero
        elem_der = {test_var: dict()}
        for var in root_vars:
            elem_der[test_var][var] = 0

        # base computation
        nodal_test_aux.compute(nodal_data, nodal_der)
        test_aux.compute(nodal_data, nodal_der, elem_data, elem_der)

        q = 0
        base = elem_data[test_var][q]
        hand_elem_der = deepcopy(elem_der)

        # compute finite difference derivatives and compute relative differences
        fd_eps = 1e-8
        fd_der = dict()
        rel_diffs = dict()
        for x in root_vars:
            nodal_data_perturbed = deepcopy(nodal_data)
            nodal_data_perturbed[x] += fd_eps
            nodal_test_aux.compute(nodal_data_perturbed, nodal_der)
            test_aux.compute(nodal_data_perturbed, nodal_der, elem_data, elem_der)
            fd_der[x] = (elem_data[test_var][q] - base) / fd_eps
            rel_diffs[x] = computeRelativeDifference(hand_elem_der[test_var][x][q], fd_der[x])

        # print results
        print("\nTest quantity:", test_var)
        for x in root_vars:
            print("\nDerivative:", x)
            print("  Hand-coded        =", hand_elem_der[test_var][x][q])
            print("  Finite difference =", fd_der[x])
            print("  Rel. difference   =", rel_diffs[x])

        # take the absolute value of the relative differences
        for x in rel_diffs:
            self.assertLessEqual(abs(rel_diffs[x]), 1e-6)
