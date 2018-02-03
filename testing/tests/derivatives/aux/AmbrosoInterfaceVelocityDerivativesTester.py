import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# interface velocity aux
test_aux = factory.createObject("AmbrosoInterfaceVelocity")

# phase-1 velocity aux
params = dict()
params["var"] = "u1"
params["other_vars"] = ["arhoA1", "arhouA1"]
params["coefs"] = [1.2, 2.2]
u1_aux = factory.createObject("TestAux", params)

# phase-2 velocity aux
params = dict()
params["var"] = "u2"
params["other_vars"] = ["arhoA2", "arhouA2"]
params["coefs"] = [1.5, 2.5]
u2_aux = factory.createObject("TestAux", params)

# beta aux
params = dict()
params["var"] = "beta"
params["other_vars"] = ["arhoA1", "arhoA2"]
params["coefs"] = [1.7, 3.2]
beta_aux = factory.createObject("TestAux", params)

other_aux = [u1_aux, u2_aux, beta_aux]
root_vars = ["arhoA1", "arhouA1", "arhoA2", "arhouA2"]


class AmbrosoInterfaceVelocityDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 5e-6)
