import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# interface pressure aux
test_aux = factory.createObject("AmbrosoVelocityRelaxationCoef")

# phase-1 temperature aux
params = dict()
params["var"] = "T1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]
params["coefs"] = [1.2, 2.2, 3.2, 4.2]
T1_aux = factory.createObject("TestAux", params)

# phase-2 temperature aux
params = dict()
params["var"] = "T2"
params["other_vars"] = ["aA1", "arhoA2", "arhouA2", "arhoEA2"]
params["coefs"] = [1.5, 2.5, 3.5, 4.5]
T2_aux = factory.createObject("TestAux", params)

# beta aux
params = dict()
params["var"] = "beta"
params["other_vars"] = ["arhoA1", "arhoA2"]
params["coefs"] = [1.7, 2.7]
beta_aux = factory.createObject("TestAux", params)

other_aux = [T1_aux, T2_aux, beta_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]


class AmbrosoVelocityRelaxationCoefDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 5e-6)
