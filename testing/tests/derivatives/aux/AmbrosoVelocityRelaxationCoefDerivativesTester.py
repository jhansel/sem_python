import unittest

from sem_python.aux.AmbrosoVelocityRelaxationCoef import AmbrosoVelocityRelaxationCoef, AmbrosoVelocityRelaxationCoefParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# interface pressure aux
params = AmbrosoVelocityRelaxationCoefParameters()
test_aux = AmbrosoVelocityRelaxationCoef(params)

# phase-1 temperature aux
params = TestAuxParameters()
params.set("var", "T1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1"])
params.set("coefs", [1.2, 2.2, 3.2, 4.2])
T1_aux = TestAux(params)

# phase-2 temperature aux
params = TestAuxParameters()
params.set("var", "T2")
params.set("other_vars", ["aA1", "arhoA2", "arhouA2", "arhoEA2"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5])
T2_aux = TestAux(params)

# beta aux
params = TestAuxParameters()
params.set("var", "beta")
params.set("other_vars", ["arhoA1", "arhoA2"])
params.set("coefs", [1.7, 2.7])
beta_aux = TestAux(params)

other_aux = [T1_aux, T2_aux, beta_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]


class AmbrosoVelocityRelaxationCoefDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 5e-6)
