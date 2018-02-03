import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# pressure relaxation aux
test_aux = factory.createObject("AmbrosoPressureRelaxationCoef", {"pressure_relaxation_time": 2.0})

# phase-1 volume fraction aux
params = dict()
params["var"] = "vf1"
params["other_vars"] = ["aA1"]
params["coefs"] = [1.3]
vf1_aux = factory.createObject("TestAux", params)

# phase-1 pressure aux
params = dict()
params["var"] = "p1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]
params["coefs"] = [1.2, 2.2, 3.2, 4.2]
p1_aux = factory.createObject("TestAux", params)

# phase-2 pressure aux
params = dict()
params["var"] = "p2"
params["other_vars"] = ["aA1", "arhoA2", "arhouA2", "arhoEA2"]
params["coefs"] = [1.5, 2.5, 3.5, 4.5]
p2_aux = factory.createObject("TestAux", params)

other_aux = [vf1_aux, p1_aux, p2_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]


class AmbrosoPressureRelaxationCoefDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 5e-6)
