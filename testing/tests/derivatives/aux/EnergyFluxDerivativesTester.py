import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# test aux
test_aux = factory.createObject("EnergyFlux", {"phase": 0})

# volume fraction aux
vf_aux = factory.createObject("VolumeFractionPhase1", {"phase": 0})

# velocity aux
params = dict()
params["var"] = "u1"
params["other_vars"] = ["arhoA1", "arhouA1"]
params["coefs"] = [1.3, 2.2]
u_aux = factory.createObject("TestAux", params)

# pressure aux
params = dict()
params["var"] = "p1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]
params["coefs"] = [1.4, 2.3, 2.1, 1.2]
p_aux = factory.createObject("TestAux", params)

other_aux = [vf_aux, u_aux, p_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]


class EnergyFluxDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(
            test_aux, other_aux, root_vars, constant_data={
                "A": 0.3
            })
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
