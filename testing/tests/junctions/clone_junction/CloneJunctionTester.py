import unittest

import sem
from InputFileModification import InputFileModification
from CSVTester import CSVTester

class CloneJunctionTester(unittest.TestCase):
  ## Tests the solution
  def testSolution(self):
    test_dir = "tests/junctions/clone_junction/"
    mods = list()
    mods.append(InputFileModification("NonlinearSolver", "verbose", False))
    mods.append(InputFileModification("Executioner", "verbose", False))
    mods.append(InputFileModification("Executioner", "end_time", 0.05))
    mods.append(InputFileModification("Output", "solution_file", "tests/junctions/clone_junction/solution_with_junction.csv"))
    mods.append(InputFileModification("Output", "plot_solution", False))
    sem.run(test_dir + "clone_junction_with_junction.in", mods)

    csv_tester = CSVTester(test_dir, "solution_with_junction.csv")
    self.assertTrue(csv_tester.filesAreEqual())
