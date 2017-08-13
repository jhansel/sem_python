import unittest
import filecmp

import sem
from InputFileModification import InputFileModification

class LaxFriedrichs2PhaseTester(unittest.TestCase):
  ## Tests the solution
  def testSolution(self):
    test_dir = "tests/stabilization/lax_friedrichs_2phase/"
    sem.run(test_dir + "lax_friedrichs_2phase.in")
    self.assertTrue(filecmp.cmp(test_dir + "solution.csv", test_dir + "gold/solution.csv"))

if __name__ == "__main__":
  mods = list()
  mods.append(InputFileModification("NonlinearSolver", "verbose", True))
  mods.append(InputFileModification("Executioner", "verbose", True))
  # mods.append(InputFileModification("Output", "save_solution", False))
  mods.append(InputFileModification("Output", "save_solution", True))
  mods.append(InputFileModification("Output", "solution_file", "solution.csv"))
  sem.run("lax_friedrichs_2phase.in", mods)
