import sem
from ParameterModification import BlockParameterModification, SubblockParameterModification
from CSVTester import CSVTester

class SolutionTester(object):
  def __init__(self, test_dir, input_file, mods=list(), solution_file_name="solution.csv"):
    self.test_dir = test_dir
    self.input_file = input_file
    self.mods = mods
    self.solution_file_name = solution_file_name

  def solutionsAreEqual(self):
    # add default modifications
    self.mods.append(BlockParameterModification("Executioner", "verbose", False))
    self.mods.append(BlockParameterModification("NonlinearSolver", "verbose", False))
    self.mods.append(BlockParameterModification("Output", "save_solution", True))
    self.mods.append(BlockParameterModification("Output", "solution_file", self.test_dir + self.solution_file_name))
    self.mods.append(BlockParameterModification("Output", "plot_solution", False))

    sem.run(self.input_file, self.mods)

    csv_tester = CSVTester(self.test_dir, self.solution_file_name)
    return csv_tester.filesAreEqual()
