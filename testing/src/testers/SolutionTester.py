import sem

from sem_python.input.InputFileModifier import InputFileModifier
from .CSVTester import CSVTester

class SolutionTester(object):
  def __init__(self, test_dir, input_file, input_file_modifier=InputFileModifier(), solution_file_name="solution.csv"):
    self.test_dir = test_dir
    self.input_file = input_file
    self.input_file_modifier = input_file_modifier
    self.solution_file_name = solution_file_name

  def solutionsAreEqual(self):
    # add default modifications
    self.input_file_modifier.modifyBlockParam("Executioner", "verbose", False)
    self.input_file_modifier.modifyBlockParam("NonlinearSolver", "verbose", False)
    self.input_file_modifier.modifySubblockParam("Output", "csv", "file_name", self.test_dir + self.solution_file_name)
    self.input_file_modifier.modifySubblockParam("Output", "csv", "save_by_mesh", False)
    self.input_file_modifier.removeSubblock("Output", "plot")

    sem.run(self.input_file, self.input_file_modifier)

    csv_tester = CSVTester(self.test_dir, self.solution_file_name)
    return csv_tester.filesAreEqual()
