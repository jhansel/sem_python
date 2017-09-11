from collections import OrderedDict

from Output import Output, OutputParameters
from file_utilities import writeCSVFile

class CSVOutputParameters(OutputParameters):
  def __init__(self):
    OutputParameters.__init__(self)
    self.registerStringParameter("file_name", "Name of solution output file", "solution.csv")
    self.registerStringListParameter("names", "List of names of quantities to save")
    self.registerIntParameter("output_precision", "Precision used in solution output file", 5)

class CSVOutput(Output):
  def __init__(self, params):
    Output.__init__(self, params)
    self.file_name = params.get("file_name")
    self.names = params.get("names")
    self.output_precision = params.get("output_precision")

  def run(self, data):
    save_data = OrderedDict()
    save_data["x"] = data["x"]
    for name in self.names:
      save_data[name] = data[name]

    writeCSVFile(save_data, self.file_name, self.output_precision)
