from collections import OrderedDict

from .Output import Output, OutputParameters
from ..utilities.file_utilities import writeCSVFile
from ..utilities.string_utilities import stripStringFromRight

class CSVOutputParameters(OutputParameters):
  def __init__(self):
    OutputParameters.__init__(self)
    self.registerStringParameter("file_name", "Name of solution output file", "solution.csv")
    self.registerStringListParameter("names", "List of names of quantities to save")
    self.registerIntParameter("output_precision", "Precision used in solution output file", 5)
    self.registerBoolParameter("save_by_mesh", "Flag to save data in separate output files for each mesh", False)

class CSVOutput(Output):
  def __init__(self, params):
    Output.__init__(self, params)
    self.file_name = params.get("file_name")
    self.names = params.get("names")
    self.output_precision = params.get("output_precision")
    self.save_by_mesh = params.get("save_by_mesh")

  def run(self, data):
    if self.save_by_mesh:
      save_data = list()
      for i_mesh in range(self.dof_handler.n_meshes):
        save_data.append(OrderedDict())

      x = self.dof_handler.separateNodalQuantityByMesh(data["x"])
      for i_mesh in range(self.dof_handler.n_meshes):
        save_data[i_mesh]["x"] = x[i_mesh]

      for name in self.names:
        data_name = self.dof_handler.separateNodalQuantityByMesh(data[name])
        for i_mesh in range(self.dof_handler.n_meshes):
          save_data[i_mesh][name] = data_name[i_mesh]

      for i_mesh in range(self.dof_handler.n_meshes):
        mesh_file_name = stripStringFromRight(self.file_name, ".csv") + "_" + str(i_mesh+1) + ".csv"
        writeCSVFile(save_data[i_mesh], mesh_file_name, self.output_precision)
    else:
      save_data = OrderedDict()
      save_data["x"] = data["x"]
      for name in self.names:
        save_data[name] = data[name]

      writeCSVFile(save_data, self.file_name, self.output_precision)
