from .Output import Output, OutputParameters

class ScreenOutputParameters(OutputParameters):
  def __init__(self):
    OutputParameters.__init__(self)
    self.registerStringListParameter("data_names", "List of names of quantities to print")

class ScreenOutput(Output):
  def __init__(self, params):
    Output.__init__(self, params)
    self.data_names = params.get("data_names")

    self.n_data_names = len(self.data_names)
    self.header_format = "\n" + "%15s" * self.n_data_names
    self.data_format = "%15g" * self.n_data_names
    self.header_tuple = tuple(self.data_names)

  def run(self, data):
    # extract data to be printed
    self.data_list = [data[name] for name in self.data_names]

    n_node = self.dof_handler.n_node

    # print solution
    print(self.header_format % self.header_tuple)
    for k in range(n_node):
      data_value_list = [self.data_list[name][k] for name in self.data_list]
      data_tuple = tuple(data_value_list)
      print(self.data_format % data_tuple)
