import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

class ModelParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    model_selection_list = ["1phase", "2phase_noninteracting", "2phase"]
    self.registerStringSelectionParameter("model", model_selection_list, "Model type")
