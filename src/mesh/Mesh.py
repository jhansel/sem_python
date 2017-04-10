import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

class MeshParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)

class Mesh(object):
  pass
