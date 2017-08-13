import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/mesh")
from Mesh import Mesh, MeshParameters

class UniformMeshParameters(MeshParameters):
  def __init__(self):
    MeshParameters.__init__(self)
    self.registerIntParameter("n_cell", "Number of cells in the mesh")
    self.registerFloatParameter("length", "Length of the spatial domain", 1.0)
    self.registerFloatParameter("x_min", "Position of left edge of domain", 0.0)

class UniformMesh(Mesh):
  def __init__(self, params):
    Mesh.__init__(self, params)
    self.n_cell = params.get("n_cell")
    self.L = params.get("length")
    self.x_min = params.get("x_min")
    h = self.L / self.n_cell
    self.h = [h for e in xrange(self.n_cell)]
    self.x = np.zeros(self.n_cell + 1)
    self.x_max = self.x_min + self.L
    self.x[0] = self.x_min
    for e in xrange(self.n_cell):
      self.x[e + 1] = self.x[e] + self.h[e]

  def getMinimumCellWidth(self):
    return self.L / self.n_cell
