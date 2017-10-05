import numpy as np

from .Mesh import Mesh, MeshParameters

class UniformMeshParameters(MeshParameters):
  def __init__(self):
    MeshParameters.__init__(self)
    self.registerIntParameter("n_cell", "Number of cells in the mesh")
    self.registerFloatParameter("length", "Length of the spatial domain", 1.0)

class UniformMesh(Mesh):
  def __init__(self, params):
    Mesh.__init__(self, params)
    self.n_cell = params.get("n_cell")
    self.n_node = self.n_cell + 1
    self.L = params.get("length")
    h = self.L / self.n_cell
    self.h = [h for e in xrange(self.n_cell)]

    self.x = np.zeros(self.n_cell + 1)
    self.x[0] = self.start[0]
    nx = self.orientation[0]
    for e in xrange(self.n_cell):
      self.x[e + 1] = self.x[e] + self.h[e] * nx

  def getMinimumCellWidth(self):
    return self.L / self.n_cell
