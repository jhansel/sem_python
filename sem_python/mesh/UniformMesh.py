import numpy as np

from .Mesh import Mesh, MeshParameters


class UniformMeshParameters(MeshParameters):

    def __init__(self, factory):
        MeshParameters.__init__(self, factory)
        self.registerIntParameter("n_cell", "Number of cells in the mesh")
        self.registerFloatParameter("length", "Length of the spatial domain", 1.0)


class UniformMesh(Mesh):

    def __init__(self, params):
        Mesh.__init__(self, params)
        self.n_cell = params.get("n_cell")
        self.L = params.get("length")
        h = self.L / self.n_cell
        self.h = [h for e in range(self.n_cell)]

        # vertex positions
        self.x_vertices = np.zeros(self.n_cell + 1)
        self.x_vertices[0] = self.start[0]
        nx = self.orientation[0]
        for e in range(self.n_cell):
            self.x_vertices[e + 1] = self.x_vertices[e] + self.h[e] * nx

        # cell center positions
        self.x_centers = np.zeros(self.n_cell)
        for e in range(self.n_cell):
            self.x_centers[e] = 0.5 * (self.x_vertices[e] + self.x_vertices[e+1])

    def getMinimumCellWidth(self):
        return self.L / self.n_cell
