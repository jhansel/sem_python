from .DoFHandler import DoFHandler, DoFHandlerParameters


class FVMDoFHandlerParameters(DoFHandlerParameters):

    def __init__(self, factory):
        DoFHandlerParameters.__init__(self, factory)


class FVMDoFHandler(DoFHandler):

    def __init__(self, params):
        DoFHandler.__init__(self, params)

    def computeNumberOfNodesInMesh(self, n_cells):
        return n_cells

    def getNodalPositionsFromMesh(self, mesh):
        return mesh.x_centers

    def getNumberOfDoFsPerCellPerVariable(self):
        return 1

    def aggregateFlux(self, F, k, left, r):
        i_min = self.i(k, 0)
        i_max = self.i(k, self.n_var - 1)
        if left:
            r[i_min:i_max + 1] -= F
        else:
            r[i_min:i_max + 1] += F
