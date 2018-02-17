from .DoFHandler import DoFHandler, DoFHandlerParameters


class FEMDoFHandlerParameters(DoFHandlerParameters):

    def __init__(self, factory):
        DoFHandlerParameters.__init__(self, factory)


class FEMDoFHandler(DoFHandler):

    def __init__(self, params):
        DoFHandler.__init__(self, params)

    def computeNumberOfNodesInMesh(self, n_cells):
        return n_cells + 1

    def getNodalPositionsFromMesh(self, mesh):
        return mesh.x

    def getNumberOfDoFsPerCellPerVariable(self):
        return 2

    ## Returns global node index for an element index and local node index
    # @param[in] e  element index
    # @param[in] k_local  local node index
    def k(self, e, k_local):
        return e + k_local + self.elem_to_mesh_index[e]

    # aggregates local cell vector into global vector
    def aggregateLocalCellVector(self, r, r_cell, e):
        i_min = self.i(self.k(e, 0), 0)
        i_max = self.i(self.k(e, 1), self.n_var - 1)
        r[i_min:i_max + 1] += r_cell

    # aggregates local node vector into global vector
    def aggregateLocalNodeVector(self, r, r_node, k):
        i_min = k * self.n_var
        i_max = (k + 1) * self.n_var - 1
        r[i_min:i_max + 1] += r_node

    # aggregates local cell matrix into global matrix
    def aggregateLocalCellMatrix(self, J, J_cell, e):
        i_min = self.i(self.k(e, 0), 0)
        i_max = self.i(self.k(e, 1), self.n_var - 1)
        J[i_min:i_max + 1, i_min:i_max + 1] += J_cell

    # aggregates local node matrix into global matrix
    def aggregateLocalNodeMatrix(self, J, J_node, k):
        i_min = k * self.n_var
        i_max = (k + 1) * self.n_var - 1
        J[i_min:i_max + 1, i_min:i_max + 1] += J_node
