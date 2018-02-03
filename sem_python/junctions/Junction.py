from ..input.Parameters import Parameters
from ..utilities.error_utilities import error


class JunctionParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerStringListParameter(
            "mesh_names", "List of names of the meshes to be connected")
        self.registerStringListParameter(
            "mesh_sides", "List of sides of the meshes to be connected")
        self.registerParameter("dof_handler", "Degree of freedom handler")
        self.registerParameter("eos_list", "List of equations of state")


## Base class for a junction between meshes
class Junction(object):

    def __init__(self, params):
        self.mesh_names = params.get("mesh_names")
        self.mesh_sides = params.get("mesh_sides")
        self.dof_handler = params.get("dof_handler")
        self.eos_list = params.get("eos_list")

        self.model_type = self.dof_handler.model_type

        # ensure that lengths of list parameters agree
        self.n_meshes = len(self.mesh_names)
        if len(self.mesh_names) != len(self.mesh_sides):
            error("The list parameters 'mesh_names' and 'mesh_sides' must have the same size.")

        # get node indices and adjacent node indices
        self.node_indices = list()
        self.adjacent_node_indices = list()
        for i, mesh_name in enumerate(self.mesh_names):
            if self.mesh_sides[i] == "left":
                self.node_indices.append(self.dof_handler.getNodeIndexFromLeft(mesh_name, 0))
                self.adjacent_node_indices.append(
                    self.dof_handler.getNodeIndexFromLeft(mesh_name, 1))
            elif self.mesh_sides[i] == "right":
                self.node_indices.append(self.dof_handler.getNodeIndexFromRight(mesh_name, 0))
                self.adjacent_node_indices.append(
                    self.dof_handler.getNodeIndexFromRight(mesh_name, 1))
            else:
                error("Side parameters must be either 'left' or 'right'.")

        # get normal vectors and areas
        self.nx = list()
        self.A = [self.dof_handler.A[k] for k in self.node_indices]
        for i in range(self.n_meshes):
            if self.mesh_sides[i] == "left":
                self.nx.append(-1.0)
            else:
                self.nx.append(1.0)

        # initialize number of constraints to zero
        self.n_constraints = 0

        self.i_constraint = None

    def setConstraintDoFIndices(self, constraint_dof_indices):
        self.i_constraint = constraint_dof_indices

        # call derived class functions that query DoF indices; it is incorrect to
        # make this call in the constructor because the DoF handler requires
        # construction of the junctions before it can know the DoF indices
        self.setDoFIndices()

    def setDoFIndices(self):
        pass

    def initializeConstraintVariables(self, U):
        pass

    def applyWeaklyToNonlinearSystem(self, U, U_old, r, J):
        pass

    def applyStronglyToNonlinearSystem(self, U, U_old, r, J):
        pass

    def applyStronglyToLinearSystemMatrix(self, A):
        pass

    def applyStronglyToLinearSystemRHSVector(self, U_old, b):
        pass
