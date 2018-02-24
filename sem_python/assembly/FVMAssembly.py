import numpy as np

from .Assembly import Assembly, AssemblyParameters
from ..base.enums import ModelType, VariableName
from ..utilities.assembly_utilities import initializeDerivativeData


class FVMAssemblyParameters(AssemblyParameters):

    def __init__(self, factory):
        AssemblyParameters.__init__(self, factory)
        self.registerParameter("numerical_flux", "Numerical flux")


class FVMAssembly(Assembly):

    def __init__(self, params):
        Assembly.__init__(self, params)
        self.numerical_flux = params.get("numerical_flux")

        self.n_faces = 0
        for mesh in self.meshes:
            self.n_faces += mesh.n_cell + 1

    def isLeftFace(self, l):


    # computes the steady-state residual and Jacobian without applying strong BC
    def assembleSteadyStateSystemWithoutConstraints(self, U):
        r = np.zeros(self.dof_handler.n_dof)
        J = np.zeros(shape=(self.dof_handler.n_dof, self.dof_handler.n_dof))

        # volumetric terms
        self.addVolumetricTerms(U, r, J)

        # BCs
        for bc in self.bcs:
            bc.applyWeakBC(U, r, J)

        # junctions
        for junction in self.junctions:
            junction.applyWeaklyToNonlinearSystem(U, r, J)

        return (r, J)

    def performTransientSetup(self):
        self.M = self.computeMassMatrix()

    def computeMassMatrix(self):
        M = np.zeros(shape=(self.dof_handler.n_dof, self.dof_handler.n_dof))

        self.addMassMatrixPhase(M, 0)
        if (self.model_type != ModelType.OnePhase):
            self.addMassMatrixPhase(M, 1)
        if (self.model_type == ModelType.TwoPhase):
            self.addMassMatrixVolumeFraction(M)

        return M

    def addMassMatrixPhase(self, M, phase):
        arhoA_index = self.dof_handler.variable_index[VariableName.ARhoA][phase]
        arhouA_index = self.dof_handler.variable_index[VariableName.ARhoUA][phase]
        arhoEA_index = self.dof_handler.variable_index[VariableName.ARhoEA][phase]

        for e in range(self.dof_handler.n_cell):
            dx = self.dof_handler.dx(e)

            i_arhoA = self.dof_handler.i(e, arhoA_index)
            i_arhouA = self.dof_handler.i(e, arhouA_index)
            i_arhoEA = self.dof_handler.i(e, arhoEA_index)

            M[i_arhoA, i_arhoA] = dx
            M[i_arhouA, i_arhouA] = dx
            M[i_arhoEA, i_arhoEA] = dx

    def addMassMatrixVolumeFraction(self, M):
        aA1_index = self.dof_handler.variable_index[VariableName.AA1][0]

        for e in range(self.dof_handler.n_cell):
            dx = self.dof_handler.dx(e)
            i_aA1 = self.dof_handler.i(e, aA1_index)
            M[i_aA1, i_aA1] = dx

    def assembleTransientSystem(self, U, U_old):
        M_dU = np.matmul(self.M, U - U_old)
        return (M_dU, self.M)

    ## Applies constraints to a nonlinear system solved with Newton's method
    # @param[in] U  implicit solution vector
    # @param[in] r  nonlinear system residual vector
    # @param[in] J  nonlinear system Jacobian matrix
    def applyConstraintsToNonlinearSystem(self, U, r, J):
        # BCs
        for bc in self.bcs:
            bc.applyStrongBCNonlinearSystem(U, r, J)

        # junctions
        for junction in self.junctions:
            junction.applyStronglyToNonlinearSystem(U, r, J)

    ## Applies constraints to a linear system matrix.
    #
    # This is separated from the corresponding RHS vector modification function
    # because the matrix needs to be modified only once; the RHS vector might
    # depend on time or the solution vector.
    #
    # @param[in] A      linear system matrix
    def applyConstraintsToLinearSystemMatrix(self, A):
        # BCs
        for bc in self.bcs:
            bc.applyStrongBCLinearSystemMatrix(A)

        # junctions
        for junction in self.junctions:
            junction.applyStronglyToLinearSystemMatrix(A)

    ## Applies strong constraints to a linear system RHS vector.
    # @param[in] U   solution, needed if Dirichlet values are solution-dependent
    # @param[in] b   linear system RHS vector
    def applyConstraintsToLinearSystemRHSVector(self, U, b):
        # BCs
        for bc in self.bcs:
            bc.applyStrongBCLinearSystemRHSVector(U, b)

        # junctions
        for junction in self.junctions:
            junction.applyStronglyToLinearSystemRHSVector(U, b)

    ## Computes the steady-state residual and Jacobian
    def addVolumetricTerms(self, U, r, J):
        # loop over faces
        for l in range(self.n_faces):
            data_left = self.getLeftSolutionData(U, l)
            data_right = self.getRightSolutionData(U, l)
            F = self.numerical_flux.compute(data_left, data_right)

            if self.isLeftFace(l):
                k_left = self.getLeftNode(l)
                self.dof_handler.aggregateFlux(F, k_left, True, r)
            elif self.isRightFace(l):
                k_right = self.getRightNode(l)
                self.dof_handler.aggregateFlux(F, k_right, False, r)
            else:
                k_left = self.getLeftNode(l)
                k_right = self.getRightNode(l)
                self.dof_handler.aggregateFlux(F, k_left, True, r)
                self.dof_handler.aggregateFlux(F, k_right, False, r)
