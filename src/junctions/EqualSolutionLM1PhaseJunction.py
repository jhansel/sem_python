from FreeBCJunction import FreeBCJunction, FreeBCJunctionParameters
from error_utilities import error

class EqualSolutionLM1PhaseJunctionParameters(FreeBCJunctionParameters):
  def __init__(self):
    FreeBCJunctionParameters.__init__(self)

## Junction that uses Lagrange multipliers to enforce solution equality between 2 meshes
class EqualSolutionLM1PhaseJunction(FreeBCJunction):
  def __init__(self, params):
    FreeBCJunction.__init__(self, params)
    if self.n_meshes != 2:
      error("Only implemented for connecting 2 meshes.")

    self.n_var = self.dof_handler.n_var
    self.n_constraints += self.n_var

  def setDoFIndices(self):
    FreeBCJunction.setDoFIndices(self)

    # create lists of DoF indices involved in constraints
    self.i1 = [self.dof_handler.i(self.node_indices[0], m) for m in xrange(self.n_var)]
    self.i2 = [self.dof_handler.i(self.node_indices[1], m) for m in xrange(self.n_var)]

  def applyWeaklyToNonlinearSystem(self, U, U_old, r, J):
    # add normal boundary fluxes
    FreeBCJunction.applyWeaklyToNonlinearSystem(self, U, U_old, r, J)

    # add contributions from Lagrange multipliers
    for m in xrange(self.n_var):
      r[self.i1[m]] += U[self.i_constraint[m]]
      J[self.i1[m]][self.i_constraint[m]] += 1

      r[self.i2[m]] -= U[self.i_constraint[m]]
      J[self.i2[m]][self.i_constraint[m]] -= 1

  def applyStronglyToNonlinearSystem(self, U, U_old, r, J):
    # apply constraints
    for m in xrange(self.n_var):
      r[self.i_constraint[m]] = U[self.i1[m]] - U[self.i2[m]]
      J[self.i_constraint[m],self.i1[m]] = 1
      J[self.i_constraint[m],self.i2[m]] = -1

  def applyStronglyToLinearSystemMatrix(self, A):
    error("Not implemented")

  def applyStronglyToLinearSystemRHSVector(self, U_old, b):
    error("Not implemented")
