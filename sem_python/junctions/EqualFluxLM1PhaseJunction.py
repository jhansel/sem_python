from .Junction1Phase import Junction1Phase, Junction1PhaseParameters
from ..base.enums import ModelType


class EqualFluxLM1PhaseJunctionParameters(Junction1PhaseParameters):

    def __init__(self, factory):
        Junction1PhaseParameters.__init__(self, factory)


## Junction that enforces that the fluxes are equal at the junction
#
#  The first mesh will have its residuals replaced strongly, and the second
#  mesh will add a boundary flux computed with its solution.
class EqualFluxLM1PhaseJunction(Junction1Phase):

    def __init__(self, params):
        Junction1Phase.__init__(self, params)

        self.n_var = self.dof_handler.n_var
        self.n_constraints += self.n_var

        self.initializeVariableVectors()

    def setDoFIndices(self):
        Junction1Phase.setDoFIndices(self)

        self.i_constraint_mass = self.i_constraint[0]
        self.i_constraint_momentum = self.i_constraint[1]
        self.i_constraint_energy = self.i_constraint[2]

    def applyWeaklyToNonlinearSystem(self, U, r, J):
        # add normal boundary fluxes
        Junction1Phase.applyWeaklyToNonlinearSystem(self, U, r, J)

        # add contributions from Lagrange multipliers
        for n in range(self.n_meshes):
            r[self.i_arhoA[n]] += U[self.i_constraint_momentum] * self.df_momentum_darhoA[n]
            r[self.i_arhoA[n]] += U[self.i_constraint_energy] * self.df_energy_darhoA[n]
            J[self.i_arhoA[n], self.i_constraint_momentum] += self.df_momentum_darhoA[n]
            J[self.i_arhoA[n], self.i_constraint_energy] += self.df_energy_darhoA[n]

            r[self.i_arhouA[n]] += U[self.i_constraint_mass] * self.df_mass_darhouA[n]
            r[self.i_arhouA[n]] += U[self.i_constraint_momentum] * self.df_momentum_darhouA[n]
            r[self.i_arhouA[n]] += U[self.i_constraint_energy] * self.df_energy_darhouA[n]
            J[self.i_arhouA[n], self.i_constraint_mass] += self.df_mass_darhouA[n]
            J[self.i_arhouA[n], self.i_constraint_momentum] += self.df_momentum_darhouA[n]
            J[self.i_arhouA[n], self.i_constraint_energy] += self.df_energy_darhouA[n]

            r[self.i_arhoEA[n]] += U[self.i_constraint_momentum] * self.df_momentum_darhoEA[n]
            r[self.i_arhoEA[n]] += U[self.i_constraint_energy] * self.df_energy_darhoEA[n]
            J[self.i_arhoEA[n], self.i_constraint_momentum] += self.df_momentum_darhoEA[n]
            J[self.i_arhoEA[n], self.i_constraint_energy] += self.df_energy_darhoEA[n]

    def applyStronglyToNonlinearSystem(self, U, r, J):
        r[self.i_constraint_mass] = sum(self.f_mass)
        r[self.i_constraint_momentum] = sum(self.f_momentum)
        r[self.i_constraint_energy] = sum(self.f_energy)

        J[self.i_constraint_mass, :] = 0
        J[self.i_constraint_momentum, :] = 0
        J[self.i_constraint_energy, :] = 0

        for n in range(self.n_meshes):
            J[self.i_constraint_mass, self.i_arhouA[n]] = self.df_mass_darhouA[n]

            J[self.i_constraint_momentum, self.i_arhoA[n]] = self.df_momentum_darhoA[n]
            J[self.i_constraint_momentum, self.i_arhouA[n]] = self.df_momentum_darhouA[n]
            J[self.i_constraint_momentum, self.i_arhoEA[n]] = self.df_momentum_darhoEA[n]

            J[self.i_constraint_energy, self.i_arhoA[n]] = self.df_energy_darhoA[n]
            J[self.i_constraint_energy, self.i_arhouA[n]] = self.df_energy_darhouA[n]
            J[self.i_constraint_energy, self.i_arhoEA[n]] = self.df_energy_darhoEA[n]

            if self.model_type == ModelType.TwoPhase:
                J[self.i_constraint_momentum, self.i_aA1[n]] = self.df_momentum_daA1[n]
                J[self.i_constraint_energy, self.i_aA1[n]] = self.df_energy_daA1[n]
