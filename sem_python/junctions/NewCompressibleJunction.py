from .Junction1Phase import Junction1Phase, Junction1PhaseParameters
from ..base.enums import ModelType
from ..utilities.error_utilities import error
from ..closures.thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeSpecificVolume, computeVelocity, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy, computeSpecificEnthalpy


class NewCompressibleJunctionParameters(Junction1PhaseParameters):

    def __init__(self, factory):
        Junction1PhaseParameters.__init__(self, factory)


## Junction that uses compressible flow assumption
class NewCompressibleJunction(Junction1Phase):

    def __init__(self, params):
        Junction1Phase.__init__(self, params)
        if self.n_meshes != 2:
            error("NewCompressibleJunction is only implemented for connecting 2 meshes.")

        self.n_constraints += 2

        self.initializeVariableVectors()

    def setDoFIndices(self):
        Junction1Phase.setDoFIndices(self)

        self.i_constraint_mass = self.i_constraint[0]
        self.i_constraint_momentum = self.i_constraint[1]

        self.i_h0_junction = self.i_constraint_mass
        self.i_s_junction = self.i_constraint_momentum

    def initializeConstraintVariables(self, U):
        # initialize the junction h0 and s with the average of all IC values for these quantities
        h0_sum = 0
        s_sum = 0
        for i in range(self.n_meshes):
            k = self.node_indices[i]

            aA1 = self.dof_handler.aA1(U, k)
            arhoA = U[self.i_arhoA[i]]
            arhouA = U[self.i_arhouA[i]]
            arhoEA = U[self.i_arhoEA[i]]

            vf, _ = computeVolumeFraction(aA1, self.A[i], self.phase, self.model_type)
            rho, _, _, _ = computeDensity(vf, arhoA, self.A[i])
            v, _ = computeSpecificVolume(rho)
            u, _, _ = computeVelocity(arhoA, arhouA)
            E, _, _ = computeSpecificTotalEnergy(arhoA, arhoEA)
            e, _, _ = computeSpecificInternalEnergy(u, E)
            p, _, _ = self.eos.p(v, e)
            s, _, _ = self.eos.s(v, e)
            h, _, _, _ = computeSpecificEnthalpy(e, p, rho)
            h0 = h + 0.5 * u**2

            h0_sum += h0
            s_sum += s

        h0_junction = h0_sum / self.n_meshes
        s_junction = s_sum / self.n_meshes

        U[self.i_h0_junction] = h0_junction
        U[self.i_s_junction] = s_junction

    def applyWeaklyToNonlinearSystem(self, U, r, J):
        self.computeFlowQuantities(U)
        self.computeFluxes(U)

        # get the junction entropy
        s_junction = U[self.i_s_junction]

        at_least_one_inlet = False
        at_least_one_outlet = False
        for i in range(self.n_meshes):
            # inlets to the junction; outlet BC
            if self.u[i] * self.nx[i] > 0:
                # flag that at least one inlet was found
                at_least_one_inlet = True

                # use the entropy from the solution
                self.addFluxes(
                    U, r, J, i, self.s[i], 0, self.ds_daA1[i], self.ds_darhoA[i],
                    self.ds_darhouA[i], self.ds_darhoEA[i])
            else:
                # flag that at least one outlet was found
                at_least_one_outlet = True

                # use the entropy from the junction
                self.addFluxes(U, r, J, i, s_junction, 1, 0, 0, 0, 0)

    def addFluxes(self, U, r, J, i, s, ds_ds_junction, ds_daA1, ds_darhoA, ds_darhouA, ds_darhoEA):
        # get the junction stagnation enthalpy
        h0_junction = U[self.i_h0_junction]

        # extract solution data
        k = self.node_indices[i]

        A = self.dof_handler.A[k]
        arhouA = U[self.i_arhouA[i]]
        arhoEA = U[self.i_arhoEA[i]]

        nx = self.nx[i]

        vf = self.vf[i]
        u = self.u[i]

        # compute the back pressure for the flux evaluations
        h = h0_junction - 0.5 * u**2
        dh_dh0_junction = 1.0
        dh_darhoA = -u * self.du_darhoA[i]
        dh_darhouA = -u * self.du_darhouA[i]

        p, dp_dh, dp_ds = self.eos.p_from_h_s(h, s)
        dp_dh0_junction = dp_dh * dh_dh0_junction
        dp_ds_junction = dp_ds * ds_ds_junction
        dp_daA1 = dp_ds * ds_daA1
        dp_darhoA = dp_dh * dh_darhoA + dp_ds * ds_darhoA
        dp_darhouA = dp_dh * dh_darhouA + dp_ds * ds_darhouA
        dp_darhoEA = dp_ds * ds_darhoEA

        # compute the fluxes
        # i_vf = self.i_aA1[i]
        i_mass = self.i_arhoA[i]
        i_momentum = self.i_arhouA[i]
        i_energy = self.i_arhoEA[i]

        r[i_mass] += arhouA * nx
        J[i_mass][i_momentum] += nx

        r[i_momentum] += (arhouA * u + vf * p * A) * nx
        J[i_momentum][self.i_h0_junction] += vf * dp_dh0_junction * A * nx
        J[i_momentum][self.i_s_junction] += vf * dp_ds_junction * A * nx
        # J[i_momentum][i_vf] += (self.dvf_daA1[i] * p + vf * dp_daA1) * A * nx
        J[i_momentum][i_mass] += (arhouA * self.du_darhoA[i] + vf * dp_darhoA * A) * nx
        J[i_momentum][i_momentum] += (arhouA * self.du_darhouA[i] + u + vf * dp_darhouA * A) * nx
        J[i_momentum][i_energy] += vf * dp_darhoEA * A * nx

        r[i_energy] += (arhoEA + vf * p * A) * u * nx
        J[i_energy][self.i_h0_junction] += vf * dp_dh0_junction * A * u * nx
        J[i_energy][self.i_s_junction] += vf * dp_ds_junction * A * u * nx
        # J[i_energy][i_vf] += (self.dvf_daA1[i] * p + vf * dp_daA1) * A * u * nx
        J[i_energy][i_mass] += (
            (arhoEA + vf * p * A) * self.du_darhoA[i] + vf * dp_darhoA * A * u) * nx
        J[i_energy][i_momentum] += (
            (arhoEA + vf * p * A) * self.du_darhouA[i] + vf * dp_darhouA * A * u) * nx
        J[i_energy][i_energy] += (1 + vf * dp_darhoEA * A) * u * nx

    def applyStronglyToNonlinearSystem(self, U, r, J):
        r[self.i_constraint_mass] = sum(self.f_mass)
        r[self.i_constraint_momentum] = sum(self.f_momentum)

        J[self.i_constraint_mass, :] = 0
        J[self.i_constraint_momentum, :] = 0

        for n in range(self.n_meshes):
            J[self.i_constraint_mass, self.i_arhouA[n]] = self.df_mass_darhouA[n]

            J[self.i_constraint_momentum, self.i_arhoA[n]] = self.df_momentum_darhoA[n]
            J[self.i_constraint_momentum, self.i_arhouA[n]] = self.df_momentum_darhouA[n]
            J[self.i_constraint_momentum, self.i_arhoEA[n]] = self.df_momentum_darhoEA[n]

            if self.model_type == ModelType.TwoPhase:
                J[self.i_constraint_momentum, self.i_aA1[n]] = self.df_momentum_daA1[n]

    def applyStronglyToLinearSystemMatrix(self, A):
        pass

    def applyStronglyToLinearSystemRHSVector(self, U, b):
        pass
