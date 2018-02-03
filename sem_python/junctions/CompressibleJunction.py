from .Junction1Phase import Junction1Phase, Junction1PhaseParameters
from ..base.enums import ModelType
from ..utilities.error_utilities import error
from ..closures.thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeSpecificVolume, computeVelocity, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy, computeSpecificEnthalpy


class CompressibleJunctionParameters(Junction1PhaseParameters):

    def __init__(self, factory):
        Junction1PhaseParameters.__init__(self, factory)
        self.registerBoolParameter(
            "use_momentum_flux_balance", "Flag to use a momentum flux balance for last equation",
            False)
        self.registerBoolParameter("use_lm", "Flag to use Lagrange multipliers", False)


## Junction that uses compressible flow assumption
class CompressibleJunction(Junction1Phase):

    def __init__(self, params):
        Junction1Phase.__init__(self, params)
        if self.n_meshes != 2:
            error("CompressibleJunction is only implemented for connecting 2 meshes.")

        # determine if momentum flux balance is to be used for last equation;
        # otherwise, stagnation pressure condition will be used
        self.use_momentum_flux_balance = params.get("use_momentum_flux_balance")

        # add constraints
        self.use_lm = params.get("use_lm")
        if self.use_lm:
            self.n_constraints += 6

        # get node indices
        self.k_left = self.node_indices[0]
        self.k_right = self.node_indices[1]
        self.k_left_adjacent = self.adjacent_node_indices[0]
        self.k_right_adjacent = self.adjacent_node_indices[1]

        self.initializeVariableVectors()

    def setDoFIndices(self):
        Junction1Phase.setDoFIndices(self)

        if (self.model_type == ModelType.TwoPhase):
            self.i_aA1L = self.dof_handler.i(self.k_left, self.aA1_index)
            self.i_aA1R = self.dof_handler.i(self.k_right, self.aA1_index)
            self.i_aA1LL = self.dof_handler.i(self.k_left_adjacent, self.aA1_index)
            self.i_aA1RR = self.dof_handler.i(self.k_right_adjacent, self.aA1_index)

        self.i_arhoAL = self.dof_handler.i(self.k_left, self.arhoA_index)
        self.i_arhouAL = self.dof_handler.i(self.k_left, self.arhouA_index)
        self.i_arhoEAL = self.dof_handler.i(self.k_left, self.arhoEA_index)

        self.i_arhoAR = self.dof_handler.i(self.k_right, self.arhoA_index)
        self.i_arhouAR = self.dof_handler.i(self.k_right, self.arhouA_index)
        self.i_arhoEAR = self.dof_handler.i(self.k_right, self.arhoEA_index)

        self.i_arhoALL = self.dof_handler.i(self.k_left_adjacent, self.arhoA_index)
        self.i_arhouALL = self.dof_handler.i(self.k_left_adjacent, self.arhouA_index)
        self.i_arhoEALL = self.dof_handler.i(self.k_left_adjacent, self.arhoEA_index)

        self.i_arhoARR = self.dof_handler.i(self.k_right_adjacent, self.arhoA_index)
        self.i_arhouARR = self.dof_handler.i(self.k_right_adjacent, self.arhouA_index)
        self.i_arhoEARR = self.dof_handler.i(self.k_right_adjacent, self.arhoEA_index)

    def applyWeaklyToNonlinearSystem(self, U_new, U_old, r, J):
        if self.use_lm:
            Junction1Phase.applyWeaklyToNonlinearSystem(self, U_new, U_old, r, J)

        n_values = 6
        L = 0
        R = 1
        L_old = 2
        R_old = 3
        LL_old = 4
        RR_old = 5

        k = [0] * n_values
        k[L] = self.k_left
        k[R] = self.k_right
        k[L_old] = self.k_left
        k[R_old] = self.k_right
        k[LL_old] = self.k_left_adjacent
        k[RR_old] = self.k_right_adjacent

        i_arhoA = [0] * n_values
        i_arhoA[L] = self.i_arhoAL
        i_arhoA[R] = self.i_arhoAR
        i_arhoA[L_old] = self.i_arhoAL
        i_arhoA[R_old] = self.i_arhoAR
        i_arhoA[LL_old] = self.i_arhoALL
        i_arhoA[RR_old] = self.i_arhoARR

        i_arhouA = [0] * n_values
        i_arhouA[L] = self.i_arhouAL
        i_arhouA[R] = self.i_arhouAR
        i_arhouA[L_old] = self.i_arhouAL
        i_arhouA[R_old] = self.i_arhouAR
        i_arhouA[LL_old] = self.i_arhouALL
        i_arhouA[RR_old] = self.i_arhouARR

        i_arhoEA = [0] * n_values
        i_arhoEA[L] = self.i_arhoEAL
        i_arhoEA[R] = self.i_arhoEAR
        i_arhoEA[L_old] = self.i_arhoEAL
        i_arhoEA[R_old] = self.i_arhoEAR
        i_arhoEA[LL_old] = self.i_arhoEALL
        i_arhoEA[RR_old] = self.i_arhoEARR

        U = [0] * n_values
        U[L] = U_new
        U[R] = U_new
        U[L_old] = U_old
        U[R_old] = U_old
        U[LL_old] = U_old
        U[RR_old] = U_old

        # initialize lists
        A = [0] * n_values

        rho = [0] * n_values
        drho_daA1 = [0] * n_values
        drho_darhoA = [0] * n_values

        u = [0] * n_values
        du_darhoA = [0] * n_values
        du_darhouA = [0] * n_values

        e = [0] * n_values
        de_darhoA = [0] * n_values
        de_darhouA = [0] * n_values
        de_darhoEA = [0] * n_values

        p = [0] * n_values
        dp_daA1 = [0] * n_values
        dp_darhoA = [0] * n_values
        dp_darhouA = [0] * n_values
        dp_darhoEA = [0] * n_values

        c = [0] * n_values
        dc_daA1 = [0] * n_values
        dc_darhoA = [0] * n_values
        dc_darhouA = [0] * n_values
        dc_darhoEA = [0] * n_values

        if not self.use_momentum_flux_balance:
            p0 = [0] * n_values
            dp0_daA1 = [0] * n_values
            dp0_darhoA = [0] * n_values
            dp0_darhouA = [0] * n_values
            dp0_darhoEA = [0] * n_values

        # loop over subscript/superscript combinations
        for i in range(n_values):
            A[i] = self.dof_handler.A[k[i]]
            aA1 = self.dof_handler.aA1(U[i], k[i])
            vf, dvf_daA1 = computeVolumeFraction(aA1, A[i], self.phase, self.model_type)
            arhoA = U[i][i_arhoA[i]]
            arhouA = U[i][i_arhouA[i]]
            arhoEA = U[i][i_arhoEA[i]]

            rho[i], drho_dvf, drho_darhoA[i], _ = computeDensity(vf, arhoA, A[i])
            drho_daA1[i] = drho_dvf * dvf_daA1

            u[i], du_darhoA[i], du_darhouA[i] = computeVelocity(arhoA, arhouA)

            v, dv_drho = computeSpecificVolume(rho[i])
            dv_daA1 = dv_drho * drho_daA1[i]
            dv_darhoA = dv_drho * drho_darhoA[i]

            E, dE_darhoA, dE_darhoEA = computeSpecificTotalEnergy(arhoA, arhoEA)

            e[i], de_du, de_dE = computeSpecificInternalEnergy(u[i], E)
            de_darhoA[i] = de_du * du_darhoA[i] + de_dE * dE_darhoA
            de_darhouA[i] = de_du * du_darhouA[i]
            de_darhoEA[i] = de_dE * dE_darhoEA

            p[i], dp_dv, dp_de = self.eos.p(v, e[i])
            dp_daA1[i] = dp_dv * dv_daA1
            dp_darhoA[i] = dp_dv * dv_darhoA + dp_de * de_darhoA[i]
            dp_darhouA[i] = dp_de * de_darhouA[i]
            dp_darhoEA[i] = dp_de * de_darhoEA[i]

            c[i], dc_dv, dc_de = self.eos.c(v, e[i])
            dc_daA1[i] = dc_dv * dv_daA1
            dc_darhoA[i] = dc_dv * dv_darhoA + dc_de * de_darhoA[i]
            dc_darhouA[i] = dc_de * de_darhouA[i]
            dc_darhoEA[i] = dc_de * de_darhoEA[i]

            if not self.use_momentum_flux_balance:
                s, ds_dv, ds_de = self.eos.s(v, e[i])
                ds_daA1 = ds_dv * dv_daA1
                ds_darhoA = ds_dv * dv_darhoA + ds_de * de_darhoA[i]
                ds_darhouA = ds_de * de_darhouA[i]
                ds_darhoEA = ds_de * de_darhoEA[i]

                h, dh_de, dh_dp, dh_drho = computeSpecificEnthalpy(e[i], p[i], rho[i])
                dh_daA1 = dh_dp * dp_daA1[i]
                dh_darhoA = dh_de * de_darhoA[i] + dh_dp * dp_darhoA[i] + dh_drho * drho_darhoA[i]
                dh_darhouA = dh_de * de_darhouA[i] + dh_dp * dp_darhouA[i]
                dh_darhoEA = dh_de * de_darhoEA[i] + dh_dp * dp_darhoEA[i]

                h0 = h + 0.5 * u[i]**2
                dh0_daA1 = dh_daA1
                dh0_darhoA = dh_darhoA + u[i] * du_darhoA[i]
                dh0_darhouA = dh_darhouA + u[i] * du_darhouA[i]
                dh0_darhoEA = dh_darhoEA

                p0[i], dp0_dh0, dp0_ds = self.eos.p_from_h_s(h0, s)
                dp0_daA1[i] = dp0_dh0 * dh0_daA1 + dp0_ds * ds_daA1
                dp0_darhoA[i] = dp0_dh0 * dh0_darhoA + dp0_ds * ds_darhoA
                dp0_darhouA[i] = dp0_dh0 * dh0_darhouA + dp0_ds * ds_darhouA
                dp0_darhoEA[i] = dp0_dh0 * dh0_darhoEA + dp0_ds * ds_darhoEA

        # compute old average quantities
        rhoL = 0.5 * (rho[L_old] + rho[LL_old])
        rhoR = 0.5 * (rho[R_old] + rho[RR_old])
        uL = 0.5 * (u[L_old] + u[LL_old])
        uR = 0.5 * (u[R_old] + u[RR_old])
        pL = 0.5 * (p[L_old] + p[LL_old])
        pR = 0.5 * (p[R_old] + p[RR_old])
        cL = 0.5 * (c[L_old] + c[LL_old])
        cR = 0.5 * (c[R_old] + c[RR_old])

        # compute residuals and Jacobians
        self.r_i1 = p[L] - pL + rhoL * (cL - uL) * (u[L] - uL)
        self.J_i1_vf1L = dp_daA1[L]
        self.J_i1_arhoAL = dp_darhoA[L] + rhoL * (cL - uL) * du_darhoA[L]
        self.J_i1_arhouAL = dp_darhouA[L] + rhoL * (cL - uL) * du_darhouA[L]
        self.J_i1_arhoEAL = dp_darhoEA[L]

        self.r_i2 = p[R] - pR - rhoR * (cR - uR) * (u[R] - uR)
        self.J_i2_vf1R = dp_daA1[R]
        self.J_i2_arhoAR = dp_darhoA[R] - rhoR * (cR - uR) * du_darhoA[R]
        self.J_i2_arhouAR = dp_darhouA[R] - rhoR * (cR - uR) * du_darhouA[R]
        self.J_i2_arhoEAR = dp_darhoEA[R]

        if u[L] >= 0:
            if u[R] < 0:
                error("Assumption violated: Both velocity conditions were true.")
            self.r_i3 = rho[L] - rhoL - (p[L] - pL) / cL**2
            self.J_i3_vf1L = drho_daA1[L] - dp_daA1[L] / cL**2
            self.J_i3_arhoAL = drho_darhoA[L] - dp_darhoA[L] / cL**2
            self.J_i3_arhouAL = -dp_darhouA[L] / cL**2
            self.J_i3_arhoEAL = -dp_darhoEA[L] / cL**2
            self.J_i3_vf1R = 0
            self.J_i3_arhoAR = 0
            self.J_i3_arhouAR = 0
            self.J_i3_arhoEAR = 0
        elif u[R] < 0:
            self.r_i3 = rho[R] - rhoR - (p[R] - pR) / cR**2
            self.J_i3_vf1R = drho_daA1[R] - dp_daA1[R] / cR**2
            self.J_i3_arhoAR = drho_darhoA[R] - dp_darhoA[R] / cR**2
            self.J_i3_arhouAR = -dp_darhouA[R] / cR**2
            self.J_i3_arhoEAR = -dp_darhoEA[R] / cR**2
            self.J_i3_vf1L = 0
            self.J_i3_arhoAL = 0
            self.J_i3_arhouAL = 0
            self.J_i3_arhoEAL = 0
        else:
            error("Assumption violated: Neither velocity condition was true.")

        self.r_i4 = rho[L] * u[L] * A[L] - rho[R] * u[R] * A[R]
        self.J_i4_vf1L = drho_daA1[L] * u[L] * A[L]
        self.J_i4_arhoAL = (drho_darhoA[L] * u[L] + rho[L] * du_darhoA[L]) * A[L]
        self.J_i4_arhouAL = rho[L] * du_darhouA[L] * A[L]
        self.J_i4_vf1R = -drho_daA1[R] * u[R] * A[R]
        self.J_i4_arhoAR = -(drho_darhoA[R] * u[R] + rho[R] * du_darhoA[R]) * A[R]
        self.J_i4_arhouAR = -rho[R] * du_darhouA[R] * A[R]

        self.r_i5 = (e[L] + p[L] / rho[L] + 0.5 * u[L]**2) * A[L] - (
            e[R] + p[R] / rho[R] + 0.5 * u[R]**2) * A[R]
        self.J_i5_vf1L = (dp_daA1[L] / rho[L] - p[L] / rho[L]**2 * drho_daA1[L]) * A[L]
        self.J_i5_arhoAL = (
            de_darhoA[L] + dp_darhoA[L] / rho[L] - p[L] / rho[L]**2 * drho_darhoA[L] +
            u[L] * du_darhoA[L]) * A[L]
        self.J_i5_arhouAL = (de_darhouA[L] + dp_darhouA[L] / rho[L] + u[L] * du_darhouA[L]) * A[L]
        self.J_i5_arhoEAL = (de_darhoEA[L] + dp_darhoEA[L] / rho[L]) * A[L]
        self.J_i5_vf1R = -(dp_daA1[R] / rho[R] - p[R] / rho[R]**2 * drho_daA1[R]) * A[R]
        self.J_i5_arhoAR = -(
            de_darhoA[R] + dp_darhoA[R] / rho[R] - p[R] / rho[R]**2 * drho_darhoA[R] +
            u[R] * du_darhoA[R]) * A[R]
        self.J_i5_arhouAR = -(de_darhouA[R] + dp_darhouA[R] / rho[R] + u[R] * du_darhouA[R]) * A[R]
        self.J_i5_arhoEAR = -(de_darhoEA[R] + dp_darhoEA[R] / rho[R]) * A[R]

        if self.use_momentum_flux_balance:
            self.r_i6 = (rho[L] * u[L]**2 + p[L]) * A[L] - (rho[R] * u[R]**2 + p[R]) * A[R]
            self.J_i6_vf1L = (drho_daA1[L] * u[L]**2 + dp_daA1[L]) * A[L]
            self.J_i6_arhoAL = (
                drho_darhoA[L] * u[L]**2 + rho[L] * 2.0 * u[L] * du_darhoA[L] + dp_darhoA[L]) * A[L]
            self.J_i6_arhouAL = (rho[L] * 2.0 * u[L] * du_darhouA[L] + dp_darhouA[L]) * A[L]
            self.J_i6_arhoEAL = dp_darhoEA[L] * A[L]
            self.J_i6_vf1R = -(drho_daA1[R] * u[R]**2 + dp_daA1[R]) * A[R]
            self.J_i6_arhoAR = -(
                drho_darhoA[R] * u[R]**2 + rho[R] * 2.0 * u[R] * du_darhoA[R] + dp_darhoA[R]) * A[R]
            self.J_i6_arhouAR = -(rho[R] * 2.0 * u[R] * du_darhouA[R] + dp_darhouA[R]) * A[R]
            self.J_i6_arhoEAR = -dp_darhoEA[R] * A[R]
        else:
            self.r_i6 = p0[L] - p0[R]
            self.J_i6_vf1L = dp0_daA1[L]
            self.J_i6_arhoAL = dp0_darhoA[L]
            self.J_i6_arhouAL = dp0_darhouA[L]
            self.J_i6_arhoEAL = dp0_darhoEA[L]
            self.J_i6_vf1R = -dp0_daA1[R]
            self.J_i6_arhoAR = -dp0_darhoA[R]
            self.J_i6_arhouAR = -dp0_darhouA[R]
            self.J_i6_arhoEAR = -dp0_darhoEA[R]

        # residual indices
        if self.use_lm:
            i1 = self.i_constraint[0]
            i2 = self.i_constraint[1]
            i3 = self.i_constraint[2]
            i4 = self.i_constraint[3]
            i5 = self.i_constraint[4]
            i6 = self.i_constraint[5]
        else:
            i1 = self.i_arhoAL
            i2 = self.i_arhouAL
            i3 = self.i_arhoEAL
            i4 = self.i_arhoAR
            i5 = self.i_arhouAR
            i6 = self.i_arhoEAR

        # create list of residual entries
        self.residual_entries = list()
        self.residual_entries.append((i1, self.r_i1))
        self.residual_entries.append((i2, self.r_i2))
        self.residual_entries.append((i3, self.r_i3))
        self.residual_entries.append((i4, self.r_i4))
        self.residual_entries.append((i5, self.r_i5))
        self.residual_entries.append((i6, self.r_i6))

        # create list of tuples of the Jacobian entries
        self.jacobian_entries = list()

        self.jacobian_entries.append((i1, self.i_arhoAL, self.J_i1_arhoAL))
        self.jacobian_entries.append((i1, self.i_arhouAL, self.J_i1_arhouAL))
        self.jacobian_entries.append((i1, self.i_arhoEAL, self.J_i1_arhoEAL))

        self.jacobian_entries.append((i2, self.i_arhoAR, self.J_i2_arhoAR))
        self.jacobian_entries.append((i2, self.i_arhouAR, self.J_i2_arhouAR))
        self.jacobian_entries.append((i2, self.i_arhoEAR, self.J_i2_arhoEAR))

        self.jacobian_entries.append((i3, self.i_arhoAL, self.J_i3_arhoAL))
        self.jacobian_entries.append((i3, self.i_arhouAL, self.J_i3_arhouAL))
        self.jacobian_entries.append((i3, self.i_arhoEAL, self.J_i3_arhoEAL))
        self.jacobian_entries.append((i3, self.i_arhoAR, self.J_i3_arhoAR))
        self.jacobian_entries.append((i3, self.i_arhouAR, self.J_i3_arhouAR))
        self.jacobian_entries.append((i3, self.i_arhoEAR, self.J_i3_arhoEAR))

        self.jacobian_entries.append((i4, self.i_arhoAL, self.J_i4_arhoAL))
        self.jacobian_entries.append((i4, self.i_arhouAL, self.J_i4_arhouAL))
        self.jacobian_entries.append((i4, self.i_arhoAR, self.J_i4_arhoAR))
        self.jacobian_entries.append((i4, self.i_arhouAR, self.J_i4_arhouAR))

        self.jacobian_entries.append((i5, self.i_arhoAL, self.J_i5_arhoAL))
        self.jacobian_entries.append((i5, self.i_arhouAL, self.J_i5_arhouAL))
        self.jacobian_entries.append((i5, self.i_arhoEAL, self.J_i5_arhoEAL))
        self.jacobian_entries.append((i5, self.i_arhoAR, self.J_i5_arhoAR))
        self.jacobian_entries.append((i5, self.i_arhouAR, self.J_i5_arhouAR))
        self.jacobian_entries.append((i5, self.i_arhoEAR, self.J_i5_arhoEAR))

        self.jacobian_entries.append((i6, self.i_arhoAL, self.J_i6_arhoAL))
        self.jacobian_entries.append((i6, self.i_arhouAL, self.J_i6_arhouAL))
        self.jacobian_entries.append((i6, self.i_arhoEAL, self.J_i6_arhoEAL))
        self.jacobian_entries.append((i6, self.i_arhoAR, self.J_i6_arhoAR))
        self.jacobian_entries.append((i6, self.i_arhouAR, self.J_i6_arhouAR))
        self.jacobian_entries.append((i6, self.i_arhoEAR, self.J_i6_arhoEAR))

        if self.model_type == ModelType.TwoPhase:
            self.jacobian_entries.append((i1, self.i_aA1L, self.J_i1_vf1L))
            self.jacobian_entries.append((i2, self.i_aA1R, self.J_i2_vf1R))
            self.jacobian_entries.append((i3, self.i_aA1L, self.J_i3_vf1L))
            self.jacobian_entries.append((i3, self.i_aA1R, self.J_i3_vf1R))
            self.jacobian_entries.append((i4, self.i_aA1L, self.J_i4_vf1L))
            self.jacobian_entries.append((i4, self.i_aA1R, self.J_i4_vf1R))
            self.jacobian_entries.append((i5, self.i_aA1L, self.J_i5_vf1L))
            self.jacobian_entries.append((i5, self.i_aA1R, self.J_i5_vf1R))
            self.jacobian_entries.append((i6, self.i_aA1L, self.J_i6_vf1L))
            self.jacobian_entries.append((i6, self.i_aA1R, self.J_i6_vf1R))

        # if using Lagrange multipliers, apply contributions to residual and Jacobian
        if self.use_lm:
            for entry in self.jacobian_entries:
                i, j, value = entry
                r[j] += U_new[i] * value
                J[j, i] += value

    def applyStronglyToNonlinearSystem(self, U_new, U_old, r, J):
        # set constraint residual entries and zero-out Jacobian rows
        for entry in self.residual_entries:
            i, value = entry
            r[i] = value
            J[i, :] = 0

        # set constraint Jacobian entries
        for entry in self.jacobian_entries:
            i, j, value = entry
            J[i, j] = value

    def applyStronglyToLinearSystemMatrix(self, A):
        pass

    def applyStronglyToLinearSystemRHSVector(self, U_old, b):
        pass
