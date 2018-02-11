from ..base.enums import ModelType
from .OnePhaseBC import OnePhaseBC, OnePhaseBCParameters
from ..closures.thermodynamic_functions import computeVolumeFraction, computeSpecificVolume, \
  computeVelocity


class InletP0T0BCParameters(OnePhaseBCParameters):

    def __init__(self, factory):
        OnePhaseBCParameters.__init__(self, factory)
        self.registerFloatParameter("p0", "Specified stagnation pressure")
        self.registerFloatParameter("T0", "Specified stagnation temperature")
        self.registerBoolParameter("strongly_enforce_mass", "Strongly enforce mass?", True)
        self.registerBoolParameter("strongly_enforce_momentum", "Strongly enforce momentum?", False)
        self.registerBoolParameter("strongly_enforce_energy", "Strongly enforce energy?", False)


class InletP0T0BC(OnePhaseBC):

    def __init__(self, params):
        OnePhaseBC.__init__(self, params)
        self.p0 = params.get("p0")
        self.T0 = params.get("T0")
        self.strongly_enforce_mass = params.get("strongly_enforce_mass")
        self.strongly_enforce_momentum = params.get("strongly_enforce_momentum")
        self.strongly_enforce_energy = params.get("strongly_enforce_energy")

    def applyWeakBC(self, U, r, J):
        self.computeQuantities(U)

        A = self.dof_handler.A[self.k]
        aA1 = self.dof_handler.aA1(U, self.k)
        vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)

        # mass
        if not self.strongly_enforce_mass:
            r[self.i_arhoA] += vf * self.rho * self.u * A * self.nx
            if self.model_type == ModelType.TwoPhase:
                J[self.i_arhoA, self.i_aA1] += dvf_daA1 * self.rho * self.u * A * self.nx
            J[self.i_arhoA, self.i_arhoA] += vf * (
                self.drho_darhoA * self.u + self.rho * self.du_darhoA) * A * self.nx
            J[self.i_arhoA, self.i_arhouA] += vf * (
                self.drho_darhouA * self.u + self.rho * self.du_darhouA) * A * self.nx

        # momentum
        if not self.strongly_enforce_momentum:
            r[self.i_arhouA] += vf * (self.rho * self.u**2 + self.p) * A * self.nx
            if (self.model_type == ModelType.TwoPhase):
                J[self.i_arhouA, self.i_aA1] += dvf_daA1 * (
                    self.rho * self.u**2 + self.p) * A * self.nx
            J[self.i_arhouA, self.i_arhoA] += vf * (
                self.drho_darhoA * self.u**2 + self.rho * 2 * self.u * self.du_darhoA +
                self.dp_darhoA) * A * self.nx
            J[self.i_arhouA, self.i_arhouA] += vf * (
                self.drho_darhouA * self.u**2 + self.rho * 2 * self.u * self.du_darhouA +
                self.dp_darhouA) * A * self.nx

        # energy
        if not self.strongly_enforce_energy:
            r[self.i_arhoEA] += self.u * vf * (self.rho * self.E + self.p) * A * self.nx
            if (self.model_type == ModelType.TwoPhase):
                J[self.i_arhoEA, self.i_aA1] += self.u * dvf_daA1 * (
                    self.rho * self.E + self.p) * A * self.nx
            J[self.i_arhoEA,self.i_arhoA] += vf * (self.du_darhoA * (self.rho * self.E + self.p) \
              + self.u * (self.drho_darhoA * self.E + self.rho * self.dE_darhoA + self.dp_darhoA)) * A * self.nx
            J[self.i_arhoEA,self.i_arhouA] += vf * (self.du_darhouA * (self.rho * self.E + self.p) \
              + self.u * (self.drho_darhouA * self.E + self.rho * self.dE_darhouA + self.dp_darhouA)) * A * self.nx

    def applyStrongBCNonlinearSystem(self, U, r, J):
        A = self.dof_handler.A[self.k]
        aA1 = self.dof_handler.aA1(U, self.k)
        vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)

        # mass
        if self.strongly_enforce_mass:
            arhoA = U[self.i_arhoA]

            arhoABC = vf * self.rho * A
            darhoABC_daA1 = dvf_daA1 * self.rho * A
            darhoABC_darhoA = vf * self.drho_darhoA * A
            darhoABC_darhouA = vf * self.drho_darhouA * A

            r[self.i_arhoA] = arhoA - arhoABC
            J[self.i_arhoA, :] = 0
            if (self.model_type == ModelType.TwoPhase):
                J[self.i_arhoA, self.i_aA1] = -darhoABC_daA1
            J[self.i_arhoA, self.i_arhoA] = 1 - darhoABC_darhoA
            J[self.i_arhoA, self.i_arhouA] = -darhoABC_darhouA

        # momentum
        if self.strongly_enforce_momentum:
            arhouA = U[self.i_arhouA]

            arhouABC = vf * (self.rho * self.u**2 + self.p) * A
            darhouABC_daA1 = dvf_daA1 * (self.rho * self.u**2 + self.p) * A
            darhouABC_darhoA = vf * (
                self.drho_darhoA * self.u**2 + self.rho * 2 * self.u * self.du_darhoA +
                self.dp_darhoA) * A
            darhouABC_darhouA = vf * (
                self.drho_darhouA * self.u**2 + self.rho * 2 * self.u * self.du_darhouA +
                self.dp_darhouA) * A

            r[self.i_arhouA] = arhouA - arhouABC
            J[self.i_arhouA, :] = 0
            if self.model_type == ModelType.TwoPhase:
                J[self.i_arhouA, self.i_aA1] = -darhouABC_daA1
            J[self.i_arhouA, self.i_arhoA] = -darhouABC_darhoA
            J[self.i_arhouA, self.i_arhouA] = -darhouABC_darhouA

        # energy
        if self.strongly_enforce_energy:
            arhoEA = U[self.i_arhoEA]

            arhoEABC = vf * self.rho * self.E * A
            darhoEABC_daA1 = dvf_daA1 * self.rho * self.E * A
            darhoEABC_darhoA = vf * (self.drho_darhoA * self.E + self.rho * self.dE_darhoA) * A
            darhoEABC_darhouA = vf * (self.drho_darhouA * self.E + self.rho * self.dE_darhouA) * A

            r[self.i_arhoEA] = arhoEA - arhoEABC
            J[self.i_arhoEA, :] = 0
            if self.model_type == ModelType.TwoPhase:
                J[self.i_arhoEA, self.i_aA1] = -darhoEABC_daA1
            J[self.i_arhoEA, self.i_arhoA] = -darhoEABC_darhoA
            J[self.i_arhoEA, self.i_arhouA] = -darhoEABC_darhouA
            J[self.i_arhoEA, self.i_arhoEA] = 1

    def applyStrongBCLinearSystemMatrix(self, A):
        pass

    def applyStrongBCLinearSystemRHSVector(self, U, b):
        pass

    ## Computes rho, u, p, and E, which are needed by the fluxes
    def computeQuantities(self, U):
        rho0, _, _ = self.eos.rho(self.p0, self.T0)
        v0, _ = computeSpecificVolume(rho0)
        h0, _, _ = self.eos.h(self.p0, self.T0)
        e0, _, _ = self.eos.e(v0, self.p0)
        s0, _, _ = self.eos.s(v0, e0)

        arhoA = U[self.i_arhoA]
        arhouA = U[self.i_arhouA]
        self.u, self.du_darhoA, self.du_darhouA = computeVelocity(arhoA, arhouA)

        h = h0 - 0.5 * self.u**2
        dh_darhoA = -self.u * self.du_darhoA
        dh_darhouA = -self.u * self.du_darhouA

        self.p, dp_dh, _ = self.eos.p_from_h_s(h, s0)
        self.dp_darhoA = dp_dh * dh_darhoA
        self.dp_darhouA = dp_dh * dh_darhouA

        self.rho, drho_dp, _ = self.eos.rho_from_p_s(self.p, s0)
        self.drho_darhoA = drho_dp * self.dp_darhoA
        self.drho_darhouA = drho_dp * self.dp_darhouA

        v, dv_drho = computeSpecificVolume(self.rho)
        dv_darhoA = dv_drho * self.drho_darhoA
        dv_darhouA = dv_drho * self.drho_darhouA

        e, de_dv, de_dp = self.eos.e(v, self.p)
        de_darhoA = de_dv * dv_darhoA + de_dp * self.dp_darhoA
        de_darhouA = de_dv * dv_darhouA + de_dp * self.dp_darhouA

        self.E = e + 0.5 * self.u**2
        self.dE_darhoA = de_darhoA + self.u * self.du_darhoA
        self.dE_darhouA = de_darhouA + self.u * self.du_darhouA
