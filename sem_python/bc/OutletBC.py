from ..base.enums import ModelType
from .OnePhaseBC import OnePhaseBC, OnePhaseBCParameters
from ..closures.thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeVelocity, computeSpecificVolume, addKineticEnergy


## Parameters class for OutletBC
class OutletBCParameters(OnePhaseBCParameters):

    def __init__(self, factory):
        OnePhaseBCParameters.__init__(self, factory)
        self.registerFloatParameter("p", "specified outlet pressure")
        self.registerBoolParameter("strongly_enforce_energy", "Strongly enforce energy?", True)


class OutletBC(OnePhaseBC):

    def __init__(self, params):
        OnePhaseBC.__init__(self, params)
        self.p = params.get("p")
        self.strongly_enforce_energy = params.get("strongly_enforce_energy")

    def computeQuantities(self, U):
        aA1 = self.dof_handler.aA1(U, self.k)
        self.vf, self.dvf_daA1 = computeVolumeFraction(aA1, self.A, self.phase, self.model_type)

        arhoA = U[self.i_arhoA]
        arhouA = U[self.i_arhouA]

        self.u, self.du_darhoA, self.du_darhouA = computeVelocity(arhoA, arhouA)

        self.rho, drho_dvf, self.drho_darhoA, _ = computeDensity(self.vf, arhoA, self.A)
        self.drho_daA1 = drho_dvf * self.dvf_daA1

        v, dv_drho = computeSpecificVolume(self.rho)
        dv_daA1 = dv_drho * self.drho_daA1
        dv_darhoA = dv_drho * self.drho_darhoA

        e, de_dv, _ = self.eos.e(v, self.p)
        de_daA1 = de_dv * dv_daA1
        de_darhoA = de_dv * dv_darhoA

        self.E, dE_de, dE_du = addKineticEnergy(e, self.u)
        self.dE_daA1 = dE_de * de_daA1
        self.dE_darhoA = dE_de * de_darhoA + dE_du * self.du_darhoA
        self.dE_darhouA = dE_du * self.du_darhouA

    def applyWeakBC(self, U, r, J):
        self.computeQuantities(U)

        arhouA = U[self.i_arhouA]

        # mass
        r[self.i_arhoA] += arhouA * self.nx
        J[self.i_arhoA, self.i_arhouA] += self.nx

        # momentum
        r[self.i_arhouA] += (arhouA * self.u + self.vf * self.p * self.A) * self.nx
        if (self.model_type == ModelType.TwoPhase):
            J[self.i_arhouA, self.i_aA1] += self.dvf_daA1 * self.p * self.A * self.nx
        J[self.i_arhouA, self.i_arhoA] += arhouA * self.du_darhoA * self.nx
        J[self.i_arhouA, self.i_arhouA] += (arhouA * self.du_darhouA + self.u) * self.nx

        # energy
        if not self.strongly_enforce_energy:
            r[self.i_arhoEA] += self.vf * self.u * (self.rho * self.E + self.p) * self.A * self.nx
            if (self.model_type == ModelType.TwoPhase):
                J[self.i_arhoEA, self.i_aA1] += (
                    self.dvf_daA1 * self.u * (self.rho * self.E + self.p) + self.vf *
                    (self.drho_daA1 * self.E + self.rho * self.dE_daA1)) * self.A * self.nx
            J[self.i_arhoEA, self.i_arhoA] += self.vf * (
                self.du_darhoA * (self.rho * self.E + self.p) + self.u *
                (self.drho_darhoA * self.E + self.rho * self.dE_darhoA)) * self.A * self.nx
            J[self.i_arhoEA, self.i_arhouA] += self.vf * (
                self.du_darhouA * (self.rho * self.E + self.p) + self.u * self.rho * self.dE_darhouA
            ) * self.A * self.nx

    def applyStrongBCNonlinearSystem(self, U, r, J):
        arhoA = U[self.i_arhoA]
        arhoEA = U[self.i_arhoEA]

        # energy
        if self.strongly_enforce_energy:
            arhoEABC = arhoA * self.E
            darhoEABC_daA1 = arhoA * self.dE_daA1
            darhoEABC_darhoA = self.E + arhoA * self.dE_darhoA
            darhoEABC_darhouA = arhoA * self.dE_darhouA

            r[self.i_arhoEA] = arhoEA - arhoEABC
            J[self.i_arhoEA, :] = 0
            if (self.model_type == ModelType.TwoPhase):
                J[self.i_arhoEA, self.i_aA1] = -darhoEABC_daA1
            J[self.i_arhoEA, self.i_arhoA] = -darhoEABC_darhoA
            J[self.i_arhoEA, self.i_arhouA] = -darhoEABC_darhouA
            J[self.i_arhoEA, self.i_arhoEA] = 1

    def applyStrongBCLinearSystemMatrix(self, A):
        A[self.i_arhoEA, :] = 0
        A[self.i_arhoEA, self.i_arhoEA] = 1

    def applyStrongBCLinearSystemRHSVector(self, U, b):
        A = self.dof_handler.A[self.k]
        aA1 = self.dof_handler.aA1(U, self.k)
        vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)
        arhoA = U[self.i_arhoA]
        arhouA = U[self.i_arhouA]

        u, du_darhoA, du_darhouA = computeVelocity(arhoA, arhouA)
        rho, drho_dvf, drho_darhoA, _ = computeDensity(vf, arhoA, A)
        v, dv_drho = computeSpecificVolume(rho)
        e, de_dv, _ = self.eos.e(v, self.p)
        arhoEA = arhoA * (e + 0.5 * u * u)

        b[self.i_arhoEA] = arhoEA
