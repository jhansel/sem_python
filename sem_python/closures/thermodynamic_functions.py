from ..base.enums import ModelType, VariableName


def computeVolumeFraction(aA1, A, phase, model_type):
    if model_type == ModelType.OnePhase:
        vf = 1
        dvf_daA1 = 0
    else:
        if phase == 0:
            vf = aA1 / A
            dvf_daA1 = 1 / A
        else:
            vf = 1 - aA1 / A
            dvf_daA1 = -1 / A
    return (vf, dvf_daA1)


def computeSpecificVolume(rho):
    v = 1.0 / rho
    dv_drho = -1.0 / rho / rho

    return (v, dv_drho)


def computeSpecificTotalEnergy(arhoA, arhoEA):
    E = arhoEA / arhoA
    dE_darhoA = -arhoEA / arhoA / arhoA
    dE_darhoEA = 1.0 / arhoA

    return (E, dE_darhoA, dE_darhoEA)


def computeSpecificInternalEnergy(u, E):
    e = E - 0.5 * u * u
    de_dE = 1.0
    de_du = -u

    return (e, de_du, de_dE)


def computeSpecificEnthalpy(e, p, rho):
    h = e + p / rho
    dh_de = 1.0
    dh_dp = 1.0 / rho
    dh_drho = -p / rho / rho

    return (h, dh_de, dh_dp, dh_drho)


## Adds kinetic energy; this can be used to compute H from h or E from e
def addKineticEnergy(e, u):
    E = e + 0.5 * u * u
    dE_de = 1
    dE_du = u

    return (E, dE_de, dE_du)


## Subtracts kinetic energy; this can be used to compute h from H or e from E
def subtractKineticEnergy(E, u):
    e = E - 0.5 * u * u
    de_dE = 1
    de_du = -u

    return (e, de_dE, de_du)


def computeVelocity(arhoA, arhouA):
    u = arhouA / arhoA
    du_darhoA = -arhouA / arhoA / arhoA
    du_darhouA = 1.0 / arhoA

    return (u, du_darhoA, du_darhouA)


def computeDensity(vf, arhoA, A):
    rho = arhoA / (vf * A)
    drho_dvf = -arhoA / (vf * vf * A)
    drho_darhoA = 1.0 / (vf * A)
    drho_dA = -arhoA / (vf * A * A)

    return (rho, drho_dvf, drho_darhoA, drho_dA)
