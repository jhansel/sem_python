from enums import ModelType, VariableName

def computeVolumeFraction(vf1, phase, model_type):
  if model_type == ModelType.OnePhase:
    vf = 0 * vf1 + 1
    dvf_dvf1 = 0 * vf1 + 1
  else:
    if phase == 0:
      vf = vf1
      dvf_dvf1 = 0 * vf1 + 1
    else:
      vf = 1 - vf1
      dvf_dvf1 = 0 * vf1 - 1
  return (vf, dvf_dvf1)

def computeSpecificVolume(rho):
  v = 1.0 / rho
  dv_drho = - 1.0 / rho / rho

  return (v, dv_drho)

def computeSpecificTotalEnergy(arhoA, arhoEA):
  E = arhoEA / arhoA
  dE_darhoA = - arhoEA / arhoA / arhoA
  dE_darhoEA = 1.0 / arhoA

  return (E, dE_darhoA, dE_darhoEA)

def computeSpecificInternalEnergy(u, E):
  e = E - 0.5 * u * u
  de_dE = 1.0
  de_du = - u

  return (e, de_du, de_dE)

def computeSpecificEnthalpy(e, p, rho):
  h = e + p / rho
  dh_de = 1.0
  dh_dp = 1.0 / rho
  dh_drho = - p / rho / rho

  return (h, dh_de, dh_dp, dh_drho)

def computeVelocity(arhoA, arhouA):
  u = arhouA / arhoA
  du_darhoA = - arhouA / arhoA / arhoA
  du_darhouA = 1.0 / arhoA

  return (u, du_darhoA, du_darhouA)

def computeDensity(vf, arhoA):
  rho = arhoA / vf
  drho_dvf = - arhoA / vf / vf
  drho_darhoA = 1.0 / vf

  return (rho, drho_dvf, drho_darhoA)
