from enums import ModelType, VariableName

def computeVolumeFraction(vf1, phase, model_type):
  if model_type == ModelType.OnePhase:
    vf = 0 * vf1 + 1
    dvf_dvf1 = float("NaN")
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

def computeSpecificTotalEnergy(arho, arhoE):
  E = arhoE / arho
  dE_darho = - arhoE / arho / arho
  dE_darhoE = 1.0 / arho

  return (E, dE_darho, dE_darhoE)

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

def computeVelocity(arho, arhou):
  u = arhou / arho
  du_darho = - arhou / arho / arho
  du_darhou = 1.0 / arho

  return (u, du_darho, du_darhou)

def computeDensity(vf, arho):
  rho = arho / vf
  drho_dvf = - arho / vf / vf
  drho_darho = 1.0 / vf

  return (rho, drho_dvf, drho_darho)
