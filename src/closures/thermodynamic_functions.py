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
