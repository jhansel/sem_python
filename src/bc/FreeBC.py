import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType

sys.path.append(base_dir + "src/bc")
from OnePhaseBC import OnePhaseBC, OnePhaseBCParameters

sys.path.append(base_dir + "src/closures")
from thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeSpecificVolume, computeVelocity, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy

class FreeBCParameters(OnePhaseBCParameters):
  def __init__(self):
    OnePhaseBCParameters.__init__(self)

class FreeBC(OnePhaseBC):
  def __init__(self, params, dof_handler, eos_map):
    OnePhaseBC.__init__(self, params, dof_handler, eos_map)

  def applyWeakBC(self, U, r, J):
    vf1 = self.dof_handler.getVolumeFraction(U, self.k)
    vf, dvf_dvf1 = computeVolumeFraction(vf1, self.phase, self.model_type)
    arho = U[self.i_arho]
    arhou = U[self.i_arhou]
    arhoE = U[self.i_arhoE]

    u, du_darho, du_darhou = computeVelocity(arho, arhou)

    rho, drho_dvf, drho_darho = computeDensity(vf, arho)
    drho_dvf1 = drho_dvf * dvf_dvf1

    v, dv_drho = computeSpecificVolume(rho)
    dv_dvf1 = dv_drho * drho_dvf1
    dv_darho = dv_drho * drho_darho

    E, dE_darho, dE_darhoE = computeSpecificTotalEnergy(arho, arhoE)

    e, de_du, de_dE = computeSpecificInternalEnergy(u, E)
    de_darho = de_du * du_darho + de_dE * dE_darho
    de_darhou = de_du * du_darhou
    de_darhoE = de_dE * dE_darhoE

    p, dp_dv, dp_de = self.eos.p(v, e)
    dp_dvf1 = dp_dv * dv_dvf1
    dp_darho = dp_dv * dv_darho + dp_de * de_darho
    dp_darhou = dp_de * de_darhou
    dp_darhoE = dp_de * de_darhoE

    # mass
    r[self.i_arho] += rhou * self.nx
    J[self.i_arho,self.i_arhou] += self.nx

    # momentum
    r[self.i_arhou] += (arhou * u + vf * p) * self.nx
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhou,self.i_vf1] += (dvf_dvf1 * p + vf * dp_dvf1) * self.nx
    J[self.i_arhou,self.i_arho] += (arhou * du_darho + vf * dp_darho) * self.nx
    J[self.i_arhou,self.i_arhou] += (arhou * du_darhou + u + vf * dp_darhou) * self.nx
    J[self.i_arhou,self.i_arhoE] += vf * dp_darhoE * self.nx

    # energy
    r[self.i_arhoE] += (arhoE + vf * p) * u * self.nx
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhoE,self.i_vf1] += (dvf_dvf1 * p + vf * dp_dvf1) * u * self.nx
    J[self.i_arhoE,self.i_arho] += ((arhoE + vf * p) * du_darho + vf * dp_darho * u) * self.nx
    J[self.i_arhoE,self.i_arhou] += ((arhoE + vf * p) * du_darhou + vf * dp_darhou * u) * self.nx
    J[self.i_arhoE,self.i_arhoE] += (1 + vf * dp_darhoE) * u * self.nx

  def applyStrongBC(self, U, r, J):
    pass
