import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType

sys.path.append(base_dir + "src/bc")
from OnePhaseBC import OnePhaseBC, OnePhaseBCParameters

sys.path.append(base_dir + "src/closures")
from thermodynamic_functions import computeDensity, computeVelocity, computeSpecificVolume

## Parameters class for OutletBC
class OutletBCParameters(OnePhaseBCParameters):
  def __init__(self):
    OnePhaseBCParameters.__init__(self)
    self.registerFloatParameter("p", "specified outlet pressure")

class OutletBC(OnePhaseBC):
  def __init__(self, params, dof_handler, eos_map):
    OnePhaseBC.__init__(self, params, dof_handler, eos_map)
    self.p = params.get("p")

  def apply(self, U, r, J):
    vf, dvf_dvf1 = self.dof_handler.getVolumeFraction(U, self.k, self.phase)
    arho = U[self.i_arho]
    arhou = U[self.i_arhou]
    arhoE_solution = U[self.i_arhoE]

    u, du_darho, du_darhou = computeVelocity(arho, arhou)

    rho, drho_dvf, drho_darho = computeDensity(vf, arho)
    drho_dvf1 = drho_dvf * dvf_dvf1

    v, dv_drho = computeSpecificVolume(rho)
    dv_dvf1 = dv_drho * drho_dvf1
    dv_darho = dv_drho * drho_darho

    e, de_dv, de_dp = self.eos.e(v, self.p)
    de_dvf1 = de_dv * dv_dvf1
    de_darho = de_dv * dv_darho

    arhoE = arho * (e + 0.5 * u * u)
    darhoE_dvf1 = arho * de_dv * dv_dvf1
    darhoE_darho = (e + 0.5 * u * u) + arho * (de_dv * dv_darho + u * du_darho)
    darhoE_darhou = arho * u * du_darhou

    # mass
    r[self.i_arho] += arhou * self.nx
    J[self.i_arho,self.i_arhou] += self.nx

    # momentum
    r[self.i_arhou] += (arhou * u + vf * self.p) * self.nx
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhou,self.i_vf1] += dvf_dvf1 * self.p * self.nx
    J[self.i_arhou,self.i_arho] += arhou * du_darho * self.nx
    J[self.i_arhou,self.i_arhou] += (arhou * du_darhou + u) * self.nx

    # energy
    r[self.i_arhoE] = arhoE_solution - arhoE
    J[self.i_arhoE,:] = 0
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhoE,self.i_vf1] = - darhoE_dvf1
    J[self.i_arhoE,self.i_arho] = - darhoE_darho
    J[self.i_arhoE,self.i_arhou] = - darhoE_darhou
    J[self.i_arhoE,self.i_arhoE] = 1
