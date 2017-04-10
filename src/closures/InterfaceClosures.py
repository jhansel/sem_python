import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

class InterfaceClosuresParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFloatParameter("chi", "Weight fraction for phase 1", 0.5)
    self.registerFloatParameter("pressure_relaxation_time", "Relaxation time for pressures")

class InterfaceClosures(object):
  def __init__(self, params):
    self.pressure_relaxation_time = params.get("pressure_relaxation_time")
    self.chi = params.get("chi") # should be in (0,1)

  def computeInterfaceVelocity(self, u1, u2, beta):
    uI = beta * u1 + (1 - beta) * u2
    duI_du1 = beta
    duI_du2 = (1 - beta)
    duI_dbeta = u1 - u2

    return (uI, duI_du1, duI_du2, duI_dbeta)

  def computeInterfacePressure(self, p1, p2, mu):
    pI = mu * p1 + (1 - mu) * p2
    dpI_dp1 = mu
    dpI_dp2 = (1 - mu)
    dpI_dmu = p1 - p2

    return (pI, dpI_dp1, dpI_dp2, dpI_dmu)

  def computeBeta(self, arho1, arho2):
    denominator = self.chi * arho1 + (1 - self.chi) * arho2
    beta = self.chi * arho1 / denominator
    dbeta_darho1 = self.chi / denominator - self.chi * arho1 / denominator / denominator * self.chi
    dbeta_darho2 = - self.chi * arho1 / denominator / denominator * (1 - self.chi)

    return (beta, dbeta_darho1, dbeta_darho2)

  def computeMu(self, T1, T2, beta):
    denominator = beta * T1 + (1 - beta) * T2
    mu = (1 - beta) * T2 / denominator
    dmu_dT1 = - (1 - beta) * T2 / denominator / denominator * beta
    dmu_dT2 = (1 - beta) / denominator - (1 - beta) * T2 / denominator / denominator * (1 - beta)
    dmu_dbeta = - T2 / denominator - (1 - beta) * T2 / denominator / denominator * (T1 - T2)

    return (mu, dmu_dT1, dmu_dT2, dmu_dbeta)

  def computeTheta(self, vf1, p1, p2):
    vf2 = 1 - vf1
    dvf2_dvf1 = 0 * vf2 - 1

    denominator = self.pressure_relaxation_time * (p1 + p2)
    theta = vf1 * vf2 / denominator
    dtheta_ddenominator = - vf1 * vf2 / denominator / denominator
    dtheta_dvf1 = (vf2 + vf1 * dvf2_dvf1) / denominator
    dtheta_dp1 = dtheta_ddenominator * self.pressure_relaxation_time
    dtheta_dp2 = dtheta_ddenominator * self.pressure_relaxation_time

    return (theta, dtheta_dvf1, dtheta_dp1, dtheta_dp2)
