from EoS import EoS, EoSParameters

class TestEoSParameters(EoSParameters):
  def __init__(self):
    EoSParameters.__init__(self)
    self.registerFloatParameter("slope_initial", "First slope in sequence", 1.0)
    self.registerFloatParameter("slope_increment", "Increment in slope in sequence", 0.1)

class TestEoS(EoS):
  def __init__(self, params):
    EoS.__init__(self)
    self.slope = params.get("slope_initial")
    self.slope_increment = params.get("slope_increment")

    self.de_dv = self.nextSlope()
    self.de_dp = self.nextSlope()
    self.dp_dv = self.nextSlope()
    self.dp_de = self.nextSlope()
    self.dT_dv = self.nextSlope()
    self.dT_de = self.nextSlope()
    self.dc_dv = self.nextSlope()
    self.dc_dp = self.nextSlope()
    self.ds_dv = self.nextSlope()
    self.ds_de = self.nextSlope()
    self.dp_dh = self.nextSlope()
    self.dp_ds = self.nextSlope()

  def nextSlope(self):
    self.slope += self.slope_increment
    return self.slope

  def rho(self, p, T):
    return p + T

  ## Computes a test property that is a linear combination of 2 properties
  # @param[in] a  first property
  # @param[in] b  second property
  # @param[in] dda  derivative of test property with respect to first property
  # @param[in] ddb  derivative of test property with respect to second property
  def testProperty(self, a, b, dda, ddb):
    c = a * dda + b * ddb
    return (c, dda, ddb)

  def e(self, v, p):
    return self.testProperty(v, p, self.de_dv, self.de_dp)

  def p(self, v, e):
    return self.testProperty(v, e, self.dp_dv, self.dp_de)

  def T(self, v, e):
    return self.testProperty(v, e, self.dT_dv, self.dT_de)

  def c(self, v, p):
    return self.testProperty(v, p, self.dc_dv, self.dc_dp)

  def s(self, v, e):
    return self.testProperty(v, e, self.ds_dv, self.ds_de)

  def p_from_h_s(self, h, s):
    return self.testProperty(h, s, self.dp_dh, self.dp_ds)
