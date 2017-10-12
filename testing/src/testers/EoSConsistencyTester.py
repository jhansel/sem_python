from sem_python.utilities.numeric_utilities import computeRelativeDifference

class EoSConsistencyTester(object):
  def __init__(self, verbose=False):
    self.verbose = verbose

  def checkConsistency(self, eos):
    p = 1e5
    T = 300

    # perform base calls
    rho, _, _ = eos.rho(p, T)
    h, _, _ = eos.h(p, T)
    v = 1.0 / rho
    e, _, _ = eos.e(v, p)
    s, _, _ = eos.s(v, e)

    # perform calls to check consistency
    rho_from_p_s, _, _ = eos.rho_from_p_s(p, s)
    p_from_v_e, _, _ = eos.p(v, e)
    T_from_v_e, _, _ = eos.T(v, e)
    p_from_h_s, _, _ = eos.p_from_h_s(h, s)
    h_from_e_p_rho = e + p / rho
    s_from_h_p, _, _ = eos.s_from_h_p(h, p)

    base_values = dict()
    base_values["rho_from_p_s"] = rho
    base_values["p_from_v_e"] = p
    base_values["T_from_v_e"] = T
    base_values["p_from_h_s"] = p
    base_values["h_from_e_p_rho"] = h
    base_values["s_from_h_p"] = s

    test_values = dict()
    test_values["rho_from_p_s"] = rho_from_p_s
    test_values["p_from_v_e"] = p_from_v_e
    test_values["T_from_v_e"] = T_from_v_e
    test_values["p_from_h_s"] = p_from_h_s
    test_values["h_from_e_p_rho"] = h_from_e_p_rho
    test_values["s_from_h_p"] = s_from_h_p

    # compute relative differences
    reldiffs = dict()
    for check in test_values:
      reldiffs[check] = computeRelativeDifference(test_values[check], base_values[check])

    if self.verbose:
      print "EoS Consistency Check:\n"
      print "%-20s%-20s%-20s%-20s" % ("Check", "Base value", "Test value", "Relative difference")
      print "=" * 80
      for check in reldiffs:
        print "%-20s%-20g%-20g%-20g" % (check, base_values[check], test_values[check], reldiffs[check])
      print ""

    # take absolute values of relative differences
    abs_reldiffs = dict()
    for check in reldiffs:
      abs_reldiffs[check] = abs(reldiffs[check])

    return abs_reldiffs
