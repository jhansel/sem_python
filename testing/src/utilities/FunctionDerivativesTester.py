from inspect import getargspec

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/utilities")
from numeric_utilities import computeRelativeDifference

class FunctionDerivativesTester(object):
  def __init__(self, in_unittest_mode=True, use_debug_mode=False):
    self.in_unittest_mode = in_unittest_mode
    self.use_debug_mode = use_debug_mode

  def checkDerivatives(self, f, n_args, fd_eps=1e-8):
    # compute hand-coded derivatives
    vals = [2.0 * (i+1) + 1.0 for i in xrange(n_args)]
    U = tuple(vals)
    f_base = f(*U)[0]
    df_dU_coded = f(*U)[1:]
    if self.use_debug_mode:
      print "U =", U
      print "f(U) =", f_base
      print "df/dU (hand-coded) =", df_dU_coded

    # print header
    if not self.in_unittest_mode:
      # get the name of the function arguments and remove "self" if applicable
      args = getargspec(f).args
      if len(args) != n_args:
        if len(args) == n_args + 1:
          args.pop(0)
        else:
          raise Exception("Unexpected number of arguments.")

      # determine the max string length of the arguments for formatting
      max_arg_length = 3
      for arg in args:
        max_arg_length = max(max_arg_length, len(arg))

      # print header for table
      print f.__name__ + "(" + ",".join(args) + ") Derivatives:\n"
      format_string = "%" + str(max_arg_length) + "s%13s%13s%13s%13s"
      print format_string % ("Arg","Coded","FD","Abs Diff","Rel Diff")
      print "=" * (max_arg_length + 52)

    # compute finite difference derivatives
    reldiffs = []
    for i in xrange(n_args):
      U_forward = list(U)
      U_forward[i] += fd_eps
      U_forward = tuple(U_forward)
      f_forward = f(*U_forward)[0]
      df_dUi_fd = (f_forward - f_base) / fd_eps
      if self.use_debug_mode:
        print "U_forward =", U_forward
        print "f(U_forward) =", f_forward
        print "df/dUi (FD) =", df_dUi_fd

      # compute absolute and relative difference
      absdiff = df_dU_coded[i] - df_dUi_fd
      reldiff = computeRelativeDifference(df_dU_coded[i], df_dUi_fd)
      reldiffs.append(reldiff)

      # print results
      if not self.in_unittest_mode:
        format_string = "%" + str(max_arg_length) + "s%13.4e%13.4e%13.4e%13.4e"
        print format_string % (args[i], df_dU_coded[i], df_dUi_fd, absdiff, reldiff)

    # return relative differences if in unittest mode
    if self.in_unittest_mode:
      return [abs(x) for x in reldiffs]
    else:
      print "\n"
