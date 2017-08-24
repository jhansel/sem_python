from numpy import sin, cos, exp, log, pi
import parser

from error_utilities import error

def stringToBool(s):
  if s == "True" or s == True:
    return True
  elif s == "False" or s == False:
    return False
  else:
    error("'" + s + "' is not a valid boolean.")
  return value

def stringToInt(s):
  try:
    value = int(s)
  except:
    error("'" + s + "' is not a valid integer.")
  return value

def stringToFloat(s):
  try:
    value = float(s)
  except:
    error("'" + s + "' is not a valid float.")
  return value

def stringToFunction(s):
  code = parser.expr(s).compile()
  def f(x):
    return eval(code)
  return f
