from numpy import sin, cos, exp, log, pi
import parser

from .error_utilities import error

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

def stringToStringList(s):
  if isinstance(s, str):
    return s.split()
  elif isinstance(s, list):
    return s
  else:
    error("'" + s + "' is not a valid string list.")

def stringToFloatList(s):
  if isinstance(s, str):
    string_list = s.split()
    return [stringToFloat(x) for x in string_list]
  elif isinstance(s, list):
    return s
  else:
    error("'" + s + "' is not a valid float list.")
