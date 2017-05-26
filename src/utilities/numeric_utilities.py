def computeRelativeDifference(a, b):
  if (abs(b) <= 1e-15):
    return a - b
  else:
    return (a - b) / b
