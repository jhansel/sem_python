import numpy as np
import sys
from termcolor import colored

def computeRelativeDifferenceMatrix(A, B):
  n = A.shape[0]
  C = np.zeros(shape=(n, n))
  for i in xrange(n):
    for j in xrange(n):
      if (abs(B[i,j]) <= 1e-15):
        C[i,j] = A[i,j] - B[i,j]
      else:
        C[i,j] = (A[i,j] - B[i,j]) / B[i,j]
  return C

def printMatrix(A, red_threshold=None, yellow_threshold=None):
  n = A.shape[0]
  for i in xrange(n):
    for j in xrange(n):
      if (red_threshold and abs(A[i,j]) >= red_threshold):
        sys.stdout.write(colored("   %11.4e" % (A[i,j]), "red"))
      elif (yellow_threshold and abs(A[i,j]) >= yellow_threshold):
        sys.stdout.write(colored("   %11.4e" % (A[i,j]), "yellow"))
      else:
        sys.stdout.write("   %11.4e" % (A[i,j]))
    sys.stdout.write("\n")
  sys.stdout.flush()
