import numpy as np
import sys
from termcolor import colored

from numeric_utilities import computeRelativeDifference

def computeRelativeDifferenceMatrix(A, B):
  n = A.shape[0]
  C = np.zeros(shape=(n, n))
  for i in xrange(n):
    for j in xrange(n):
      C[i,j] = computeRelativeDifference(A[i,j], B[i,j])
  return C

def printMatrix(A):
  n = A.shape[0]
  for i in xrange(n):
    for j in xrange(n):
      if (abs(A[i,j]) > 0):
        sys.stdout.write(colored("   %11.4e" % (A[i,j]), "blue"))
      else:
        sys.stdout.write("   %11.4e" % (A[i,j]))
    sys.stdout.write("\n")
  sys.stdout.flush()

def printMatrixDifference(A, red_threshold=None, yellow_threshold=None):
  n = A.shape[0]
  for i in xrange(n):
    for j in xrange(n):
      if (red_threshold and abs(A[i,j]) >= red_threshold):
        sys.stdout.write(colored("   %11.4e" % (A[i,j]), "red"))
      elif (yellow_threshold and abs(A[i,j]) >= yellow_threshold):
        sys.stdout.write(colored("   %11.4e" % (A[i,j]), "yellow"))
      elif (abs(A[i,j]) > 0):
        sys.stdout.write(colored("   %11.4e" % (A[i,j]), "green"))
      else:
        sys.stdout.write("   %11.4e" % (A[i,j]))
    sys.stdout.write("\n")
  sys.stdout.flush()
