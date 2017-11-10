"""
This script prints a set of rows from a matrix.

Usage:
  print_matrix_rows.py -h | --help
  print_matrix_rows.py

Options:
  -h --help  Display help information.
"""

from copy import deepcopy
from docopt import docopt
import numpy as np
from numpy.linalg import matrix_rank
from sys import exit

docopt(__doc__)

A = np.loadtxt("jacobian.csv", delimiter=",")
i_list = [63, 64, 65]

n = A.shape[0]

# print out entries of rows
j_min_nonzero = n
j_max_nonzero = 0
for i in i_list:
  # find minimum j that is nonzero
  row_is_zero = True
  for j in range(n):
    if abs(A[i,j]) > 1e-14:
      row_is_zero = False
      j_min_nonzero = min(j_min_nonzero, j)
      j_max_nonzero = max(j_max_nonzero, j)
  if row_is_zero:
    print(("Row %d is a zero row." % (i)))

print(("Nonzero range: j in (%d, %d)" % (j_min_nonzero, j_max_nonzero)))

print("Nonzeros for each selected row:")
n_entries = j_max_nonzero - j_min_nonzero + 1
line_format = "%10.2e" * n_entries
for i in i_list:
  print((line_format % tuple(A[i,j_min_nonzero:j_max_nonzero+1])))
