"""
This script returns a list of a MINIMUM set of linearly dependent rows of a
matrix. Currently, the script looks for a file called jacobian.csv in the
current directory, which has the matrix in CSV format.

Usage:
  analyze_linear_dependencies.py -h | --help
  analyze_linear_dependencies.py

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

n = A.shape[0]
rank = matrix_rank(A)

if rank == n:
    print("Matrix is of full rank.")
    exit()
else:
    print(("Matrix is rank-deficient: Rank = ", rank, ", Size = ", n))

i_dep = []
not_finished = True
while not_finished:
    # Each iteration adds an entry to i_dep.

    # The stopping criteria is that the rows in i_dep are linearly dependent, by themselves.
    if len(i_dep) > 0:
        n_rows = len(i_dep)
        rank = matrix_rank(A[i_dep, :])
        if rank < n_rows:
            print(("Minimum set of linearly dependent rows:", i_dep))
            break

    i_test = deepcopy(i_dep)
    for i in range(n):
        if i not in i_test:
            i_test.append(i)
            n_rows = len(i_test)
            rank = matrix_rank(A[i_test, :])
            if rank < n_rows:
                i_dep.append(i)
                break

# print out entries of rows
j_min_nonzero = n
j_max_nonzero = 0
for k in i_dep:
    # find minimum j that is nonzero
    row_is_zero = True
    for j in range(n):
        if abs(A[k, j]) > 1e-14:
            row_is_zero = False
            j_min_nonzero = min(j_min_nonzero, j)
            j_max_nonzero = max(j_max_nonzero, j)
    if row_is_zero:
        print(("Row %d is a zero row." % (k)))

print(("Nonzero range: j in (%d, %d)" % (j_min_nonzero, j_max_nonzero)))

print("Nonzeros for each linearly dependent row:")
n_entries = j_max_nonzero - j_min_nonzero + 1
line_format = "%10.2e" * n_entries
for k in i_dep:
    print((line_format % tuple(A[k, j_min_nonzero:j_max_nonzero + 1])))
