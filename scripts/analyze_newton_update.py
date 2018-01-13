"""
This script analyzes an entry of a Newton update vector to determine why it
attained an unphysical value.

Usage:
  analyze_newton_update.py -h | --help
  analyze_newton_update.py

Options:
  -h --help  Display help information.
"""

from docopt import docopt
import numpy as np

docopt(__doc__)

i = 1

J = np.loadtxt("jacobian.csv", delimiter=",")
r = np.loadtxt("residual.csv", delimiter=",")
n = J.shape[0]

J_inv = np.linalg.inv(J)
dU = J_inv.dot(-r)

print(("norm(J) = %e" % (np.linalg.cond(J))))

print("")

print(("%4s: %12s %12s %12s" % ("i", "Jinv[i,j]", "-r[j]", "-Jinvr[i,j]")))
print((44 * "="))
for j in range(n):
    print(("%4d: %12.2e %12.2e %12.2e" % (j, J_inv[i, j], -r[j], J_inv[i, j] * -r[j])))
