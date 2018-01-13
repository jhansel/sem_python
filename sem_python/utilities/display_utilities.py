import numpy as np
import sys
from termcolor import colored

from .numeric_utilities import computeRelativeDifference


def computeRelativeDifferenceMatrix(A, B):
    n = A.shape[0]
    C = np.zeros(shape=(n, n))
    for i in range(n):
        for j in range(n):
            C[i, j] = computeRelativeDifference(A[i, j], B[i, j])
    return C


def printMatrix(A):
    n = A.shape[0]
    for i in range(n):
        for j in range(n):
            # underline the diagonal
            if i == j:
                attrs = ["underline"]
            else:
                attrs = list()

            # color the nonzero values
            if (abs(A[i, j]) > 0):
                color = "blue"
            else:
                color = "white"

            sys.stdout.write(colored("   %11.4e" % (A[i, j]), color, attrs=attrs))

        sys.stdout.write("\n")
    sys.stdout.flush()


def printMatrixDifference(A, red_threshold=None, yellow_threshold=None):
    n = A.shape[0]
    for i in range(n):
        for j in range(n):
            # underline the diagonal
            if i == j:
                attrs = ["underline"]
            else:
                attrs = list()

            # color by match
            if (red_threshold and abs(A[i, j]) >= red_threshold):
                color = "red"
            elif (yellow_threshold and abs(A[i, j]) >= yellow_threshold):
                color = "yellow"
            elif (abs(A[i, j]) > 0):
                color = "green"
            else:
                color = "white"

            sys.stdout.write(colored("   %11.4e" % (A[i, j]), color, attrs=attrs))

        sys.stdout.write("\n")
    sys.stdout.flush()


def printRelativeMatrixDifference(
        A, A_diff, red_threshold=None, yellow_threshold=None, abs_threshold=1e-6):
    n = A.shape[0]
    for i in range(n):
        for j in range(n):
            # underline the diagonal
            if i == j:
                attrs = ["underline"]
            else:
                attrs = list()

            # color by match
            not_near_zero = abs(A_diff[i, j]) > abs_threshold
            if (not_near_zero and red_threshold and abs(A[i, j]) >= red_threshold):
                color = "red"
            elif (not_near_zero and yellow_threshold and abs(A[i, j]) >= yellow_threshold):
                color = "yellow"
            elif (abs(A[i, j]) > 0):
                color = "green"
            else:
                color = "white"

            sys.stdout.write(colored("   %11.4e" % (A[i, j]), color, attrs=attrs))

        sys.stdout.write("\n")
    sys.stdout.flush()


def printDoFVector(U, dof_handler):
    n_var = dof_handler.n_var
    header_items = ("i", ) + tuple(dof_handler.variable_names)
    header_format = "%4s" + " %12s" * n_var
    entry_format = "%4i" + " %12.3e" * n_var
    print(header_format % header_items)
    for k in range(dof_handler.n_node):
        items_k = list()
        for m in range(n_var):
            i_k_m = dof_handler.i(k, m)
            items_k.append(U[i_k_m])
        entry_items = (k, ) + tuple(items_k)
        print(entry_format % entry_items)
    print("")


def printMatrixSparsity(A):
    n = A.shape[0]

    print("")
    row_is_zero = n * [True]
    for i in range(n):
        line = ""
        for j in range(n):
            if abs(A[i, j]) > 1e-14:
                line += "X"
                row_is_zero[j] = False
            else:
                line += " "
        print(line)

    column_is_zero = n * [True]
    for j in range(n):
        for i in range(n):
            if abs(A[i, j]) > 1e-14:
                column_is_zero[j] = False

    print("")
    no_zero_rows = True
    no_zero_columns = True
    for i in range(n):
        if row_is_zero[i]:
            print("WARNING: Row ", i, " is zero.")
            no_zero_rows = False
        if column_is_zero[i]:
            print("WARNING: Column ", i, " is zero.")
            no_zero_columns = False
    if no_zero_rows:
        print("There were no zero rows.")
    if no_zero_columns:
        print("There were no zero columns.")
    print("")
