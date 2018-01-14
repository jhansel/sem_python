import numpy as np


## Initializes derivative data to zero for each solution variable
# @param[in] names  list of aux quantity names
# @param[in] n  number of values for each entry
def initializeDerivativeData(names, n):
    variable_names = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]
    der = dict()
    for name in names:
        der[name] = dict()
        for var in variable_names:
            der[name][var] = np.zeros(n)
            der[name]["grad_" + var] = np.zeros(n)
    return der
