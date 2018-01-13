def computeRelativeDifference(a, b):
    a_is_zero = abs(a) <= 1e-15
    b_is_zero = abs(b) <= 1e-15
    if a_is_zero and b_is_zero:
        return 0
    elif a_is_zero:
        return 1
    elif b_is_zero:
        return 1
    else:
        return (a - b) / b
