import numpy as np

def find_nearest_index(arr, value):
    """For a given value, the function finds the nearest value
    in the array and returns its index."""
    arr = np.array(arr)
    index = (abs(arr-value)).argmin()
    return index

def normal_equation(x, y):
    """
    X: matrix of features, one sample per row (without bias unit)
    y: values (continuous) corresponding to rows (samples) in X
    """
    num_samples = y.size
    x = np.insert(x, 0, np.ones(num_samples), axis=1)

    return pinv(x)*y

