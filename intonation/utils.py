from numpy import array, abs

def find_nearest_index(arr, value):
    """For a given value, the function finds the nearest value
    in the array and returns its index."""
    arr = array(arr)
    index = (abs(arr-value)).argmin()
    return index