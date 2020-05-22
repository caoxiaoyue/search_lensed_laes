import numpy as np

def find_blocks(value, a):
    
    """
    a: the searching array
    value: the value
    searching value a in array a;
    """
    # Create an array that is 1 where a is `value`, 
    # and pad each end with an extra 0.
    
    isvalue = np.concatenate(([0], np.equal(a, value).view(np.int8), [0]))
    absdiff = np.abs(np.diff(isvalue))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges
