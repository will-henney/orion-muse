import numpy as np

def fuzz(data, variance, nfuzz=1):
    '''Return `nfuzz` fuzzed versions of `data` in which each pixel is
replaced by a value taken from a Gaussian distribution with mean equal
to `data` and sigma equal to sqrt of `variance` (which must have same
shape as `data`).  Returns an array of shape `(nfuzz,) + shape(data)` in
which each (hyper) plane is a distinct fuzzed version of the original

    '''
    newshape = (nfuzz,) + data.shape
    return np.random.normal(data, np.sqrt(variance), newshape)
