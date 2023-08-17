'''
This file implements some helper functions and tests for the ccaModels.py file
'''

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
import logging, sys

def find_zero_columns(P):
    return np.where(~P.any(axis=0))[0]

def remove_indices(P,indices):
    return np.delete(np.delete(P,indices,axis=1),indices,axis=0)

def replace_indices_with_zeros(v,indices):
    if len(indices) == 0:
        return v
    newV = np.zeros(v.size+len(indices))
    offset = 0
    for i in range(len(v)):
        if i in indices:
            offset += 1
        newV[i+offset] = v[i]
    return newV

def well_defined_pi(pi):
    if np.any(pi<-1e-5):
        print(f"Error in computing the stationnary distribution. Negative values in the distribution.")
        sys.exit()
    if abs(np.sum(pi) - 1) > 1e-5:
        print(f"Error in computing the stationnary distribution. Sum of the distribution is {np.sum(pi)}")
    return True