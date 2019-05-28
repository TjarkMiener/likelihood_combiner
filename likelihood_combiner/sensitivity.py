import numpy as np

def compute_sensitivity(sigmav, tstable, confidence_level = 2.71):
    limits = []
    for imass in np.arange(0,tstable.shape[0],1):
        min_ind = list(tstable[imass]).index(np.min(tstable[imass]))
        up_pos = min_ind
        for ind in np.arange(min_ind,0,-1):
            if tstable[imass][ind] > confidence_level:
                up_pos = ind+1
                break
        if up_pos != min_ind:
            limit = (confidence_level - tstable[imass][up_pos]) * (sigmav[up_pos-1] - sigmav[up_pos]) / (tstable[imass][up_pos-1] - tstable[imass][up_pos]) + sigmav[up_pos]
        limits.append(limit)
    return np.array(limits)
