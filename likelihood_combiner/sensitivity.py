import numpy as np
from scipy.interpolate import CubicSpline

def compute_sensitivity(sigmav, tstable, confidence_level = 2.71):
    # We need to reverse the sigmav array, because x values in CubicSpline() must be in strictly increasing order.
    sigmav = sigmav[::-1]
    limits = []
    for imass in np.arange(0,tstable.shape[0],1):
        # Reverse the tstable[imass] as well (see above)
        tstable[imass] = tstable[imass][::-1]
        min_ind = list(tstable[imass]).index(np.min(tstable[imass]))
        cSpline = CubicSpline(sigmav, tstable[imass])
        if min_ind != 0:
            x = np.linspace(sigmav[min_ind-1], sigmav[min_ind+1], num=25, endpoint=True)
            y = cSpline(x)
            min_sigmav = x[list(y).index(np.min(y))]
        else:
            min_sigmav = sigmav[min_ind]

        for ind in np.arange(min_ind,sigmav.shape[0]):
            if tstable[imass][ind] > confidence_level:
                x = np.linspace(sigmav[ind-1], sigmav[ind], num=25, endpoint=True)
                y = cSpline(x)
                for i,y_val in enumerate(y):
                    if y_val > confidence_level:
                        limit = (confidence_level - y[i]) * (x[i-1] - x[i]) / (y[i-1] - y[i]) + x[i]
                        limit -= min_sigmav
                        break
                break
        limits.append(limit)
    return np.array(limits)
