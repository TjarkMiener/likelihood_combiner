import numpy as np
from scipy.interpolate import interp1d

__all__ = [
    'compute_sensitivity'
]

def compute_sensitivity(sigmav, ts_dict, confidence_level = 2.71):
    limits = {}
    sensitivities  = {}
    for key, tstable in ts_dict.items():
        if np.all(tstable==0) or tstable[-1] < confidence_level:
            limits[key] = np.nan
            sensitivities[key] = np.nan
        else:
            min_val = np.min(tstable)
            min_ind = list(tstable).index(np.min(tstable))
            if min_val > 0.0:
                ts = tstable - min_val
            elif min_val < 0.0:
                ts = tstable + np.abs(min_val)
            else:
                ts = tstable
            cub_interpolation = interp1d(sigmav, ts, kind='cubic')
            if min_ind == 0:
                min_sigmav = sigmav[min_ind]
            elif min_ind == sigmav.shape[0]-1:
                limits[key] = np.nan
                sensitivities[key] = np.nan
                continue
            else:
                x = np.linspace(sigmav[min_ind-1], sigmav[min_ind+1], num=25, endpoint=True)
                y = cub_interpolation(x)
                min_sigmav = x[list(y).index(np.min(y))]
            for ind in np.arange(min_ind,sigmav.shape[0]):
                if ts[ind] > confidence_level:
                    x = np.linspace(sigmav[ind-1], sigmav[ind], num=25, endpoint=True)
                    y = cub_interpolation(x)
                    for i,y_val in enumerate(y):
                        if y_val > confidence_level:
                            limit = (confidence_level - y[i]) * (x[i-1] - x[i]) / (y[i-1] - y[i]) + x[i]
                            sensitivty = limit - min_sigmav
                            break
                    break
                else:
                    limit = np.nan
                    sensitivty = np.nan
            limits[key] = limit
            sensitivities[key] = sensitivty
    return limits, sensitivities
