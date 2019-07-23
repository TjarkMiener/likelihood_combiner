import numpy as np
from scipy.interpolate import CubicSpline

def compute_sensitivity(sigmav, ts_dict, confidence_level = 2.71):
    # We need to reverse the sigmav array, because x values in CubicSpline() must be in strictly increasing order.
    limits = {}
    for key, tstable in ts_dict.items():
        if np.all(tstable==0):
            limits[key] = np.nan
        else:
            min_ind = list(tstable).index(np.min(tstable))
            cSpline = CubicSpline(sigmav, tstable)
            if min_ind != 0:
                x = np.linspace(sigmav[min_ind-1], sigmav[min_ind+1], num=25, endpoint=True)
                y = cSpline(x)
                min_sigmav = x[list(y).index(np.min(y))]
            else:
                min_sigmav = sigmav[min_ind]

            for ind in np.arange(min_ind,sigmav.shape[0]):
                if tstable[ind] > confidence_level:
                    x = np.linspace(sigmav[ind-1], sigmav[ind], num=25, endpoint=True)
                    y = cSpline(x)
                    for i,y_val in enumerate(y):
                        if y_val > confidence_level:
                            limit = (confidence_level - y[i]) * (x[i-1] - x[i]) / (y[i-1] - y[i]) + x[i]
                            limit -= min_sigmav
                            break
                    break
            limits[key] = limit
    return limits


def compute_Jnuisance(sigmav, ts_dict, DlogJ_dict):
    ts_dict_Jnuisance = {}
    for key, tstable in ts_dict.items():
        if np.all(tstable==0):
            ts_dict_Jnuisance[key] = tstable
        else:
            source = key.split("_")
            DlogJ = DlogJ_dict[source[0]]
            maxdev = 6. * DlogJ
            profiling_steps = 10000
            l = np.linspace(-maxdev, maxdev, num=profiling_steps, endpoint=True)
            lLkl = -2.*np.log((np.exp(-np.power(l,2)/(2.*np.power(DlogJ,2)))/(np.sqrt(2.*np.pi)*DlogJ*np.log(10))))
            cSpline = CubicSpline(sigmav, tstable)
            min_ind = np.min(np.where(tstable == np.min(tstable)))
            ts_val = []
            for sv in sigmav:
                g = sv/np.power(10,l)
                gSpline = cSpline(g)
                gLkl = np.where(((g>1e-28) & (g<1e-18)), gSpline,  np.nan)
                totLkl = gLkl + lLkl
                ts_val.append(np.nanmin(totLkl))
            ts_val= np.array(ts_val)
            dis = ts_val[min_ind] - tstable[min_ind]
            ts_dict_Jnuisance[key] = ts_val - dis
    return ts_dict_Jnuisance
