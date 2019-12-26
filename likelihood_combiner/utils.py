import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from likelihood_combiner.reader import LklComReader,JFactor_Reader

def compute_sensitivity(sigmav, ts_dict, confidence_level = 2.71):
    limits = {}
    sensitivities  = {}
    for key, tstable in ts_dict.items():
        if np.all(tstable==0):
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
            lin_interpolation = interp1d(sigmav, tstable, kind='linear', fill_value='extrapolate')
            min_ind = np.min(np.where(tstable == np.min(tstable)))
            ts_val = []
            for sv in sigmav:
                g = sv/np.power(10,l)
                gSpline = lin_interpolation(g)
                gLkl = np.where(((g>1e-28) & (g<1e-18)), gSpline,  np.nan)
                totLkl = gLkl + lLkl
                ts_val.append(np.nanmin(totLkl))
            ts_val= np.array(ts_val)
            dis = ts_val[min_ind] - tstable[min_ind]
            ts_dict_Jnuisance[key] = ts_val - dis
    return ts_dict_Jnuisance

def plot_sigmavULs(config):
    
    try:
        hdf5file = config['Data']['hdf5_dataset']
        if hdf5file is None:
            raise KeyError
    except KeyError:
        hdf5file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/lklcom.h5"))
    
    try:
        output_dir = config['Output']['output_directory']
        if output_dir is None:
            raise KeyError
    except KeyError:
        output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output/"))
    # Create output directory if it doesn't exist already
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    channels = config['Configuration']['channels']
    sources = config['Configuration']['sources']
    collaborations = config['Configuration']['collaborations']

    # Corresponding LaTex notation
    channels_LaTex = {'bb':'b\\bar{b}', 'tautau':'\\tau^{+}\\tau^{-}', 'mumu':'\mu^{+}\mu^{-}', 'WW':'W^{+}W^{-}', 'gammagamma':'\gamma\gamma', 'ZZ':'ZZ', 'ee':'e^{+}e^{-}'}
    
    reader = LklComReader()
    sigmavULs, sigmavULs_Jnuisance, massvals = reader.read_sigmavULs(hdf5file,channels)
    del reader

    for channel in channels:
        for ul_dict in [sigmavULs,sigmavULs_Jnuisance]:
            if ul_dict is sigmavULs_Jnuisance and not config['Data']['J_nuisance']:
                continue
            fig, ax = plt.subplots()
            ax.plot(massvals[channel+"_Combination"], ul_dict[channel+"_Combination"], label='Combined limit')
            ax.set_xscale('log')
            ax.set_xbound(lower=np.min(massvals[channel+"_Combination"]),upper=np.max(massvals[channel+"_Combination"]))
            ax.set_xlabel(r'$m_{\chi} \: [GeV]$')
            ax.set_yscale('log')
            ax.set_ylabel(r'$95\%$ CL $\langle\sigma v\rangle^{UL} \, [cm^{3}/s]$')
            if ul_dict is sigmavULs:
                ax.set_title(r'$\langle\sigma v\rangle$ ULs vs mass - J fixed')
            else:
                ax.set_title(r'$\langle\sigma v\rangle$ ULs vs mass - J as nuisance')

            ax.text(0.2, 0.85, 'All dSphs', fontsize=18,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            ax.text(0.85, 0.1, r'$\chi\chi \to {}$'.format(channels_LaTex[str(channel)]), fontsize=15,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
            # Shrink current axis by 20%
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        
            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8)
            ax.grid(b=True,which='both',color='grey', linestyle='--', linewidth=0.25)
            if ul_dict is sigmavULs:
                plt.savefig('{}lklcom_{}.pdf'.format(output_dir,channel))
                print("Saved plot in {}lklcom_{}.pdf".format(output_dir,channel))
                plt.close()
            else:
                plt.savefig('{}lklcom_{}_Jnuisance.pdf'.format(output_dir,channel))
                print("Saved plot in {}lklcom_{}_Jnuisance.pdf".format(output_dir,channel))
                plt.close()
            ax.clear()

    return
    
