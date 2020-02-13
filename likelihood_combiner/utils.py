import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.interpolate import interp1d

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


def compute_Jnuisance(sigmav, ts_dict, DlogJ_dict, sigmav_min=1e-28, sigmav_max=1e-18):
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
                gLkl = np.where(((g>sigmav_min) & (g<sigmav_max)), gSpline,  np.nan)
                totLkl = gLkl + lLkl
                ts_val.append(np.nanmin(totLkl))
            ts_val= np.array(ts_val)
            dis = ts_val[min_ind] - tstable[min_ind]
            ts_dict_Jnuisance[key] = ts_val - dis
    return ts_dict_Jnuisance

def plot_sigmavULs(hdf5file, output_dir, config):
        
    channels = config['Configuration']['channels']
    # Corresponding LaTex notation
    channels_LaTex = {'bb':'b\\bar{b}', 'tautau':'\\tau^{+}\\tau^{-}', 'mumu':'\mu^{+}\mu^{-}', 'tt':'t\\bar{t}', 'WW':'W^{+}W^{-}', 'gammagamma':'\gamma\gamma', 'hh':'H^{0}H^{0}', 'ZZ':'Z^{0}Z^{0}', 'ee':'e^{+}e^{-}'}
    
    for channel in channels:
        sigmavULs = pd.read_hdf(hdf5file, key='{}/sigmavULs'.format(channel))
        sigmavULs_Jnuisance = {}
        if config['Data']['J_nuisance']:
            sigmavULs_Jnuisance = pd.read_hdf(hdf5file, key='{}/sigmavULs_Jnuisance'.format(channel))
        for ul_dict in [sigmavULs,sigmavULs_Jnuisance]:
            if not dict(ul_dict):
                continue 
            masses = np.squeeze(ul_dict[['masses']].to_numpy())
            data = np.squeeze(ul_dict[['data']].to_numpy())
            if len(data)==0:
                continue
            y_low = 1e-26
            y_up = np.nanmax(data) * 2
        
            fig, ax = plt.subplots()
            ax.plot(masses,data,label='Combined limit',c='k')
            if config['Data']['cl_bands']:
                simulations = config['Data']['simulations']
                # Calculate the median (null hypothesis)
                null_hypothesis, sv_plus1, sv_minus1, sv_plus2, sv_minus2 = [],[],[],[],[]
                nh_index = np.int(simulations/2)
                svp1_index = simulations-np.int(0.5+16*simulations/100)
                svm1_index = np.int(16*simulations/100)
                svp2_index = simulations-np.int(0.5+2.5*simulations/100)
                svm2_index = np.int(2.5*simulations/100)
                for mass in np.arange(len(np.squeeze(ul_dict[['masses']].to_numpy(dtype=np.float32)))):
                    uls_mass = np.sort([np.squeeze(ul_dict[['simu{}'.format(simulation)]].to_numpy(dtype=np.float32))[mass] for simulation in np.arange(1,simulations)])
                    null_hypothesis.append(uls_mass[nh_index])
                    sv_plus1.append(uls_mass[svp1_index])
                    sv_minus1.append(uls_mass[svm1_index])
                    sv_plus2.append(uls_mass[svp2_index])
                    sv_minus2.append(uls_mass[svm2_index])

                y_up = np.nanmax(sv_plus2) * 2
                ax.plot(masses,null_hypothesis,label=r'$ H_{0} $ median',c='k',linewidth=0.75,linestyle='--')
                ax.fill_between(masses,sv_plus1,sv_minus1,color='green', alpha=0.5, linewidth=0)
                ax.fill_between(masses,sv_plus1,sv_plus2,color='yellow', alpha=0.5, linewidth=0)
                ax.fill_between(masses,sv_minus1,sv_minus2,color='yellow', alpha=0.5, linewidth=0)
            
                # Creating dummy data, which is needed for beautiful legend
                dummy_val = np.ones(len(masses))*1e-32
                plt.plot(masses,dummy_val,label='$ H_{0} \, 68\% $ containment',c='green', alpha=0.5, linewidth=6)
                plt.plot(masses,dummy_val,label=r'$ H_{0} \, 95\% $ containment',c='yellow', alpha=0.5, linewidth=6)

                # Plot the thermal relic, which was taken from Steigman G., Dasgupta B, and Beacom J. F., 
                # Precise relic WIMP abundance and its impact onsearches for dark matter annihilation, 
                # Phys.Rev. D86(2012) 023506, [arXiv:1204.3622]
                thermal_relic_mass = [1.00e-01, 1.78e-01, 3.16e-01, 5.62e-01, 1.00e+00, 1.78e+00, 3.16e+00, 5.62e+00, 1.00e+01, 1.78e+01, 3.16e+01, 5.62e+01, 1.00e+02, 1.78e+02, 3.16e+02, 5.62e+02, 1.00e+03, 1.78e+03,3.16e+03, 5.62e+03, 1.00e+04, 1.00e+05]
                thermal_relic_sigmav = [4.8e-26, 4.9e-26, 5.1e-26, 5.0e-26, 4.7e-26, 4.5e-26, 3.9e-26, 2.8e-26, 2.5e-26, 2.3e-26, 2.2e-26, 2.2e-26, 2.2e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26,2.4e-26, 2.4e-26]
                plt.plot(thermal_relic_mass,thermal_relic_sigmav,label=r'Thermal relic $\langle\sigma v\rangle$',c='r',linewidth=1.5,linestyle='--')
            
            ax.set_xscale('log')
            ax.set_xbound(lower=masses[0],upper=masses[-1])
            ax.set_xlabel(r'$m_{\chi} \: [GeV]$')
            ax.set_yscale('log')
            ax.set_ybound(lower=y_low,upper=y_up)
            ax.set_ylabel(r'$\langle\sigma v\rangle \, [cm^{3}/s]$')
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
    
def plot_sigmavULs_collaborations(hdf5file, output_dir, config):

    channels = config['Configuration']['channels']
    collaborations = config['Configuration']['collaborations']
    # Corresponding LaTex notation
    channels_LaTex = {'bb':'b\\bar{b}', 'tautau':'\\tau^{+}\\tau^{-}', 'mumu':'\mu^{+}\mu^{-}', 'tt':'t\\bar{t}', 'WW':'W^{+}W^{-}', 'gammagamma':'\gamma\gamma', 'ZZ':'Z^{0}Z^{0}', 'ee':'e^{+}e^{-}'}
    
    for channel in channels:
        for table_name in ['sigmavULs', 'sigmavULs_Jnuisance']:
            sigmavULs = pd.read_hdf(hdf5file, key='{}/{}'.format(channel,table_name))
            masses = np.squeeze(sigmavULs[['masses']].to_numpy())
            data = np.squeeze(sigmavULs[['data'.format(channel)]].to_numpy())
            if len(data)==0:
                continue
            
            fig, ax = plt.subplots()
            ax.plot(masses,data,label='Combined limit',c='k')
            for collaboration in collaborations:
                sigmavULs_col = pd.read_hdf(hdf5file, key='{}/{}/{}'.format(channel,collaboration,table_name))
                masses_col = np.squeeze(sigmavULs_col[['masses']].to_numpy())
                data_col = np.squeeze(sigmavULs_col[['data']].to_numpy())
                ax.plot(masses_col,data_col,label='{} limit'.format(collaboration))
                
            ax.set_xscale('log')
            ax.set_xbound(lower=masses[0],upper=masses[-1])
            ax.set_xlabel(r'$m_{\chi} \: [GeV]$')
            ax.set_yscale('log')
            ax.set_ylabel(r'$\langle\sigma v\rangle^{UL} \, [cm^{3}/s]$')
            if table_name == 'sigmavULs':
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
            if table_name == 'sigmavULs':
                plt.savefig('{}lklcom_{}_collaborations.pdf'.format(output_dir,channel))
                print("Saved plot in {}lklcom_{}_collaborations.pdf".format(output_dir,channel))
                plt.close()
            else:
                plt.savefig('{}lklcom_{}_Jnuisance_collaborations.pdf'.format(output_dir,channel))
                print("Saved plot in {}lklcom_{}_Jnuisance_collaborations.pdf".format(output_dir,channel))
                plt.close()
            ax.clear()

    return
