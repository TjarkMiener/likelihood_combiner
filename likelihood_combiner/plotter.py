import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

__all__ = [
    'plot_sigmavULs',
    'plot_sigmavULs_collaborations',
]

def plot_sigmavULs(hdf5file, output_dir, config, channel=None):

    if channel is None:
        channels = config['Configuration']['channels']
    else:
        channels = [channel]

    # Corresponding LaTex notation
    channels_LaTex = {
        "bb":"b\\bar{b}",
        "tautau":"\\tau^{+}\\tau^{-}",
        "mumu":"\mu^{+}\mu^{-}",
        "tt":"t\\bar{t}",
        "WW":"W^{+}W^{-}",
        "gammagamma":"\gamma\gamma",
        "hh":"hh",
        "ZZ":"ZZ",
        "ee":"e^{+}e^{-}"
    }

    for channel in channels:
        print("LklCom results:")
        tables = ['sigmavULs']
        if config['Data']['j_nuisance']:
            tables.append('sigmavULs_Jnuisance')
        for table in tables:
            sigmavULs = pd.read_hdf(hdf5file, key='{}/{}'.format(channel,table))
            masses = np.squeeze(sigmavULs[['masses']].to_numpy())
            data = np.squeeze(sigmavULs[['data'.format(channel)]].to_numpy())

            y_low = 1e-26
            y_up = np.nanmax(data) * 2

            fig, ax = plt.subplots()
            print("  observational_limits: {}".format(data.tolist()))
            ax.plot(masses,data,label='Combined limit',c='k')
            if config['Data']['cl_bands']:
                simulations = len(sigmavULs.columns)-2
                # Calculate the median (null hypothesis)
                simu_tab = sigmavULs.drop(['data','masses'], axis=1).to_numpy()
                null_hypothesis = np.mean(simu_tab, axis=1)
                sv_plus1 = np.percentile(simu_tab, 84.14, axis=1)
                sv_plus2 = np.percentile(simu_tab, 97.725, axis=1)
                sv_minus1 = np.percentile(simu_tab, 15.87, axis=1)
                sv_minus2 = np.percentile(simu_tab, 2.275, axis=1)

                print("  null_hypothesis: {}".format(null_hypothesis))
                print("  cl_minus1: {}".format(sv_minus1))
                print("  cl_minus2: {}".format(sv_minus2))
                print("  cl_plus1: {}".format(sv_plus1))
                print("  cl_plus2: {}".format(sv_plus2))

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
            if table == 'sigmavULs':
                ax.set_title(r'$\langle\sigma v\rangle$ ULs vs mass - J fixed')
            else:
                ax.set_title(r'$\langle\sigma v\rangle$ ULs vs mass - J as nuisance')

            if len(config['Configuration']['sources']) == 1:
                ax.text(0.2, 0.85, config['Configuration']['sources'][0], fontsize=18,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            else:
                ax.text(0.2, 0.85, 'All dSphs', fontsize=18,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

            ax.text(0.85, 0.1, r'$\chi\chi \to {}$'.format(channels_LaTex[str(channel)]), fontsize=15,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
            # Shrink current axis by 20%
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8)
            ax.grid(b=True,which='both',color='grey', linestyle='--', linewidth=0.25)
            if table == 'sigmavULs':
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
    channels_LaTex = {
        "bb":"b\\bar{b}",
        "tautau":"\\tau^{+}\\tau^{-}",
        "mumu":"\mu^{+}\mu^{-}",
        "tt":"t\\bar{t}",
        "WW":"W^{+}W^{-}",
        "gammagamma":"\gamma\gamma",
        "hh":"hh",
        "ZZ":"ZZ",
        "ee":"e^{+}e^{-}"
    }

    for channel in channels:
        tables = ['sigmavULs']
        if config['Data']['j_nuisance']:
            tables.append('sigmavULs_Jnuisance')
        for table in tables:
            sigmavULs = pd.read_hdf(hdf5file, key='{}/{}'.format(channel,table))
            masses = np.squeeze(sigmavULs[['masses']].to_numpy())
            data = np.squeeze(sigmavULs[['data'.format(channel)]].to_numpy())

            fig, ax = plt.subplots()
            ax.plot(masses,data,label='Combined limit',c='k')
            for collaboration in collaborations:
                sigmavULs_col = pd.read_hdf(hdf5file, key='{}/{}/{}'.format(channel,collaboration,table))
                masses_col = np.squeeze(sigmavULs_col[['masses']].to_numpy())
                data_col = np.squeeze(sigmavULs_col[['data']].to_numpy())
                ax.plot(masses_col,data_col,label='{} limit'.format(collaboration))

            ax.set_xscale('log')
            ax.set_xbound(lower=masses[0],upper=masses[-1])
            ax.set_xlabel(r'$m_{\chi} \: [GeV]$')
            ax.set_yscale('log')
            ax.set_ylabel(r'$\langle\sigma v\rangle^{UL} \, [cm^{3}/s]$')
            if table == 'sigmavULs':
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
            if table == 'sigmavULs':
                plt.savefig('{}lklcom_{}_collaborations.pdf'.format(output_dir,channel))
                print("Saved plot in {}lklcom_{}_collaborations.pdf".format(output_dir,channel))
                plt.close()
            else:
                plt.savefig('{}lklcom_{}_Jnuisance_collaborations.pdf'.format(output_dir,channel))
                print("Saved plot in {}lklcom_{}_Jnuisance_collaborations.pdf".format(output_dir,channel))
                plt.close()
            ax.clear()

    return
