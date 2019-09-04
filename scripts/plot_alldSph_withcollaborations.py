import argparse
import tables
import numpy as np
import os
import yaml
import matplotlib.pyplot as plt

from likelihood_combiner.gloryduck import gloryduckInfo
from likelihood_combiner.reader import gloryduckReader

def plot_alldSph_withcollaborations(config):
    
    try:
        data_dir = config['Data']['data_directory']
        if data_dir is None:
            raise KeyError
    except KeyError:
        data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/"))
    
    try:
        hdf5file = config['Data']['hdf5_dataset']
        if hdf5file is None:
            raise KeyError
    except KeyError:
        hdf5file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/gloryduck_all.h5"))
    
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
    collaborations = config['Configuration']['collaborations']
    
    # Get information from the gloryduck class
    gloryduck = gloryduckInfo()
    if channels is None:
        channels = np.array(gloryduck.channels)
    if collaborations is None:
        collaborations = np.array(gloryduck.collaborations)
    channels_LaTex = gloryduck.channels_LaTex

    gd_reader = gloryduckReader()
    sigmavULs, sigmavULs_Jnuisance, massvals = gd_reader.read_gloryduck_sigmavULs_collaborations(data_dir, hdf5file, channels, collaborations)
    del gd_reader
  
    for channel in channels:
        for ul_dict in [sigmavULs,sigmavULs_Jnuisance]:
            if ul_dict is sigmavULs_Jnuisance and not config['Data']['J_nuisance']:
                continue
            fig, ax = plt.subplots()
            for collaboration in collaborations:
                if channel+"_"+collaboration in ul_dict:
                    ax.plot(massvals[channel+"_"+collaboration],ul_dict[channel+"_"+collaboration],label='{} limit'.format(collaboration))
            ax.plot(massvals[channel+"_Combination_all"], ul_dict[channel+"_Combination_all"], label='Combined limit')
            ax.set_xscale('log')
            ax.set_xbound(lower=np.min(massvals[channel+"_Combination_all"]),upper=np.max(massvals[channel+"_Combination_all"]))
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
                plt.savefig('{}/{}_alldSph_withcollaborations.pdf'.format(output_dir,channel))
                print("Saved plot in {}/{}_alldSph_withcollaborations.pdf".format(output_dir,channel))
                plt.close()
            else:
                plt.savefig('{}/{}_alldSph_withcollaborations_Jnuisance.pdf'.format(output_dir,channel))
                print("Saved plot in {}/{}_alldSph_withcollaborations_Jnuisance.pdf".format(output_dir,channel))
                plt.close()
            ax.clear()



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
            description=("Combining likelihoods from different experiments."))
    parser.add_argument(
            'config_file',
            help="path to YAML configuration file with combining options")
    args = parser.parse_args()
                            
    with open(args.config_file, 'r') as config_file:
        config = yaml.load(config_file)

    plot_alldSph_withcollaborations(config)
