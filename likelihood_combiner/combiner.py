import argparse
import numpy as np
import os
import yaml
import matplotlib.pyplot as plt

from likelihood_combiner.reader import gloryduckReader,JFactor_Reader
from likelihood_combiner.writer import gloryduckWriter
from likelihood_combiner.gloryduck import gloryduckInfo
from likelihood_combiner.utils import compute_sensitivity,compute_Jnuisance

def run_combiner(config):
    
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
        hdf5file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/gloryduck_dataset.h5"))

    try:
        output_dir = config['Output']['output_directory']
        if output_dir is None:
            raise KeyError
    except KeyError:
        output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output/"))
    # Create output directory if it doesn't exist already
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if config['Data']['overwrite_hdf5file']:
        gd_writer = gloryduckWriter()
        gd_writer.convert_txts2hdf5(hdf5file,data_dir)
        del gd_writer

    channels = config['Configuration']['channels']
    sources = config['Configuration']['sources']
    collaborations = config['Configuration']['collaborations']
    
    # Get information from the gloryduck class
    gloryduck = gloryduckInfo()
    if channels is None:
        channels = np.array(gloryduck.channels)
    if sources is None:
        sources = np.array(gloryduck.sources)
    if collaborations is None:
        collaborations = np.array(gloryduck.collaborations)
    channels_LaTex = gloryduck.channels_LaTex

    gd_reader = gloryduckReader()
    tstables, massvals = gd_reader.read_gloryduck_tstables(hdf5file,channels,sources,collaborations)
    del gd_reader

    try:
        JFactor_file = config['Data']['JFactor_table']
        if JFactor_file is None:
            raise KeyError
    except KeyError:
        JFactor_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/Jfactor_Geringer-SamethTable.txt"))

    if config['Data']['J_nuisance']:
        Jfactor_reader = JFactor_Reader()
        sources_logJ,sources_DlogJ = Jfactor_reader.read_JFactor(JFactor_file,sources)
        del Jfactor_reader
    
    print("Combining limits for the selected configuration:")
    for channel in channels:
        sigmav = None
        mass_axis = []
        # Combine ts values for each dSphs
        combined_sources_ts = {}
        for source in sources:
            tstable_ref = None
            for key,mass,tstable in zip(massvals.keys(),massvals.values(),tstables.values()):
                key_split = key.split("_")
                if channel == key_split[0] and source == key_split[1]:
                    # Checking that all ranges (mass and sigmav) and the J-Factor (first element of the mass arrays) are equal
                    if tstable_ref is None:
                        tstable_ref = sigmav = tstable[0]
                        exponent = (np.abs(np.floor(np.log10(np.abs(sigmav))).astype(int))+3).astype(int)
                        for i,e in enumerate(exponent):
                            sigmav[i] = np.around(sigmav[i],decimals=e)
                        tstable_key_ref = key
                    else:
                        if (sigmav!=tstable_ref).any():
                            raise ValueError("The sigma values have to be equal! Discrepancy in '{}.txt' and '{}.txt'".format(key,tstable_key_ref))
                    for i,m in enumerate(mass[1:]):
                        if m not in mass_axis:
                            mass_axis.append(m)
                        if source+"_"+str(m) in combined_sources_ts:
                            combined_sources_ts[source+"_"+str(m)] += tstable[i+1][::-1]
                        else:
                            combined_sources_ts[source+"_"+str(m)] = tstable[i+1][::-1]
        sigmav = sigmav[::-1]
        if config['Data']['J_nuisance']:
            combined_sources_ts_Jnuisance = compute_Jnuisance(sigmav, combined_sources_ts, sources_DlogJ)
            combined_sources_dict = compute_sensitivity(sigmav, combined_sources_ts_Jnuisance)
        else:
            combined_sources_dict = compute_sensitivity(sigmav, combined_sources_ts)
        combined_sources = []
        combined_sources_masses = []
        combined_sources_limits = []
        for source in sources:
            source_mass = []
            source_limit = []
            for mass,limit in combined_sources_dict.items():
                if source in mass:
                    source_mass.append(float(mass[mass.find("_")+1:]))
                    source_limit.append(limit)
            combined_sources.append(source)
            combined_sources_masses.append(source_mass)
            combined_sources_limits.append(source_limit)
        combined_sources = np.array(combined_sources)
        combined_sources_masses = np.array(combined_sources_masses)
        combined_sources_limits = np.array(combined_sources_limits)
    
        # Combine ts values for each collaboration
        combined_collaborations_ts = {}
        for collaboration in collaborations:
            for key,mass,tstable in zip(massvals.keys(),massvals.values(),tstables.values()):
                if channel in key and collaboration in key:
                    for i,m in enumerate(mass[1:]):
                        if m not in mass_axis:
                            mass_axis.append(m)
                        if collaboration+"_"+str(m) in combined_collaborations_ts:
                            combined_collaborations_ts[collaboration+"_"+str(m)] += tstable[i+1][::-1]
                        else:
                            combined_collaborations_ts[collaboration+"_"+str(m)] = tstable[i+1][::-1]

        combined_collaborations_dict = compute_sensitivity(sigmav, combined_collaborations_ts)
        combined_collaborations = []
        combined_collaborations_masses = []
        combined_collaborations_limits = []
        for collaboration in collaborations:
            collaboration_mass = []
            collaboration_limit = []
            for mass,limit in combined_collaborations_dict.items():
                if collaboration in mass:
                    collaboration_mass.append(float(mass[mass.find("_")+1:]))
                    collaboration_limit.append(limit)
            combined_collaborations.append(collaboration)
            combined_collaborations_masses.append(collaboration_mass)
            combined_collaborations_limits.append(collaboration_limit)
        combined_collaborations = np.array(combined_collaborations)
        combined_collaborations_masses = np.array(combined_collaborations_masses)
        combined_collaborations_limits = np.array(combined_collaborations_limits)

        # Combine ts values for all dSphs
        combined_all_dSphs_ts = {}
        for m in mass_axis:
            for comb_source,tstable in combined_collaborations_ts.items():
                if str(m) in comb_source:
                    if str(m) in combined_all_dSphs_ts:
                        combined_all_dSphs_ts[str(m)] += tstable
                    else:
                        combined_all_dSphs_ts[str(m)] = tstable
        combined_all_dSphs_dict = compute_sensitivity(sigmav, combined_all_dSphs_ts)
        combined_all_dSphs_masses = [ float(m) for m in combined_all_dSphs_dict ]
        combined_all_dSphs_limits = [ l for l in combined_all_dSphs_dict.values() ]

        # Plotting
        # Combined with sources
        fig, ax = plt.subplots()
        for i in np.arange(combined_sources_masses.shape[0]):
            ax.plot(combined_sources_masses[i],combined_sources_limits[i],label='{} limit'.format(combined_sources[i]))
        ax.plot(combined_all_dSphs_masses, combined_all_dSphs_limits, label='Combined limit')
        ax.set_xscale('log')
        ax.set_xbound(lower=np.min(combined_all_dSphs_masses),upper=np.max(combined_all_dSphs_masses))
        ax.set_xlabel(r'$m_{\chi} \: [GeV]$')
        ax.set_yscale('log')
        ax.set_ylabel(r'$95\%$ CL $\langle\sigma v\rangle^{UL} \, [cm^{3}/s]$')
        ax.set_title(r'$\langle\sigma v\rangle$ ULs vs mass')
        ax.text(0.4, 0.85, 'All dSphs', fontsize=18,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.85, 0.1, r'$\chi\chi \to {}$'.format(channels_LaTex[str(channel)]), fontsize=15,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8)
        #fig.tight_layout()
        #ax.legend(loc='upper right')
        ax.grid(b=True,which='both',color='grey', linestyle='--', linewidth=0.25)
        plt.savefig('{}/{}_alldSph_withsources.pdf'.format(output_dir,channel))
        print("Saved plot in {}/{}_alldSph_withsources.pdf".format(output_dir,channel))
        ax.clear()

        fig, ax = plt.subplots()
        for i in np.arange(combined_collaborations_masses.shape[0]):
            ax.plot(combined_collaborations_masses[i],combined_collaborations_limits[i],label='{} limit'.format(combined_collaborations[i]))
        ax.plot(combined_all_dSphs_masses, combined_all_dSphs_limits, label='Combined limit')
        ax.set_xscale('log')
        ax.set_xbound(lower=np.min(combined_all_dSphs_masses),upper=np.max(combined_all_dSphs_masses))
        ax.set_xlabel(r'$m_{\chi} \: [GeV]$')
        ax.set_yscale('log')
        ax.set_ylabel(r'$95\%$ CL $\langle\sigma v\rangle^{UL} \, [cm^{3}/s]$')
        ax.set_title(r'$\langle\sigma v\rangle$ ULs vs mass')
        ax.text(0.4, 0.85, 'All dSphs', fontsize=18,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.85, 0.1, r'$\chi\chi \to {}$'.format(channels_LaTex[str(channel)]), fontsize=15,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8)
        #fig.tight_layout()
        #ax.legend(loc='upper right')
        ax.grid(b=True,which='both',color='grey', linestyle='--', linewidth=0.25)
        plt.savefig('{}/{}_alldSph_withcollaborations.pdf'.format(output_dir,channel))
        print("Saved plot in {}/{}_alldSph_withcollaborations.pdf".format(output_dir,channel))
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

    run_combiner(config)
