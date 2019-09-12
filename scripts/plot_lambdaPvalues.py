import argparse
import tables
import numpy as np
import os
import yaml
import matplotlib.pyplot as plt

from likelihood_combiner.gloryduck import gloryduckInfo
from likelihood_combiner.reader import gloryduckReader

def plot_lambdaPvalues(config):
    
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
        output_dir = config['Output']['output_directory'] + "lambdaP/"
        if output_dir is None:
            raise KeyError
    except KeyError:
        output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output/lambdaP/"))
    # Create output directory if it doesn't exist already
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    channels = config['Configuration']['channels']
    sources = config['Configuration']['sources']
    collaborations = config['Configuration']['collaborations']
    
    # Get information from the gloryduck class
    gloryduck = gloryduckInfo()
    if channels is None:
        channels = gloryduck.channels
    channels = np.array(channels)
    if sources is None:
        sources = gloryduck.sources
    sources = np.array(sources)
    if collaborations is None:
        collaborations = gloryduck.collaborations
    collaborations = np.array(collaborations)
    channels_LaTex = gloryduck.channels_LaTex

    gd_reader = gloryduckReader()
    tstables, massvals = gd_reader.read_gloryduck_tstables(hdf5file,channels,sources,collaborations)
    del gd_reader

    for key,mass,tstable in zip(massvals.keys(),massvals.values(),tstables.values()):
        for i in np.arange(1,mass.shape[0]):
            if np.all(tstable[i]==0):
                continue
            
            fig, ax = plt.subplots()
            ax.plot(tstable[0],tstable[i],label=r'$\lambda_{p}$')
            ax.plot(tstable[0],2.71*np.ones(tstable[i].shape[0]),'r--',label=r'$95\%$ CL')

            ax.set_xscale('log')
            ax.set_xbound(lower=np.min(tstable[0]),upper=np.max(tstable[0]))
            ax.set_ybound(lower=-2.5,upper=5.)
            ax.legend()
            ax.set_xlabel(r'$\langle \sigma v\rangle \, [cm^{3}/s]$')
            ax.set_ylabel(r'$\lambda_{p}$')
            ax.set_title(r'$\lambda_{p}$ vs $\langle \sigma v\rangle$')
            
            # Put a legend to the right of the current axis
            ax.grid(b=True,which='both',color='grey', linestyle='--', linewidth=0.25)
            plt.savefig('{}{}_{}.pdf'.format(output_dir,key,mass[i]))
            print("Saved plot in {}{}_{}.pdf".format(output_dir,key,mass[i]))
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

    plot_lambdaPvalues(config)
