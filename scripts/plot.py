import argparse
import numpy as np
import pandas as pd
import os
import yaml

from likelihood_combiner.utils import *

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
            description=("Combining likelihoods from different experiments."))
    parser.add_argument(
            'config_file',
            help="path to YAML configuration file with combining options")

    parser.add_argument(
            '--channel')
        
    args = parser.parse_args()
                                     
    with open(args.config_file, 'r') as config_file:
        config = yaml.safe_load(config_file)
    channel = args.channel

    try:
        output_dir = config['Output']['output_directory']
        if output_dir is None:
            raise KeyError
    except KeyError:
        output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output/"))
    if not os.path.exists(output_dir):
         raise ValueError("Error 404: Directory '{}' not found. ".format(directory))
    files = np.array([x for x in os.listdir(output_dir) if x.endswith(".h5")])
    
    try:
        h5file = config['Output']['hdf5_dataset']
        if h5file is None:
            raise KeyError
    except KeyError:
        h5file = "lklcom.h5"
    
    # Confidence level bands
    config['Data']['cl_bands'] = np.int(config['Data']['simulations']) > 0
    # Merge h5 files
    if h5file not in files:
        svUL, svUL_Jnuisance = {}, {}
        data_file = h5file.replace('.h5','') + "_{}_data.h5".format(channel)
        data = pd.HDFStore(output_dir+data_file, 'r')
        if '/masses' in data.keys():
            svUL['masses'] = data['masses'][0]
            svUL['data'] = data['sigmavULs'][0]
            if config['Data']['j_nuisance']:
               svUL_Jnuisance['masses'] = data['masses'][0]
               svUL_Jnuisance['data'] = data['sigmavULs_Jnuisance'][0]
        simu_files = files[(data_file!=files)]
        for file in simu_files:
            data = pd.HDFStore(output_dir+file, 'r')
            run = file.replace('.h5','').split("_")[-1]
            if '/sigmavULs' in data.keys():
                svUL[run] = data['sigmavULs'][0]
            if '/sigmavULs_Jnuisance' in data.keys():
                    svUL_Jnuisance[run] = data['sigmavULs_Jnuisance'][0]
        svUL = pd.DataFrame(data=svUL)
        # Write the panda DataFrames into the hdf5 file
        pd.DataFrame(data=svUL).to_hdf(output_dir+h5file, key='{}/sigmavULs'.format(channel), mode='a')
        if config['Data']['j_nuisance']:
            # Write the panda DataFrames into the hdf5 file
            pd.DataFrame(data=svUL_Jnuisance).to_hdf(output_dir+h5file, key='{}/sigmavULs_Jnuisance'.format(channel), mode='a')
 
    # Plot the sigmav upper limits    
    plot_sigmavULs(output_dir+h5file, output_dir, config, channel)

    
