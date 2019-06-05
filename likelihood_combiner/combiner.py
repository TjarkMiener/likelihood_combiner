import argparse
import numpy as np
import os
import yaml
import matplotlib.pyplot as plt

from likelihood_combiner.reader import gloryduckReader
from likelihood_combiner.writer import gloryduckWriter
from likelihood_combiner.gloryduck import gloryduckInfo
from likelihood_combiner.sensitivity import compute_sensitivity

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
        writer = gloryduckWriter()
        writer.convert_txts2hdf5(hdf5file,data_dir)
        del writer

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

    reader = gloryduckReader()
    tstables, massvals = reader.read_gloryduck_tstables(hdf5file,channels,sources,collaborations)
    del reader

    for channel in channels:
        for source in sources:
            mass_ref = None
            tstable_ref = None
            combined_ts = None
            for key,mass,tstable in zip(massvals.keys(),massvals.values(),tstables.values()):
                if channel in key and source in key:
                    # Checking that all ranges (mass and sigmav) and the J-Factor (first element of the mass arrays) are equal
                    if mass_ref is None:
                        mass_ref = mass
                        mass_key_ref = key
                    else:
                        if mass[0] != mass_ref[0]:
                            print("Warning: Discrepancy in the J-Factor value of '{}.txt' and '{}.txt'! Make sure each collabration took the right J-Factor value from the GS table.".format(key,mass_key_ref))
                        if (mass[1:]!=mass_ref[1:]).any():
                            raise ValueError("The mass values have to be equal! Discrepancy in '{}.txt' and '{}.txt'".format(key,mass_key_ref))
            
                    if tstable_ref is None:
                        tstable_ref = tstable[0]
                        tstable_key_ref = key
                    else:
                        if (tstable[0]!=tstable_ref).any():
                            raise ValueError("The sigma values have to be equal! Discrepancy in '{}.txt' and '{}.txt'".format(key,tstable_key_ref))
                    sigmav = tstable[0]
                    if combined_ts is None:
                        combined_ts = np.zeros((tstable.shape[0]-1, tstable.shape[1]))
                    for i in np.arange(tstable.shape[0]-1):
                        # Combining the likelihood values
                        combined_ts[i] += tstable[i+1]

                    # Calculating the limits
                    limits = compute_sensitivity(sigmav, tstable[1:])
                    print(limits)

            combined_limit = compute_sensitivity(sigmav, combined_ts)
            print(combined_limit)

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
