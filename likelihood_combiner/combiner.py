import argparse
import numpy as np
import os
import yaml

from likelihood_combiner.reader import gloryduckReader
from likelihood_combiner.writer import gloryduckWriter
from likelihood_combiner.gloryduck import gloryduckInfo

def run_combiner(config):
    
    try:
        data_dir = config['Data']['data_directory']
        if data_dir is None:
            raise KeyError
    except KeyError:
        data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data"))

    try:
        hdf5file = config['Data']['hdf5_dataset']
        if hdf5file is None:
            raise KeyError
    except KeyError:
        hdf5file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/gloryduck_dataset.h5"))

    writer = gloryduckWriter()

    writer.convert_txts2hdf5(hdf5file,data_dir)

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

    for channel in channels:
        for source in sources:
            # Checking that all ranges (mass and sigmav) and the J-Factor (first element of the mass arrays) are equal
            mass_ref = None
            sigmav_ref = None
            for key,mass,sigmav in zip(massvals.keys(),massvals.values(),tstables.values()):
                if channel in key and source in key:
                    if mass_ref is None:
                        mass_ref = mass
                        mass_key_ref = key
                    else:
                        if mass[0] != mass_ref[0]:
                            raise ValueError("The J-Factor value have to be equal! Discrepancy in '{}.txt' and '{}.txt'".format(key,mass_key_ref))
                        if (mass[1:]!=mass_ref[1:]).any():
                            raise ValueError("The mass values have to be equal! Discrepancy in '{}.txt' and '{}.txt'".format(key,mass_key_ref))
            
                    if sigmav_ref is None:
                        sigmav_ref = sigmav[0]
                        sigmav_key_ref = key
                    else:
                        if (sigmav[0]!=sigmav_ref).any():
                            raise ValueError("The sigma values have to be equal! Discrepancy in '{}.txt' and '{}.txt'".format(key,sigmav_key_ref))

    del writer
    del reader

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
