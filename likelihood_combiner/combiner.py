import argparse
import numpy as np
import os
import yaml

from likelihood_combiner.reader import gloryduckReader
from likelihood_combiner.writer import gloryduckWriter

def run_combiner(config):
    
    print(config)
    
    try:
        data_dir = config['Data']['data_directory']
        if data_dir is None:
            raise KeyError
    except KeyError:
        data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data"))
    # Check if the data directory is existing
    if not os.path.exists(data_dir):
            raise ValueError("'{}' is not existing!".format(data_dir))
    try:
        hdf5file = config['Data']['hdf5_dataset']
        if hdf5file is None:
            raise KeyError
    except KeyError:
        hdf5file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/gloryduck_dataset.h5"))
    # Check if the data directory is existing
    if not os.path.exists(hdf5file):
        raise ValueError("'{}' is not existing!".format(hdf5file))

    writer = gloryduckWriter()

    writer.convert_txts2hdf5(hdf5file,data_dir)

    channels = config['Configuration']['channels']
    sources = config['Configuration']['sources']
    collaborations = config['Configuration']['collaborations']

    reader = gloryduckReader()

    tstables, massvals = reader.read_gloryduck_tstables(hdf5file,channels,sources,collaborations)

    print(tstables)
    print(massvals)

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
