import argparse
import numpy as np
from multiprocessing import Process,Manager
import pandas as pd
import os
import yaml

from likelihood_combiner.combiner import *
from likelihood_combiner.utils import *


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description=("Combining likelihoods from different experiments."))
    parser.add_argument(
            'config_file',
            help="path to YAML configuration file with combining options")

    args = parser.parse_args()

    with open(args.config_file, 'r') as config_file:
        config = yaml.safe_load(config_file)

    try:
        output_dir = config['Output']['output_directory']
        if output_dir is None:
            raise KeyError
    except KeyError:
        output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output/"))
    # Create output directory if it doesn't exist already
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    try:
        hdf5filename = config['Output']['hdf5_dataset']
        if hdf5filename is None:
            raise KeyError
    except KeyError:
        hdf5filename = "lklcom.h5"

    # Check if hdf5file exist
    for file in os.listdir(output_dir):
        if file == hdf5filename:
            raise KeyError("'{}{}' already exist.".format(output_dir, hdf5filename))

    # Set up the path to the hdf5 file
    hdf5file = output_dir + hdf5filename

    # Combining limits using real observational data sets
    for channel in config['Configuration']['channels']:
        # Create a multiprocessing.Manager dict to share memory between the parallel processes
        manager = Manager()
        sigmavULs = manager.dict()
        sigmavULs_Jnuisance = manager.dict()
        print("Observations for the '{}' channel:".format(channel))
        combiner(config, channel, sigmavULs, sigmavULs_Jnuisance)
        print("Combined limits for observational data!")
        # Confidence level bands
        config['Data']['cl_bands'] = np.int(config['Data']['simulations']) > 0
        if config['Data']['cl_bands']:
            # Calculate the confidence level bands
            # Set up the hardware settings for the parallel processing
            try:
                cpu_counts = os.cpu_count() if config['Hardware']['cpu_counts'] == 'all' else config['Hardware']['cpu_counts']
                if cpu_counts is None:
                    raise KeyError
            except KeyError:
                cpu_counts = 1
            simulations = np.int(config['Data']['simulations'])
            if cpu_counts > simulations:
                cpu_counts = simulations

            # Set up all processes
            simulation_counter = manager.Value("i", 0)
            print("Combining {} simulations:".format(simulations))
            progress_bar(simulation_counter.value, simulations)
            parallel_simulations = np.array_split(np.arange(1,simulations+1), cpu_counts)
            jobs = []
            for parallel_simulation in parallel_simulations:
                process = Process(target=combiner, args=(config, channel, sigmavULs, sigmavULs_Jnuisance, simulation_counter, parallel_simulation))
                jobs.append(process)
            try:
                # Start all parallel processes
                for j in jobs:
                    j.start()
                # Wait for all processes to complete
                for j in jobs:
                    j.join()

            except KeyboardInterrupt:
                print("Caught keyboard interrupt, killing all processes...")
                for j in jobs:
                    j.terminate()

        # Convert multiprocessing.managers.DictProxy to python 'dict'
        sigmavULs = dict(sigmavULs)
        sigmavULs_Jnuisance= dict(sigmavULs_Jnuisance)

        svUL = {'masses': sigmavULs['{}_masses'.format(channel)]}
        for key, value in dict(sigmavULs).items():
            if channel in key: svUL[key.replace('{}_'.format(channel),'')] = value
        svUL = pd.DataFrame(data=svUL)
        # Write the panda DataFrames into the hdf5 file
        svUL.to_hdf(hdf5file, key='{}/sigmavULs'.format(channel), mode='a')
        if config['Data']['j_nuisance']:
            svUL_Jnuisance = {'masses': sigmavULs_Jnuisance['{}_masses'.format(channel)]}
            for key, value in dict(sigmavULs_Jnuisance).items():
                if channel in key: svUL_Jnuisance[key.replace('{}_'.format(channel),'')] = value
            svUL_Jnuisance = pd.DataFrame(data=svUL_Jnuisance)
            # Write the panda DataFrames into the hdf5 file
            svUL_Jnuisance.to_hdf(hdf5file, key='{}/sigmavULs_Jnuisance'.format(channel), mode='a')

        # Plot the sigmav upper limits
        plot_sigmavULs(hdf5file, output_dir, config, channel)

    if config['Output']['collaboration_plot'] and len(config['Configuration']['collaborations']) > 1:
        collaborations = config['Configuration']['collaborations']
        for collaboration in collaborations:
            config['Configuration']['collaborations'] = [collaboration]
            sigmavULs = manager.dict()
            sigmavULs_Jnuisance = manager.dict()
            for channel in config['Configuration']['channels']:
                combiner(config, channel, sigmavULs, sigmavULs_Jnuisance)
                print("Combined {} limits (observational data)".format(collaboration))
                # Convert multiprocessing.managers.DictProxy to python 'dict'
                sigmavULs = dict(sigmavULs)
                sigmavULs_Jnuisance= dict(sigmavULs_Jnuisance)
                svUL = {'masses': sigmavULs['{}_masses'.format(channel)]}
                svUL['data'] = sigmavULs['{}_data'.format(channel)]
                svUL = pd.DataFrame(data=svUL)
                # Write the panda DataFrames into the hdf5 file
                svUL.to_hdf(hdf5file, key='{}/{}/sigmavULs'.format(channel,collaboration), mode='a')
                if config['Data']['j_nuisance']:
                    svUL_Jnuisance = {'masses': sigmavULs_Jnuisance['{}_masses'.format(channel)]}
                    svUL_Jnuisance['data'] = sigmavULs_Jnuisance['{}_data'.format(channel)]
                    svUL_Jnuisance = pd.DataFrame(data=svUL_Jnuisance)
                    # Write the panda DataFrames into the hdf5 file
                    svUL_Jnuisance.to_hdf(hdf5file, key='{}/{}/sigmavULs_Jnuisance'.format(channel,collaboration), mode='a')

        config['Configuration']['collaborations'] = collaborations
        plot_sigmavULs_collaborations(hdf5file, output_dir, config)
