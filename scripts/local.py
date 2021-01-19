import argparse
import numpy as np
from multiprocessing import Process, Manager
import pandas as pd
import os
import yaml

import likelihood_combiner as lklcom
from likelihood_combiner.combiner import combiner
from likelihood_combiner.utils import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description=("Combining likelihoods from different experiments."))
    parser.add_argument(
            'config_file',
            help="path to YAML configuration file with combining options")
    parser.add_argument(
            '--input',
            help="path to input file or directory")
    parser.add_argument(
            '--output',
            help="path to output file or directory")

    args = parser.parse_args()
    
    with open(args.config_file, 'r') as config_file:
        config = yaml.safe_load(config_file)
    
    # Initializing of the LklCom jfactor class
    if config['Data']['buildin_j_factors'] == "GeringerSameth":
        LklCom_jfactor_class = lklcom.jfactor.GeringerSameth(sources=config['Configuration']['sources'],
                                                    collaborations=config['Configuration']['collaborations'],
                                                    combination_data=args.input,
                                                    jnuisance=config['Data']['j_nuisance'])
    elif config['Data']['buildin_j_factors'] == "Bonnivard":
        LklCom_jfactor_class = lklcom.jfactor.Bonnivard(sources=config['Configuration']['sources'],
                                                    collaborations=config['Configuration']['collaborations'],
                                                    combination_data=args.input,
                                                    jnuisance=config['Data']['j_nuisance'])
    else:
        LklCom_jfactor_class = lklcom.jfactor.Custom(logJ=config['Data']['custom_logJ'],
                                                    DlogJ=config['Data']['custom_DlogJ'],
                                                    jnuisance=config['Data']['j_nuisance'])

    # Constructing in the the sigmav range and spacing
    sigmav_range = get_sigmav_range()
    # Overwriting from the config file, if provided.
    if "sigmav_min" in config['Data'] and "sigmav_max" in config['Data'] and "sigmav_nPoints" in config['Data']:
        sigmav_range = get_sigmav_range(config['Data']['sigmav_min'], config['Data']['sigmav_max'], config['Data']['sigmav_nPoints'], 3)
    
    for channel in config['Configuration']['channels']:
        print("\nChannel '{}'".format(channel))
        
        # Initializing of the LklCom reader class
        if args.input.endswith(".h5") or args.input.endswith(".hdf5"):
            LklCom_reader_class = lklcom.reader.LklCom_hdf5(channel=channel,
                                                            LklCom_jfactor_class=LklCom_jfactor_class)
        if os.path.isdir(args.input):
            LklCom_reader_class = lklcom.reader.LklCom_txtdir(channel=channel,
                                                            LklCom_jfactor_class=LklCom_jfactor_class)
        
        # Set up the hardware settings for the parallel processing
        try:
            cpu_counts = os.cpu_count() if config['Hardware']['cpu_counts'] == 'all' or config['Hardware']['cpu_counts'] > os.cpu_count() else config['Hardware']['cpu_counts']
            if cpu_counts is None:
                raise KeyError
        except KeyError:
            cpu_counts = 1
        if 'simulations' in config['Data']:
            simulations = np.int(config['Data']['simulations']+1)
        else:
            simulations = 1
        if cpu_counts > simulations:
            cpu_counts = simulations

        # Set up all processes
        if simulations <= 1:
            combiner(sigmav_range=sigmav_range,
                    LklCom_reader_class=LklCom_reader_class,
                    output=args.output)
        else:
            # Create a multiprocessing.Manager dict to share memory between the parallel processes
            manager = Manager()
            sigmavULs = manager.dict()
            sigmavULs_Jnuisance = manager.dict()
            
            simulation_counter = manager.Value("i", 0)
            progress_bar(simulation_counter.value, simulations)
            parallel_simulations = np.array_split(np.arange(0,simulations), cpu_counts)
            jobs = []
            for parallel_simulation in parallel_simulations:
                process = Process(target=combiner,
                                    args=(sigmav_range,
                                        LklCom_reader_class,
                                        args.output,
                                        sigmavULs,
                                        sigmavULs_Jnuisance,
                                        simulation_counter,
                                        simulations,
                                        parallel_simulation))
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
            svUL.to_hdf(args.output, key='{}/sigmavULs'.format(channel), mode='a')
            if config['Data']['j_nuisance']:
                svUL_Jnuisance = {'masses': sigmavULs_Jnuisance['{}_masses'.format(channel)]}
                for key, value in dict(sigmavULs_Jnuisance).items():
                    if channel in key: svUL_Jnuisance[key.replace('{}_'.format(channel),'')] = value
                svUL_Jnuisance = pd.DataFrame(data=svUL_Jnuisance)
                # Write the panda DataFrames into the hdf5 file
                svUL_Jnuisance.to_hdf(args.output, key='{}/sigmavULs_Jnuisance'.format(channel), mode='a')
