import argparse
import numpy as np
from multiprocessing import Process,Manager
import pandas as pd
from scipy.interpolate import interp1d
import os
import yaml

from likelihood_combiner.reader import LklComReader
from likelihood_combiner.utils import *

def run_combiner(config, channel, sigmavULs, sigmavULs_Jnuisance, simulation_counter=None, simulations=[-1]):

    try:
        data_dir = config['Data']['data_directory']
        if data_dir is None:
            raise KeyError
    except KeyError:
        data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/"))


    sources = config['Configuration']['sources']
    collaborations = config['Configuration']['collaborations']
    jnuisance = config['Data']['J_nuisance']

    # Initializing the LklComReader
    reader = LklComReader(channel, sources, collaborations)
    
    if jnuisance:
        try:
            JFactor_file = config['Data']['JFactor_table']
            if JFactor_file is None:
                raise KeyError
        except KeyError:
            JFactor_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/Jfactor_Geringer-SamethTable.txt"))
            
        sources_logJ,sources_DlogJ = reader.read_JFactor(JFactor_file)

    # Reading in the the sigmav range and spacing. Round of the third digit to avoid an interpolation.
    sigmavMin = -(np.abs(np.floor(np.log10(np.abs(float(config['Data']['sigmaV_min'])))).astype(int))).astype(int)
    sigmavMax = -(np.abs(np.floor(np.log10(np.abs(float(config['Data']['sigmaV_max'])))).astype(int))).astype(int)
    sigmavNPoints = int(config['Data']['sigmaV_nPoints'])
    sigmav = np.logspace(sigmavMin, sigmavMax, sigmavNPoints, endpoint=True, dtype=np.float32)
    exponent = (np.abs(np.floor(np.log10(np.abs(sigmav))).astype(int))+3).astype(int)
    for i,e in enumerate(exponent):
        sigmav[i] = np.around(sigmav[i],decimals=e)

    for simulation in simulations:
    
        tstables = reader.read_tstables(data_dir, simulation)


        combined_ts = {}
        combined_ts_Jnuisance = {}
        mass_axis = []
        for source in sources:
            # Combine ts values for each dSphs
            combined_source_ts = {}
            combined_source_ts_Jnuisance = {}
            for collaboration in collaborations:
                key = "{}_{}_{}".format(channel,source,collaboration)
                if key+'_ts' in tstables:
                    tstable = tstables[key+'_ts']
                    masses = tstables[key+'_masses']
                    sigmavFile = tstable[0]
                    exponent = (np.abs(np.floor(np.log10(np.abs(sigmavFile))).astype(int))+3).astype(int)
                    for i,e in enumerate(exponent):
                        sigmavFile[i] = np.around(sigmavFile[i],decimals=e)
                    sigmavFile = sigmavFile[::-1]
                    source_ts_dict = {}
                    for i,m in enumerate(masses[1:]):
                        if m not in mass_axis:
                            mass_axis.append(m)
                        if not np.array_equal(sigmav,sigmavFile):
                            lin_interpolation = interp1d(sigmavFile, tstable[i+1][::-1], kind='linear', fill_value='extrapolate')
                            source_ts_dict[source+"_"+str(m)] = lin_interpolation(sigmav)
                        else:
                            source_ts_dict[source+"_"+str(m)] = tstable[i+1][::-1]
                        if source+"_"+str(m) in combined_source_ts:
                            combined_source_ts[source+"_"+str(m)] += source_ts_dict[source+"_"+str(m)]
                        else:
                            combined_source_ts[source+"_"+str(m)] = source_ts_dict[source+"_"+str(m)]
        
            if jnuisance:
                combined_source_ts_Jnuisance = compute_Jnuisance(sigmav, combined_source_ts, sources_DlogJ)
                combined_source_limits_Jnuisance,combined_source_sensitivity_Jnuisance = compute_sensitivity(sigmav, combined_source_ts_Jnuisance)
            for m in mass_axis:
                if str(m) in combined_ts:
                    if source+"_"+str(m) in combined_source_ts:
                        combined_ts[str(m)] += combined_source_ts[source+"_"+str(m)]
                        if jnuisance:
                            combined_ts_Jnuisance[str(m)] += combined_source_ts_Jnuisance[source+"_"+str(m)]
                else:
                    if source+"_"+str(m) in combined_source_ts:
                        combined_ts[str(m)] = combined_source_ts[source+"_"+str(m)]
                        if jnuisance:
                            combined_ts_Jnuisance[str(m)] = combined_source_ts_Jnuisance[source+"_"+str(m)]
            
        combined_sources_limits,combined_sources_sensitivity = compute_sensitivity(sigmav, combined_ts)
        if jnuisance:
            combined_sources_limits_Jnuisance,combined_sources_sensitivity_Jnuisance = compute_sensitivity(sigmav, combined_ts_Jnuisance)
                
        # Filling the data of the txt file into the table of the hdf5 file.
        sorted_mass = sorted(mass_axis)
        svUL = []
        svUL_Jnuisance = []
        for m in sorted_mass:
            svUL.append(combined_sources_limits[str(m)])
            if jnuisance:
                svUL_Jnuisance.append(combined_sources_limits_Jnuisance[str(m)])
            else:
                svUL_Jnuisance.append(np.nan)

        if simulation >= 0:
            h5_key = "{}_simu{}".format(channel, simulation)
        else:
            h5_key = "{}_data".format(channel)
            sigmavULs["{}_masses".format(channel)] = sorted_mass
            sigmavULs_Jnuisance["{}_masses".format(channel)] = sorted_mass

        sigmavULs[h5_key] = svUL
        sigmavULs_Jnuisance[h5_key] = svUL_Jnuisance
            
        if simulation >= 0:
            simulation_counter.value += 1
            progress_bar(simulation_counter.value,np.int(config['Data']['simulations']))

    del reader
    return

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

    # Create a multiprocessing.Manager dict to share memory between the parallel processes
    manager = Manager()
    sigmavULs = manager.dict()
    sigmavULs_Jnuisance = manager.dict()
    # Combining limits using real observational data sets
    for channel in config['Configuration']['channels']:
        print("Observations for the '{}' channel:".format(channel))
        run_combiner(config, channel, sigmavULs, sigmavULs_Jnuisance)
        print("Combined limits for observational data!")
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
                process = Process(target=run_combiner, args=(config, channel, sigmavULs, sigmavULs_Jnuisance, simulation_counter, parallel_simulation))
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
    for channel in config['Configuration']['channels']:
        svUL = {'masses': sigmavULs['{}_masses'.format(channel)]}
        for key, value in dict(sigmavULs).items():
            if channel in key: svUL[key.replace('{}_'.format(channel),'')] = value
        svUL = pd.DataFrame(data=svUL)
        # Write the panda DataFrames into the hdf5 file
        svUL.to_hdf(hdf5file, key='{}/sigmavULs'.format(channel), mode='a')
        if config['Data']['J_nuisance']:
            svUL_Jnuisance = {'masses': sigmavULs_Jnuisance['{}_masses'.format(channel)]}
            for key, value in dict(sigmavULs_Jnuisance).items():
                if channel in key: svUL_Jnuisance[key.replace('{}_'.format(channel),'')] = value
            svUL_Jnuisance = pd.DataFrame(data=svUL_Jnuisance)
            # Write the panda DataFrames into the hdf5 file
            svUL_Jnuisance.to_hdf(hdf5file, key='{}/sigmavULs_Jnuisance'.format(channel), mode='a')

    # Plot the sigmav upper limits    
    plot_sigmavULs(hdf5file, output_dir, config)
    
    if config['Output']['collaboration_plot'] and len(config['Configuration']['collaborations']) > 1:
        collaborations = config['Configuration']['collaborations']
        for collaboration in collaborations:
            config['Configuration']['collaborations'] = [collaboration]
            sigmavULs = manager.dict()
            sigmavULs_Jnuisance = manager.dict()
            run_combiner(config, sigmavULs, sigmavULs_Jnuisance)
            print("Combined {} limits (observational data)".format(collaboration))
            # Convert multiprocessing.managers.DictProxy to python 'dict'
            sigmavULs = dict(sigmavULs)
            sigmavULs_Jnuisance= dict(sigmavULs_Jnuisance)
            for channel in config['Configuration']['channels']:
                svUL = {'masses': sigmavULs['{}_masses'.format(channel)]}
                svUL['data'] = sigmavULs['{}_data'.format(channel)]
                svUL = pd.DataFrame(data=svUL)
                # Write the panda DataFrames into the hdf5 file
                svUL.to_hdf(hdf5file, key='{}/{}/sigmavULs'.format(channel,collaboration), mode='a')
                if config['Data']['J_nuisance']:
                    svUL_Jnuisance = {'masses': sigmavULs_Jnuisance['{}_masses'.format(channel)]}
                    svUL_Jnuisance['data'] = sigmavULs_Jnuisance['{}_data'.format(channel)]
                    svUL_Jnuisance = pd.DataFrame(data=svUL_Jnuisance)
                    # Write the panda DataFrames into the hdf5 file
                    svUL_Jnuisance.to_hdf(hdf5file, key='{}/{}/sigmavULs_Jnuisance'.format(channel,collaboration), mode='a')
                    
        config['Configuration']['collaborations'] = collaborations
        plot_sigmavULs_collaborations(hdf5file, output_dir, config)
