import argparse
import numpy as np
from multiprocessing import Process,Manager
import pandas as pd
from scipy.interpolate import interp1d
import os
import yaml

from likelihood_combiner.reader import LklComReader,JFactor_Reader
from likelihood_combiner.writer import LklComWriter
from likelihood_combiner.utils import compute_sensitivity,compute_Jnuisance,plot_sigmavULs

def run_combiner(config, hdf5file, sigmavUL, sigmavUL_Jnuisance, simulations=[-1]):

    try:
        data_dir = config['Data']['data_directory']
        if data_dir is None:
            raise KeyError
    except KeyError:
        data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/"))
    
    try:
        simu_dir = config['Data']['simu_directory']
        if simu_dir is None:
            raise KeyError
    except KeyError:
        simu_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../simu/"))

    channels = config['Configuration']['channels']
    sources = config['Configuration']['sources']
    collaborations = config['Configuration']['collaborations']
    jnuisance = config['Data']['J_nuisance']

    try:
        JFactor_file = config['Data']['JFactor_table']
        if JFactor_file is None:
            raise KeyError
    except KeyError:
        JFactor_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/Jfactor_Geringer-SamethTable.txt"))

    if jnuisance:
        Jfactor_reader = JFactor_Reader()
        sources_logJ,sources_DlogJ = Jfactor_reader.read_JFactor(JFactor_file,sources)
        del Jfactor_reader

    # Reading in the the sigmav range and spacing. Round of the third digit to avoid an interpolation.
    sigmavMin = -(np.abs(np.floor(np.log10(np.abs(float(config['Data']['sigmaV_min'])))).astype(int))).astype(int)
    sigmavMax = -(np.abs(np.floor(np.log10(np.abs(float(config['Data']['sigmaV_max'])))).astype(int))).astype(int)
    sigmavNPoints = int(config['Data']['sigmaV_nPoints'])
    sigmav = np.logspace(sigmavMin, sigmavMax, sigmavNPoints, endpoint=True, dtype=np.float32)
    exponent = (np.abs(np.floor(np.log10(np.abs(sigmav))).astype(int))+3).astype(int)
    for i,e in enumerate(exponent):
        sigmav[i] = np.around(sigmav[i],decimals=e)
        
    reader = LklComReader()
    if simulations == [-1]:
        writer = LklComWriter()
        writer.convert_txts2hdf5(hdf5file, data_dir, channels, sources, collaborations)
        del writer
        tstables, massvals = reader.read_tstables(hdf5file, channels, sources, collaborations)

    for simulation in simulations:
        if simulation >= 0:
            tstables, massvals = reader.read_simutstables(simu_dir, simulation, channels, sources, collaborations)
        for channel in channels:
            combined_channel_ts = {}
            combined_channel_ts_Jnuisance = {}
            mass_axis = []
            for source in sources:
                # Combine ts values for each dSphs
                combined_source_ts = {}
                combined_source_ts_Jnuisance = {}
                for collaboration in collaborations:
                    key = "{}_{}_{}".format(channel,source,collaboration)
                    if key in massvals:
                        tstable = tstables[key]
                        masses = massvals[key]
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
                    if str(m) in combined_channel_ts:
                        if source+"_"+str(m) in combined_source_ts:
                            combined_channel_ts[str(m)] += combined_source_ts[source+"_"+str(m)]
                            if jnuisance:
                                combined_channel_ts_Jnuisance[str(m)] += combined_source_ts_Jnuisance[source+"_"+str(m)]
                    else:
                        if source+"_"+str(m) in combined_source_ts:
                            combined_channel_ts[str(m)] = combined_source_ts[source+"_"+str(m)]
                            if jnuisance:
                                combined_channel_ts_Jnuisance[str(m)] = combined_source_ts_Jnuisance[source+"_"+str(m)]
            
            combined_sources_limits,combined_sources_sensitivity = compute_sensitivity(sigmav, combined_channel_ts)
            if jnuisance:
                combined_sources_limits_Jnuisance,combined_sources_sensitivity_Jnuisance = compute_sensitivity(sigmav, combined_channel_ts_Jnuisance)
                
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
                sigmavUL['masses'] = sorted_mass
                sigmavUL_Jnuisance['masses'] = sorted_mass

            sigmavUL[h5_key] = svUL
            sigmavUL_Jnuisance[h5_key] = svUL_Jnuisance
            
        if simulation >= 0:
            print("Combined MC limits (simulation {})".format(simulation))
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
        config = yaml.load(config_file)

    try:
        hdf5file = config['Data']['hdf5_dataset']
        if hdf5file is None:
            raise KeyError
    except KeyError:
        hdf5file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/lklcom.h5"))
      
    
    # Create a multiprocessing.Manager dict to share memory between the parallel processes
    manager = Manager()
    sigmavUL = manager.dict()
    sigmavUL_Jnuisance = manager.dict()
    # Combining limits using real observational data sets
    run_combiner(config, hdf5file, sigmavUL, sigmavUL_Jnuisance)
    print("Combined limits (observational data)")
    if config['Data']['cl_bands']:
        # Calculate the convendance level bands
        # Set up the hardware settings for the parallel processing
        try:
            cpu_counts = os.cpu_count() if config['Hardware']['cpu_counts'] == 'all' else config['Hardware']['cpu_counts']
            if cpu_counts is None:
                raise KeyError
        except KeyError:
            cpu_counts = 1
        simulations = config['Data']['simulations']
        if cpu_counts > simulations:
            cpu_counts = simulations
                
        # Set up all processes
        parallel_simulations = np.array_split(np.arange(simulations), cpu_counts)
        jobs = []
        for parallel_simulation in parallel_simulations:
            process = Process(target=run_combiner, args=(config, hdf5file, sigmavUL, sigmavUL_Jnuisance, parallel_simulation))
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
                    
    # Convert multiprocessing.managers.<>Dict to python 'dict'
    sigmavUL = dict(sigmavUL)
    sigmavUL_Jnuisance= dict(sigmavUL_Jnuisance)
    for channel in config['Configuration']['channels']:
        svUL = {'masses': sigmavUL['masses']}
        for key, value in dict(sigmavUL).items():
            if channel in key: svUL[key.replace('{}_'.format(channel),'')] = value
        svUL = pd.DataFrame(data=svUL)
        # Write the panda DataFrames into the hdf5 file
        svUL.to_hdf(hdf5file, key='{}/Combination/sigmavUL'.format(channel), mode='a')
        if config['Data']['J_nuisance']:
            svUL_Jnuisance = {'masses': sigmavUL['masses']}
            for key, value in dict(sigmavUL_Jnuisance).items():
                if channel in key: svUL_Jnuisance[key.replace('{}_'.format(channel),'')] = value
            svUL_Jnuisance = pd.DataFrame(data=svUL_Jnuisance)
            # Write the panda DataFrames into the hdf5 file
            svUL_Jnuisance.to_hdf(hdf5file, key='{}/Combination/sigmavUL_Jnuisance'.format(channel), mode='a')
    
    try:
        output_dir = config['Output']['output_directory']
        if output_dir is None:
            raise KeyError
    except KeyError:
        output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output/"))
    # Create output directory if it doesn't exist already
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    plot_sigmavULs(hdf5file, output_dir, config)
