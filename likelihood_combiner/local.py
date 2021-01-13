"""
local.py
========
Function to automatically set up multiprocessing (only Linux supported)
"""

import argparse
import numpy as np
from multiprocessing import Process, Manager, freeze_support
import tables
import pandas as pd
import os
import yaml

import likelihood_combiner as lklcom
from likelihood_combiner.combiner import combiner
from likelihood_combiner.utils import *

__all__ = [
    'run_local_on_linux'
]    

def run_local_on_linux(settings,
                       input,
                       output):
    """
    Run LklCom on a local linux machine. Support multiprocessing.
    This function only works for linux users, because MacOS or Windows don't allow you to set up multiprocessing this way.
    See: https://www.pythonforthelab.com/blog/differences-between-multiprocessing-windows-and-linux/  

    Parameters
    ----------
    settings: dict
        settings of the combination, following the minimum skeleton (see example config file for more documentation):
        {'Hardware' : {'cpu_counts': `int`},
        'Data' : {'buildin_j_factors': `string`, 'j_nuisance': `boolean`, 'simulations': `int`},
        'Configuration' : {'channels': `numpy.ndarray of type string`, 'sources': `numpy.ndarray of type string`, 'collaborations': `dictionary`}}
    input: str
        path to the input file or directory
    output: str
        path to the output hdf5 file
    """

    # Initializing of the LklCom jfactor class
    if 'buildin_j_factors' not in settings['Data']:
        settings['Data']['buildin_j_factors'] = "Custom"
    
    if settings['Data']['buildin_j_factors'] == "GeringerSameth":
        LklCom_jfactor_class = lklcom.jfactor.GeringerSameth(sources=settings['Configuration']['sources'],
                                                    collaborations=settings['Configuration']['collaborations'],
                                                    combination_data=input,
                                                    jnuisance=settings['Data']['j_nuisance'])
    elif settings['Data']['buildin_j_factors'] == "Bonnivard":
        LklCom_jfactor_class = lklcom.jfactor.Bonnivard(sources=settings['Configuration']['sources'],
                                                    collaborations=settings['Configuration']['collaborations'],
                                                    combination_data=input,
                                                    jnuisance=settings['Data']['j_nuisance'])
    else:
        LklCom_jfactor_class = lklcom.jfactor.Custom(logJ=settings['Data']['custom_logJ'],
                                                    DlogJ=settings['Data']['custom_DlogJ'],
                                                    jnuisance=settings['Data']['j_nuisance'])

    # Constructing in the the sigmav range and spacing
    sigmav_range = get_sigmav_range()
    # Overwriting from the settings, if provided.
    if "sigmav_min" in settings['Data'] and "sigmav_max" in settings['Data'] and "sigmav_nPoints" in settings['Data']:
        sigmav_range = get_sigmav_range(settings['Data']['sigmav_min'], settings['Data']['sigmav_max'], settings['Data']['sigmav_nPoints'], 3)
    
    for channel in settings['Configuration']['channels']:
        print("\nChannel '{}'".format(channel))
        
        # Initializing of the LklCom reader class
        if input.endswith(".h5") or input.endswith(".hdf5"):
            LklCom_reader_class = lklcom.reader.LklCom_hdf5(channel=channel,
                                                            LklCom_jfactor_class=LklCom_jfactor_class)
        if os.path.isdir(input):
            LklCom_reader_class = lklcom.reader.LklCom_txtdir(channel=channel,
                                                            LklCom_jfactor_class=LklCom_jfactor_class)
        
        # Set up the hardware settings for the parallel processing
        try:
            cpu_counts = os.cpu_count() if settings['Hardware']['cpu_counts'] == 'all' or settings['Hardware']['cpu_counts'] > os.cpu_count() else settings['Hardware']['cpu_counts']
            if cpu_counts is None:
                raise KeyError
        except KeyError:
            cpu_counts = 1
        if 'simulations' in settings['Data']:
            simulations = np.int(settings['Data']['simulations']+1)
        else:
            simulations = 1
        if cpu_counts > simulations:
            cpu_counts = simulations

        # Set up all processes
        if simulations <= 1:
            combiner(sigmav_range=sigmav_range,
                    LklCom_reader_class=LklCom_reader_class,
                    output=output)
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
                                        output,
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
            svUL.to_hdf(output, key='{}/sigmavULs'.format(channel), mode='a')
            if settings['Data']['j_nuisance']:
                svUL_Jnuisance = {'masses': sigmavULs_Jnuisance['{}_masses'.format(channel)]}
                for key, value in dict(sigmavULs_Jnuisance).items():
                    if channel in key: svUL_Jnuisance[key.replace('{}_'.format(channel),'')] = value
                svUL_Jnuisance = pd.DataFrame(data=svUL_Jnuisance)
                # Write the panda DataFrames into the hdf5 file
                svUL_Jnuisance.to_hdf(output, key='{}/sigmavULs_Jnuisance'.format(channel), mode='a')
    return
