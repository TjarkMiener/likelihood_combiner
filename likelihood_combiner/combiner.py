"""
combiner.py
===========
Function to combine data read by a LklCom reader class
"""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import os

import likelihood_combiner as lklcom
from likelihood_combiner.utils import *

__all__ = [
    'combiner'
]


def combiner(sigmav_range,
            LklCom_reader_class,
            output,
            sigmavULs=None,
            sigmavULs_Jnuisance=None,
            simulation_counter=None,
            total_number_of_simulation=None,
            simulations=[0]):
    """
    Combine data read by one of the LklCom reader class for a given sigmav range.

    Parameters
    ----------
    sigmav_range: `numpy.ndarray of type numpy.float32`
        sigmav range (ascending).
    LklCom_reader_class: `likelihood_combiner.reader.LklCom`
        class of the lklcom to handle the reading of the likelihood or ts tables.        
    output: path
        path to the lklcom results directory or hdf5 file.
    sigmavULs: `multiprocessing.managers.DictProxy`
        shared multiprocessing dict of the sigmav upper limits.
    sigmavULs_Jnuisance: `multiprocessing.managers.DictProxy`
        shared multiprocessing dict of the sigmav upper limits including the J-Factors as a nuisance parameter.
    simulation_counter: `multiprocessing.Manager.Value`
        current status of the simulation progress (needed for the progress bar with multiprocessing).
    total_number_of_simulation: int
        total number of the simulation considered (needed for the progress bar with multiprocessing).
    simulations: `list of type int`
        list of the numbers of the simulations.
    """
            
    channel = LklCom_reader_class.get_channel()
    LklCom_jfactor_class = LklCom_reader_class.get_LklCom_jfactor_class()
    jnuisance = LklCom_jfactor_class.get_jnuisance()
    logJ = LklCom_jfactor_class.get_DlogJ_comb() if jnuisance else LklCom_jfactor_class.get_logJ()
    
    for simulation in simulations:
    
        tstables = LklCom_reader_class(simulation)
        
        combined_ts, combined_ts_Jnuisance = {}, {}
        total_masses = []
        for source in logJ:
            combined_source_ts, combined_source_ts_Jnuisance = {}, {}
            # Combine ts values for the particular dSph
            for collaboration in logJ[source].keys():
                table = tstables[source][collaboration]
                
                sigmav_range_from_file = table["sigmav_range"]
                sigmav_range_from_file = round_sigmav_range(sigmav_range_from_file,LklCom_reader_class.get_sigmav_precision())

                masses = table["masses"]
                ts_values_per_mass = table["ts_values"]
                for mass, ts_values in zip(masses, ts_values_per_mass):
                    source_mass = "{}_{}".format(source, mass)

                    if mass not in total_masses:
                        total_masses.append(mass)
                    
                    if source_mass not in combined_source_ts:
                        combined_source_ts[source_mass] = np.zeros(len(sigmav_range))
                        if jnuisance:
                            combined_source_ts_Jnuisance[source_mass] = np.zeros(len(sigmav_range))

                    if not np.all(ts_values == ts_values[0]):
                        if not np.array_equal(sigmav_range, sigmav_range_from_file):
                            lin_interpolation = interp1d(sigmav_range_from_file, ts_values, kind='linear', fill_value='extrapolate')
                            ts_values = lin_interpolation(sigmav_range)
                        combined_source_ts[source_mass] += ts_values
                        if jnuisance:
                            combined_source_ts_Jnuisance[source_mass] += ts_values
                            if logJ[source][collaboration] != 0.0:
                                combined_source_ts_Jnuisance[source_mass] = LklCom_jfactor_class.compute_Jnuisance(sigmav_range, combined_source_ts_Jnuisance[source_mass], np.float32(logJ[source][collaboration]))
                    else:
                        if jnuisance and logJ[source][collaboration] != 0.0 and collaboration == list(logJ[source].keys())[-1]:
                            combined_source_ts_Jnuisance[source_mass] = LklCom_jfactor_class.compute_Jnuisance(sigmav_range, combined_source_ts_Jnuisance[source_mass], np.float32(logJ[source][collaboration]))


            # Combine ts values for all dSphs
            for mass in total_masses:
                source_mass = "{}_{}".format(source, mass)
                if str(mass) not in combined_ts:
                    combined_ts[str(mass)] = np.zeros(len(sigmav_range))
                    if jnuisance:
                        combined_ts_Jnuisance[str(mass)] = np.zeros(len(sigmav_range))
                if source_mass in combined_source_ts:
                    if not np.all(combined_source_ts[source_mass] == combined_source_ts[source_mass][0]):
                        combined_ts[str(mass)] += combined_source_ts[source_mass]
                        if jnuisance:
                            combined_ts_Jnuisance[str(mass)] += combined_source_ts_Jnuisance[source_mass]

        # Compute the sigmav ULs per mass
        combined_sources_limits, combined_sources_sensitivity = lklcom.sensitivity.compute_sensitivity(sigmav_range, combined_ts)
        if jnuisance:
            combined_sources_limits_Jnuisance, combined_sources_sensitivity_Jnuisance = lklcom.sensitivity.compute_sensitivity(sigmav_range, combined_ts_Jnuisance)

        # Output handling
        total_masses = sorted(total_masses)
        svUL = []
        svUL_Jnuisance = []
        for mass in total_masses:
            svUL.append(combined_sources_limits[str(mass)])
            if jnuisance:
                svUL_Jnuisance.append(combined_sources_limits_Jnuisance[str(mass)])
            else:
                svUL_Jnuisance.append(np.nan)

        if sigmavULs is None and sigmavULs_Jnuisance is None:
            if simulation > 0:
                output += "/{}_simu{}.h5".format(channel, simulation)
            else:
                output += "/{}_data.h5".format(channel)
            # Dumping the upper limits in the h5 file
            pd.DataFrame(data=total_masses).to_hdf(output, key='/masses', mode='w')
            pd.DataFrame(data=svUL).to_hdf(output, key='/sigmavULs', mode='a')
            if jnuisance:
                pd.DataFrame(data=svUL_Jnuisance).to_hdf(output, key='/sigmavULs_Jnuisance', mode='a')

        else:
            if simulation > 0:
                h5_key = "{}_simu{}".format(channel, simulation)
            else:
                h5_key = "{}_data".format(channel)
                sigmavULs["{}_masses".format(channel)] = total_masses
                sigmavULs_Jnuisance["{}_masses".format(channel)] = total_masses

            sigmavULs[h5_key] = svUL
            sigmavULs_Jnuisance[h5_key] = svUL_Jnuisance

            simulation_counter.value += 1
            progress_bar(simulation_counter.value,np.int(total_number_of_simulation))

    return

