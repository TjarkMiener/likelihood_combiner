"""
reader.py
=========
Collection of classes to read likelihood or ts tables.
"""

import numpy as np
from scipy.interpolate import interp1d
import tables
import os

__all__ = [
    'LklCom',
    'LklCom_hdf5',
    'LklCom_txtdir',
    'LklCom_custom'
]

class LklCom:
    """
    Abstract class for reading likelihood or ts tables.
    """

    def __init__(self,
                 LklCom_jfactor_class,
                 channel=None,
                 combination_data=None,
                 simulations=None,
                 sigmav_precision=3):
        """
        Parameters
        ----------
        LklCom_jfactor_class : `likelihood_combiner.jfactor.JFactor`
            class of the lklcom to handle the J-Factor.
        channel: str
            name of the channel.
        combination_data : path
            path to the dataset, which will be used in the combination.
        simulations: `list of int`
            list of simulations, which will be used in the combination (only for subclass LklCom_custom).
        sigmav_precision: int
            precision of the returning sigmav range.
        """


        self.LklCom_jfactor_class = LklCom_jfactor_class
        if channel is None:
            channel = LklCom_jfactor_class.get_channel()
        self.channel = channel
        self.logJ = LklCom_jfactor_class.logJ
        self.combination_info = self._construct_combination_info()
        if combination_data is None:
            combination_data = LklCom_jfactor_class.combination_data
        self.combination_data = combination_data
        self.simulations= simulations
        self.sigmav_precision = sigmav_precision

    def get_channel(self):
        return self.channel

    def get_LklCom_jfactor_class(self):
        return self.LklCom_jfactor_class

    def get_logJ(self):
        return self.logJ

    def get_combination_info(self):
        return self.combination_info

    def get_combination_data(self):
        return self.combination_data
    
    def get_simulations(self):
        return self.simulations
    
    def get_sigmav_precision(self):
        return self.sigmav_precision
    
    def _construct_combination_info(self):
        combination_info = []
        for source in self.logJ:
            for collaboration in self.logJ[source]:
                if self.logJ[source][collaboration]:
                    combination_info.append((source, collaboration))
        return combination_info

class LklCom_hdf5(LklCom):
    """
    Subclass to read likelihood or ts tables from the lklcom hdf5 file.
    """

    def __init__(self,
                LklCom_jfactor_class,
                channel=None,
                combination_data=None,
                sigmav_precision=3):
        """
        Parameters
        ----------
        LklCom_jfactor_class : `likelihood_combiner.jfactor.JFactor`
            class of the lklcom to handle the J-Factor.
        channel: str
            name of the channel.
        combination_data : path
            path to the dataset, which will be used in the combination.
        sigmav_precision: int
            presicion of the returning sigmav range.
        """

        super().__init__(LklCom_jfactor_class=LklCom_jfactor_class, channel=channel, combination_data=combination_data, sigmav_precision=sigmav_precision)

    def __call__(self, simulation=0):
        """
        Read the likelihood or tstables from the lklcom hdf5 file.

        Parameters
        ----------
        simulation: `int`
            number of the simulation.
            Default: -1 (= data sample)

        Returns
        -------
        tstables: `dict`
            nested dictionaries (name of the source -> name of the collaboration):
                keys: `string`
                    "masses" with values `numpy.ndarray of type numpy.float32`:
                        DM masses.
                    "sigmav_range" with values `numpy.ndarray of type numpy.float32`:
                        sigmav range (ascending).
                    "ts_values" with values `numpy.ndarray of type numpy.float32`:
                        likelihood or ts values (ascending).
        """

        # Opening the hdf5 file.
        h5 = tables.open_file(self.combination_data, "r")
        
        # Initializing of the nested dictionary.
        tstables = {}
        for source in self.logJ:
            tstables[source] = {}
        for (source, collaboration) in self.combination_info:
            # Initializing of the returning dictionary.
            tstables[source][collaboration] = {}
            
            # Getting the table from the lklcom hdf5 file.
            table_name = "data" if simulation == 0 else "simu_{}".format(simulation)
            table = eval("h5.root.{}.{}.{}.{}".format(collaboration, source, self.channel, table_name))

            # Reading the data from the table.
            masses = np.array(table.cols._f_col('masses'))
            ts_values = np.array(table.cols._f_col('ts_values'))
            # The first element of the first row correpond to the logJ-Factor and to the sigmav range.
            logJ_from_table = np.float32(masses[0])
            # Delete the logJ-factor of the
            masses = np.array(masses[1:], dtype=np.float32)
            # The first element of the first row correpond to the logJ-Factor.
            sigmav_range = np.float32(ts_values[0])
            # Delete the logJ-factor of the
            ts_values = np.array(ts_values[1:], dtype=np.float32)

            # Delete the first entry, so only TS values are stored in ts_val.
            # Re-scale the sigmav values, if the logJ is not matching with the total logJ-Factor.
            if self.logJ[source][collaboration] != logJ_from_table:
                sigmav_range *= np.power(10.0,logJ_from_table)/np.power(10.0,self.logJ[source][collaboration])

            # Store the arrays in the dictionary.
            tstables[source][collaboration]["masses"] = masses
            tstables[source][collaboration]["sigmav_range"] = sigmav_range
            tstables[source][collaboration]["ts_values"] = ts_values

        # Closing hdf5 file.
        h5.close()
        
        return tstables

class LklCom_txtdir(LklCom):
    """
    Subclass to read likelihood or ts tables from the gLike txt files.
    """

    def __init__(self,
                LklCom_jfactor_class,
                channel=None,
                combination_data=None,
                sigmav_precision=3):
        """
        Parameters
        ----------
        LklCom_jfactor_class : `likelihood_combiner.jfactor.JFactor`
            class of the lklcom to handle the J-Factor.
        channel: str
            name of the channel.
        combination_data : path
            path to the dataset, which will be used in the combination.
        sigmav_precision: int
            presicion of the returning sigmav range.
        """

        super().__init__(LklCom_jfactor_class=LklCom_jfactor_class, channel=channel, combination_data=combination_data, sigmav_precision=sigmav_precision)
            
    def __call__(self, simulation=0):
        """
        Read the likelihood or tstables from the gLike txt files.

        Parameters
        ----------
        simulation: `int`
            number of the simulation.
            Default: -1 (= data sample)

        Returns
        -------
        tstables: `dict`
            nested dictionaries (name of the source -> name of the collaboration):
                keys: `string`
                    "masses" with values `numpy.ndarray of type numpy.float32`:
                        DM masses.
                    "sigmav_range" with values `numpy.ndarray of type numpy.float32`:
                        sigmav range (ascending).
                    "ts_values" with values `numpy.ndarray of type numpy.float32`:
                        likelihood or ts values (ascending).
        """
        
        # Initializing of the nested dictionary.
        tstables = {}
        for source in self.logJ:
            tstables[source] = {}
        for (source, collaboration) in self.combination_info:
            # Initializing of the returning dictionary.
            tstables[source][collaboration] = {}
            
            # Opening the txt files.
            if simulation == 0:
                ts_file = open("{}/{}_{}_{}.txt".format(self.combination_data, self.channel, source, collaboration), "r")
            else:
                ts_file = open("{}/{}_{}_{}_{}.txt".format(self.combination_data, self.channel, source, collaboration, simulation), "r")
                
            # Going through the table in the txt file and storing the entries in a 2D array.
            table = np.array([[i for i in line.split()] for line in ts_file], dtype=np.float32).T
            # The first element of the first row correpond to the logJ-Factor.
            logJ_from_file = np.float32(table[0][0])
            # Storing the inverted sigmav range.
            sigmav_range = np.array(table[0][1:], dtype=np.float32)[::-1]
            # The first entry of each row correponds to the mass (or log J-Factor).
            # Detect the first entry and store it in a separate array.
            masses = np.array([table[i][0] for i in np.arange(1,table.shape[0],1)], dtype=np.float32)
            # Delete the first entry, so only TS values are stored in ts_val.
            ts_values = []
            for row in table[1:]:
                ts_values.append(np.array(row[1:], dtype=np.float32)[::-1])
            ts_values = np.array(ts_values, dtype=np.float32)
            # Re-scale the sigmav values, if the logJ is not matching with the total logJ-Factor.
            if self.logJ[source][collaboration] != logJ_from_file:
                sigmav_range *= np.power(10.0,logJ_from_file)/np.power(10.0,self.logJ[source][collaboration])

            # Store the arrays in the dictionary.
            tstables[source][collaboration]["masses"] = masses
            tstables[source][collaboration]["sigmav_range"] = sigmav_range
            tstables[source][collaboration]["ts_values"] = ts_values

        return tstables


class LklCom_custom(LklCom):
    """
    Subclass to read likelihood or ts tables from the custom tstables dict. Only recommended to use for by passing lklcom or gLike data formats. 
    """

    def __init__(self,
                LklCom_jfactor_class,
                channel,
                combination_data,
                simulations,
                sigmav_precision=3):
        """
        Parameters
        ----------
        LklCom_jfactor_class : `likelihood_combiner.jfactor.JFactor`
            class of the lklcom to handle the J-Factor.
        channel: str
            name of the channel.
        combination_data : path
            path to the dataset, which will be used in the combination.
        simulations: `list of int`
            list of simulations, which will be used in the combination (only for subclass LklCom_custom).
        sigmav_precision: int
            presicion of the returning sigmav range.
        """

        super().__init__(LklCom_jfactor_class=LklCom_jfactor_class, channel=channel, combination_data=combination_data, simulations=simulations, sigmav_precision=sigmav_precision)
    
        self.tstables = {}
        for simulation in self.simulations:
            self.tstables[simulation] = combination_data[simulation]
        self.sigmav_precision = sigmav_precision

    def __call__(self, simulation=0):
        """
        Read the likelihood or tstables from the custom tstables dict.

        Parameters
        ----------
        simulation: `int`
            number of the simulation.
            Default: -1 (= data sample)

        Returns
        -------
        tstables: `dict`
            nested dictionaries (name of the source -> name of the collaboration):
                keys: `string`
                    "masses" with values `numpy.ndarray of type numpy.float32`:
                        DM masses.
                    "sigmav_range" with values `numpy.ndarray of type numpy.float32`:
                        sigmav range (ascending).
                    "ts_values" with values `numpy.ndarray of type numpy.float32`:
                        likelihood or ts values (ascending).
        """
        return self.tstables[simulation]
