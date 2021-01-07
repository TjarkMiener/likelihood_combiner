import numpy as np
from scipy.interpolate import interp1d
import tables
import os

__all__ = [
    'LklCom',
]

class LklCom:

    def __init__(self,
                 channel,
                 lklcom_JFactor,
                 lklcom_file=None,
                 data_directory=None
                 ):

        self.channel = channel
        self.lklcom_logJ = lklcom_JFactor.logJ
        self.combination_info = self._construct_combination_info()

    def read_tstables_from_lklcom_hdf5(self,
                                       lklcom_file,
                                       simulation=-1):
        """
        Read the likelihood or ts tables from the lklcom hdf5 file.
        Parameters
        ----------
        lklcom_file: `string`
            path to the lklcom hdf5 file.
        simulation: `int`
            number of the simulation.
            Default: -1 (= data sample)
        Returns
        -------
        ts_tables: `dict`
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
        h5 = tables.open_file(lklcom_file, "r")
        
        # Initializing of the nested dictionary.
        ts_tables = {}
        for source in self.lklcom_logJ:
            ts_tables[source] = {}
        for (source, collaboration) in self.combination_info:
            # Initializing of the returning dictionary.
            ts_tables[source][collaboration] = {}
            
            # Getting the table from the lklcom hdf5 file.
            table_name = "data" if simulation == -1 else "simu_{}".format(simulation)
            table = eval("h5.root.{}.{}.{}.{}".format(collaboration, source, self.channel, table_name))

            # Reading the data from the table.
            masses = np.array(table.cols._f_col('masses'))
            ts_values = np.array(table.cols._f_col('ts_values'))
            # The first element of the first row correpond to the logJ-Factor and to the sigmav range.
            logJ_from_file = np.float32(masses[0])
            # Delete the logJ-factor of the
            masses = np.array(masses[1:], dtype=np.float32)
            # The first element of the first row correpond to the logJ-Factor.
            sigmav_range = np.float32(ts_values[0])
            # Delete the logJ-factor of the
            ts_values = np.array(ts_values[1:], dtype=np.float32)

            # Delete the first entry, so only TS values are stored in ts_val.
            # Re-scale the sigmav values, if the logJ is not matching with the total logJ-Factor.
            if self.lklcom_logJ[source][collaboration] != logJ_from_file:
                sigmav_range *= np.power(10.0,logJ_from_file)/np.power(10.0,self.lklcom_logJ[source][collaboration])

            # Store the arrays in the dictionary.
            ts_tables[source][collaboration]["masses"] = masses
            ts_tables[source][collaboration]["sigmav_range"] = sigmav_range
            ts_tables[source][collaboration]["ts_values"] = ts_values
        return ts_tables


    def read_tstables_from_gLike_txt(self,
                                     data_directory,
                                     simulation=-1):
        """
        Read the likelihood or ts tables from the gLike txt files.
        Parameters
        ----------
        lklcom_file: `string`
            path to the lklcom hdf5 file.
        simulation: `int`
            number of the simulation.
            Default: -1 (= data sample)
        Returns
        -------
        ts_tables: `dict`
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
        ts_tables = {}
        for source in self.lklcom_logJ:
            ts_tables[source] = {}
        for (source, collaboration) in self.combination_info:
            # Initializing of the returning dictionary.
            ts_tables[source][collaboration] = {}
            
            # Opening the txt files.
            if simulation == -1:
                ts_file = open("{}/{}_{}_{}.txt".format(data_directory, self.channel, source, collaboration), "r")
            else:
                ts_file = open("{}/{}_{}_{}_{}.txt".format(data_directory, self.channel, source, collaboration, simulation), "r")
                
            # Going through the table in the txt file and storing the entries in a 2D array.
            values = np.array([[i for i in line.split()] for line in ts_file], dtype=np.float32).T

            # The first element of the first row correpond to the logJ-Factor. 
            logJ_from_file = np.float32(values[0][0])

            # The first entry of each row correponds to the mass.
            # Detect the first entry and store it in a separate array.
            masses = np.array([values[i][0] for i in np.arange(values.shape[0])], dtype=np.float32)
                
            # Delete the first entry, so only TS values are stored in ts_val.
            # Re-scale the sigmav values, if the logJ is not matching with the total logJ-Factor.
            ts_val = []
            for val in values:
                ts_val.append(val[1:])
            ts_val = np.array(ts_val[1:], dtype=np.float32)
            if self.lklcom_logJ[source][collaboration] != logJ_from_file:
                sigmav_range *= np.power(10.0,logJ_from_file)/np.power(10.0,self.lklcom_logJ[source][collaboration])

            # Store the arrays in the dictionary.
            ts_tables[source][collaboration]["masses"] = masses
            ts_tables[source][collaboration]["sigmav_range"] = sigmav_range
            ts_tables[source][collaboration]["ts_values"] = ts_values
        return tstables

    def _construct_combination_info(self):
        combination_info = []
        for source in self.lklcom_logJ:
            for collaboration in self.lklcom_logJ[source]:
                if self.lklcom_logJ[source][collaboration]:
                    combination_info.append((source, collaboration))
        return combination_info

'''
class JFactor:

    def __init__(self,
                 sources=None,
                 collaborations=None,
                 build_in_JFactors=None,
                 resource=None,
                 lklcom_file=None,
                 precision=2,
                 logJ=None,
                 DlogJ=None
                 ):

        
        self.sources = sources
        self.collaborations = collaborations
        self.lklcom_file = lklcom_file

        # The arguments logJ and DlogJ should be only used to hardcode the log J-factor
        # and it's uncertainty.
        if logJ is None:
            logJ = {}
        else:
            self.sources = logJ.keys()
        self.logJ = logJ
        if DlogJ is None:
            DlogJ = {}
        self.DlogJ = DlogJ

        # Build-in JFactor table
        if build_in_JFactors not in ["GeringerSameth", "Bonnivard"]:
            build_in_JFactors = "Custom"
        self.build_in_JFactors = build_in_JFactors
    
        if self.build_in_JFactors != "Custom":
            self.angular_separations, self.logJ_profile, self.DlogJ_profile = self._construct_jprofile(build_in_JFactors, resource)
            self.combination_info = self._construct_combination_info_from_hdf5()
            self.logJ, self.DlogJ = self._get_jfactors(precision=precision)
        
        self.DlogJ_comb = self._compute_DlogJ_for_combination()


    def _compute_DlogJ_for_combination(self):
        # Compute the actual uncertainties, which is introduced in the combination.
        # Therefore, we have to sort the DlogJ[source] in descending order and compute the
        # actual uncertainty using DlogJ_diff = (DlogJ^2 - DlogJ_next^2)^1/2. For the last
        # element DlogJ_diff is set to (the smallest) DlogJ.
        # After this computation, DlogJ_comb[source] holds the J Factor uncertainty
        # for the combined analysis.
        DlogJ_comb = {}
        for source in self.sources:
            DlogJ_comb[source] = {}
            DlogJ[source] = dict(sorted(DlogJ[source].items(), key=lambda x: x[1], reverse=True))
            prev_collaboration = None
            for collaboration in DlogJ[source]:
                if prev_collaboration:
                    DlogJ_diff = 0.0
                    if DlogJ_comb[source][prev_collaboration] != DlogJ[source][collaboration]:
                        DlogJ_diff = np.sqrt(np.power(DlogJ_comb[source][prev_collaboration],2) - np.power(DlogJ[source][collaboration],2))
                    DlogJ_comb[source][prev_collaboration] = DlogJ_diff
                prev_collaboration = collaboration
                prev_DlogJ = DlogJ_comb[source][prev_collaboration] = DlogJ[source][collaboration]
        return DlogJ_comb
    
    
    def _get_jfactors(self, precision=2):
        logJ, DlogJ = {}, {}
        for source in self.sources:
            logJ[source] = {}
            DlogJ[source] = {}
            for collaboration in self.collaborations:
                # Check which collaboration observed the particular dSph
                if (source, collaboration) in self.combination_info:
                    angular_separation = self.collaborations[collaboration]
                    logJ_profile_interpolation = interp1d(self.angular_separations[source], self.logJ_profile[source], kind='linear', fill_value='extrapolate')
                    DlogJ_profile_interpolation = interp1d(self.angular_separations[source], self.DlogJ_profile[source], kind='linear', fill_value='extrapolate')
                    logJ[source][collaboration] =  np.around(logJ_profile_interpolation(angular_separation), precision)
                    DlogJ[source][collaboration] =  np.around(DlogJ_profile_interpolation(angular_separation), precision)
        return logJ, DlogJ

    def _construct_jprofile(self, build_in_JFactors, resource):
        # Create dictionary to map the selected build-in
        # JFactor table to the corresponding read function.
        construct_jprofile_functions = {
            'GeringerSameth': self._construct_GeringerSameth_JFactor,
            'Bonnivard': self._construct_Bonnivard_JFactor
        }
        return construct_jprofile_functions[build_in_JFactors](resource)
        
    def _construct_GeringerSameth_JFactor(self, resource):
        angular_separations, logJ_profiles, DlogJ_profiles = {}, {}, {}
        for source in self.sources:
            # Opening txt file.
            file = open(resource, "r")
            angular_separation = []
            logJ = []
            DlogJ = []
            for i,line in enumerate(file):
                if i > 41:
                    # The source in the GS file for this line
                    source_GS = line[:18].replace(" ", "")
                    line = line[18:].split()
                    if source == source_GS:
                        angular_separation.append(np.float32(line[0]))
                        logJ.append(np.float32(line[3]))
                        DlogJ.append(np.float32(line[3]) - np.float32(line[2]))
            file.close()
            angular_separations[source] = angular_separation
            logJ_profiles[source] = logJ
            DlogJ_profiles[source] = DlogJ
        return angular_separations, logJ_profiles, DlogJ_profiles
                    
    def _construct_Bonnivard_JFactor(self, resource):
        # Credits: V.Poireu
        angular_separations, logJ_profiles, DlogJ_profiles = {}, {}, {}
        for source in self.sources:
            file = source + "_Jalphaint_cls.output"
            angular_separation, J_array, J_low_array = np.genfromtxt(resource + file, skip_header=4, usecols = (0,1,3), unpack=True)
            # The conversion factor is 1 GeV2.cm−5 = 2.25 × 10−7 Msun2.kpc−5
            conversion_factor = 2.25e-7
            logJ = np.log10(J_array / conversion_factor)
            DlogJ = logJ - np.log10(J_low_array / conversion_factor)
            angular_separations[source] = angular_separation
            logJ_profiles[source] = logJ
            DlogJ_profiles[source] = DlogJ
        return angular_separations, logJ_profiles, DlogJ_profiles

    def _construct_combination_info_from_hdf5(self):
        combination_info = []
        if self.lklcom_file is None:
            # Writing all possible combinations.
            for source in self.sources:
                for collaboration in self.collaborations:
                    combination_info.append((source, collaboration))
        else:
            # Writing only valid combinations.
            # Opening the hdf5 file.
            h5 = tables.open_file(self.lklcom_file, "r")
            # Looping over the likelihood or ts tables and append valid combinations.
            for h5_groups in h5.walk_groups("/"):
                for table in h5.list_nodes(h5_groups, classname='Table'):
                    table_info = table._v_pathname.split("/")
                    if (table_info[2], table_info[1]) not in combination_info:
                        combination_info.append((table_info[2], table_info[1]))
        return combination_info
        
'''
