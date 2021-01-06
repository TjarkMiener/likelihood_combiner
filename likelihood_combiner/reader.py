import numpy as np
from scipy.interpolate import interp1d
import os
import pandas as pd

__all__ = [
    'LklComReader'
]

class LklComReader:
    def __init__(self,
                 channel,
                 sources,
                 collaborations,
                 data_directory,
                 j_factors):

        self.channel = channel
        self.sources = sources
        self.collaborations = collaborations
        self.data_directory = data_directory
        self.files_in_comb = self._construct_files_in_combination()
        self.j_factors = j_factors
        self.angular_separations, self.logJ_profile, self.DlogJ_profile = self._construct_jprofile()
        
        
    def read_tstables(self, logJ, simulation=-1):
        tstables = {}
        for file in self.files_in_comb:
            file_key = file.replace('.txt','')
            file_info = file_key.split("_")
            
            # Opening the txt files.
            if simulation == -1:
                ts_file = open("{}/{}".format(self.data_directory+self.channel,file),"r")
            else:
                ts_file = open("{}/simulations/{}_{}.txt".format(self.data_directory+self.channel,file_key,simulation),"r")
                
            # Going through the table in the txt file and storing the entries in a 2D array.
            values = np.array([[i for i in line.split()] for line in ts_file], dtype=np.float32).T

            # The first element of the first row correpond to the logJ-Factor. 
            logJ_file = np.float32(values[0][0])

            # The first entry of each row correponds to the mass.
            # Detect the first entry and store it in a separate array.
            mass = np.array([values[i][0] for i in np.arange(values.shape[0])], dtype=np.float32)
                
            # Delete the first entry, so only TS values are stored in ts_val.
            # Re-scale the sigmav values, if the logJ is not matching with the total logJ-Factor.
            ts_val = []
            for val in values:
                ts_val.append(val[1:])
            ts_val = np.array(ts_val, dtype=np.float32)
            if logJ[file_info[1]][file_info[2]] != logJ_file:
                ts_val[0] *= np.power(10.0,logJ_file)/np.power(10.0,logJ[file_info[1]][file_info[2]])

            # Store the arrays in the dictionary.
            tstables[file_key+'_ts'] = ts_val
            tstables[file_key+'_masses'] = mass
        return tstables
    
    def read_jfactor(self):
        logJ, DlogJ = {}, {}
        for source in self.sources:
            logJ[source] = {}
            DlogJ[source] = {}
            for collaboration in self.collaborations:
                # Check which collaboration observed the particular dSph
                if "{}_{}_{}.txt".format(self.channel,source,collaboration) in self.files_in_comb:
                    angular_separation = self.collaborations[collaboration]
                    logJ_profile_interpolation = interp1d(self.angular_separations[source], self.logJ_profile[source], kind='linear', fill_value='extrapolate')
                    DlogJ_profile_interpolation = interp1d(self.angular_separations[source], self.DlogJ_profile[source], kind='linear', fill_value='extrapolate')
                    logJ[source][collaboration] =  np.around(logJ_profile_interpolation(angular_separation), 2)
                    DlogJ[source][collaboration] =  np.around(DlogJ_profile_interpolation(angular_separation),2)
        return logJ, DlogJ

    def _construct_files_in_combination(self):
        files_in_comb = []
        files = np.array([x for x in os.listdir(self.data_directory+self.channel) if x.endswith(".txt")])
        for file in files:
            # Parsing the file names and checking validation.
            file_info = file.replace('.txt','').split("_")
            if file_info[0] == self.channel and file_info[1] in self.sources and file_info[2] in self.collaborations.keys():
                files_in_comb.append(file)
        return files_in_comb
        
    def _construct_jprofile(self):
        # Create dictionary to map the selected build-in
        # JFactor table to the corresponding read function.
        construct_jprofile_functions = {
            'GeringerSameth': self._construct_GeringerSameth_JFactor,
            'Bonnivard': self._construct_Bonnivard_JFactor
        }
        return construct_jprofile_functions[self.j_factors](self.data_directory)
        
    def _construct_GeringerSameth_JFactor(self, data_dir):
        JFactor_file = data_dir + "/jfactor/GeringerSameth/table.txt"
        angular_separations, logJ_profiles, DlogJ_profiles = {}, {}, {}
        for source in self.sources:
            # Opening txt file.
            file = open(JFactor_file, "r")
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
                    
    def _construct_Bonnivard_JFactor(self, data_dir):
        # Credits: V.Poireu
        JFactor_dir = data_dir + "/jfactor/Bonnivard/"
        angular_separations, logJ_profiles, DlogJ_profiles = {}, {}, {}
        for source in self.sources:
            file = source + "_Jalphaint_cls.output"
            angular_separation, J_array, J_low_array = np.genfromtxt(JFactor_dir + file, skip_header=4, usecols = (0,1,3), unpack=True)
            # The conversion factor is 1 GeV2.cm−5 = 2.25 × 10−7 Msun2.kpc−5
            conversion_factor = 2.25e-7
            logJ = np.log10(J_array / conversion_factor)
            DlogJ = logJ - np.log10(J_low_array / conversion_factor)
            angular_separations[source] = angular_separation
            logJ_profiles[source] = logJ
            DlogJ_profiles[source] = DlogJ
        return angular_separations, logJ_profiles, DlogJ_profiles
