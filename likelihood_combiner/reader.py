import numpy as np
import os
import pandas as pd

class LklComReader:
    def __init__(self,
                 channel,
                 sources,
                 collaborations):

        self.channel = channel
        self.sources = sources
        self.collaborations = collaborations
    
    def read_tstables(self, path2txts, logJ, simulation=-1):
    
        tstables = {}

        if simulation == -1:
            files = np.array([x for x in os.listdir(path2txts+self.channel) if x.endswith(".txt")])
        else:
            files = np.array([x for x in os.listdir(path2txts+self.channel+"/simulations/") if x.endswith(".txt")])
        for file in files:
            # Parsing the file names and checking validation.
            if simulation == -1:
                file_key = file.replace('.txt','')
                file_info = file_key.split("_")
            else:
                file_key = file.replace('_{}.txt'.format(simulation),'')
                file_info = file_key.split("_")
            if file_info[0] != self.channel or file_info[1] not in self.sources or file_info[2] not in self.collaborations.keys():
                continue
                    
            # Printing the files, which are included in the combination.
            if simulation == -1:
                print(file_info)
                
            # Opening the txt files.
            if simulation == -1:
                ts_file = open("{}/{}".format(path2txts+self.channel,file),"r")
            else:
                ts_file = open("{}/simulations/{}".format(path2txts+self.channel,file),"r")
                
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
    
    def read_JFactor(self, JFactor_file):
        logJ, DlogJ = {}, {}
        for source in self.sources:
            source_logJ, source_DlogJ = {}, {}
            for collaboration in self.collaborations:
                angular_separation = self.collaborations[collaboration]
                # Opening txt file.
                file = open(JFactor_file, "r")
                for i,line in enumerate(file):
                    if i > 41:
                        # The source in the GS file for this line
                        source_GS = line[:18].replace(" ", "")
                        line = line[18:].split()
                        if source == source_GS and float(line[0]) == angular_separation:
                            source_logJ[collaboration] = np.float32(line[3])
                            source_DlogJ[collaboration] = np.float32(line[3]) - np.float32(line[2])
                            file.close()
                            break
                if collaboration in source_logJ and collaboration in source_DlogJ:
                    logJ[source] = source_logJ
                    DlogJ[source] = source_DlogJ
                else:
                    raise ValueError("'{}' J-Factor not in '{}' for angular separation {} deg.".format(source, JFactor_file, angular_separation))
        return logJ, DlogJ
