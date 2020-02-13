import numpy as np
import os
import pandas as pd

class LklComReader:
    def __init__(self,
                 channels,
                 sources,
                 collaborations,
                 angular_separation=None):

        self.channels = channels
        self.sources = sources
        self.collaborations = collaborations
        
        # Define the angular separation from dwarf center:
        if angular_separation is None:
            angular_separation = 0.53086117
        self.angular_separation = angular_separation
    
    def read_tstables(self, path2txts, simulation=-1):
    
        tstables = {}
        for channel in self.channels:
            if simulation == -1:
                files = np.array([x for x in os.listdir(path2txts+channel) if x.endswith(".txt")])
            else:
                files = np.array([x for x in os.listdir(path2txts+channel+"/simulations/") if x.endswith(".txt")])
            for file in files:
                # Parsing the file names and checking validation.
                if simulation == -1:
                    file_key = file.replace('.txt','')
                    file_info = file_key.split("_")
                else:
                    file_key = file.replace('_{}.txt'.format(simulation),'')
                    file_info = file_key.split("_")
                if file_info[0] != channel or file_info[1] not in self.sources or file_info[2] not in self.collaborations:
                    continue
                    
                # Printing the files, which are included in the combination
                if simulation == -1:
                    print(file_info)
                
                # Opening the txt files.
                if simulation == -1:
                    ts_file = open("{}/{}".format(path2txts+channel,file),"r")
                else:
                    ts_file = open("{}/simulations/{}".format(path2txts+channel,file),"r")
                
                # Going through the table in the txt file and storing the entries in a 2D array.
                values = np.array([[i for i in line.split()] for line in ts_file], dtype=np.float32).T
                    
                # The first entry of each row correponds to the mass.
                # Detect the first entry and store it in a separate array.
                mass = np.array([values[i][0] for i in np.arange(0,values.shape[0],1)], dtype=np.float32)
                
                # Delete the first entry, so only TS values are stored in ts_val.
                ts_val = []
                for val in values:
                    ts_val.append(val[1:])
                
                # Store the arrays in the dictionary.
                tstables[file_key+'_ts'] = np.array(ts_val, dtype=np.float32)
                tstables[file_key+'_masses'] = mass
        return tstables
    
    def read_JFactor(self, JFactor_file):
        sources_logJ = {}
        sources_DlogJ = {}
        for source in self.sources:
            # Opening txt file.
            file = open(JFactor_file, "r")
            for i,line in enumerate(file):
                if i > 41:
                    # The source in the GS file for this line
                    source_GS = line[:18].replace(" ", "")
                    line = line[18:].split()
                    if source == source_GS and float(line[0]) == self.angular_separation:
                        sources_logJ[source] = float(line[3])
                        sources_DlogJ[source] = float(line[3]) - float(line[2])
                        file.close()
                        break
        return sources_logJ,sources_DlogJ
