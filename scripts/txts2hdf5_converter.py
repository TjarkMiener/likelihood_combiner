#!/usr/bin/python'
import tables
from tables import *
import h5py
import numpy as np
import os

######################################################
#              Modify information here:              #
######################################################

# List of the collaborations:
collaborations = ['MAGIC','HESS','FermiLAT','VERITAS','HAWC']

# List of sources:
sources = ['Booetes1', 'Carina', 'ComaBerenices', 'Draco', 'Fornax', 'UrsaMajorII', 'UrsaMinor', 'Sagittarius', 'Sculptor', 'Segue1', 'Willman1']

# List of annihilation channels:
channels = ['bb', 'tautau', 'mumu', 'WW', 'gammagamma', 'ZZ','ee']

#######################################################
#                     Thank you!                      #
#######################################################

# Opening the hdf5 file.
h5file_name = "gloryduck_dataset.h5"
h5file = open_file(h5file_name, mode="w", title="Test file")

# Reading the .txt files in the folder.
files = np.array([x for x in os.listdir() if x.endswith(".txt")])
print("The files written in '{}' ({}):".format(h5file.title,h5file_name))

for counter,file in enumerate(files):
    
    # Parsing the file name and checking validation with the predefined information.
    file_info = file.replace('.txt','').split("_")
    if file_info[0] not in channels:
        raise ValueError("'{}' is not a valid channel!".format(file_info[0]))
    if file_info[1] not in sources:
        raise ValueError("'{}' is not a valid source!".format(file_info[1]))
    if file_info[2] not in collaborations:
        raise ValueError("'{}' is not a valid collaboration!".format(file_info[2]))
    
    # Creating a new node in the hdf5 file, if it isn't already existing.
    if "/{}".format(file_info[1]) not in h5file:
        h5file.create_group(h5file.root,file_info[1],"Further information about the source {}".format(file_info[1]))
    if "/{}/{}".format(file_info[1],file_info[0]) not in h5file:
        h5file.create_group(eval("h5file.root.{}".format(file_info[1])), file_info[0], "Further information about the annihilation channel {}".format(file_info[0]))
    if "/{}/{}/{}".format(file_info[1],file_info[0],file_info[2]) not in h5file:
        h5file.create_group(eval("h5file.root.{}.{}".format(file_info[1],file_info[0])), file_info[2], "Further information about the collaboration {}".format(file_info[2]))

    # Opening the txt files.
    ts_file = open(file,"r")
    print("    {}) '{}'".format(counter+1,file))
        
    # Going through the table in the txt file and storing the entries in a 2D array.
    ts_val = np.array([[i for i in line.split()] for line in ts_file])
            
    # The first entry of each row correponds to the mass.
    # Detect the first entry and store it in a separate array.
    mass = np.array([ts_val[i][0] for i in np.arange(0,ts_val.shape[0],1)])
    
    # Closing the txt files.
    ts_file.close()
        
    # Creating the table structure for the hdf5 file.
    sigmav_shape = (ts_val.shape[1]-1,)
    columns_dict={"mass":tables.Float32Col(),
                  "ts":tables.Float32Col(shape=sigmav_shape)}
    description = type("description", (tables.IsDescription,), columns_dict)
        
    # Creating the table mass vs sigmav for each source and for each channel.
    table_name = "sigmavVsMass_{}_{}_{}".format(file_info[1],file_info[0],file_info[2])
    table = h5file.create_table(eval("h5file.root.{}.{}.{}".format(file_info[1],file_info[0],file_info[2])),table_name,description,"Table of the {} collaboration for the source {} with the annihilation channel {}.".format(file_info[2],file_info[1],file_info[0]))
        
    # Filling the data of the txt file into the table of the hdf5 file.
    for i, mass_val in enumerate(mass):
        table = eval("h5file.root.{}.{}.{}.{}".format(file_info[1],file_info[0],file_info[2],table_name))
        row = table.row
        row['mass']  = mass_val
        # The first element of the ts array (mass value) will be ignored, since it's in the mass column.
        row['ts'] = ts_val[i][1:]
        row.append()
        table.flush()

# Closing hdf5 file.
h5file.close()
