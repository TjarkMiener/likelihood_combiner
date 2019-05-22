#!/usr/bin/python'
import tables
from tables import *
import h5py
import numpy as np

######################################################
#        Please insert your information here:        #
######################################################

# Name of your collaboration:
experiment = 'MAGIC'

# List of sources your collaboration has observed:
#source_str = ['Booetes1', 'Carina', 'ComaBerenices', 'Draco', 'Fornax', 'UrsaMajorII', 'UrsaMinor', 'Sagittarius', 'Sculptor', 'Segue1', 'Willman1']
source_str = ['Booetes1', 'Carina', 'ComaBerenices', 'Draco', 'Fornax', 'UrsaMajorII', 'UrsaMinor', 'Sagittarius', 'Sculptor', 'Segue1', 'Willman1']

# List of annihilation channels you want to provide:
channels_str = ['bb', 'tautau', 'mumu', 'WW', 'gammagamma', 'pi0pi0', 'pi0gamma']

#######################################################
#                     Thank you!                      #
#######################################################

# Opening the hdf5 file
h5file = open_file("glory_duck_{}.h5".format(experiment), mode="w", title="Test file")

# Creating the structure of the hdf5 file
h5file.create_group(h5file.root, experiment, "Further information about {}".format(experiment))
for i in source_str:
    h5file.create_group(eval("h5file.root.{}".format(experiment)), i, "Further information about  the source {}".format(i))
    for j in channels_str:
        h5file.create_group(eval("h5file.root.{}.{}".format(experiment,i)), j, "Further information about the annihilation channel {}".format(j))

# Looping over the structure of the hdf5 file
for source in source_str:
    for channel in channels_str:
        
        # Getting data from the txt files
        minus2logL_file = open("MAGIC1/minus2logL_{}_{}_{}.txt".format(experiment,source,channel),"r")
        ts_file = open("MAGIC1/ts_{}_{}_{}.txt".format(experiment,source,channel),"r")
        
        # Going through the table in the txt file and storing the entries in a 2D array
        minus2logL_val = np.array([[i for i in line.split()] for line in minus2logL_file])
        ts_val = np.array([[i for i in line.split()] for line in ts_file])
        
        # The first entry of each row correponds to the mass.
        # Detect the first entry and store it in a separate array.
        mass = np.array([minus2logL_val[i][0] for i in np.arange(0,minus2logL_val.shape[0],1)])
        
        # Closing the txt files
        minus2logL_file.close()
        ts_file.close()
        
        # Creating the table structure for the hdf5 file
        sigmav_shape = (minus2logL_val.shape[1]-1,)
        columns_dict={"mass":tables.Float32Col(),
                      "minus2logL":tables.Float32Col(shape=sigmav_shape),
                      "ts":tables.Float32Col(shape=sigmav_shape)}
        description = type("description", (tables.IsDescription,), columns_dict)
        
        # Creating the table mass vs sigmav for each source and for each channel
        table_name = "sigmavVSmass_{}_{}_{}".format(experiment,source,channel)
        table = h5file.create_table(eval("h5file.root.{}.{}.{}".format(experiment,source,channel)),table_name,description,"Table of the {} collaboration for the source {} with the annihilation channel {}.".format(experiment,source,channel))
        
        # Filling the data of the txt file into the table of the hdf5 file
        for i, mass_val in enumerate(mass):
            table = eval("h5file.root.{}.{}.{}.{}".format(experiment,source,channel,table_name))
            row = table.row
            row['mass']  = mass_val
            # The first element of the minus2logL and ts array (mass value) will be ignored, since it's in the mass column
            row['minus2logL'] = minus2logL_val[i][1:]
            row['ts'] = ts_val[i][1:]
            row.append()
            table.flush()

# Closing hdf5 file
h5file.close()

