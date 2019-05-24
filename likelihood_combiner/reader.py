#!/usr/bin/python'

import tables
import numpy as np
import os

class gloryduckReader:
    def __init__(self):
        """Constructor"""
    
    def read_gloryduck_tstable(self,hdf5file, source=None, channel=None, collaboration=None):
        # Opening hdf5 file.
        h5 = tables.open_file(hdf5file, 'r')
        
        channels = np.array(list(h5.keys()))
        print(channels)
        #for row in f.root.MAGIC.iterrows():
        
        for i,channel in enumerate(channels):
            print(i)
            print(channel)
            print("channels:")
            print(h5.root.channels[i])
        # Closing hdf5 file.
        h5.close()
        return channels

class cirelliReader:
    def __init__(self):
        """Constructor"""

    def read_dNdE(self,hdf5file):
        """Read the dN/dE from the crielli table in hdf5 format"""
        return
