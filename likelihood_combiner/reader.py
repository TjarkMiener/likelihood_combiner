import tables
import numpy as np
import os

from likelihood_combiner.gloryduck import gloryduckInfo

class gloryduckReader:
    def __init__(self):
        """Constructor"""
        # Get information from the gloryduck class
        gloryduck = gloryduckInfo()
        # List of valid collaborations:
        self.gd_collaborations = gloryduck.collaborations
        # List of valid sources:
        self.gd_sources = gloryduck.sources
        # List of valid annihilation channels:
        self.gd_channels = gloryduck.channels
    
    def read_gloryduck_tstables(self, hdf5file, channels=None, sources=None, collaborations=None):
        
        # Opening hdf5 file.
        h5 = tables.open_file(hdf5file, 'r')
        
        if channels is None:
            channels = self.gd_channels
        channels = np.array(channels)
        if sources is None:
            sources = self.gd_sources
        sources = np.array(sources)
        if collaborations is None:
            collaborations = self.gd_collaborations
        collaborations = np.array(collaborations)

        print("The tables read from '{}' ({}):".format(h5.title,hdf5file))
        counter = 1
        tstables = {}
        massvals = {}
        for channel in channels:
            if channel not in self.gd_channels:
                raise ValueError("'{}' is not a valid channel!".format(channel))
            for source in sources:
                if source not in self.gd_sources:
                    raise ValueError("'{}' is not a valid source!".format(source))
                for collaboration in collaborations:
                    if collaboration not in self.gd_collaborations:
                        raise ValueError("'{}' is not a valid collaboration!".format(collaboration))

                    if "/{}/{}/{}".format(channel,source,collaboration) in h5:
                        print("    {}) ('{}', '{}', '{}')".format(counter,channel,source,collaboration))
                        sigmavVsMassTable = eval("h5.root.{}.{}.{}.sigmavVsMass".format(channel,source,collaboration))
                        table = []
                        mass = []
                        for x in sigmavVsMassTable.iterrows():
                            table.append(x['ts'])
                            mass.append(x['mass'])
                        table_info = "{}_{}_{}".format(channel,source,collaboration)
                        tstables[table_info] = np.array(table)
                        massvals[table_info] = np.array(mass)
                        counter+=1
         
        # Closing hdf5 file.
        h5.close()
        return tstables, massvals


class JFactor_Reader:
    def __init__(self):
        """Constructor"""
        # Get information from the gloryduck class
        gloryduck = gloryduckInfo()
        # List of valid collaborations:
        self.gd_collaborations = gloryduck.collaborations
        # List of valid sources:
        self.gd_sources = gloryduck.sources
        
        # Define the angular separation from dwarf center:
        self.angular_separation = 0.53086117
    
    def read_JFactor(self, JFactor_file, sources=None):
        if sources is None:
            sources = self.gd_sources
        sources = np.array(sources)
        
        sources_logJ = {}
        sources_DlogJ = {}
        angle = self.angular_separation
        for source in sources:
            # Opening txt file.
            file = open(JFactor_file, "r")
            for i,line in enumerate(file):
                if i > 41:
                    # The source in the GS file for this line
                    source_GS = line[:18].replace(" ", "")
                    line = line[18:].split()
                    if source == source_GS and float(line[0]) == angle:
                        sources_logJ[source] = float(line[3])
                        sources_DlogJ[source] = float(line[3]) - float(line[2])
                        file.close()
                        break
        return sources_logJ,sources_DlogJ

class cirelliReader:
    def __init__(self):
        """Constructor"""

    def read_dNdE(self,hdf5file):
        """Read the dN/dE from the crielli table in hdf5 format"""
        return
