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

    def read_gloryduck_sigmavULs(self, hdf5file, channels=None, sources=None, collaborations=None):
    
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
        collaborations.append('Combination')
        collaborations = np.array(collaborations)

        print("The sigmav ULs read from '{}' ({}):".format(h5.title,hdf5file))
        counter = 1
        sigmavULs = {}
        sigmavULs_Jnuisance = {}
        massvals = {}
        for channel in channels:
            if channel not in self.gd_channels:
                raise ValueError("'{}' is not a valid channel!".format(channel))
            
            if "/{}/Combination_all".format(channel) in h5:
                print("    {}) ('{}', 'Combination_all')".format(counter,channel))
                ULsigmavVsMass = eval("h5.root.{}.Combination_all.ULsigmavVsMass".format(channel))
                sigmavUL = []
                sigmavUL_Jnuisance = []
                mass = []
                for x in ULsigmavVsMass.iterrows():
                    sigmavUL.append(x['sigmav_UL'])
                    sigmavUL_Jnuisance.append(x['sigmav_UL_Jnuisance'])
                    mass.append(x['mass'])
                table_info = "{}_Combination_all".format(channel)
                sigmavULs[table_info] = np.array(sigmavUL)
                sigmavULs_Jnuisance[table_info] = np.array(sigmavUL_Jnuisance)
                massvals[table_info] = np.array(mass)
                counter+=1
            for source in sources:
                if source not in self.gd_sources:
                    raise ValueError("'{}' is not a valid source!".format(source))
                for collaboration in collaborations:
                    if collaboration not in self.gd_collaborations and collaboration != 'Combination':
                        raise ValueError("'{}' is not a valid collaboration!".format(collaboration))
                    
                    if "/{}/{}/{}".format(channel,source,collaboration) in h5:
                        print("    {}) ('{}', '{}', '{}')".format(counter,channel,source,collaboration))
                        ULsigmavVsMass = eval("h5.root.{}.{}.{}.ULsigmavVsMass".format(channel,source,collaboration))
                        sigmavUL = []
                        sigmavUL_Jnuisance = []
                        mass = []
                        for x in ULsigmavVsMass.iterrows():
                            sigmavUL.append(x['sigmav_UL'])
                            sigmavUL_Jnuisance.append(x['sigmav_UL_Jnuisance'])
                            mass.append(x['mass'])
                        table_info = "{}_{}_{}".format(channel,source,collaboration)
                        sigmavULs[table_info] = np.array(sigmavUL)
                        sigmavULs_Jnuisance[table_info] = np.array(sigmavUL_Jnuisance)
                        massvals[table_info] = np.array(mass)
                        counter+=1
        # Closing hdf5 file.
        h5.close()
        return sigmavULs, sigmavULs_Jnuisance, massvals

    def read_gloryduck_sigmavULs_collaborations(self, data_dir, hdf5file, channels=None, collaborations=None):
        
        if channels is None:
            channels = self.gd_channels
        channels = np.array(channels)
        if collaborations is None:
            collaborations = self.gd_collaborations
        collaborations = np.array(collaborations)
   
        sigmav_UL = {}
        sigmav_UL_Jnuisance = {}
        massvals = {}
        for channel in channels:
            for collaboration in collaborations:
                h5 = tables.open_file("{}gloryduck_{}.h5".format(data_dir,collaboration), 'r')
                if "/{}/Combination_all".format(channel) in h5:
                    ULsigmavVsMassTable = eval("h5.root.{}.Combination_all.ULsigmavVsMass".format(channel))
                    ul = []
                    ul_Jnuisance = []
                    mass = []
                    for x in ULsigmavVsMassTable.iterrows():
                        ul.append(x['sigmav_UL'])
                        ul_Jnuisance.append(x['sigmav_UL_Jnuisance'])
                        mass.append(x['mass'])
                    table_info = "{}_{}".format(channel,collaboration)
                    sigmav_UL[table_info] = np.array(ul)
                    sigmav_UL_Jnuisance[table_info] = np.array(ul_Jnuisance)
                    massvals[table_info] = np.array(mass)
            
                # Closing hdf5 file
                h5.close()
            h5 = tables.open_file(hdf5file, 'r')
            if "/{}/Combination_all".format(channel) in h5:
                ULsigmavVsMassTable = eval("h5.root.{}.Combination_all.ULsigmavVsMass".format(channel))
                ul = []
                ul_Jnuisance = []
                mass = []
                for x in ULsigmavVsMassTable.iterrows():
                    ul.append(x['sigmav_UL'])
                    ul_Jnuisance.append(x['sigmav_UL_Jnuisance'])
                    mass.append(x['mass'])
                table_info = "{}_Combination_all".format(channel)
                sigmav_UL[table_info] = np.array(ul)
                sigmav_UL_Jnuisance[table_info] = np.array(ul_Jnuisance)
                massvals[table_info] = np.array(mass)
        return sigmav_UL, sigmav_UL_Jnuisance, massvals

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
