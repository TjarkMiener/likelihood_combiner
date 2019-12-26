import tables
import numpy as np
import os

class LklComReader:
    def __init__(self):
        """Constructor"""
    
    def read_tstables(self, hdf5file, channels, sources, collaborations):
        
        # Opening hdf5 file.
        h5 = tables.open_file(hdf5file, 'r')
        
        print("The tables read from '{}' ({}):".format(h5.title,hdf5file))
        counter = 1
        tstables = {}
        massvals = {}
        for channel in channels:
            for source in sources:
                for collaboration in collaborations:
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

    def read_sigmavULs(self, hdf5file, channels):
    
        # Opening hdf5 file.
        h5 = tables.open_file(hdf5file, 'r')
    
        print("The sigmav ULs read from '{}' ({}):".format(h5.title,hdf5file))
        counter = 1
        sigmavULs = {}
        sigmavULs_Jnuisance = {}
        massvals = {}
        for channel in channels:
            if "/{}/Combination".format(channel) in h5:
                print("    {}) ('{}', 'Combination')".format(counter,channel))
                ULsigmavVsMass = eval("h5.root.{}.Combination.ULsigmavVsMass".format(channel))
                sigmavUL = []
                sigmavUL_Jnuisance = []
                mass = []
                for x in ULsigmavVsMass.iterrows():
                    sigmavUL.append(x['sigmav_UL'])
                    sigmavUL_Jnuisance.append(x['sigmav_UL_Jnuisance'])
                    mass.append(x['mass'])
                table_info = "{}_Combination".format(channel)
                sigmavULs[table_info] = np.array(sigmavUL)
                sigmavULs_Jnuisance[table_info] = np.array(sigmavUL_Jnuisance)
                massvals[table_info] = np.array(mass)
                counter+=1
        # Closing hdf5 file.
        h5.close()
        return sigmavULs, sigmavULs_Jnuisance, massvals

    def read_sigmavULs_collaborations(self, data_dir, hdf5file, channels, collaborations):
        sigmav_UL = {}
        sigmav_UL_Jnuisance = {}
        massvals = {}
        for channel in channels:
            for collaboration in collaborations:
                h5 = tables.open_file("{}lklcom_{}.h5".format(data_dir,collaboration), 'r')
                if "/{}/Combination".format(channel) in h5:
                    ULsigmavVsMassTable = eval("h5.root.{}.Combination.ULsigmavVsMass".format(channel))
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
            if "/{}/Combination".format(channel) in h5:
                ULsigmavVsMassTable = eval("h5.root.{}.Combination.ULsigmavVsMass".format(channel))
                ul = []
                ul_Jnuisance = []
                mass = []
                for x in ULsigmavVsMassTable.iterrows():
                    ul.append(x['sigmav_UL'])
                    ul_Jnuisance.append(x['sigmav_UL_Jnuisance'])
                    mass.append(x['mass'])
                table_info = "{}_Combination".format(channel)
                sigmav_UL[table_info] = np.array(ul)
                sigmav_UL_Jnuisance[table_info] = np.array(ul_Jnuisance)
                massvals[table_info] = np.array(mass)
        return sigmav_UL, sigmav_UL_Jnuisance, massvals

class JFactor_Reader:
    def __init__(self):
        """Constructor"""
        # Define the angular separation from dwarf center:
        self.angular_separation = 0.53086117
    
    def read_JFactor(self, JFactor_file, sources):
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
