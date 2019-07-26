import argparse
import tables
import numpy as np
import os
import yaml
import matplotlib.pyplot as plt

from likelihood_combiner.reader import gloryduckReader,JFactor_Reader
from likelihood_combiner.writer import gloryduckWriter
from likelihood_combiner.gloryduck import gloryduckInfo
from likelihood_combiner.utils import compute_sensitivity,compute_Jnuisance,plot_sigmavULs

def run_combiner(config):
    
    try:
        data_dir = config['Data']['data_directory']
        if data_dir is None:
            raise KeyError
    except KeyError:
        data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/"))

    try:
        hdf5file = config['Data']['hdf5_dataset']
        if hdf5file is None:
            raise KeyError
    except KeyError:
        hdf5file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/gloryduck_dataset.h5"))

    try:
        output_dir = config['Output']['output_directory']
        if output_dir is None:
            raise KeyError
    except KeyError:
        output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output/"))
    # Create output directory if it doesn't exist already
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    gd_writer = gloryduckWriter()
    gd_writer.convert_txts2hdf5(hdf5file,data_dir)
    del gd_writer

    channels = config['Configuration']['channels']
    sources = config['Configuration']['sources']
    collaborations = config['Configuration']['collaborations']
    
    # Get information from the gloryduck class
    gloryduck = gloryduckInfo()
    if channels is None:
        channels = np.array(gloryduck.channels)
    if sources is None:
        sources = np.array(gloryduck.sources)
    if collaborations is None:
        collaborations = np.array(gloryduck.collaborations)
    channels_LaTex = gloryduck.channels_LaTex

    gd_reader = gloryduckReader()
    tstables, massvals = gd_reader.read_gloryduck_tstables(hdf5file,channels,sources,collaborations)
    del gd_reader

    try:
        JFactor_file = config['Data']['JFactor_table']
        if JFactor_file is None:
            raise KeyError
    except KeyError:
        JFactor_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/Jfactor_Geringer-SamethTable.txt"))

    if config['Data']['J_nuisance']:
        Jfactor_reader = JFactor_Reader()
        sources_logJ,sources_DlogJ = Jfactor_reader.read_JFactor(JFactor_file,sources)
        del Jfactor_reader
    
    print("Combining limits for the selected configuration...")
    for channel in channels:
        sigmav = None
        mass_axis = []
        combined_sources = []
        combined_sources_ts = {}
        combined_sources_ts_Jnuisance = {}
        for source in sources:
            # Combine ts values for each dSphs
            combined_source_ts = {}
            combined_source_ts_Jnuisance = {}
            combined_collaborations = []
            tstable_ref = None
            for key,mass,tstable in zip(massvals.keys(),massvals.values(),tstables.values()):
                key_split = key.split("_")
                if channel == key_split[0] and source == key_split[1]:
                    # Checking that all ranges (mass and sigmav) and the J-Factor (first element of the mass arrays) are equal
                    collaboration = key_split[2]
                    combined_collaborations.append(collaboration)
                    if source not in combined_sources:
                        combined_sources.append(source)
                    if tstable_ref is None:
                        tstable_ref = sigmav = tstable[0]
                        exponent = (np.abs(np.floor(np.log10(np.abs(sigmav))).astype(int))+3).astype(int)
                        for i,e in enumerate(exponent):
                            sigmav[i] = np.around(sigmav[i],decimals=e)
                        tstable_key_ref = key
                    else:
                        if (sigmav!=tstable_ref).any():
                            raise ValueError("The sigma values have to be equal! Discrepancy in '{}.txt' and '{}.txt'".format(key,tstable_key_ref))
                    sigmav_rev = sigmav[::-1]
                    source_ts_dict = {}
                    for i,m in enumerate(mass[1:]):
                        if m not in mass_axis:
                            mass_axis.append(m)
                        source_ts_dict[source+"_"+str(m)] = tstable[i+1][::-1]
                        if source+"_"+str(m) in combined_source_ts:
                            combined_source_ts[source+"_"+str(m)] += source_ts_dict[source+"_"+str(m)]
                        else:
                            combined_source_ts[source+"_"+str(m)] = source_ts_dict[source+"_"+str(m)]
                        if config['Data']['J_nuisance']:
                            source_ts_dict_Jnuisance = compute_Jnuisance(sigmav_rev, source_ts_dict, sources_DlogJ)
                            if source+"_"+str(m) in combined_source_ts_Jnuisance:
                                combined_source_ts_Jnuisance[source+"_"+str(m)] += source_ts_dict_Jnuisance[source+"_"+str(m)]
                            else:
                                combined_source_ts_Jnuisance[source+"_"+str(m)] = source_ts_dict_Jnuisance[source+"_"+str(m)]
                    
                    h5 = tables.open_file(hdf5file, 'a')
                    columns_dict={"mass":tables.Float32Col(),
                                  "sigmav_UL":tables.Float32Col(),
                                  "sigmav_UL_Jnuisance":tables.Float32Col()}
                    description = type("description", (tables.IsDescription,), columns_dict)
                    table_name = "ULsigmavVsMass"
                    table = h5.create_table(eval("h5.root.{}.{}.{}".format(channel,source,collaboration)),table_name,description,"95% CL sigmav UL vs mass of the {} collaboration for the source {} with the annihilation channel {}.".format(collaboration,source,channel))

                    sensitivity = compute_sensitivity(sigmav_rev, source_ts_dict)
                    if config['Data']['J_nuisance']:
                        sensitivity_Jnuisance = compute_sensitivity(sigmav_rev, source_ts_dict_Jnuisance)

                    # Filling the sigmav UL into the table of the hdf5 file.
                    for key in sensitivity.keys():
                        key_split = key.split("_")
                        table = eval("h5.root.{}.{}.{}.{}".format(channel,source,collaboration,table_name))
                        row = table.row
                        row['mass']  = key_split[1]
                            # The first element of the ts array (mass value) will be ignored, since it's in the mass column.
                        row['sigmav_UL'] = sensitivity[key]
                        if config['Data']['J_nuisance']:
                            row['sigmav_UL_Jnuisance'] = sensitivity_Jnuisance[key]
                        else:
                            row['sigmav_UL_Jnuisance'] = np.nan

                        row.append()
                        table.flush()
                    # Closing hdf5 file.
                    h5.close()
        
            if not ','.join(combined_collaborations):
                continue

            h5 = tables.open_file(hdf5file, 'a')
            if "/{}/{}/Combination".format(channel,source) not in h5:
                h5.create_group(eval("h5.root.{}.{}".format(channel,source)), "Combination", "Further information about the combination of collaborations ({}).".format(','.join(combined_collaborations)))

            combined_source_sensitivity = compute_sensitivity(sigmav_rev, combined_source_ts)
            if config['Data']['J_nuisance']:
                combined_source_sensitivity_Jnuisance = compute_sensitivity(sigmav_rev, combined_source_ts_Jnuisance)
            
            # Creating the table structure for the hdf5 file.
            sigmav_shape = (sigmav_rev.shape[0],)
            columns_dict={"mass":tables.Float32Col(),
                          "ts":tables.Float32Col(shape=sigmav_shape),
                          "ts_Jnuisance":tables.Float32Col(shape=sigmav_shape)}
            description = type("description", (tables.IsDescription,), columns_dict)
                            
            # Creating the table mass vs sigmav for each source and for each channel.
            table_name = "sigmavVsMass"
            table = h5.create_table(eval("h5.root.{}.{}.Combination".format(channel,source)),table_name,description,"Table of the combination of collaborations ({}) for the source {} with the annihilation channel {}.".format(','.join(combined_collaborations),source,channel))
                                
            # Filling the data of the txt file into the table of the hdf5 file.
            for key in combined_source_sensitivity.keys():
                key_split = key.split("_")
                table = eval("h5.root.{}.{}.Combination.{}".format(channel,source,table_name))
                row = table.row
                row['mass']  = key_split[1]
                # The first element of the ts array (mass value) will be ignored, since it's in the mass column.
                row['ts'] = combined_source_ts[key]
                if config['Data']['J_nuisance']:
                    row['ts_Jnuisance'] = combined_source_ts_Jnuisance[key]
                else:
                    row['ts_Jnuisance'] = np.nan

                row.append()
                table.flush()
            # Closing hdf5 file.
            h5.close()
            h5 = tables.open_file(hdf5file, 'a')
            columns_dict={"mass":tables.Float32Col(),
                          "sigmav_UL":tables.Float32Col(),
                          "sigmav_UL_Jnuisance":tables.Float32Col()}
            description = type("description", (tables.IsDescription,), columns_dict)
            table_name = "ULsigmavVsMass"
            table = h5.create_table(eval("h5.root.{}.{}.Combination".format(channel,source)),table_name,description,"95% CL sigmav UL vs mass of the combination of collaborations ({}) for the source {} with the annihilation channel {}.".format(','.join(combined_collaborations),source,channel))
            # Filling the sigmav UL into the table of the hdf5 file.
            for key in combined_source_sensitivity.keys():
                key_split = key.split("_")
                table = eval("h5.root.{}.{}.Combination.{}".format(channel,source,table_name))
                row = table.row
                row['mass']  = key_split[1]
                # The first element of the ts array (mass value) will be ignored, since it's in the mass column.
                row['sigmav_UL'] = combined_source_sensitivity[key]
                if config['Data']['J_nuisance']:
                    row['sigmav_UL_Jnuisance'] = combined_source_sensitivity_Jnuisance[key]
                else:
                    row['sigmav_UL_Jnuisance'] = np.nan
                
                row.append()
                table.flush()
            # Closing hdf5 file.
            h5.close()

            for m in mass_axis:
                if str(m) in combined_sources_ts:
                    combined_sources_ts[str(m)] += combined_source_ts[source+"_"+str(m)]
                    if config['Data']['J_nuisance']:
                        combined_sources_ts_Jnuisance[str(m)] += combined_source_ts_Jnuisance[source+"_"+str(m)]
                else:
                    combined_sources_ts[str(m)] = combined_source_ts[source+"_"+str(m)]
                    if config['Data']['J_nuisance']:
                        combined_sources_ts_Jnuisance[str(m)] = combined_source_ts_Jnuisance[source+"_"+str(m)]
                            
        h5 = tables.open_file(hdf5file, 'a')
        if "/{}/Combination_ALL".format(channel) not in h5:
            h5.create_group(eval("h5.root.{}".format(channel)), "Combination_ALL", "Further information about the combination of sources ({}).".format(','.join(combined_sources)))
        
        columns_dict={"mass":tables.Float32Col(),
                      "ts":tables.Float32Col(shape=sigmav_shape),
                      "ts_Jnuisance":tables.Float32Col(shape=sigmav_shape)}
        description = type("description", (tables.IsDescription,), columns_dict)

        # Creating the table mass vs sigmav for each source and for each channel.
        table_name = "sigmavVsMass"
        table = h5.create_table(eval("h5.root.{}.Combination_ALL".format(channel)),table_name,description,"Table of the combination of collaborations ({}) for the sources ({}) with the annihilation channel {}.".format(','.join(collaborations),','.join(combined_sources),channel))
        
        # Filling the data of the txt file into the table of the hdf5 file.
        for key in combined_sources_ts.keys():
            table = eval("h5.root.{}.Combination_ALL.{}".format(channel,table_name))
            row = table.row
            row['mass']  = key
            # The first element of the ts array (mass value) will be ignored, since it's in the mass column.
            row['ts'] = combined_sources_ts[key]
            if config['Data']['J_nuisance']:
                row['ts_Jnuisance'] = combined_sources_ts_Jnuisance[key]
            else:
                row['ts_Jnuisance'] = np.nan
                
            row.append()
            table.flush()

        combined_sources_sensitivity = compute_sensitivity(sigmav_rev, combined_sources_ts)
        if config['Data']['J_nuisance']:
            combined_sources_sensitivity_Jnuisance = compute_sensitivity(sigmav_rev, combined_sources_ts_Jnuisance)

        columns_dict={"mass":tables.Float32Col(),
                      "sigmav_UL":tables.Float32Col(),
                      "sigmav_UL_Jnuisance":tables.Float32Col()}
        description = type("description", (tables.IsDescription,), columns_dict)
        table_name = "ULsigmavVsMass"
        table = h5.create_table(eval("h5.root.{}.Combination_ALL".format(channel)),table_name,description,"95% CL sigmav UL vs mass of the combination of collaborations ({}) for the sources ({}) with the annihilation channel {}.".format(','.join(collaborations),','.join(combined_sources),channel))
        # Filling the sigmav UL into the table of the hdf5 file.
        for key in combined_sources_sensitivity.keys():
            table = eval("h5.root.{}.Combination_ALL.{}".format(channel,table_name))
            row = table.row
            row['mass']  = key
            # The first element of the ts array (mass value) will be ignored, since it's in the mass column.
            row['sigmav_UL'] = combined_sources_sensitivity[key]
            if config['Data']['J_nuisance']:
                row['sigmav_UL_Jnuisance'] = combined_sources_sensitivity_Jnuisance[key]
            else:
                row['sigmav_UL_Jnuisance'] = np.nan
                
            row.append()
            table.flush()
        # Closing hdf5 file.
        h5.close()

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
            description=("Combining likelihoods from different experiments."))
    parser.add_argument(
            'config_file',
            help="path to YAML configuration file with combining options")
        
    args = parser.parse_args()
                                     
    with open(args.config_file, 'r') as config_file:
        config = yaml.load(config_file)

    if not config['Output']['only_plotting']:
        run_combiner(config)

    plot_sigmavULs(config)
