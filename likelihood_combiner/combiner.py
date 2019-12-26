import argparse
import tables
import numpy as np
from scipy.interpolate import interp1d
import os
import yaml

from likelihood_combiner.reader import LklComReader,JFactor_Reader
from likelihood_combiner.writer import LklComWriter
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
        hdf5file = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/lklcom.h5"))

    try:
        output_dir = config['Output']['output_directory']
        if output_dir is None:
            raise KeyError
    except KeyError:
        output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output/"))
    # Create output directory if it doesn't exist already
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    channels = config['Configuration']['channels']
    sources = config['Configuration']['sources']
    collaborations = config['Configuration']['collaborations']
    
    writer = LklComWriter()
    writer.convert_txts2hdf5(hdf5file,data_dir,channels,sources,collaborations)
    del writer

    reader = LklComReader()
    tstables, massvals = reader.read_tstables(hdf5file,channels,sources,collaborations)
    del reader

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

    # Reading in the the sgimav range and spacing. Round of the third digit to avoid an interpolation.
    sigmavMin = -(np.abs(np.floor(np.log10(np.abs(float(config['Data']['sigmavMin'])))).astype(int))).astype(int)
    sigmavMax = -(np.abs(np.floor(np.log10(np.abs(float(config['Data']['sigmavMax'])))).astype(int))).astype(int)
    sigmavNPoints = int(config['Data']['sigmavNPoints'])
    sigmav = np.logspace(sigmavMin, sigmavMax, sigmavNPoints, endpoint=True, dtype=np.float32)
    exponent = (np.abs(np.floor(np.log10(np.abs(sigmav))).astype(int))+3).astype(int)
    for i,e in enumerate(exponent):
        sigmav[i] = np.around(sigmav[i],decimals=e)

    print("Combining limits for the selected configuration...")
    for channel in channels:
        print("Channel: '{}'".format(channel))
        combined_channel_ts = {}
        combined_channel_ts_Jnuisance = {}
        mass_axis = []
        for source in sources:
            for key in massvals.keys():
                if channel+"_"+source in key:
                    print("  - {}".format(source))
                    break
            # Combine ts values for each dSphs
            combined_source_ts = {}
            combined_source_ts_Jnuisance = {}
            for collaboration in collaborations:
                key = "{}_{}_{}".format(channel,source,collaboration)
                if key in massvals:
                    tstable = tstables[key]
                    masses = massvals[key]
                    sigmavFile = tstable[0]
                    exponent = (np.abs(np.floor(np.log10(np.abs(sigmavFile))).astype(int))+3).astype(int)
                    for i,e in enumerate(exponent):
                        sigmavFile[i] = np.around(sigmavFile[i],decimals=e)
                    sigmavFile = sigmavFile[::-1]
                    source_ts_dict = {}
                    for i,m in enumerate(masses[1:]):
                        if m not in mass_axis:
                            mass_axis.append(m)
                        if not np.array_equal(sigmav,sigmavFile):
                            lin_interpolation = interp1d(sigmavFile, tstable[i+1][::-1], kind='linear', fill_value='extrapolate')
                            source_ts_dict[source+"_"+str(m)] = lin_interpolation(sigmav)
                        else:
                            source_ts_dict[source+"_"+str(m)] = tstable[i+1][::-1]
                        if source+"_"+str(m) in combined_source_ts:
                            combined_source_ts[source+"_"+str(m)] += source_ts_dict[source+"_"+str(m)]
                        else:
                            combined_source_ts[source+"_"+str(m)] = source_ts_dict[source+"_"+str(m)]
        
            if config['Data']['J_nuisance']:
                combined_source_ts_Jnuisance = compute_Jnuisance(sigmav, combined_source_ts, sources_DlogJ)
                combined_source_limits_Jnuisance,combined_source_sensitivity_Jnuisance = compute_sensitivity(sigmav, combined_source_ts_Jnuisance)
            for m in mass_axis:
                if str(m) in combined_channel_ts:
                    if source+"_"+str(m) in combined_source_ts:
                        combined_channel_ts[str(m)] += combined_source_ts[source+"_"+str(m)]
                        if config['Data']['J_nuisance']:
                           combined_channel_ts_Jnuisance[str(m)] += combined_source_ts_Jnuisance[source+"_"+str(m)]
                else:
                    if source+"_"+str(m) in combined_source_ts:
                        combined_channel_ts[str(m)] = combined_source_ts[source+"_"+str(m)]
                        if config['Data']['J_nuisance']:
                            combined_channel_ts_Jnuisance[str(m)] = combined_source_ts_Jnuisance[source+"_"+str(m)]
   
        h5 = tables.open_file(hdf5file, 'a')
        if "/{}/Combination".format(channel) not in h5:
            h5.create_group(eval("h5.root.{}".format(channel)), "Combination", "Further information about the combination of sources ({}).".format(','.join(sources)))
        
        sigmav_shape = (sigmav.shape[0],)
        columns_dict={"mass":tables.Float32Col(),
                      "ts":tables.Float32Col(shape=sigmav_shape),
                      "ts_Jnuisance":tables.Float32Col(shape=sigmav_shape)}
        description = type("description", (tables.IsDescription,), columns_dict)

        # Creating the table mass vs sigmav for each source and for each channel.
        table_name = "sigmavVsMass"
        table = h5.create_table(eval("h5.root.{}.Combination".format(channel)),table_name,description,"Table of the combination of collaborations ({}) for the sources ({}) with the annihilation channel {}.".format(','.join(collaborations),','.join(sources),channel))
        
        # Filling the data of the txt file into the table of the hdf5 file.
        sorted_mass = sorted(mass_axis)
        for m in sorted_mass:
            table = eval("h5.root.{}.Combination.{}".format(channel,table_name))
            row = table.row
            row['mass']  = str(m)
            # The first element of the ts array (mass value) will be ignored, since it's in the mass column.
            row['ts'] = combined_channel_ts[str(m)]
            if config['Data']['J_nuisance']:
                row['ts_Jnuisance'] = combined_channel_ts_Jnuisance[str(m)]
            else:
                row['ts_Jnuisance'] = np.nan
                
            row.append()
            table.flush()

        combined_sources_limits,combined_sources_sensitivity = compute_sensitivity(sigmav, combined_channel_ts)
        if config['Data']['J_nuisance']:
            combined_sources_limits_Jnuisance,combined_sources_sensitivity_Jnuisance = compute_sensitivity(sigmav, combined_channel_ts_Jnuisance)

        columns_dict={"mass":tables.Float32Col(),
                      "sigmav_UL":tables.Float32Col(),
                      "sigmav_UL_Jnuisance":tables.Float32Col()}
        description = type("description", (tables.IsDescription,), columns_dict)
        table_name = "ULsigmavVsMass"
        table = h5.create_table(eval("h5.root.{}.Combination".format(channel)),table_name,description,"95% CL sigmav UL vs mass of the combination of collaborations ({}) for the sources ({}) with the annihilation channel {}.".format(','.join(collaborations),','.join(sources),channel))
        # Filling the sigmav UL into the table of the hdf5 file.
        for m in sorted_mass:
            table = eval("h5.root.{}.Combination.{}".format(channel,table_name))
            row = table.row
            row['mass']  = str(m)
            # The first element of the ts array (mass value) will be ignored, since it's in the mass column.
            row['sigmav_UL'] = combined_sources_limits[str(m)]
            if config['Data']['J_nuisance']:
                row['sigmav_UL_Jnuisance'] = combined_sources_limits_Jnuisance[str(m)]
            else:
                row['sigmav_UL_Jnuisance'] = np.nan
                
            row.append()
            table.flush()
        # Closing hdf5 file.
        h5.close()
    return

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
