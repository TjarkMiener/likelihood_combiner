import tables
import numpy as np
import os

__all__ = [
    'write_to_lklcom',
    'gLike_to_lklcom',
    'lklcom_to_gLike'
]

def write_to_lklcom(collaboration,
                    source,
                    channel,
                    logJ,
                    sigmav_range,
                    lkl_dict,
                    output_file,
                    mode="a",
                    simulation=-1):
    """
    Write/append a given likelihood table into lklcom hdf5 file.
    Parameters
    ----------
    collaboration: `string`
        name of the collaboration.
    source: `string`
        name of the source.
    channel: `string`
        name of the channel.
    logJ: `numpy.float32`
        value of the log J-Factor.
    sigmav_range: `numpy.ndarray of type numpy.float32`
        sigmav range (ascending).
    lkl_dict: `dict`
        dictionary with items:
            keys: `string`
                DM mass.
            values: `numpy.ndarray of type numpy.float32`
                likelihood or ts values (ascending).
    output_file: `string`
        path to the lklcom hdf5 file.
    mode: `string`
        mode to open the lklcom hdf5 file.
        Valid (recommended) option:
            "w": Write; a new file is created (an existing file
                with the same name would be deleted).
            "a": Append; an existing file is opened for reading and
                writing, and if the file does not exist it is created.
            "r+": It is similar to ‘a’, but the file must already exist.
        Default: "a"
    simulation: `int`
        number of the simulation.
        Default: -1 (= data sample)
    Returns
    -------
    
    """

    # Opening the hdf5 file.
    h5 = tables.open_file(output_file, mode=mode, title="LklCom database")
    # Creating a new node in the hdf5 file, if it isn't already existing.
    if "/{}".format(collaboration) not in h5:
        h5.create_group(h5.root, collaboration, "Further information about the collaboration {}".format(collaboration))
    if "/{}/{}".format(collaboration, source) not in h5:
         h5.create_group(eval("h5.root.{}".format(collaboration)), source, "Further information about the source {}".format(source))
    if "/{}/{}/{}".format(collaboration, source, channel) not in h5:
        h5.create_group(eval("h5.root.{}.{}".format(collaboration, source)), channel, "Further information about the annihilation channel {}".format(channel))
    
    # Creating the table structure for the hdf5 file.
    sigmav_shape = (sigmav_range.shape[0],)
    columns_dict={"masses":tables.Float32Col(),
                  "ts_values":tables.Float32Col(shape=sigmav_shape)}
    description = type("description", (tables.IsDescription,), columns_dict)
    
    # Creating the table mass vs sigmav for each source and for each channel.
    table_name = "data"
    if simulation != -1:
        table_name = "simu_{}".format(simulation)
    table = h5.create_table(eval("h5.root.{}.{}.{}".format(collaboration, source, channel)),table_name,description,"Table of the {} collaboration for the source {} with the annihilation channel {}.".format(collaboration, source, channel))
    
    # Filling the data into the table of the hdf5 file.
    # First row (row 0) contains the log J-Factor and the sigmav range (ascending).
    table = eval("h5.root.{}.{}.{}.{}".format(collaboration, source, channel, table_name))
    row = table.row
    row['masses'] = np.float32(logJ)
    row['ts_values'] = np.array(sigmav_range, dtype=np.float32)
    row.append()
    table.flush()
    # Other row (row >0) contains the mass and the likelihood or ts values (ascending).
    for mass in lkl_dict:
        table = eval("h5.root.{}.{}.{}.{}".format(collaboration, source, channel, table_name))
        row = table.row
        row['masses'] = np.float32(mass)
        row['ts_values'] = np.array(lkl_dict[mass], dtype=np.float32)
        row.append()
        table.flush()

    # Closing hdf5 file.
    h5.close()
    return


def gLike_to_lklcom(input_dir,
                    output_file,
                    mode="w"):
    """
    Translate gLike txt files into lklcom hdf5 file.
    Parameters
    ----------
    input_dir: `string`
        path to the input directory, which holds txt files in gLike format.
    output_file: `string`
        path to the lklcom hdf5 file.
    mode: `string`
        mode to open the lklcom hdf5 file.
        Valid (recommended) option:
            "w": Write; a new file is created (an existing file
                with the same name would be deleted).
            "a": Append; an existing file is opened for reading and
                writing, and if the file does not exist it is created.
            "r+": It is similar to "a", but the file must already exist.
        Default: "w"
    Returns
    -------
    
    """
    
    # Deleting the output file, when the mode "w" (write) is selected.
    # Overwriting the mode to "a" (append), since the lklcom hdf5 file is opened
    # inside the function write_to_lklcom().
    if mode == "w":
        mode = "a"
        if os.path.isfile(output_file):
            os.remove(output_file)

    # Getting the txt files of the input directory.
    files = np.array([x for x in os.listdir(input_dir) if x.endswith(".txt")])
    # Looping over the files and store the likelihood or ts tables into the lklcom hdf5 file.
    for counter, file in enumerate(files):
    
        # Parsing the file name.
        file_info = file.replace('.txt','').split("_")
        # Getting the number of the simulation.
        simulation=-1
        if len(file_info) == 4:
            simulation=file_info[3]
            
        # Opening the txt file.
        txt_file = open("{}/{}".format(input_dir, file), "r")
        
        # Going through the table in the txt file and storing the entries in a 2D array.
        table = np.array([[i for i in line.split()] for line in txt_file]).T
        # Storing the log J-Factor.
        logJ = np.float32(table[0][0])
        # Storing the inverted sigmav range.
        sigmav_range = np.array(table[0][1:], dtype=np.float32)[::-1]
        # The first entry of each row correponds to the mass (or log J-Factor).
        # Detect the first entry and store it in a separate array.
        masses = np.array([table[i][0] for i in np.arange(0,table.shape[0],1)], dtype=np.float32)
        # Constructing the likelihood or ts dictionary with the DM mass `string` as keys
        # and the inverted likelihood or ts values `numpy.ndarray of type numpy.float32` as values.
        lkl_dict = {}
        for mass, ts_values in zip(masses[1:], table[1:]):
            lkl_dict[mass] = np.array(ts_values[1:], dtype=np.float32)[::-1]
            
        # Closing the txt files.
        txt_file.close()
        
        # Writing the likelihood table into the lklcom hdf5 file.
        write_to_lklcom(collaboration=file_info[2],
                        source=file_info[1],
                        channel=file_info[0],
                        logJ=logJ,
                        sigmav_range=sigmav_range,
                        lkl_dict=lkl_dict,
                        output_file=output_file,
                        mode=mode,
                        simulation=simulation)
    return


def lklcom_to_gLike(input_file,
                    output_dir,
                    reduce=True):
    """
    Translate the lklcom hdf5 file into gLike txt files.
    Parameters
    ----------
    input_file: `string`
        path to the lklcom hdf5 input file.
    output_dir: `string`
        path to the output directory.
    reduce: `boolean`
        flag, if the txt files should be reduced/compressed.
        Default: True
    Returns
    -------

    """

    # Opening the hdf5 file.
    h5 = tables.open_file(input_file, "r")
    
    # Looping over the likelihood or ts tables and store them into the gLike txt files.
    for h5_groups in h5.walk_groups("/"):
        for table in h5.list_nodes(h5_groups, classname='Table'):
            
            # Constructing the filename from the lklcom hdf5 file.
            table_info = table._v_pathname.split("/")
            filename = "{}_{}_{}".format(table_info[3], table_info[2], table_info[1])
            filename += ".txt" if table_info[4] == "data" else "_{}.txt".format(table_info[4].split("_")[-1])
            # Opening the gLike txt file to dump the likelihood or ts table.
            gLike_file = open("{}/{}".format(output_dir, filename), "w+")

            # Writing the first line, which contains the J-Factor and the different masses
            # Detecting the kinematic limit to skip the masses below that limit
            kinematic_limit = 0
            masses = np.array(table.cols._f_col('masses'))
            ts_values = np.array(table.cols._f_col('ts_values')).T
                    
            for i, mass in enumerate(masses):

                if i == 0:
                    # log J-Factor
                    gLike_file.write("{:.2f} ".format(mass))
                else:
                    # Detecting which channel correspond to which file
                    # Channel have to be in the filename!
                    if 'WW' == table_info[3] and mass < 80.39:
                        kinematic_limit = i
                        continue
                    if 'ZZ' == table_info[3] and mass < 91.19:
                        kinematic_limit = i
                        continue
                    if 'tt' == table_info[3] and mass < 173.1:
                        kinematic_limit = i
                        continue
                    # the bb limit, because Cirelli et al. doesn't provide 5 GeV
                    if 'bb' == table_info[3] and mass < 5.5:
                        kinematic_limit = i
                        continue
                    gLike_file.write("{:.0f} ".format(mass))
            gLike_file.write("\n")

            # Reducing the TS values (as agreed within the GloryDuck project).
            for line in ts_values[::-1]:
                for i, value in enumerate(line):
                    if reduce:
                        # Fist element of each line is the <sigma v> value
                        if i == 0:
                            gLike_file.write("{:.3e} ".format(value))
                        else:
                            # Skipping the masses below the kinematic limit
                            if i <= kinematic_limit:
                                continue
                            # Different strength of reducing depending on ROI
                            if value > 100:
                                gLike_file.write("{:.0f} ".format(value))
                            elif value > 1:
                                gLike_file.write("{:.3f} ".format(value))
                            elif value > 0.01:
                                gLike_file.write("{:.3e} ".format(value))
                            elif value > 1e-5:
                                gLike_file.write("{:.2e} ".format(value))
                            else:
                                gLike_file.write("0 ")
                    else:
                        gLike_file.write("{} ".format(value))

                gLike_file.write("\n")
            # Closing gLike txt file.
            gLike_file.close()
        
    # Closing hdf5 file.
    h5.close()
    return
