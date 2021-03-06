"""
io.py
=====
Functions to translate between data formats
"""

import argparse
import tables
import tarfile
import numpy as np
import os
import pandas as pd

__all__ = [
    'write_to_lklcom',
    'gLike_to_lklcom',
    'lklcom_to_gLike',
    'gLikeLimits_to_lklcomLimits',
    'merge_to_lklcom',
    'gloryduck_reducer'
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
    Write/append a given likelihood table into lklcom hdf5 file. lkl_dict is a dictionary
    with the DM mass as keys (`str`) and likelihood or ts values (ascending) as values
    (`numpy.ndarray of type numpy.float32`).

    Parameters
    ----------
    collaboration: str
        name of the collaboration.
    source: str
        name of the source.
    channel: str
        name of the channel.
    logJ: `numpy.float32`
        value of the log J-Factor.
    sigmav_range: `numpy.ndarray of type numpy.float32`
        sigmav range (ascending).
    lkl_dict: dict
        likelihood data as dictionary with the DM mass as keys (str) and likelihood or ts values (ascending) as values (`numpy.ndarray of type numpy.float32`).
    output_file: path
        path to the lklcom hdf5 file.
    mode: str
        mode to open the lklcom hdf5 file.
    simulation: int
        number of the simulation.
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


def gLike_to_lklcom(input_dir,
                    output_file,
                    mode="w"):
    """
    Translate gLike txt files into lklcom hdf5 file.

    Parameters
    ----------
    input_dir: path
        path to the input directory, which holds txt files in gLike format.
    output_file: path
        path to the lklcom hdf5 file.
    mode: str
        mode to open the lklcom hdf5 file.
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
        file_info[2] = file_info[2].replace('-','')
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


def _gLike_to_lklcom():
    """
    Interface with gLike_to_lklcom() from the command line. 
    """

    parser = argparse.ArgumentParser(
            description=("Translate gLike txt files into lklcom hdf5 file."))
    parser.add_argument(
            '--input', '-i',
            help="path to input directory")
    parser.add_argument(
            '--output', '-o',
            help="path to output file")
    parser.add_argument(
            '--mode', '-m',
            default="w",
            help="mode to open the lklcom hdf5 file")

    args = parser.parse_args()

    gLike_to_lklcom(input_dir=args.input, output_file=args.output, mode=args.mode)


def lklcom_to_gLike(input_file,
                    output_dir,
                    reduce=True):
    """
    Translate the lklcom hdf5 file into gLike txt files.

    Parameters
    ----------
    input_file: path
        path to the lklcom hdf5 input file.
    output_dir: path
        path to the output directory.
    reduce: bool
        flag, if the txt files should be reduced/compressed.
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
                if reduce:
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
                else:
                    gLike_file.write("{} ".format(mass))
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
                            elif value < -1:
                                gLike_file.write("{:.3f} ".format(value))
                            elif value < -0.01:
                                gLike_file.write("{:.3e} ".format(value))
                            elif value < -1e-5:
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

def _lklcom_to_gLike():
    """
    Interface with lklcom_to_gLike() from the command line. 
    """

    parser = argparse.ArgumentParser(
            description=("Translate the lklcom hdf5 file into gLike txt files."))
    parser.add_argument(
            '--input', '-i',
            help="path to input file")
    parser.add_argument(
            '--output', '-o',
            help="path to output directory")
    parser.add_argument(
            '--reduce', '-r',
            default=True,
            action=argparse.BooleanOptionalAction,
            help="flag, if the txt files should be reduced/compressed")

    args = parser.parse_args()

    lklcom_to_gLike(input_file=args.input, output_dir=args.output, reduce=args.reduce)


def gLikeLimits_to_lklcomLimits(input_dir,
                                output_file):
    """
    Translate gLike limits in txt files into lklcom results hdf5 file.

    Parameters
    ----------
    input_dir: path
        path to the input directory, which holds txt files of the results of gLike or any other framework.
    output_file: path
        path to the lklcom results hdf5 file.
    """

    # Getting the txt files of the input directory.
    files = np.array([x for x in os.listdir(input_dir) if x.endswith(".txt")])
    channels = np.unique([file.split("_")[0] for file in files])
    # Looping over the files and store the likelihood or ts tables into the lklcom hdf5 file.
    svUL = {}
    for channel in channels:
        svUL[channel] = {}
    for file in files:
        # Parsing the file name.
        file_info = file.replace('.txt','').split("_")
        # Getting the number of the simulation.
        simulation=-1
        if file_info[-1].isnumeric():
            simulation=file_info[-1]

        # Opening the txt file.
        txt_file = open("{}/{}".format(input_dir, file), "r")

        # Going through the table in the txt file and storing the entries in a 2D array.
        table = np.array([[i for i in line.split()] for line in txt_file], dtype=np.float32)

        # Dumping the upper limits

        if simulation != -1:
            col_name = "simu_{}".format(simulation)
        else:
            col_name = "data"
            svUL[file_info[0]]['masses'] = table[0]
        svUL[file_info[0]][col_name] = table[1]

    for channel in channels:
         pd.DataFrame(data=svUL[channel]).to_hdf(output_file, key='{}/{}'.format(channel, file_info[1]), mode="a")

def _gLikeLimits_to_lklcomLimits():
    """
    Interface with gLikeLimits_to_lklcomLimits() from the command line. 
    """

    parser = argparse.ArgumentParser(
            description=("Translate gLike limits in txt files into lklcom results hdf5 file."))
    parser.add_argument(
            '--input', '-i',
            help="path to input directory")
    parser.add_argument(
            '--output', '-o',
            help="path to output directory")

    args = parser.parse_args()

    gLikeLimits_to_lklcomLimits(input_dir=args.input, output_file=args.output)


def merge_to_lklcom(input_dir,
                    output_file):
    """
    Merge single lklcom hdf5 file produced by the cluster to the lklcom results hdf5 file.

    Parameters
    ----------
    input_dir: path
        path to the input directory, which holds txt files of the results of gLike or any other framework.
    output_file: path
        path to the lklcom results hdf5 file.
    """

    # Getting the h5 files of the input directory.
    files = np.array([x for x in os.listdir(input_dir) if x.endswith(".h5") or x.endswith(".hdf5")])
    channels = np.unique([file.split("_")[0] for file in files])

    j_nuisance = False
    svUL, svUL_Jnuisance = {}, {}
    for channel in channels:
        svUL[channel] = {}
        svUL_Jnuisance[channel] = {}
    for file in files:
        # Parsing the file name.
        file_info = file.replace('.hdf5','').replace('.h5','').split("_")

        # Opening the hdf5 file and reading the data
        data = pd.HDFStore("{}/{}".format(input_dir, file), 'r')

        if '/masses' in data.keys():
            svUL[file_info[0]]['masses'] = svUL_Jnuisance[file_info[0]]['masses'] = data['masses'][0]
        if '/sigmavULs' in data.keys():
            svUL[file_info[0]][file_info[1]] = data['sigmavULs'][0]
        if '/sigmavULs_Jnuisance' in data.keys():
            j_nuisance = True
            svUL_Jnuisance[file_info[0]][file_info[1]] = data['sigmavULs_Jnuisance'][0]

    # Write the panda DataFrames into the hdf5 file
    for channel in channels:
        pd.DataFrame(data=svUL[channel]).to_hdf(output_file, key='{}/sigmavULs'.format(channel), mode='a')
        if j_nuisance:
            pd.DataFrame(data=svUL_Jnuisance[channel]).to_hdf(output_file, key='{}/sigmavULs_Jnuisance'.format(channel), mode='a')
       
def _merge_to_lklcom():
    """
    Interface with merge_to_lklcom() from the command line. 
    """

    parser = argparse.ArgumentParser(
            description=("Merge single lklcom hdf5 file produced by the cluster to the lklcom results hdf5 file."))
    parser.add_argument(
            '--input', '-i',
            help="path to input file or directory")
    parser.add_argument(
            '--output', '-o',
            help="path to output directory")

    args = parser.parse_args()

    merge_to_lklcom(input_dir=args.input, output_file=args.output)

def gloryduck_reducer(input_dir,
                     reduce=True,
                     shift=True):
    """
    Reduce and/or shift the GloryDuck txt files.
    Important note: The channel have to be in the filename!

    Parameters
    ----------
    input_dir: path
        path to input directory.
    reduce: bool
        flag, if the txt files should be reduced/compressed.
    shift: bool
        flag, if the txt files should be shifted (TS values minimum to zero).
    """

    # Checking if the main directory exits
    if not os.path.exists(input_dir):
        raise ValueError("Error 404: Directory '{}' not found. ".format(input_dir))

    # Creating a new directory (for the reduce files) inside the main directory
    reduce_directory = "{}/reduced/".format(input_dir)

    if not os.path.exists(reduce_directory):
        os.makedirs(reduce_directory)

    # Collacting all txt files
    files = np.array([x for x in os.listdir(input_dir) if x.endswith(".txt")])

    # Opening the tar gz file
    tar_gz_filename = "reduced_files_GD.tar.gz"
    tar_gz_file = tarfile.open(tar_gz_filename, "w:gz")

    # Looping over the files to reduce them
    file_counter = 0.0
    nb_of_files = len(files)
    for f in files:
        # Showing progress bar
        file_counter += 1.0
        percentage = file_counter / nb_of_files * 100.0
        print(f + " " + "{:.1f}".format(percentage) + "%")

        # Opening the original txt files
        filename = "{}/{}".format(input_dir, f)
        file = open(filename, "r")

        # Creating and opening the reduced txt files
        reduce_filename = "{}/{}".format(reduce_directory, f)
        reduce_filename = reduce_filename.replace("//","/")
        reduce_file = open(reduce_filename, "w+")

        # Going through the table in the txt file and storing the entries in a 2D array.
        values = np.array([[i for i in line.split()] for line in file]).astype(np.float)
        masses = np.squeeze(values[:1])

        # Writing the first line, which contains the J-Factor and the different masses
        # Detecting the kinematic limit to skip the masses below that limit

        kinematic_limit = 0
        for i, mass in enumerate(masses):
            if reduce:
                if i == 0:
                    # log J-Factor
                    reduce_file.write("{:.2f} ".format(mass))
                else:
                    # Detecting which channel correspond to which file
                    # Channel have to be in the filename!
                    if 'WW' in f and mass < 80.39:
                        kinematic_limit = i
                        continue
                    if 'ZZ' in f and mass < 91.19:
                        kinematic_limit = i
                        continue
                    if 'tt' in f and mass < 173.1:
                        kinematic_limit = i
                        continue
                    # the bb limit, because Cirelli et al. doesn't provide 5 GeV
                    if 'bb' in f and mass < 5.5:
                        kinematic_limit = i
                        continue
                    reduce_file.write("{:.0f} ".format(mass))
            else:
                reduce_file.write("{} ".format(mass))
        reduce_file.write("\n")

        # Shifting the TS values minimum to zero
        vals = values[1:]
        if shift:
            vals = vals.T
            for i in np.arange(1,len(vals)):
                vals[i] -= np.min(vals[i])
            vals = vals.T

        # Reducing the TS values as discussed in the GD call
        for line in vals:
            for i, value in enumerate(line):
                if reduce:
                    # Fist element of each line is the <sigma v> value
                    if i == 0:
                        reduce_file.write("{:.3e} ".format(value))
                    else:
                        # Skipping the masses below the kinematic limit
                        if i <= kinematic_limit:
                            continue
                        # Different strength of reducing depending on ROI
                        if value > 100:
                            reduce_file.write("{:.0f} ".format(value))
                        elif value > 1:
                            reduce_file.write("{:.3f} ".format(value))
                        elif value > 0.01:
                            reduce_file.write("{:.3e} ".format(value))
                        elif value > 1e-5:
                            reduce_file.write("{:.2e} ".format(value))
                        elif value < -1:
                            reduce_file.write("{:.3f} ".format(value))
                        elif value < -0.01:
                            reduce_file.write("{:.3e} ".format(value))
                        elif value < -1e-5:
                            reduce_file.write("{:.2e} ".format(value))
                        else:
                            reduce_file.write("0 ")
                else:
                    reduce_file.write("{} ".format(value))

            reduce_file.write("\n")

        file.close()
        reduce_file.close()

        # Adding the final reduced file to the tar gz file
        tar_gz_file.add(reduce_filename)

    # Closing the tar gz file
    tar_gz_file.close()

    # Printing
    print("---------------------------------")
    print("Wrote the zipped tar file {}".format(tar_gz_filename))

def _gloryduck_reducer():
    """
    Interface with gloryduck_reducer() from the command line. 
    """

    parser = argparse.ArgumentParser(
            description=("Reduce the GloryDuck txt files."))
    parser.add_argument(
            '--input', '-i',
            help="path to input directory")
    parser.add_argument(
            '--reduce', '-r',
            default=True,
            action=argparse.BooleanOptionalAction,
            help="flag, if txt files should be reduced/compressed")
    parser.add_argument(
            '--shift', '-s',
            default=True,
            action=argparse.BooleanOptionalAction,
            help="flag, if txt files should be shifted; TS values minimum to zero")


    args = parser.parse_args()

    gloryduck_reducer(input_dir=args.input, reduce=args.reduce, shift=args.shift)

