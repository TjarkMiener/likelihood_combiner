import tables
import numpy as np
import os
    

def gLike2gloryduck(input_dir, output_file, mode="a"):

    # Opening the hdf5 file.
    h5 = tables.open_file(output_file, mode=mode, title="LklCom database")
    print("The files are written in '{}' ({}):".format(h5.title, output_file))
    # Getting the txt files of the input directory.
    files = np.array([x for x in os.listdir(input_dir) if x.endswith(".txt")])
    
    counter = 1
    for file in files:
        # Parsing the file name.
        file_info = file.replace('.txt','').split("_")
        
        # Creating a new node in the hdf5 file, if it isn't already existing.
        if "/{}".format(file_info[2]) not in h5:
            h5.create_group(h5.root,file_info[2],"Further information about the collaboration {}".format(file_info[2]))
        if "/{}/{}".format(file_info[2],file_info[1]) not in h5:
             h5.create_group(eval("h5.root.{}".format(file_info[2])), file_info[1], "Further information about the source {}".format(file_info[1]))
        if "/{}/{}/{}".format(file_info[2],file_info[1],file_info[0]) not in h5:
            h5.create_group(eval("h5.root.{}.{}".format(file_info[2],file_info[1])), file_info[0], "Further information about the annihilation channel {}".format(file_info[0]))

        # Opening the txt file.
        ts_file = open("{}/{}".format(input_dir, file), "r")
        print("    {}) '{}'".format(counter, file))
        counter += 1
        
        # Going through the table in the txt file and storing the entries in a 2D array.
        ts_val = np.array([[i for i in line.split()] for line in ts_file]).T
            
        # The first entry of each row correponds to the mass.
        # Detect the first entry and store it in a separate array.
        mass = np.array([ts_val[i][0] for i in np.arange(0,ts_val.shape[0],1)])
    
        # Closing the txt files.
        ts_file.close()
        
        # Creating the table structure for the hdf5 file.
        sigmav_shape = (ts_val.shape[1]-1,)
        columns_dict={"masses":tables.Float32Col(),
                      "ts_values":tables.Float32Col(shape=sigmav_shape)}
        description = type("description", (tables.IsDescription,), columns_dict)
        
        # Creating the table mass vs sigmav for each source and for each channel.
        table_name = "data"
        if len(file_info) == 4:
            table_name = "simu_{}".format(file_info[3])
        table = h5.create_table(eval("h5.root.{}.{}.{}".format(file_info[2],file_info[1],file_info[0])),table_name,description,"Table of the {} collaboration for the source {} with the annihilation channel {}.".format(file_info[2],file_info[1],file_info[0]))
        
        # Filling the data of the txt file into the table of the hdf5 file.
        for i, mass_val in enumerate(mass):
            table = eval("h5.root.{}.{}.{}.{}".format(file_info[2],file_info[1],file_info[0],table_name))
            row = table.row
            row['masses']  = mass_val
            # The first element of the ts array (mass value) will be ignored, since it's in the mass column.
            row['ts_values'] = ts_val[i][1:]
            row.append()
            table.flush()

    # Closing hdf5 file.
    h5.close()
    return


def gloryduck2gLike(input_file, output_dir, reduce=True):

    # Opening the hdf5 file.
    h5 = tables.open_file(input_file, "r")
    
    print("The files are written in '{}':".format(output_dir))
    counter = 1

    for h5_groups in h5.walk_groups("/"):
        for table in h5.list_nodes(h5_groups, classname='Table'):
            
            table_info = table._v_pathname.split("/")
            filename = "{}_{}_{}".format(table_info[3], table_info[2], table_info[1])
            filename += ".txt" if table_info[4] == "data" else "_{}.txt".format(table_info[4].split("_")[-1])
            gLike_file = open("{}/{}".format(output_dir, filename), "w+")
            print("    {}) '{}'".format(counter, filename))
            counter += 1
            # Writing the first line, which contains the J-Factor and the different masses
            # Detecting the kinematic limit to skip the masses below that limit
            kinematic_limit = 0
            masses = np.array(table.cols._f_col('masses'))
            ts_values = np.array(table.cols._f_col('ts_values')).T
                    
            for i, mass in enumerate(masses):

                if i == 0:
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

            # Reducing the TS values as discussed in the GD call
            for line in ts_values:
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
                                gLike_file.write("0  ")
                    else:
                        gLike_file.write("{} ".format(value))

                gLike_file.write("\n")
            # Closing gLike txt file.
            gLike_file.close()
        
    # Closing hdf5 file.
    h5.close()
    return
