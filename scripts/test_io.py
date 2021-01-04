from likelihood_combiner.io import *

lklcom_dir = "/Users/tmiener/deeplearning/likelihood_combiner/"

input_dir = lklcom_dir+"data/"
output_dir = lklcom_dir+ "test/"
input_file = output_file = lklcom_dir+"test_io.hdf5"

gLike2gloryduck(input_dir, output_file, mode="w")
gloryduck2gLike(input_file, output_dir)

