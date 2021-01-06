import likelihood_combiner as lklcom

lklcom_dir = "/Users/tmiener/deeplearning/likelihood_combiner/"

input_dir = lklcom_dir+"data/"
output_dir = lklcom_dir+ "test/"
input_file = output_file = lklcom_dir+"test_io.hdf5"

lklcom.io.gLike_to_lklcom(input_dir, output_file)
lklcom.io.lklcom_to_gLike(input_file, output_dir)

