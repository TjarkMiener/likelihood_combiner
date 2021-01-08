import likelihood_combiner as lklcom

lklcom_dir = "/Users/tmiener/deeplearning/likelihood_combiner/resources/"

input_dir = lklcom_dir+"mock_txt_data/"
output_dir = lklcom_dir+ "test/"
input_file = output_file = lklcom_dir+"mock_data.hdf5"

lklcom.io.gLike_to_lklcom(input_dir, output_file)
#lklcom.io.lklcom_to_gLike(input_file, output_dir)

