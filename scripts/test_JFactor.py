import likelihood_combiner as lklcom

lklcom_dir = "/Users/tmiener/deeplearning/likelihood_combiner/"
resource_GS = lklcom_dir + "/resources/GeringerSameth/intJ_cf.txt"
resource_BV = lklcom_dir + "/resources/Bonnivard/"
lklcom_file = lklcom_dir + "/resources/test_io.hdf5"
output_dir = lklcom_dir+ "/test/"

geringer_jfactor = lklcom.jfactor.GeringerSameth(
                                        sources=["Segue1", "Draco"],
                                        collaborations={"IACT": 0.12, "MAGIC": 0.1, "HAWC":2.6})
print(geringer_jfactor)
print(geringer_jfactor.get_logJ())
print(geringer_jfactor.get_DlogJ())
print(geringer_jfactor.get_DlogJ_comb())

geringer_jfactor = lklcom.jfactor.GeringerSameth(resource=resource_GS,
                                        sources=["Segue1", "Draco"],
                                        collaborations={"IACT": 2.6},
                                        combination_data=lklcom_file)
print(geringer_jfactor)
print(geringer_jfactor.get_logJ())
print(geringer_jfactor.get_DlogJ())
print(geringer_jfactor.get_DlogJ_comb())

lklcom_reader = lklcom.reader.LklCom_hdf5(channel="bb",
                                        LklCom_jfactor_class=geringer_jfactor)

tstables = lklcom_reader.read_tstables()
print(tstables)


bonnivard_jfactor = lklcom.jfactor.Bonnivard(
                                        sources=["UrsaMajorII"],
                                        collaborations={"IACT": 0.1},
                                        combination_data=output_dir,
                                        precision=4)
print(bonnivard_jfactor)
print(bonnivard_jfactor.get_logJ())
print(bonnivard_jfactor.get_DlogJ())
print(bonnivard_jfactor.get_DlogJ_comb())

lklcom_reader = lklcom.reader.LklCom_txtdir(channel="bb",
                                            LklCom_jfactor_class=bonnivard_jfactor)

tstables = lklcom_reader.read_tstables()
print(tstables)


custom_logJ = {"Segue1": {"IACT": 19.25}}
custom_DlogJ = {"Segue1": {"IACT": 0.25}}
custom_jfactor = lklcom.jfactor.Custom(logJ=custom_logJ,
                                       DlogJ=custom_DlogJ)
print(custom_jfactor)
print(custom_jfactor.get_logJ())
print(custom_jfactor.get_DlogJ())
print(custom_jfactor.get_DlogJ_comb())
