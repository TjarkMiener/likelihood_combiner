#!/usr/bin/python'

import numpy as np
import os
from likelihood_combiner.reader import gloryduckReader
from likelihood_combiner.writer import gloryduckWriter

hdf5file = os.path.join(os.path.dirname(__file__), "../data/gloryduck_dataset.h5")
path2txts = os.path.join(os.path.dirname(__file__), "../data/")

writer = gloryduckWriter()

writer.convert_txts2hdf5(hdf5file,path2txts)


del writer

