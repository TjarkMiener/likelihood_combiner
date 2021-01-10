import argparse
import yaml
import os

import likelihood_combiner as lklcom
from likelihood_combiner.combiner import combiner
from likelihood_combiner.utils import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description=("Combining likelihoods from different experiments."))
    parser.add_argument(
            'config_file',
            help="path to YAML configuration file with combining options")
    parser.add_argument(
            '--channel',
            default='bb')
    parser.add_argument(
            '--simulation',
            default=0,
            type=int,
            help="number of the simulation")
    args = parser.parse_args()

    with open(args.config_file, 'r') as config_file:
        config = yaml.safe_load(config_file)
    
    try:
        input = config['Data']['input']
        if input is None:
            raise KeyError
    except KeyError:
        input = os.path.abspath(os.path.join(os.path.dirname(__file__), "../input/mock_data.hdf5"))

    try:
        output_dir = config['Output']['directory']
        if output_dir is None:
            raise KeyError
    except KeyError:
        output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output/"))
    # Create output directory if it doesn't exist already
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Initializing of the LklCom jfactor class
    if config['Data']['buildin_j_factors'] == "GeringerSameth":
        LklCom_jfactor_class = lklcom.jfactor.GeringerSameth(sources=config['Configuration']['sources'],
                                                    collaborations=config['Configuration']['collaborations'],
                                                    combination_data=input,
                                                    jnuisance=config['Data']['j_nuisance'])
    elif config['Data']['buildin_j_factors'] == "Bonnivard":
        LklCom_jfactor_class = lklcom.jfactor.Bonnivard(sources=config['Configuration']['sources'],
                                                    collaborations=config['Configuration']['collaborations'],
                                                    combination_data=input,
                                                    jnuisance=config['Data']['j_nuisance'])
    else:
        LklCom_jfactor_class = lklcom.jfactor.Custom(logJ=config['Data']['custom_logJ'],
                                                    DlogJ=config['Data']['custom_DlogJ'],
                                                    jnuisance=config['Data']['j_nuisance'])

    # Initializing of the LklCom reader class
    if input.endswith(".h5") or input.endswith(".hdf5"):
        LklCom_reader_class = lklcom.reader.LklCom_hdf5(channel=args.channel,
                                                        LklCom_jfactor_class=LklCom_jfactor_class)
    if os.path.isdir(input):
        LklCom_reader_class = lklcom.reader.LklCom_txtdir(channel=args.channel,
                                                        LklCom_jfactor_class=LklCom_jfactor_class)

    # Constructing in the the sigmav range and spacing
    sigmav_range = get_sigmav_range()
    # Overwriting from the config file, if provided.
    if "sigmav_min" in config['Data'] and "sigmav_max" in config['Data'] and "sigmav_nPoints" in config['Data']:
        sigmav_range = get_sigmav_range(config['Data']['sigmav_min'], config['Data']['sigmav_max'], config['Data']['sigmav_nPoints'], 3)

    combiner(sigmav_range=sigmav_range,
            LklCom_reader_class=LklCom_reader_class,
            output=output_dir)
