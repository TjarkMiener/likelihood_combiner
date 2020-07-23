import argparse
import yaml

from likelihood_combiner.combiner import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description=("Combining likelihoods from different experiments."))
    parser.add_argument(
            'config_file',
            help="path to YAML configuration file with combining options")
    parser.add_argument(
            '--channel')
    parser.add_argument(
            '--simulation',
            default=-1,
            type=int,
            help="number of the simulation")
    args = parser.parse_args()

    with open(args.config_file, 'r') as config_file:
        config = yaml.safe_load(config_file)

    combiner(config=config,
             channel=args.channel,
             sigmavULs=None,
             sigmavULs_Jnuisance=None,
             simulation_counter=None,
             simulations=[args.simulation])

