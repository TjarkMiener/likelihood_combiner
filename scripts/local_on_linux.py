import argparse
import yaml

import likelihood_combiner as lklcom

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description=("Combining likelihoods from different experiments."))
    parser.add_argument(
            'config_file',
            help="path to YAML configuration file with combining options")
    parser.add_argument(
            '--input',
            help="path to input file or directory")
    parser.add_argument(
            '--output',
            help="path to output file or directory")

    args = parser.parse_args()
    
    with open(args.config_file, 'r') as config_file:
        config = yaml.safe_load(config_file)
    
    lklcom.local.run_local_on_linux(settings=config, input=args.input, output=args.output)
