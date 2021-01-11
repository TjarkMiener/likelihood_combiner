import likelihood_combiner as lklcom


settings = {
        	'Hardware': {'cpu_counts': 8},
                'Data': {'buildin_j_factors': 'GeringerSameth', 'j_nuisance': True, 'simulations': 30},
                'Configuration': {'sources': ['Segue1', 'UrsaMajorII'], 'collaborations': {'IACT': 2.6}}
                #'Configuration': {'channels': ['bb'], 'sources': ['Segue1', 'UrsaMajorII'], 'collaborations': {'IACT': 2.6}}
           }
input = "/Users/tmiener/deeplearning/likelihood_combiner/input/mock_data.hdf5"
#output = "/Users/tmiener/deeplearning/likelihood_combiner/output/output3.hdf5"
#lklcom.local.run_local_on_linux(settings, input, output)

output = "/Users/tmiener/deeplearning/likelihood_combiner/output/"
lklcom.cluster.run_cluster(settings, "bb", input, output, simulation=25)
