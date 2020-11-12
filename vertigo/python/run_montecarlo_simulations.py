#!/usr/bin/env python2

import subprocess
import rospkg

rospack = rospkg.RosPack()
package_folder = rospack.get_path('vertigo_ros')


# list dataset
dataset_folder = package_folder + '/vertigo/datasets/'

# all the datasets should be listed here
datasets_to_process = {'intel' :    {'file': 'intel/originalDataset/intel.g2o', 
                                     'inliers': 100}, # this must be changed
                       'manhattan': {'file': 'manhattan/originalDataset/Olson/manhattanOlson3500.g2o',
                                     'inliers': 100} # this must be changed
                    # here we add more datasets in the same way
                      }

# outlier percentage
outlier_percentages = [10, 20, 30, 40, 50]

# number of monte carlo runs
N = 10


# Iterate over all the datasets
for dataset_name in datasets_to_process:
    print('Running dataset [' + dataset_name) + ']'

    dataset = datasets_to_process[dataset_name]

    # Iterate ovr all the percentages
    for percentage in outlier_percentages:
        print('  Experiments with ' + str(percentage) + '% of outliers')

        # Compute number of outliers for the dataset
        num_inliers = dataset['inliers']
        
        # Compute total number of inliers+outliers
        num_total = int(num_inliers*(1 + percentage/100.0))
        
        # Compute number of outliers to add
        num_outliers = num_total - num_inliers

        # Monte carlo runs
        for n in range(N):
            print('      Running Monte Carlo simulation: ' + str(n))

            # Create name of output file
            dataset_file = dataset_name + '_dataset_' + 'outliers' + str(percentage) + 'perc' + '_run' + str(n) + '.g2o'



            # Create new dataset with num_outliers
            # Prepare command to create datasets
            cmd = 'rosrun vertigo_ros generateDataset.py' + \
                                ' -i ' +  dataset_folder + dataset['file'] + \
                                ' -s ' + \
                                ' -n ' + str(num_outliers) + \
                                ' -o ' + dataset_file
            #print(cmd)
            # Run the command
            subprocess.call(cmd, shell=True)
            
            # Adaptive algorithm
            # Prepare output result
            results_adaptive_file = dataset_name + '_results_adaptive_' + 'run' + str(n) + '_outliers' + str(percentage) + 'perc' + '.csv'

            # Run adaptive algorithm
            cmd = 'rosrun vertigo_ros robustISAM2-2d' + \
                                ' -i ' + dataset_file + \
                                ' --adaptive ' + \
                                ' -o ' + results_adaptive_file
            #print(cmd)
            # Run adaptive algorithm
            subprocess.call(cmd, shell=True)
            

            
