##########################################################
## ccAFv2:  Simulate AMI for cell cycle classifiers     ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier, Samantha O'Connor          ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

#docker run -it -v '/home/soconnor/old_home/ccNN/ccAFv2:/files' cplaisier/ccnn

from sklearn.metrics import adjusted_mutual_info_score as ami
import random
import numpy as np
import copy
import pandas as pd
import os

savedir = 'compare_classifiers'
newpath = savedir+'/ami_simulation'
if not os.path.exists(newpath):
    os.makedirs(newpath)

############
## Guides ##
############
guides = {}
guides['ccSeurat_cyclone'] = {1:[1,1], 2:[2,2], 3:[3,3], 'k':3}
guides['peco_tricycle'] = {1:[1,1], 2:[2,2], 3:[3,4], 'k':4}
guides['SchwabeCC'] = {1:[1,1], 2:[2,2], 3:[3,5], 'k':5}
guides['reCAT'] = {1:[1,2], 2:[3,3], 3:[4,6], 'k':6}
guides['ccAFv2'] = {1:[1,3], 2:[4,4], 3:[5,7], 'k':7}
guides['ccAF'] = {1:[1,4], 2:[5,5], 3:[6,8], 'k':8}


###############
## Functions ##
###############

# Randomly sample a list of 'c' number cell cycle classifications for 'k' number of clusters
def sample_ccc(c, k_ref):
    return [random.randint(1,k_ref) for i in range(c)]

# Use sampled cell cycle classifications as a guide to simulate larger number of clusers 'k_test'
# Note: We supply guides to constrain the space!!! This could be done by adding a parameter which takes a dictionary as input.
def expand_ccc(ccc, k_ref, k_test, guide):
    return [random.randrange(guide[i][0], guide[i][1]+1) for i in ccc]

# Randomize a proportion of a cell cycle classifier list
def shuffle_ccc_prop(ccc, prop):
    if prop==0:
        return ccc
    else:
        rand_ccc = sample_ccc(len(ccc), max(ccc))
        # Sample out indices to swap
        sampled = random.sample(range(len(ccc)), int(len(ccc)*prop))
        ccc_shuffled = copy.deepcopy(ccc)
        for i in sampled:
            ccc_shuffled[i] = rand_ccc[i]
        return ccc_shuffled


## Example of how use shuffling of simulated cell cycle classifications
ccc_k3 = sample_ccc(1200, 3)
ccc_ex7 = expand_ccc(ccc_k3, 3, 7, guides['ccAFv2'])
ccc_k7 = sample_ccc(1200, 7)
print('Random: ',ami(ccc_k3, ccc_k7)) # Worst case
print('Expanded: ',ami(ccc_k3, ccc_ex7)) # Best case
'''for s1 in [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1]:
    print('Expanded Shuffled(',s1,'): ',np.mean(shuffled), np.std(shuffled))
    print('Expanded Shuffled(',s1,'): ',ami(ccc_k3, shuffle_ccc_prop(ccc_ex7, s1)))

# Example of how use shuffling of simulated cell cycle classifications to get range with mean and standard deviation
print('Random: ',ami(ccc_k3, ccc_k7)) # Worst case
print('Expanded: ',ami(ccc_k3, ccc_ex7)) # Best case
for s1 in [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1]:
    shuffled = [ami(ccc_k3, shuffle_ccc_prop(ccc_ex7, s1)) for i in range(100)]
    print('Expanded Shuffled(',s1,'): ',np.mean(shuffled), np.std(shuffled))

# Upping the ante, resimulate ccc_ex7 each time. (Very similar answer to 57).
for s1 in [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1]:
    shuffled = [ami(ccc_k3, shuffle_ccc_prop(expand_ccc(ccc_k3, max(ccc_k3), 7, guide_ccAFv2), s1)) for i in range(100)]
    print('Expanded Shuffled(',s1,'): ',np.mean(shuffled), np.std(shuffled))
'''

### Use this!!!
# New ccc_ref each iteration, new ccc_expanded each iteration, new_random shuffling of a proportion. This is the one.
c = 2962
k_ref = 3
results_mean = {}
results_std = {}
for guide_curr in ['ccSeurat_cyclone', 'peco_tricycle', 'SchwabeCC', 'reCAT', 'ccAFv2', 'ccAF']:
    results_mean[guide_curr] = {}
    results_std[guide_curr] = {}
    for s1 in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1]:
        shuffled = []
        for i in range(100):
            ccc_ref = sample_ccc(c, k_ref)
            shuffled.append(ami(ccc_ref, shuffle_ccc_prop(expand_ccc(ccc_ref, max(ccc_ref), guides[guide_curr]['k'], guides[guide_curr]), s1)))
        results_mean[guide_curr][s1] = np.mean(shuffled)
        results_std[guide_curr][s1] = np.std(shuffled)


df_results_mean = pd.DataFrame(results_mean)
df_results_std = pd.DataFrame(results_std)
df_results_mean.to_csv(newpath+'/results_mean.csv')
df_results_mean.to_csv(newpath+'/results_stddev.csv')
