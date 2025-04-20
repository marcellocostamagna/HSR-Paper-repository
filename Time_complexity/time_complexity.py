import numpy as np  
import os 
import sys
from hsr.pre_processing import *
from hsr.pca_transform import * 
from hsr.fingerprint import *
from hsr.similarity import *
from hsr.utils import *
sys.path.append(os.path.abspath('../')) 
from perturbations import *
import time

# np.set_printoptions(precision=4, suppress=True)

PROTON_NEUTRON_FEATURES = {
    'protons' : extract_proton_number,
    'delta_neutrons' : extract_neutron_difference_from_common_isotope,
}

cwd = os.getcwd()

molecule = load_molecules_from_sdf(f'{cwd}/molecule.sdf', removeHs=False, sanitize=False)[0]

num_experiments = 1000
timings = {'3D': [], '4D': [], '5D': [], '6D': []}
features = {'3D': None, '4D': PROTON_FEATURES, '5D': PROTON_NEUTRON_FEATURES, '6D': DEFAULT_FEATURES}

# Perform the experiment 1000 times
for _ in range(num_experiments):
    temp_timings = {'3D': [], '4D': [], '5D': [], '6D': []}
    
    for feature in features:
        start = time.time()
        fingerprint = generate_fingerprint_from_molecule(molecule, features[feature])
        end = time.time()
        temp_timings[feature].append(end - start)
    
    # Store results
    for key in timings:
        timings[key].append(np.mean(temp_timings[key]))

# Compute final averages
final_timings = {key: np.mean(timings[key]) for key in timings}

print(f'Final averaged results over {num_experiments} experiments:')
for feature in ['3D', '4D', '5D', '6D']:
    print(f'\n{feature} timings:')
    print(f'Average time: {round(final_timings[feature], 5)}')
    if feature != '3D':
        print(f'Increase: {round((final_timings[feature]/final_timings["3D"] - 1) * 100, 2)}%')
