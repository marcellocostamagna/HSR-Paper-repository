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

# np.set_printoptions(precision=4, suppress=True)

cwd = os.getcwd()

molecule_1 = load_molecules_from_sdf(f'{cwd}/molecule_1.sdf', removeHs=False, sanitize=False)[0]

molecule_2 = load_molecules_from_sdf(f'{cwd}/molecule_2.sdf', removeHs=False, sanitize=False)[0]

molecule_3 = load_molecules_from_sdf(f'{cwd}/molecule_3.sdf', removeHs=False, sanitize=False)[0]

molecule_4 = load_molecules_from_sdf(f'{cwd}/molecule_4.sdf', removeHs=False, sanitize=False)[0]

molecule_5 = load_molecules_from_sdf(f'{cwd}/molecule_5.sdf', removeHs=False, sanitize=False)[0]    
   
print(f'HSR similarity scores: \n')
print(f'1 vs 2: {compute_similarity(molecule_1, molecule_2, features=DEFAULT_FEATURES)}\n')

print(f'1 vs 3: {compute_similarity(molecule_1, molecule_3, features=DEFAULT_FEATURES)}\n')

print(f'1 vs 4: {compute_similarity(molecule_1, molecule_4, features=DEFAULT_FEATURES)}\n')

print(f'1 vs 5: {compute_similarity(molecule_1, molecule_5, features=DEFAULT_FEATURES)}\n')


    