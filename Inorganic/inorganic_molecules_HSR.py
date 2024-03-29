import numpy as np  
import sys
import os 
from hsr.pre_processing import *
from hsr.pca_transform import * 
from hsr.fingerprint import *
from hsr.similarity import *
from hsr.utils import *
sys.path.append(os.path.abspath('../'))
from perturbations import *

cwd = os.getcwd()

print(f'\nHSR similarity of inorganic compounds: \n')
sorted_files = sorted(os.listdir(f'{cwd}/molecules'), key=lambda x: int(x.split('-')[0]))
for file in sorted_files:
    if file.endswith('.sdf'):
        molecules = load_molecules_from_sdf(f'{cwd}/molecules/{file}', removeHs=False, sanitize=False)
        
    fingerprints = [generate_fingerprint_from_molecule(molecule, DEFAULT_FEATURES, scaling='matrix', chirality=False) for molecule in molecules]

    # COMPARE MOLECULES
    # Compute similarity between all pairs of fingerprints
    n_molecules = len(fingerprints)
    for i in range(n_molecules):
        for j in range(i+1, n_molecules):
            similarity = compute_similarity_score(fingerprints[i], fingerprints[j])
            # print file name and the similarity score
            print(f'{file[:-4]}: {similarity:.4f} \n')