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

init_dir = 'starting_molecules'

# working directory
dir = 'enantiomeric_pairs'

# Loop through each XYZ file
sorted_files = sorted(os.listdir(init_dir), key=lambda x: int(x.split('-')[0]))

# For each file, reflect the molecule and save both together
for file in sorted_files:

    # Read the sdf content
    molecule = load_molecules_from_sdf(f'{init_dir}/{file}', removeHs=False, sanitize=False)[0]
    
    with Chem.SDWriter(f'{cwd}/{dir}/{file}') as writer:
        writer.write(molecule)
        
        # Reflect the molecule
        reflected_molecule = reflect_molecule_coordinate(molecule, coordinate='x')
        writer.write(reflected_molecule)
    writer.close()
    
    