import numpy as np  
import sys
import os 
from hsr.pre_processing import *
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSR, GetUSRCAT
sys.path.append(os.path.abspath('../'))
from usr import *
from perturbations import *
import hsr
cwd = os.getcwd()

print(f'\nUSR, USRCAT and HSR similarity of structural analogues')

# Load molecules from SDF files
co_3 = load_molecules_from_sdf(f'{cwd}/Co_3.sdf', removeHs=False, sanitize=False)[0]
co_2 = load_molecules_from_sdf(f'{cwd}/Co_2.sdf', removeHs=False, sanitize=False)[0]
ir_3 = load_molecules_from_sdf(f'{cwd}/Ir_3.sdf', removeHs=False, sanitize=False)[0]

mols = [co_3, co_2, ir_3]

pairs = {'Co_2 vs Co_3': (co_2, co_3),
        'Co_2 vs Ir_3': (co_2, ir_3),
        'Co_3 vs Ir_3': (co_3, ir_3)}

# Compute HSR, USR and USRCAT similarity between all pairs of molecules
for pair_name, (mol1, mol2) in pairs.items():
    print(f"\n{'-'*20}")
    print(f"\033[1m{pair_name}\033[0m")
    
    # Compute USR similarity 
    usrs = [GetUSR(mol) for mol in [mol1, mol2]]
    similarity = GetUSRScore(usrs[0], usrs[1])
    print(f"USR: {similarity:.4f}")

    # Compute USRCAT similarity
    usrcats = [GetUSRCAT(mol) for mol in [mol1, mol2]]
    similarity = GetUSRScore(usrcats[0], usrcats[1])
    print(f"USRCAT: {similarity:.4f}")
    
    # Compute HSR similarity
    hsr_similarity = hsr.similarity.compute_similarity(mol1, mol2)
    print(f"HSR: {hsr_similarity:.4f}")
    


