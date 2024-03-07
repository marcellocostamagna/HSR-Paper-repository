import numpy as np  
from hsr.pre_processing import *
from hsr.pca_transform import * 
from hsr.fingerprint import *
from hsr.similarity import *
from hsr.utils import *
import os 

np.set_printoptions(precision=4, suppress=True)

cwd = os.getcwd()

# List of molecules from SDF file
molecules = load_molecules_from_sdf(f'{cwd}/simple_case.sdf', removeHs=False, sanitize=False)

# MOLECULES
molecule_A = molecules[1]
molecule_B = molecules[0]


# PRE-PROCESSING
molecule_A_matrix = molecule_to_ndarray(molecule_A, features=EXAMPLE_FEATURES)
molecule_B_matrix = molecule_to_ndarray(molecule_B, features=EXAMPLE_FEATURES)

print(f'Molecules initial matrix: \nmolecule_A:\n {molecule_A_matrix} \nmolecule_B\n {molecule_B_matrix}')


# TODO: covariance, eigenvalues,eigenvectors and transformed data
print(f'Molecule_A:\n')
molecule_A_tranformed_data = compute_pca_using_covariance(molecule_A_matrix, print_steps=True)

print(f'Molecule_B:\n')
molecule_B_tranformed_data = compute_pca_using_covariance(molecule_B_matrix, print_steps=True)


# FINGERPRINTS
molecule_A_scaling = compute_scaling_matrix(molecule_A_tranformed_data)
molecule_A_fingerprint = generate_fingerprint_from_transformed_data(molecule_A_tranformed_data, molecule_A_scaling)
molecule_B_scaling = compute_scaling_matrix(molecule_B_tranformed_data)
molecule_B_fingerprint = generate_fingerprint_from_transformed_data(molecule_B_tranformed_data, molecule_B_scaling)

print(f'Fingerprints:\n \nmolecule_A:\n {molecule_A_fingerprint}\n \nmolecule_B\n {molecule_B_fingerprint}\n')

# SIMILARITY
similarity = compute_similarity_score(molecule_A_fingerprint, molecule_B_fingerprint)
print(f'Similarity: {similarity:.4f}')
