import numpy as np
import os
from hsr.pre_processing import *
from hsr.pca_transform import *
import pymol
from pymol import cmd
from rdkit.Geometry import Point3D
import pymol
from pymol import cmd

cwd = os.getcwd()

np.set_printoptions(precision=4, suppress=True)

molecule = load_molecules_from_sdf(f'{cwd}/original_molecule.sdf', removeHs=False, sanitize=True)[0]
molecule_array = molecule_to_ndarray(molecule, features=None)

# Transformed molecule with PCs orientation adjustment
transformed_molecule, PCs = compute_pca_using_covariance(molecule_array, return_axes=True)

# Simple PCA without PCs orientation adjustment
covariance_matrix = np.cov(molecule_array, rowvar=False, ddof=0,) 
eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
eigenvalues, eigenvectors = eigenvalues[::-1], eigenvectors[:, ::-1]

# PCA_molecule = np.dot(molecule_array, eigenvectors)
PCA_molecule = transformed_molecule
# print(transformed_molecule)

# Check if PCs and eigenvectors have sign changes
sign_changes = np.sign(np.sum(PCs * eigenvectors, axis=0))

new_molecule = Chem.Mol(molecule)

conf = new_molecule.GetConformer()
for i, pos in enumerate(PCA_molecule):
    conf.SetAtomPosition(i, Point3D(pos[0], pos[1], pos[2]))

# get max value for each coordinate (x,y and z) and the corresponding coordinates
PC1_max = np.max(PCA_molecule[:, 0])
PC1_max_coords = PCA_molecule[np.argmax(PCA_molecule[:, 0])].tolist()
PC2_max = np.max(PCA_molecule[:, 1])
PC2_max_coords = PCA_molecule[np.argmax(PCA_molecule[:, 1])].tolist()
PC3_max = np.max(PCA_molecule[:, 2])
PC3_max_coords = PCA_molecule[np.argmax(PCA_molecule[:, 2])].tolist()  
  
# print(f'PC1 max: {PC1_max}')
# print(f'PC2 max: {PC2_max}')
# print(f'PC3 max: {PC3_max}')

    
# save the molecule to an SDF file
# writer = Chem.SDWriter(f'{cwd}/transformed_molecule.sdf')
writer = Chem.SDWriter(f'{cwd}/PCA_molecule.sdf')
writer.write(new_molecule)
writer.close()

# Initialize PyMOL
pymol.finish_launching()

# Load the molecule
# molecule_path = f'{cwd}/transformed_molecule.sdf'
molecule_path = f'{cwd}/PCA_molecule.sdf'
cmd.load(molecule_path, 'molecule')

# create_axes()
from pymol.cgo import *
from pymol import cmd

w = 0.06 # cylinder width 
x = sign_changes[0]*(PC1_max+0.5) # cylinder length
y = sign_changes[1]*(PC2_max+0.5) # cylinder length
z = sign_changes[2]*(PC3_max+0.5) # cylinder length
h = 0.25 # cone hight
d = w * 1.618 # cone base diameter

obj = [CYLINDER, 0.0, 0.0, 0.0,   x, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
       CYLINDER, 0.0, 0.0, 0.0, 0.0,   y, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
       CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   z, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
       CONE,   x, 0.0, 0.0, x+sign_changes[0]*h, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 
       CONE, 0.0,   y, 0.0, 0.0, y+sign_changes[1]*h, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 
       CONE, 0.0, 0.0,   z, 0.0, 0.0, z+sign_changes[2]*h, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0,]

x_offset = x + (h + 0.1)*sign_changes[0]  
y_offset = y + (h + 0.1)*sign_changes[1]
z_offset = z + (h + 0.1)*sign_changes[2]
cmd.pseudoatom('label_x', pos=[x_offset, 0, 0], label='PC1', color=[1.0, 0.0, 0.0])
cmd.pseudoatom('label_y', pos=[0, y_offset, 0], label='PC2', color=[0.0, 1.0, 0.0])
cmd.pseudoatom('label_z', pos=[0, 0, z_offset], label='PC3', color=[0.0, 0.0, 1.0])

# cmd.pseudoatom('label_x', pos=[x_offset, 0, 0], color=[1.0, 0.0, 0.0])
# cmd.pseudoatom('label_y', pos=[0, y_offset, 0], color=[0.0, 1.0, 0.0])
# cmd.pseudoatom('label_z', pos=[0, 0, z_offset], color=[0.0, 0.0, 1.0])

# Create pseudatoms at each point of interest
# ctd = [0,0,0]
# ref_1 = [PC1_max, 0, 0]
# ref_2 = [0, PC2_max, 0]
# ref_3 = [0, 0, PC3_max]

# cmd.pseudoatom('reference_point_1', pos=ctd, color='black', label='Ref_1', )
# cmd.pseudoatom('reference_point_2', pos=ref_1, color='red', label='Ref_2')
# cmd.pseudoatom('reference_point_3', pos=ref_2, color='green', label='Ref_3')
# cmd.pseudoatom('reference_point_4', pos=ref_3, color='blue', label='Ref_4')

# cmd.pseudoatom('reference_point_1', pos=ctd, color='orange', label=' ', )
# cmd.pseudoatom('reference_point_2', pos=ref_1, color='lightpink', label=' ',)
# cmd.pseudoatom('reference_point_3', pos=ref_2, color='br5', label=' ',)
# cmd.pseudoatom('reference_point_4', pos=ref_3, color='brightorange', label=' ', )

# Markers settings
# cmd.show('spheres', 'reference_point_1')
# cmd.show('spheres', 'reference_point_2')
# cmd.show('spheres', 'reference_point_3')
# cmd.show('spheres', 'reference_point_4')

# cmd.set('sphere_scale', 0.1, 'reference_point_1')
# cmd.set('sphere_scale', 0.1, 'reference_point_2')
# cmd.set('sphere_scale', 0.1, 'reference_point_3')
# cmd.set('sphere_scale', 0.1, 'reference_point_4')

# def create_dotted_line(start, end, color=[1.0, 1.0, 1.0], segments=10):
#     """
#     Creates a series of short cylinders between two points to simulate a dotted line.
    
#     Parameters:
#     - start: The starting point of the line (list or tuple of x, y, z).
#     - end: The ending point of the line (list or tuple of x, y, z).
#     - color: The color of the line (list or tuple of r, g, b).
#     - segments: The number of segments (int).
#     """
#     obj = []
#     for i in range(segments):
#         fraction = i / float(segments)
#         next_fraction = (i + 1) / float(segments)
#         # Interpolate between start and end points
#         intermediate_point = [start[j] + (end[j] - start[j]) * fraction for j in range(3)]
#         next_intermediate_point = [start[j] + (end[j] - start[j]) * next_fraction for j in range(3)]
#         if i % 2 == 0:  
#             obj.extend([
#                 CYLINDER,
#                 intermediate_point[0], intermediate_point[1], intermediate_point[2],
#                 next_intermediate_point[0], next_intermediate_point[1], next_intermediate_point[2],
#                 0.03,  # Width of the cylinder
#                 color[0], color[1], color[2],
#                 color[0], color[1], color[2]
#             ])
#     return obj

# # Example usage
# dotted_line_obj = create_dotted_line(ref_1, PC1_max_coords, color=[1.0, 0.75, 0.87], segments=20)  # Red dotted line
# cmd.load_cgo(dotted_line_obj, 'dotted_line_1')

# # Repeat for other reference points as needed
# dotted_line_obj_2 = create_dotted_line(ref_2, PC2_max_coords, color=[0.6, 0.1, 0.5], segments=20)  # Green dotted line
# cmd.load_cgo(dotted_line_obj_2, 'dotted_line_2')

# dotted_line_obj_3 = create_dotted_line(ref_3, PC3_max_coords, color=[1.0, 0.7, 0.2], segments=20)  # Blue dotted line
# cmd.load_cgo(dotted_line_obj_3, 'dotted_line_3')


cmd.load_cgo(obj, 'axes')
