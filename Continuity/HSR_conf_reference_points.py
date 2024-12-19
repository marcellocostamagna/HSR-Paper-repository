import numpy as np
import os
from hsr.pre_processing import *
from hsr.pca_transform import *
from rdkit.Geometry import Point3D
import pymol
from pymol import cmd
from pymol.cgo import *
from rdkit import Chem

cwd = os.getcwd()

np.set_printoptions(precision=4, suppress=True)


### FUNCTIONS ###

def original_array(molecule, removeHs=False):
    """
    Extract the true original (uncentered) coordinates of a molecule.
    """
    coordinates = []
    for atom in molecule.GetAtoms():
        # Skip hydrogens if removeHs is True
        if removeHs and atom.GetAtomicNum() == 1:
            continue
        position = molecule.GetConformer().GetAtomPosition(atom.GetIdx())
        coordinates.append([position.x, position.y, position.z])
    return np.array(coordinates)

def process_molecule(file_path):
    """
    Reads, processes a molecule file, and returns the original molecule, transformed molecule, PCs, and transformed coordinates.
    """
    molecule = read_mol_from_file(file_path, removeHs=False, sanitize=False)
    molecule_array = molecule_to_ndarray(molecule, features=None)
    transformed_molecule, PCs = compute_pca_using_covariance(molecule_array, return_axes=True)
    new_molecule = Chem.Mol(molecule)
    conf = new_molecule.GetConformer()
    for i, pos in enumerate(transformed_molecule):
        conf.SetAtomPosition(i, Point3D(pos[0], pos[1], pos[2]))
    return molecule, new_molecule, PCs, transformed_molecule, molecule_array

def rotate_and_translate(points, PCs, translation_vector):
    """
    Applies rotation and translation to a set of points.
    """
    rotated = np.dot(points, PCs.T)
    translated = rotated - translation_vector
    return translated

def create_dotted_line(start, end, color=[1.0, 1.0, 1.0], segments=10):
    """
    Creates a series of short cylinders between two points to simulate a dotted line.
    """
    obj = []
    for i in range(segments):
        fraction = i / float(segments)
        next_fraction = (i + 1) / float(segments)
        intermediate_point = [start[j] + (end[j] - start[j]) * fraction for j in range(3)]
        next_intermediate_point = [start[j] + (end[j] - start[j]) * next_fraction for j in range(3)]
        if i % 2 == 0:  
            obj.extend([
                CYLINDER,
                intermediate_point[0], intermediate_point[1], intermediate_point[2],
                next_intermediate_point[0], next_intermediate_point[1], next_intermediate_point[2],
                0.03,  # Width of the cylinder
                color[0], color[1], color[2],
                color[0], color[1], color[2]
            ])
    return obj

def create_rotated_axes_with_dotted_lines(PCs, name_prefix, molecule, base_color, translation_vector, transformed_molecule):
    """
    Creates rotated axes for visualization, applies rotation and translation,
    generates reference points, and adds dotted lines to max atoms.
    """
    # Compute axis ends based on the longest atom projection in the transformed space
    PC1_max = np.max(transformed_molecule[:, 0])
    PC2_max = np.max(transformed_molecule[:, 1])
    PC3_max = np.max(transformed_molecule[:, 2])

    # Find the coordinates of the atoms with maximum values
    PC1_max_coords = transformed_molecule[np.argmax(transformed_molecule[:, 0])].tolist()
    PC2_max_coords = transformed_molecule[np.argmax(transformed_molecule[:, 1])].tolist()
    PC3_max_coords = transformed_molecule[np.argmax(transformed_molecule[:, 2])].tolist()

    ref1 = np.array([PC1_max, 0, 0])  # End of PC1
    ref2 = np.array([0, PC2_max, 0])  # End of PC2
    ref3 = np.array([0, 0, PC3_max])  # End of PC3
    origin = np.array([0, 0, 0])

    # Apply rotation and translation to axes and reference points
    ref1_transformed = rotate_and_translate(ref1, PCs, translation_vector)
    ref2_transformed = rotate_and_translate(ref2, PCs, translation_vector)
    ref3_transformed = rotate_and_translate(ref3, PCs, translation_vector)
    origin_transformed = rotate_and_translate(origin, PCs, translation_vector)

    PC1_max_transformed = rotate_and_translate(PC1_max_coords, PCs, translation_vector)
    PC2_max_transformed = rotate_and_translate(PC2_max_coords, PCs, translation_vector)
    PC3_max_transformed = rotate_and_translate(PC3_max_coords, PCs, translation_vector)

    # Extend the cylinder lengths slightly beyond the reference points
    cylinder_extension = 0.5  # Extend cylinders by this factor
    cone_extension = 0.25  # Additional length for cones

    # Compute extended positions for the cylinders and cones
    ref1_extended = rotate_and_translate(ref1 + np.array([cylinder_extension, 0, 0]), PCs, translation_vector)
    ref2_extended = rotate_and_translate(ref2 + np.array([0, cylinder_extension, 0]), PCs, translation_vector)
    ref3_extended = rotate_and_translate(ref3 + np.array([0, 0, cylinder_extension]), PCs, translation_vector)

    # Compute cone tips
    cone_tip1 = rotate_and_translate(ref1 + np.array([cylinder_extension + cone_extension, 0, 0]), PCs, translation_vector)
    cone_tip2 = rotate_and_translate(ref2 + np.array([0, cylinder_extension + cone_extension, 0]), PCs, translation_vector)
    cone_tip3 = rotate_and_translate(ref3 + np.array([0, 0, cylinder_extension + cone_extension]), PCs, translation_vector)

    w = 0.06  # cylinder width
    d = w * 1.618  # cone base diameter

    # Construct the CGO object
    obj = [
        CYLINDER, *origin_transformed.tolist(), *ref1_extended.tolist(), w, *base_color, *base_color,
        CYLINDER, *origin_transformed.tolist(), *ref2_extended.tolist(), w, *base_color, *base_color,
        CYLINDER, *origin_transformed.tolist(), *ref3_extended.tolist(), w, *base_color, *base_color,
        CONE, *ref1_extended.tolist(), *cone_tip1.tolist(),
              d, 0.0, *base_color, *base_color, 1.0, 1.0,
        CONE, *ref2_extended.tolist(), *cone_tip2.tolist(),
              d, 0.0, *base_color, *base_color, 1.0, 1.0,
        CONE, *ref3_extended.tolist(), *cone_tip3.tolist(),
              d, 0.0, *base_color, *base_color, 1.0, 1.0,
    ]
    cmd.load_cgo(obj, f'{name_prefix}_axes')

    # Create reference points as pseudoatoms
    if base_color == [0.1, 0.6, 0.6]:
        psa_color = 'deepteal'
    elif base_color == [0.698, 0.13, 0.13]:
        psa_color = 'firebrick'
    
    # cmd.pseudoatom(f'{name_prefix}_origin', pos=origin_transformed.tolist(), color=psa_color, label='Origin')
    # cmd.pseudoatom(f'{name_prefix}_ref1', pos=ref1_transformed.tolist(), color=psa_color, label='PC1')
    # cmd.pseudoatom(f'{name_prefix}_ref2', pos=ref2_transformed.tolist(), color=psa_color, label='PC2')
    # cmd.pseudoatom(f'{name_prefix}_ref3', pos=ref3_transformed.tolist(), color=psa_color, label='PC3')

    cmd.pseudoatom(f'{name_prefix}_origin', pos=origin_transformed.tolist(), color=psa_color, label=' ')
    cmd.pseudoatom(f'{name_prefix}_ref1', pos=ref1_transformed.tolist(), color=psa_color, label=' ')
    cmd.pseudoatom(f'{name_prefix}_ref2', pos=ref2_transformed.tolist(), color=psa_color, label=' ')
    cmd.pseudoatom(f'{name_prefix}_ref3', pos=ref3_transformed.tolist(), color=psa_color, label=' ')

    # Add dotted lines
    dotted_line1 = create_dotted_line(ref1_transformed.tolist(), PC1_max_transformed.tolist(), color=base_color, segments=20)
    dotted_line2 = create_dotted_line(ref2_transformed.tolist(), PC2_max_transformed.tolist(), color=base_color, segments=20)
    dotted_line3 = create_dotted_line(ref3_transformed.tolist(), PC3_max_transformed.tolist(), color=base_color, segments=20)

    cmd.load_cgo(dotted_line1, f'{name_prefix}_dotted_line1')
    cmd.load_cgo(dotted_line2, f'{name_prefix}_dotted_line2')
    cmd.load_cgo(dotted_line3, f'{name_prefix}_dotted_line3')

    # Show spheres for reference points
    cmd.show('spheres', f'{name_prefix}_origin')
    cmd.show('spheres', f'{name_prefix}_ref1')
    cmd.show('spheres', f'{name_prefix}_ref2')
    cmd.show('spheres', f'{name_prefix}_ref3')

    cmd.set('sphere_scale', 0.12, f'{name_prefix}_origin')
    cmd.set('sphere_scale', 0.12, f'{name_prefix}_ref1')
    cmd.set('sphere_scale', 0.12, f'{name_prefix}_ref2')
    cmd.set('sphere_scale', 0.12, f'{name_prefix}_ref3')

### MAIN ###

# Choose what confermers to process: A or B

conformer = 'A' # 'B'

# Process both molecules
if conformer == 'A':
    # Conformer A
    file1 = f'{cwd}/conformers/conformersA_H.sdf'
    file2 = f'{cwd}/conformers/conformersA_Met.sdf'
elif conformer == 'B':
    # Conformer B
    file1 = f'{cwd}/conformers/conformersB_H.sdf'
    file2 = f'{cwd}/conformers/conformersB_Met.sdf'

# Define colors for each molecule
molecule1_color = [0.1, 0.6, 0.6]  # Deepteal
molecule2_color = [0.698, 0.13, 0.13]  # Firebrick

original_molecule1, molecule1, PCs1, transformed_molecule1, _ = process_molecule(file1)
original_molecule2, molecule2, PCs2, transformed_molecule2, _ = process_molecule(file2)

# Get the original coordinates of the molecules
original_array1 = original_array(original_molecule1)
original_array2 = original_array(original_molecule2)

# Compute the centers of the original molecules
center1 = np.mean(original_array1, axis=0)
center2 = np.mean(original_array2, axis=0)

# Rotate molecules and their PCs back to their original positions
rotated_molecule1 = rotate_and_translate(transformed_molecule1, PCs1, np.zeros(3))
rotated_molecule2 = rotate_and_translate(transformed_molecule2, PCs2, center1 - center2)

# Save rotated and translated molecules back to SDF files
rotated_mol1 = Chem.Mol(molecule1)
rotated_conf1 = rotated_mol1.GetConformer()
for i, pos in enumerate(rotated_molecule1):
    rotated_conf1.SetAtomPosition(i, Point3D(pos[0], pos[1], pos[2]))
writer1 = Chem.SDWriter(f'{cwd}/rotated_molecule1.sdf')
writer1.write(rotated_mol1)
writer1.close()

rotated_mol2 = Chem.Mol(molecule2)
rotated_conf2 = rotated_mol2.GetConformer()
for i, pos in enumerate(rotated_molecule2):
    rotated_conf2.SetAtomPosition(i, Point3D(pos[0], pos[1], pos[2]))
writer2 = Chem.SDWriter(f'{cwd}/rotated_molecule2.sdf')
writer2.write(rotated_mol2)
writer2.close()

# Initialize PyMOL
pymol.finish_launching()

# Load rotated and translated molecules into PyMOL
cmd.load(f'{cwd}/rotated_molecule1.sdf', 'molecule1')
cmd.load(f'{cwd}/rotated_molecule2.sdf', 'molecule2')

# Set molecule colors
cmd.color('deepteal', 'molecule1')
cmd.color('firebrick', 'molecule2')

# Create rotated axes with reference points and dotted lines
create_rotated_axes_with_dotted_lines(PCs1, 'molecule1', original_molecule1, molecule1_color, np.zeros(3), transformed_molecule1)
create_rotated_axes_with_dotted_lines(PCs2, 'molecule2', original_molecule2, molecule2_color, center1 - center2, transformed_molecule2)

# Zoom to view
cmd.zoom()
