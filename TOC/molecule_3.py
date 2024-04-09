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

molecule = load_molecules_from_sdf(f'{cwd}/molecule_3.sdf', removeHs=False, sanitize=False)[0]

# cyclopentadienyl carbons: 2,3,6,8,10
# Get coordinates of teh center of the cyclopentadienyl ligand
penta_coords = []
for atom in molecule.GetAtoms():
    if atom.GetSymbol() == 'C':
        if atom.GetIdx() in [1, 2, 5, 7, 9]:
            position = molecule.GetConformer().GetAtomPosition(atom.GetIdx())
            penta_coords.append([position.x, position.y, position.z])
penta_coords = np.array(penta_coords)
penta_ctd = (np.mean(penta_coords, axis=0)).tolist()

# Manganese coordinates
for atom in molecule.GetAtoms():
    if atom.GetSymbol() == 'Mn':
        position = molecule.GetConformer().GetAtomPosition(atom.GetIdx())
        Mn_coords = [position.x, position.y, position.z]

# Initialize PyMOL
pymol.finish_launching()

# Load the molecule
molecule_path = f'{cwd}/molecule_3.sdf'
cmd.load(molecule_path, 'molecule')

# create_axes()
from pymol.cgo import *
from pymol import cmd


cmd.pseudoatom('penta_ctd', pos=penta_ctd, color='red')
cmd.pseudoatom('Mn', pos=Mn_coords, color='yellow')


def create_dotted_line(start, end, color=[1.0, 1.0, 1.0], segments=10):
    """
    Creates a series of short cylinders between two points to simulate a dotted line.
    
    Parameters:
    - start: The starting point of the line (list or tuple of x, y, z).
    - end: The ending point of the line (list or tuple of x, y, z).
    - color: The color of the line (list or tuple of r, g, b).
    - segments: The number of segments (int).
    """
    obj = []
    for i in range(segments):
        fraction = i / float(segments)
        next_fraction = (i + 1) / float(segments)
        # Interpolate between start and end points
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

# Example usage
dotted_line_obj = create_dotted_line(Mn_coords, penta_ctd, color=[0.5, 0.5, 0.5], segments=20)
cmd.load_cgo(dotted_line_obj, 'dotted_line_1')



