from rdkit import Chem
import numpy as np
import os
from hsr.pre_processing import *
from hsr.similarity import *
import matplotlib.pyplot as plt

def generate_all_single_deuterium_variants_and_compute_similarity(molecule, output_path):
    hydrogens = [atom for atom in molecule.GetAtoms() if atom.GetSymbol() == "H"]
    similarity_values = {}

    for i, hydrogen in enumerate(hydrogens):
        modified_molecule = Chem.Mol(molecule)  
        hydrogen_atom = modified_molecule.GetAtomWithIdx(hydrogen.GetIdx())
        hydrogen_atom.SetIsotope(2) 

        similarity = compute_similarity(molecule, modified_molecule)
        
        # Tagging the original hydrogen in the molecule
        molecule.GetAtomWithIdx(hydrogen.GetIdx()).SetProp("hydrogenID", f"H_{i+1}")
        # get atom coordinates
        pos = molecule.GetConformer().GetAtomPosition(hydrogen.GetIdx())
        
        molecule.SetProp(f"H_{i+1}", f"{pos.x},{pos.y},{pos.z}")
        
        # Store the similarity score with this tag
        similarity_values[f"H_{i+1}"] = similarity

    # Save the original molecule with hydrogen tags to an SDF file
    writer = Chem.SDWriter(output_path)
    writer.write(molecule)
    writer.close()

    return similarity_values

# Load your molecule here
cwd = os.getcwd()
sdf_file_path = os.path.join(cwd, f'{cwd}/original_molecule.sdf')
molecules = Chem.SDMolSupplier(sdf_file_path, removeHs=False, sanitize=False)
original_molecule = molecules[0]

# Assuming the output path for the tagged molecule
sdf_output_path = os.path.join(cwd, f'{cwd}/tagged_molecule.sdf')

# Generate deuterium variants, compute similarities, and tag hydrogens
similarity_values = generate_all_single_deuterium_variants_and_compute_similarity(original_molecule, sdf_output_path)

import pymol
from pymol import cmd
from pymol.cgo import *

# Initialize PyMOL
pymol.finish_launching()

def plot_similarity_values_and_get_colors(similarity_values_dict):
    tags, similarity_values = zip(*similarity_values_dict.items())
        
    min_val, max_val = min(similarity_values), max(similarity_values)
    normalized_values = [(val - min_val) / (max_val - min_val) for val in similarity_values]
    
    cmap = plt.get_cmap('viridis')
    
    fig, ax = plt.subplots()
    x_values = range(1, len(similarity_values) + 1)
    
    scatter = ax.scatter(x_values, similarity_values, c=normalized_values, cmap=cmap)
    for i, tag in enumerate(tags):
        ax.text(x_values[i], similarity_values[i], tag, fontsize=9, ha='right')
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Similarity')
    cbar_ticks = [0, 0.25, 0.5, 0.75, 1]  
    interp_vals = [min_val + (max_val - min_val) * t for t in cbar_ticks]
    cbar_tick_labels = [f"{val:.4f}" for val in interp_vals]
    print(cbar_tick_labels)
    cbar.set_ticks(cbar_ticks)  
    cbar.set_ticklabels(cbar_tick_labels) 
    
    ax.set_xlabel('Index')
    ax.set_ylabel('Similarity Value')
    ax.set_ylim([min_val - (max_val - min_val) * 0.05, max_val + (max_val - min_val) * 0.05])
    
    plt.tight_layout()
    plt.savefig(f'{os.getcwd()}/similarity_plot_with_tags.svg', format='svg')
    
    # Saving colorbar in a separate figure 
    fig, ax = plt.subplots(figsize=(2, 6 ))
    fig.subplots_adjust(right=0.5)
    cbar = fig.colorbar(scatter, cax=ax)
    # cbar.set_label('Similarity')
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1])
    cbar.set_ticklabels(cbar_tick_labels, fontsize=10, font='Arial')
    plt.tight_layout()
    plt.savefig(f'{os.getcwd()}/colorbar_only.svg', format='svg')
    plt.show()

    tag_to_color = {tag: cmap(norm)[:3] for tag, norm in zip(tags, normalized_values)}
    
    return tag_to_color

# Load the tagged molecule in PyMOL
cmd.load(sdf_output_path, 'tagged_molecule')

# Load the same molecule in RDKit
rdkit_mol = Chem.SDMolSupplier(sdf_output_path)[0]

coord_to_tag = {}
for tag in similarity_values.keys():
    pos = [float(x) for x in rdkit_mol.GetProp(tag).split(',')]
    coord_to_tag[tag] = pos

# MAKING SPHERES FOR THE HYDROGENS
colors = plot_similarity_values_and_get_colors(similarity_values)

# create a list of spheres with the correct coordinates and colors, and radius 1
d = 0.3
obj = []
for tag, color in zip(coord_to_tag.keys(), colors.values()):
    pos = coord_to_tag[tag]
    obj.extend([COLOR, color[0], color[1], color[2],  SPHERE, pos[0], pos[1], pos[2], d])

# show the spheres
cmd.load_cgo(obj, 'gradients')

# modify the transparency of the spheres objects
cmd.set("cgo_transparency", 0.0 , "gradients")
