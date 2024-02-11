import os
import time
import numpy as np
import multiprocessing
from rdkit import Chem
from hsr.pre_processing import *
from hsr.pca_transform import * 
from hsr.fingerprint import *
from hsr.similarity import *
from hsr.utils import *
from oddt import toolkit

MAX_CORES = 4

def read_molecules_from_file(file_path):
    mols = []
    for mol in toolkit.readfile(os.path.splitext(file_path)[1][1:], file_path):
        mols.append(mol)
    return mols

def read_molecules_from_sdf(sdf_file):
    supplier = Chem.SDMolSupplier(sdf_file, removeHs=False, sanitize=False)
    mols = [mol for mol in supplier if mol]
    return mols

def calculate_enrichment_factor(y_true, y_scores, percentage):
    n_molecules = len(y_true)
    n_actives = sum(y_true)
    f_actives = n_actives / n_molecules
    
    n_top_molecules = int(n_molecules * percentage)
    top_indices = np.argsort(y_scores)[-n_top_molecules:]
    n_top_actives = sum(y_true[i] for i in top_indices)
    
    expected_actives = percentage * n_molecules * f_actives
    enrichment_factor = n_top_actives / expected_actives
    return enrichment_factor

def compute_fingerprints(molecules, method):
    if method == "pseudo_usr":
        # return {mol: get_nd_fingerprint(mol, features=None, scaling_method='matrix') for mol in molecules}
        return {mol: generate_fingerprint_from_molecule(mol, features=None, scaling='matrix', removeHs=True) for mol in molecules}
    elif method == "pseudo_usr_cat":
        return {mol: get_pseudo_usrcat_fingerprint(mol) for mol in molecules}
    elif method == "pseudo_electroshape":
        # return {mol: get_nd_fingerprint(mol, features=PSEUDO_ELECTROSHAPE_FEATURES, scaling_method='matrix') for mol in molecules}
        return {mol: generate_fingerprint_from_molecule(mol, features=PSEUDO_ELECTROSHAPE_FEATURES, scaling='matrix', removeHs=True) for mol in molecules}

# PSEUDO_USRCAT PARAMETERS & FUNCTIONS
USRCAT_SMARTS = {'hydrophobic' : "[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]",        
                  'aromatic' : "[a]",                                           
                  'acceotor' :"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]",   
                  'donor' : "[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]",                             
                }

def get_pseudo_usrcat_fingerprint(mol):
    mol_3d = molecule_to_ndarray(mol, features=None, removeHs=True)
    mol_3d_pca = compute_pca_using_covariance(mol_3d)
    
    pseudo_usrcat_fingerprint = []
    scaling_matrix = compute_scaling_matrix(mol_3d_pca)
    pseudo_usrcat_fingerprint.append(generate_fingerprint_from_transformed_data(mol_3d_pca, scaling_matrix))
    
     # Create a mapping from original indices to mol_3d indices
    non_hydrogen_atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    index_mapping = {original_idx: new_idx for new_idx, original_idx in enumerate(non_hydrogen_atom_indices)}

    for smarts in USRCAT_SMARTS.values():
        # Collect the atoms indexes that match the SMARTS pattern in the query molecule removing the hydrogens
        query_atoms_matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
        
        if not query_atoms_matches:
            pseudo_usrcat_fingerprint.append(np.zeros(12))
            continue
        # Translate original indices to mol_3d indices
        query_atoms = [index_mapping[idx] for match in query_atoms_matches for idx in match if idx in index_mapping]
        try:
            sub_molecule_3d = mol_3d_pca[query_atoms]
        except:
            print('ERROR')
            
        scaling_matrix = compute_scaling_matrix(sub_molecule_3d)
        pseudo_usrcat_fingerprint.append(generate_fingerprint_from_transformed_data(sub_molecule_3d, scaling_matrix))
    return pseudo_usrcat_fingerprint

def get_pseudo_usrcat_similarity(query_usrcat_fingerprint, target_usrcat_fingerprint):
    distances = []
    n_components = len(query_usrcat_fingerprint)
    for i in range(len(query_usrcat_fingerprint)):
        distances.append(calculate_manhattan_distance(query_usrcat_fingerprint[i], target_usrcat_fingerprint[i]))
    final_distance= np.sum(distances)
    return calculate_similarity_from_distance(final_distance, n_components)

# PSEUDO_ELECTROSHAPE PARAMETERS & FUNCTIONS
def partial_charges(mol, mol_charges):
    # Set atomic property partial charge from pybel molecule
    charges = mol_charges.atom_dict['charge']
    for atom in mol.GetAtoms():
        # Set partial charge = 0 if hydrogen
        if atom.GetAtomicNum() == 1:
            atom.SetProp('partial_charge', str(0.0))
        else:
            atom.SetProp('partial_charge', str(charges[atom.GetIdx()]))   
    return mol

def get_partial_charges(atom):
    try:
        partial_charge = float(atom.GetProp('partial_charge'))
    except KeyError:
        print('property not found')
        return 0.0
    # Handle the case where the partial charge is NaN or Inf with np.nan_to_num
    partial_charge = np.nan_to_num(partial_charge)
    result = partial_charge * 25
    # Handle possible overflows of maximum value allowed for float
    if result > np.finfo(np.float32).max:
        result = np.finfo(np.float32).max
    return result

PSEUDO_ELECTROSHAPE_FEATURES = { 'partial_charge' : get_partial_charges }

def process_folder(args):
    folder, root_directory, method, enrichment_factors = args
    print(f"\nProcessing folder: {folder}")
    
    folder_path = os.path.join(root_directory, folder)
    folder_enrichments = {k: [] for k in enrichment_factors.keys()}
    results = {}
    
    if os.path.isdir(folder_path):
        actives_file = os.path.join(folder_path, "actives_final.sdf")
        decoys_file = os.path.join(folder_path, "decoys_final.sdf")
        
        actives = read_molecules_from_sdf(actives_file)
        decoys = read_molecules_from_sdf(decoys_file)
        
        # Set partial charges as in the oddt implementation (Pybel molecules)
        if method == 'pseudo_electroshape':
            actives_charges = read_molecules_from_file(actives_file)
            decoys_charges = read_molecules_from_file(decoys_file)
            actives = [partial_charges(mol, mol_charge) for mol, mol_charge in zip(actives, actives_charges)]
            decoys = [partial_charges(mol, mol_charge) for mol, mol_charge in zip(decoys, decoys_charges)]

        all_mols = actives + decoys

        # Pre-compute fingerprints
        fingerprints = compute_fingerprints(all_mols, method)

        for query_mol in actives:
            other_mols = list(all_mols)  # create a new list
            other_mols.remove(query_mol)  # remove the query molecule from the list

            query_fp = fingerprints[query_mol]
            if method == 'pseudo_usr_cat':
                y_scores = [get_pseudo_usrcat_similarity(query_fp, fingerprints[mol]) for mol in other_mols]
            else:
                y_scores = [compute_similarity_score(query_fp, fingerprints[mol]) for mol in other_mols]

            y_true = [1 if mol in actives else 0 for mol in other_mols]
            
            for percentage in enrichment_factors.keys():
                ef = calculate_enrichment_factor(y_true, y_scores, percentage)
                folder_enrichments[percentage].append(ef)

        # Return the enrichment factors for this folder
        avg_folder_enrichments = {percentage: np.mean(values) for percentage, values in folder_enrichments.items()}
        results[folder] = avg_folder_enrichments
        
        # Print the average enrichment factors for this folder
        print(f"\nAverage Enrichment Factors for folder {folder}:")
        for percentage, ef in avg_folder_enrichments.items():
            print(f"Enrichment Factor at {percentage*100}%: {ef}")

    return results

if __name__ == "__main__":
    print(f'CWD: {os.getcwd()}')
    root_directory = f"{os.getcwd()}/validation/all"
    methods =  [ 'pseudo_electroshape']#['pseudo_usr', 'pseudo_usr_cat', 'pseudo_electroshape']
    enrichment_factors = {0.0025: [], 0.005: [], 0.01: [], 0.02: [], 0.03: [], 0.05: []}
    
    overall_results = {}

    for method in methods:
        print(f"\nRunning benchmark for method: {method}")
        start_time_total = time.time()

        folders = sorted(os.listdir(root_directory))
        args_list = [(folder, root_directory, method, enrichment_factors) for folder in folders]
        
        with multiprocessing.Pool(processes=MAX_CORES) as pool:
            results_list = pool.map(process_folder, args_list)

        # Combine all results
        results = {k: v for res in results_list for k, v in res.items()}
        
        # Calculate average enrichments
        avg_enrichments = {percentage: np.mean([res[percentage] for res in results.values()]) for percentage in enrichment_factors.keys()}
        overall_results[method] = avg_enrichments

        end_time_total = time.time()
        print(f"\nTotal processing time for {method}: {end_time_total - start_time_total:.2f} seconds")

    # Print the overall enrichment factors for the different methods
    print("\nOverall Enrichment Factors for Different Methods:")
    for method, results in overall_results.items():
        print(f"\nResults for {method}:")
        for percentage, ef in results.items():
            print(f"Enrichment Factor at {percentage*100}%: {ef}")

