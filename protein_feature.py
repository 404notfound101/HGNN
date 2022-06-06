from Bio.PDB import PDBParser
from scipy.spatial import cKDTree
import numpy as np
import os

def generate_prot_graph(PDB_dir, save_dir, threshold):
    parser = PDBParser()
    for file in os.listdir(PDB_dir):
        protein = parser.get_structure('input', os.path.join(PDB_dir, file))
        chain = protein.get_chains[0]
        residues = []
        for res in chain:
            try:
                atom_a = res['CA'].get_coord()
                residues.append(atom_a)
            except KeyError:
                continue
        residues = np.array(residues)
        tree = cKDTree(residues)
        prot_edge = tree.sparse_distance_matrix(tree, threshold,p=2.0)
        prot_edge = prot_edge.toarray()
        prot_mask = np.where(prot_edge>0, 1.0, 0.0)
        edge_dict = {}
        edge_dict['edge_feature'] = prot_edge
        edge_dict['edge_mask'] = prot_mask
        np.save(os.path.join(save_dir, file), edge_dict)
        
