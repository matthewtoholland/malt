"""
Module with the methods to generate TRIPOS Mol2 blocks, and save them to a file
"""

import rdkit
from .classes import Molecule

def molecule_block(path_to_pdb_file):
    """
    Generates the TRIPOS Mol2 block for a given molecule, returned as a string
    """
    mol = Molecule(path_to_pdb_file)
    block = mol.molecule_block() + mol.atom_block() + mol.bond_block() + '\n'

    return block

def write_mol2(input, filename):
    """
    Writes molecule to filename.mol2 file, input can either be path to .pdb file, or a string of Mol2 blocks
    """
    if input[-4:] == '.pdb':
        block = molecule_block(input)
    else: 
        block = input
    
    if filename[-4:] != '.mol2':
        filename += '.mol2'

    with open(filename, 'w') as file:
        file.write(block)

    return None

