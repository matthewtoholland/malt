"""
Module with the methods to generate TRIPOS Mol2 blocks, and save them to a file
"""

from .classes import Molecule

def molecule_block(*args, **kwargs):
    """
    Generates the TRIPOS Mol2 block for a given molecule, returned as a string
    """
    mol = Molecule(*args, **kwargs)
    block = mol.molecule_block() + mol.atom_block() + mol.bond_block() + '\n'

    return block

def string2mol2(filename, string):
    """
    Writes molecule to filename.mol2 file, input is a string of Mol2 blocks
    """
    block = string
    
    if filename[-4:] != '.mol2':
        filename += '.mol2'

    with open(filename, 'w') as file:
        file.write(block)

    return None

def molecule2mol2(filename, *args, **kwargs):
    """"
    Creates a Mol2 block out of molecule files, and then saves to filename.mol2
    """
    block = molecule_block(*args, **kwargs)

    if filename[-4:] != '.mol2':
        filename += '.mol2'

    with open(filename, 'w') as file:
        file.write(block)
    
    return None

def mols2mol2(filename, num_of_molecules, *args, **kwargs):
    """
    loops through 1 to num_of_molecules and writes them sequentially to a mol2 file
    """
    