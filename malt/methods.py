"""
Module with the methods to generate TRIPOS Mol2 blocks, and save them to a file
"""
import os
import oddt
from oddt import shape

from malt.classes import Molecule

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


def mols2mol2(filename, num_of_molecules, **kwargs):
    """
    loops through 1 to num_of_molecules and writes them sequentially to a mol2 file
    args:
        path_to_pdb = path to directory of pdb files
        path_to_xyz = path to directory of xyz files
        path_to_charges = file containing charges for each atom of each molecule
        CalculateCharges = Bool - False if providing a file in path_to_charges
    """
    for key, value in kwargs.items():
        if key =='path_to_pdb':
            pdb_path = value
        elif key == 'path_to_xyz':
            xyz_path = value
        elif key == 'path_to_charges':
            charge_path = value
        elif key == 'CalculateCharges':
            Calculate = value
    
    for mol in range(1, num_of_molecules +1):
        pdb = os.path.join(pdb_path, f'S{mol}.pdb')
        xyz = os.path.join(xyz_path, f's{mol}.xyz')


        block = molecule_block(pdb, xyz, charge_path, CalculateCharges=Calculate)
        
        with open(filename, 'a+') as file:
            file.write(block)

    return None

def electroshape(filename, num_of_molecules, **kwargs):
    """
    Calculates electroshape, as implemented in oddt, for molecules from 1 to num_of_molecules, and saves to csv file
    """
    for key, value in kwargs.items():
        if key =='path_to_pdb':
            pdb_path = value
        elif key == 'path_to_xyz':
            xyz_path = value
        elif key == 'path_to_charges':
            charge_path = value
        elif key == 'CalculateCharges':
            Calculate = value

    for mol in range(1, num_of_molecules +1):
        pdb = os.path.join(pdb_path, f'S{mol}.pdb')
        xyz = os.path.join(xyz_path, f's{mol}.xyz')

        molecule = Molecule(pdb, xyz, charge_path, CalculateCharges=Calculate)
        molecule.create_dict()
        elecshape = list(shape.electroshape(molecule))

        line = f'S{mol}, {elecshape} \n'

        with open(filename, 'a+') as file:
            file.write(line)

    return None
