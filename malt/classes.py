"""
A module for generating key properties of molecules needed to generate molecular file formats, specifically aimed at 
Tripos .mol2 files.

Currently this implementation only works if the xyz files have the same atom indices as the pdb files.

"""

import rdkit
import os
import csv
import numpy as np
import pandas as pd
import pkg_resources
from csv import reader 
from oddt import shape
from rdkit import Chem
from rdkit.Chem import rdPartialCharges
from .xyztomol import mol_from_xyz

#Set path to package data files referenced within the code
VEHICLE_ELECTROSHAPE = pkg_resources.resource_filename('malt', 'Data/dft_electroshape.csv')

#Some user-defined variables - at the moment these are hard-coded, I may find a better way of doing it later
NO_SUBSTRUCTS = 1
MOLECULE_TYPE = 'SMALL'

#Dictionaries that map RDKit properties to those expected in TRIPOS Mol2 Files
bond_types = {
    'AROMATIC': 'ar',
    'SINGLE': '1',
    'DOUBLE': '2',
    'TRIPLE': '3',
}

#The SP3D entry is a patch for a bug in the xyz implementation of S169, which does not register as aromatic.
atom_types = {
    'SP': '1',
    'SP2': '2',
    'SP3': '3',
    'AROMATIC': 'ar',
    'S': '', 
    'SP3D': '2'
}

class Molecule:

    def __init__(self, *args, CalculateCharges=True):
        self.path_to_xyz = None
        self._mol = None
        self._pdb_mol = None
        self._xyz_mol = None
        self.charges = None
        self.name = None
        self.num_bonds = None
        self.num_atoms = None
        self.index = None
        self._CalculateCharges = CalculateCharges
        self.atom_dict = None
        self.atoms = []
        self.smiles = None


        #Intialise rdkit mol objects for the input files
        if len(args) > 0:
            for arg in args:
                if arg.endswith('.pdb'):
                    path_to_pdb = arg
                    self._pdb_mol = Chem.MolFromPDBFile(path_to_pdb, removeHs=False)
                    self.name = os.path.basename(path_to_pdb)[:-4]
                    self.index = int(self.name[1:])
                elif arg.endswith('.xyz'):
                    self.path_to_xyz = arg
                    self._xyz_mol = mol_from_xyz(self.path_to_xyz)
                elif CalculateCharges == False:
                    path_to_charges = arg

        #If xyz file is provided, preferentially use the information from this over that of a pdb file
        if self._xyz_mol != None:
            self._mol = self._xyz_mol
        else:
            self._mol = self._pdb_mol 

        #Set various attributes
        self.num_bonds = len(self._mol.GetBonds())
        self.num_atoms = len(self._mol.GetAtoms())

        #Set charges
        if CalculateCharges == True:
            rdPartialCharges.ComputeGasteigerCharges(self._mol)
            gasteiger_charges = []
            for atom in self._mol.GetAtoms():
                charge = float(atom.GetProp('_GasteigerCharge'))
                gasteiger_charges.append(charge)
            self.charges = gasteiger_charges 
        else:
            self.get_external_charges(path_to_charges)

        #Set smiles
        self.smiles = Chem.MolToSmiles(self._mol)
                

    def elements(self):
        """
        Returns a list of all the elements present in the input molecule, in increasing atomic number order with H last
        """
        element_dict = {}
        for atom in self._mol.GetAtoms():
            if atom.GetSymbol() not in element_dict:
                element_dict[atom.GetSymbol()] = atom.GetAtomicNum()

        #Now sort dictionary into TRIPOS format - increasing atomic no. with H at the end
        element_dict = dict(sorted(element_dict.items(), key = lambda item: item[1]))
        elements = [key for key in element_dict]
        if 'H' in elements:
            elements.remove('H')
            elements.append('H')
        return elements

    def elements_by_index(self):
        """
        Returns a dictionary with element symbols as keys, and the atom ID's of those elements in the molecule as values
        """
        mol_by_elements = {}
        for atom in self._mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol not in mol_by_elements:
                mol_by_elements[symbol] = []
            mol_by_elements[symbol].append(atom.GetIdx())
        return mol_by_elements

    def get_external_charges(self, filename):
        """
        Retrieves a list of partial charges for the molecule from the list of DFT VEHICLE Charges
        args:
            filepath - path to file of partial charges
        returns
            charges - list of atomic partial charges for the molecule
        """
        with open(filename, 'r') as file:
            csv_reader = csv.reader(file)
            all_charges = list(csv_reader)
        
        self.charges = all_charges[self.index -1]
        self.charges = [float(charge) for charge in self.charges]

        return None

    @property
    def coords(self):
        #Set up coordinates. If xyz file is present then these coordinates are used over .pdb
        if self._xyz_mol != None:
            return self._xyz_mol.GetConformer(0).GetPositions()
        else:
            return self._mol.GetConformer(0).GetPositions()

    def molecule_block(self):
        """
        Calculates and returns, in the correct format, the '@<TRIPOS>MOLECULE' block for the instance of the molecule
        The number of features and sets are hard-coded to 0
        """
        if self._CalculateCharges == True:
            charge_type = 'GASTEIGER'
        else:
            charge_type = 'DFT'
        
        molecule_block = (f'@<TRIPOS>MOLECULE\n'
                            f'{self.name}\n'
                            f'{self.num_atoms} {self.num_bonds} {NO_SUBSTRUCTS} 0 0\n'
                            f'{MOLECULE_TYPE}\n'
                            f'{charge_type}\n')

        return molecule_block 

    def atom_block(self):
        """
        Computes and returns, in the correct format, the '@<TRIPOS>ATOM block for the instance of the molecule. This calculates partial charges using rdkit's
        implementation of the Gasteiger partial charges. Atomic coordinates are extracted from the pdb file.
        """
        tripos_atom = pd.DataFrame(columns=['rdkit_index', 'atom_name', 'x_coord', 'y_coord',
                                            'z_coord', 'sybyl', 'substruct', 'substruct_name',
                                             'partial_charge', 'atom_symbol', 'atom_index_label'])

        mol_by_elements = self.elements_by_index()
        substruct_name = self.name

        for atom in self._mol.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()

            charge = self.charges[idx]

            #Get co-ordinates for each atom
            x_coord = self.coords[idx][0]
            y_coord = self.coords[idx][1]
            z_coord = self.coords[idx][2]

            atom_index_label = mol_by_elements[symbol].index(idx)+1
            atom_name = f'{symbol}{atom_index_label}'

            #Generate the sybyl code for each atom - the symbol, and the atom type e.g. aromatic, sp3 etc.
            if symbol == 'C' or symbol == 'N':
                if atom.GetIsAromatic():
                    sybyl = f'{symbol}.ar'
                else:
                    sybyl = f'{symbol}.{atom_types[str(atom.GetHybridization())]}'
            elif symbol == 'H':
                sybyl = 'H'
            else:
                sybyl = f'{symbol}.{atom_types[str(atom.GetHybridization())]}'

            #Append atom to dataframe of atomic information
            tripos_atom = tripos_atom.append({'rdkit_index': idx, 'atom_name': atom_name, 'x_coord': x_coord,
                                             'y_coord': y_coord, 'z_coord': z_coord,
                                            'sybyl': sybyl, 'substruct': NO_SUBSTRUCTS, 'substruct_name': substruct_name, 
                                            'partial_charge': "%.3f" % charge, 
                                            'atom_symbol': symbol, 'atom_index_label': atom_index_label}, ignore_index=True)
        
        #Sort dataframe by sybyl (elements first, followed by atom index label
        tripos_atom['atom_symbol'] = pd.Categorical(tripos_atom['atom_symbol'], self.elements())
        tripos_atom = tripos_atom.sort_values(by=['atom_symbol', 'atom_index_label'])
        tripos_atom.index = np.arange(1, len(tripos_atom)+1)
        
        #TRIPOS atom ID and RDKit index are not the same, need to generate a mapping from one to t'other
        index_mapper = pd.DataFrame(data=np.arange(1, len(tripos_atom)+1), index=tripos_atom['rdkit_index'])
        global index_lookup
        index_lookup = index_mapper.to_dict()[0]

        #Generate final dataframe, and return as a string
        atom_block_df = tripos_atom[['atom_name', 'x_coord',
                                  'y_coord', 'z_coord', 'sybyl', 'partial_charge']]
        
        atom_block = '@<TRIPOS>ATOM\n' + atom_block_df.to_string(header=False) + '\n'

        return atom_block 

    def bond_block(self):
        tripos_bond = pd.DataFrame(
            columns=['begin_atom_rdkit', 'end_atom_rdkit', 'bond_type'])
    
        for index, bond in enumerate(self._mol.GetBonds()):
            beginning = bond.GetBeginAtom().GetIdx()
            end = bond.GetEndAtom().GetIdx()
            
            #Get bond type and convert to TRIPOS bond type from rdkit bond type
            bond_type = str(bond.GetBondType())
            bond_type = bond_types[bond_type]

            tripos_bond = tripos_bond.append(
                {'begin_atom_rdkit': beginning, 'end_atom_rdkit': end, 'bond_type': bond_type}, ignore_index=True)
            
        #Convert RDKit indices to TRIPOS atom ID numbers so there is consistency between ATOM block and BOND block
        tripos_bond['begin'] = tripos_bond.apply(lambda row: index_lookup[row.begin_atom_rdkit], axis=1)
        tripos_bond['end'] = tripos_bond.apply(lambda row: index_lookup[row.end_atom_rdkit], axis=1)
        tripos_bond['bond_type'] = pd.Categorical(tripos_bond['bond_type'], list(bond_types.values()))
        tripos_bond = tripos_bond.sort_values(by=['bond_type', 'begin'])
        tripos_bond.index = np.arange(1, len(tripos_bond)+1)
        tripos_bond = tripos_bond[['begin', 'end', 'bond_type']]

        bond_block = '@<TRIPOS>BOND\n' + tripos_bond.to_string(header=False) + '\n'

        return bond_block

    def create_dict(self):
        """
        Creates dictionary of atom information - required for calculating electroshape as implemented in oddt
        """
        self.create_atom_list()
        atoms = self.atoms
        atom_dtype = [('name', object), ('coords', np.float32, 3), ('radius', np.float32), ('charge', np.float32)]
        atom_dict = np.empty(len(atoms), dtype=atom_dtype)
        for atom in self._mol.GetAtoms():
            idx = atom.GetIdx()
            name = atom.GetSymbol() + str(idx)
            radius = Chem.GetPeriodicTable().GetRcovalent(atom.GetAtomicNum())
            charge = self.charges[idx]
            all_coords = self.coords[idx]
            coords = [all_coords[0], all_coords[1], all_coords[2]]
            atom_dict[idx] = (name, coords, radius, charge)
        self.atom_dict = atom_dict

    def create_atom_list(self):
        """
        Creates list of atom information required for calculating electroshape as implemented in oddt
        """
        atoms = []
        for atom in self._mol.GetAtoms():
            idx = atom.GetIdx()
            name = atom.GetSymbol() + str(idx)
            charge = self.charges[idx]
            xcoor = self.coords[idx][0]
            ycoor = self.coords[idx][1]
            zcoor = self.coords[idx][2]
            radius = Chem.GetPeriodicTable().GetRcovalent(atom.GetAtomicNum())
            atom = [name, xcoor, ycoor, zcoor, radius, charge]
            atoms.append(atom)
        self.atoms = atoms

    def electroshape(self):
        """
        Calculates the electroshape vector for the molecule
        """
        self.create_dict()
        elecshape = list(shape.electroshape(self))

        return elecshape

    def filter_electroshape(self, n=100):
        """
        returns the top n (default is 100) most similar molecules in the VEHICLE database by electroshape
        """

        with open(VEHICLE_ELECTROSHAPE, 'r') as file:
            csv_reader = reader(file)
            vehicle_data = []
            for row in csv_reader:
                molecule_index = row[0]
                smiles = row[1]
                electroshape_vector = row[2:]
                electroshape_vector = [float(num) for num in electroshape_vector]
                mol_data = [molecule_index, smiles, electroshape_vector]
                vehicle_data.append(mol_data)

        electroframe = pd.DataFrame(vehicle_data, columns=['Index', 'smiles', 'electro_vector'])
        electroframe = electroframe.set_index('Index')
        electroframe.insert(2, 'similarity', 0.0)

        query = np.array(self.electroshape())
        
        for idx in electroframe.index:
            electroframe.at[idx, 'similarity'] = shape.usr_similarity(query, np.array(electroframe.at[idx, 'electro_vector']))

        electroframe = electroframe.sort_values(by=['similarity'], ascending=False, ignore_index=False)
        electroframe = electroframe.round({'similarity': 3})

        electroframe = electroframe.head(n)

        electroframe.to_csv(f'{self.name}_electroshape_similarity', columns=['smiles', 'similarity'], index=True)

        return None
