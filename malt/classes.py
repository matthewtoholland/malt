"""
A module for generating key properties of molecules needed to generate molecular file formats, specifically aimed at Tripos .mol2 files

"""

import rdkit
import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdqueries, rdPartialCharges
from biopandas.pdb import PandasPdb

#Some user-defined variables - at the moment these are hard-coded, I may find a better way of doing it later
NO_SUBSTRUCTS = 1
MOLECULE_TYPE = 'SMALL'
PARTIAL_CHARGES = 'GASTEIGER'

#Dictionaries that map RDKit properties to those expected in TRIPOS Mol2 Files
bond_types = {
    'AROMATIC': 'ar',
    'SINGLE': '1',
    'DOUBLE': '2',
    'TRIPLE': '3',
}
atom_types = {
    'SP': '1',
    'SP2': '2',
    'SP3': '3',
    'AROMATIC': 'ar',
    'S': ''
}

class Molecule(Chem.Mol):

    def __init__(self, path_to_pdb):
        self.path_to_pdb = path_to_pdb
        mol = Chem.MolFromPDBFile(path_to_pdb, removeHs=False)
        super().__init__(mol)

    def no_bonds(self):
        no_bonds = len(self.GetBonds())
        return no_bonds

    def no_atoms(self):
        no_atoms = len(self.GetAtoms())
        return no_atoms

    def name(self):
        """
        Extracts a rudimentary name from the molecule pdb file

        Returns: molecule_name
        Type: string
        """
        filename = os.path.basename(self.path_to_pdb)
        molecule_name = filename[:-4]
        return molecule_name

    def elements(self):
        """
        Returns a list of all the elements present in the input molecule, in increasing atomic number order with H last
        """
        element_dict = {}
        for atom in self.GetAtoms():
            if atom.GetSymbol() not in element_dict:
                element_dict[atom.GetSymbol()] = atom.GetAtomicNum()
        #now sort dictionary into TRIPOS format - increasing atomic no. with H at the end
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
        for atom in self.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol not in mol_by_elements:
                mol_by_elements[symbol] = []
            mol_by_elements[symbol].append(atom.GetIdx())
        return mol_by_elements


    def pdb_data(self):
        """
        Returns all of the information in the .pdb file in a Pandas dataframe. Useful for extracting atomic co-ordinates
        """
        ppdb = PandasPdb()
        ppdb.read_pdb(self.path_to_pdb)
        #For various reasons, all of the atoms in the VEHICLE Pdb datafile are classed as HETATMs
        pdb_data = ppdb.df['HETATM']
        pdb_data.set_index('atom_number')
        return pdb_data

    def get_partial_charges(func):
        """
        Computes and embeds Gasteiger partial charges for molecule, as implemented in RDKit rdPartialCharges
        """
        def wrapper(self, *args, **kwargs):
            rdPartialCharges.ComputeGasteigerCharges(self)
            result = func(self, *args, **kwargs)
            return result
        return wrapper

    def molecule_block(self):
        """
        Calculates and returns, in the correct format, the '@<TRIPOS>MOLECULE' block for the instance of the molecule
        The number of features and sets are hard-coded to 0
        """
        molecule_block = f'@<TRIPOS>MOLECULE\n{self.name()}\n{self.no_atoms()} {self.no_bonds()} {NO_SUBSTRUCTS} 0 0\n{MOLECULE_TYPE}\n{PARTIAL_CHARGES}\n'
        return molecule_block

    @get_partial_charges
    def atom_block(self):
        """
        Computes and returns, in the correct format, the '@<TRIPOS>ATOM block for the instance of the molecule
        """
        tripos_atom = pd.DataFrame(columns=['rdkit_index', 'atom_name', 'x_coord', 'y_coord',
                                            'z_coord', 'sybyl', 'substruct', 'substruct_name', 'partial_charge', 'atom_symbol', 'atom_index_label'])

        mol_by_elements = self.elements_by_index()
        substruct_name = self.name()
        pdb_data = self.pdb_data()
        #rdPartialCharges.ComputeGasteigerCharges(self)
        for atom in self.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()

            charge = float(atom.GetProp('_GasteigerCharge'))

            #Get co-ordinates for each atom
            x_coord = pdb_data['x_coord'][idx]
            y_coord = pdb_data['y_coord'][idx]
            z_coord = pdb_data['z_coord'][idx]

            atom_index_label = mol_by_elements[symbol].index(idx)+1
            atom_name = f'{symbol}{atom_index_label}'

            #Generate the sybyl code for each atom - the symbol, and the atom type e.g. aromatic, sp3 etc.
            if symbol == 'C' or symbol == 'N':
                if atom.GetIsAromatic():
                    sybyl = f'{symbol}.ar'
            elif symbol == 'H':
                sybyl = 'H'
            else:
                sybyl = f'{symbol}.{atom_types[str(atom.GetHybridization())]}'

            #Append atom to dataframe of atomic information
            tripos_atom = tripos_atom.append({'rdkit_index': idx, 'atom_name': atom_name, 'x_coord': x_coord, 'y_coord': y_coord, 'z_coord': z_coord,
                                            'sybyl': sybyl, 'substruct': NO_SUBSTRUCTS, 'substruct_name': substruct_name, 'partial_charge': "%.3f" % charge, 
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
        
        for index, bond in enumerate(self.GetBonds()):
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
