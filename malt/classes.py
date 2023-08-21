"""
A module for generating key properties of molecules needed to generate
molecular file formats, specifically aimed at Tripos .mol2 files.

Currently this implementation only works if the xyz files have the same atom
indices as the pdb files.

"""

import os
import numpy as np
import pandas as pd
import pkg_resources
from csv import reader
from rdkit import Chem
from rdkit.Chem import rdPartialCharges, AllChem
from malt.xyztomol import mol_from_xyz
from malt.mol2tomol import vehicle_mol2, mol2


# Some user-defined variables - at the moment these are hard-coded,
# I may find a better way of doing it later
VEHICLE_MOL2 = pkg_resources.resource_filename('malt', 'Data/vehicle_dft.mol2')


NO_SUBSTRUCTS = 1
MOLECULE_TYPE = 'SMALL'

# Dictionaries that map RDKit properties to those expected in TRIPOS Mol2 Files
bond_types = {
    'AROMATIC': 'ar',
    'SINGLE': '1',
    'DOUBLE': '2',
    'TRIPLE': '3',
}

# The SP3D entry is a patch for a bug in the xyz implementation of S169,
# which does not register as aromatic.
atom_types = {
    'SP': '1',
    'SP2': '2',
    'SP3': '3',
    'AROMATIC': 'ar',
    'S': '', 
    'SP3D': '2'
}


class Molecule:

    def __init__(self, *args, CalculateCharges=True, smiles=None, name=None):
        self.path_to_xyz = None
        self._mol = None
        self._mol2 = None
        self._pdb_mol = None
        self._xyz_mol = None
        self.charges = None
        self.name = name
        self.num_bonds = None
        self.num_atoms = None
        self.index = None
        self._CalculateCharges = CalculateCharges
        self.atom_dict = None
        self.atoms = []
        self.smiles = None
        self.s_flag = False

        # Initialise rdkit mol objects for the input files
        if len(args) > 0:
            for arg in args:
                if arg is None:
                    continue
                # Initialise xyz details
                    self._pdb_mol = Chem.MolFromPDBFile(path_to_pdb, removeHs=False)
                    if self.name == None:
                        self.name = os.path.basename(path_to_pdb)[:-4]
                        self.index = int(self.name[1:])

                elif arg.endswith('.xyz'):
                    self.path_to_xyz = arg
                    self._xyz_mol = mol_from_xyz(self.path_to_xyz)

        elif smiles is not None:
            if arg.endswith('.mol2'):
                path_to_mol2 = arg
                print('chicken')
                self._mol2 = mol2(path_to_mol2)
                self._mol = Chem.MolFromMol2File(path_to_mol2, sanitize=False, removeHs=False)
                if self.name == None:
                    self.name = path_to_mol2[:-5]
            elif CalculateCharges == False:
                self._path_to_charges = arg

            #Molecule is already within VEHICLe database
            elif arg[0] == 'S' or 's':
                self._mol2 = vehicle_mol2(arg)
                self._mol =  self._mol2.mol_from_mol2()
                self.s_flag = True
                if self.name == None:
                    self.name = arg

            self.smiles = smiles
            self._mol = Chem.MolFromSmiles(smiles)
            self._mol = Chem.AddHs(self._mol)
            AllChem.EmbedMolecule(self._mol)
            AllChem.MMFFOptimizeMolecule(self._mol)

        # If xyz file is provided, preferentially use the information from
        # this over that of a pdb file
        if self._xyz_mol is not None:
            self._mol = self._xyz_mol
        elif self._pdb_mol is not None:
            self._mol = self._pdb_mol

            # Set various attributes
        self.num_bonds = len(self._mol.GetBonds())
        self.num_atoms = len(self._mol.GetAtoms())

        # Set charges
        self.set_charges()

        # Set smiles
        if self.smiles is None:
            self.smiles = Chem.MolToSmiles(self._mol)

    def set_charges(self):
        """
        Set the self.charge attribute - this depends on the charge information provided
        """

        if self.charges is not None:
            return

        elif self._mol2 is not None:
            self.charges = self._mol2.get_charges()
            return
        elif self._CalculateCharges is True:
            rdPartialCharges.ComputeGasteigerCharges(self._mol)
            gasteiger_charges = []
            for atom in self._mol.GetAtoms():
                charge = float(atom.GetProp('_GasteigerCharge'))
                gasteiger_charges.append(charge)
            self.charges = gasteiger_charges
            return
        else:
            self.get_external_charges(self._path_to_charges)
            return

    def elements(self):
        """
        Returns a list of all the elements present in the input molecule,
        in increasing atomic number order with H last
        """
        element_dict = {}
        for atom in self._mol.GetAtoms():
            if atom.GetSymbol() not in element_dict:
                element_dict[atom.GetSymbol()] = atom.GetAtomicNum()

        # Now sort dictionary into TRIPOS format - increasing atomic no.
        # with H at the end
        element_dict = dict(
            sorted(element_dict.items(), key=lambda item: item[1]))
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
        Retrieves a list of partial charges for the molecule from the list
        of DFT VEHICLE Charges args: filepath - path to file of partial
        charges returns charges - list of atomic partial charges for the
        molecule
        """
        with open(filename, 'r') as file:
            csv_reader = reader(file)
            all_charges = list(csv_reader)

        if self.index is None:
            index = 1
        else:
            index = self.index

        self.charges = all_charges[index - 1]
        self.charges = [float(charge) for charge in self.charges]

        return None

    @property
    def coords(self):
        # Set up coordinates. If xyz file is present then these coordinates
        # are used over .pdb
        if self._xyz_mol is not None:
            return self._xyz_mol.GetConformer(0).GetPositions()
        else:
            return self._mol.GetConformer(0).GetPositions()

    def molecule_block(self):
        """
        Calculates and returns, in the correct format,
        the '@<TRIPOS>MOLECULE' block for the instance of the molecule The
        number of features and sets are hard-coded to 0
        """
        if self.s_flag is True:
            molecule_block = self._mol2.get_molecule_block()

            return molecule_block

        if self._CalculateCharges is True:
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
        Computes and returns, in the correct format, the '@<TRIPOS>ATOM
        block for the instance of the molecule. This calculates partial
        charges using rdkit's implementation of the Gasteiger partial
        charges. Atomic coordinates are extracted from the pdb file.
        """
        if self.s_flag is True:
            atom_block = self._mol2.get_atom_block()

            return atom_block

        tripos_atom = pd.DataFrame(
            columns=['rdkit_index', 'atom_name', 'x_coord', 'y_coord',
                     'z_coord', 'sybyl', 'substruct', 'substruct_name',
                     'partial_charge', 'atom_symbol', 'atom_index_label'])

        mol_by_elements = self.elements_by_index()
        substruct_name = self.name

        for atom in self._mol.GetAtoms():

            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            charge = self.charges[idx]

            # Get co-ordinates for each atom
            x_coord = self.coords[idx][0]
            y_coord = self.coords[idx][1]
            z_coord = self.coords[idx][2]

            atom_index_label = mol_by_elements[symbol].index(idx) + 1
            atom_name = f'{symbol}{atom_index_label}'

            # Generate the sybyl code for each atom - the symbol, and the
            # atom type e.g. aromatic, sp3 etc.
            if symbol == 'C' or symbol == 'N':
                if atom.GetIsAromatic():
                    sybyl = f'{symbol}.ar'
                else:
                    sybyl = f'{symbol}.{atom_types[str(atom.GetHybridization())]}'
            elif symbol == 'H':
                sybyl = 'H'
            else:
                sybyl = f'{symbol}.{atom_types[str(atom.GetHybridization())]}'

            # Append atom to dataframe of atomic information
            tripos_atom = tripos_atom.append(
                {'rdkit_index': idx, 'atom_name': atom_name,
                 'x_coord': x_coord,
                 'y_coord': y_coord, 'z_coord': z_coord,
                 'sybyl': sybyl, 'substruct': NO_SUBSTRUCTS,
                 'substruct_name': substruct_name,
                 'partial_charge': "%.3f" % charge,
                 'atom_symbol': symbol, 'atom_index_label': atom_index_label},
                ignore_index=True)

        # Sort dataframe by sybyl (elements first, followed by atom index label
        tripos_atom['atom_symbol'] = pd.Categorical(tripos_atom['atom_symbol'],
                                                    self.elements())
        tripos_atom = tripos_atom.sort_values(
            by=['atom_symbol', 'atom_index_label'])
        tripos_atom.index = np.arange(1, len(tripos_atom) + 1)

        # TRIPOS atom ID and RDKit index are not the same, need to generate
        # a mapping from one to t'other
        index_mapper = pd.DataFrame(data=np.arange(1, len(tripos_atom) + 1),
                                    index=tripos_atom['rdkit_index'])
        global index_lookup
        index_lookup = index_mapper.to_dict()[0]

        # Generate final dataframe, and return as a string
        atom_block_df = tripos_atom[['atom_name', 'x_coord',
                                     'y_coord', 'z_coord', 'sybyl',
                                     'partial_charge']]

        atom_block = '@<TRIPOS>ATOM\n' + atom_block_df.to_string(
            header=False) + '\n'

        return atom_block

    def bond_block(self):

        if self.s_flag is True:
            bond_block = self._mol2.get_bond_block()

            return bond_block

        tripos_bond = pd.DataFrame(
            columns=['begin_atom_rdkit', 'end_atom_rdkit', 'bond_type'])

        for index, bond in enumerate(self._mol.GetBonds()):
            beginning = bond.GetBeginAtom().GetIdx()
            end = bond.GetEndAtom().GetIdx()

            # Get bond type and convert to TRIPOS bond type from rdkit bond
            # type
            bond_type = str(bond.GetBondType())
            bond_type = bond_types[bond_type]

            tripos_bond = tripos_bond.append(
                {'begin_atom_rdkit': beginning, 'end_atom_rdkit': end,
                 'bond_type': bond_type}, ignore_index=True)

        # Convert RDKit indices to TRIPOS atom ID numbers so there is
        # consistency between ATOM block and BOND block
        tripos_bond['begin'] = tripos_bond.apply(
            lambda row: index_lookup[row.begin_atom_rdkit], axis=1)
        tripos_bond['end'] = tripos_bond.apply(
            lambda row: index_lookup[row.end_atom_rdkit], axis=1)
        tripos_bond['bond_type'] = pd.Categorical(tripos_bond['bond_type'],
                                                  list(bond_types.values()))
        tripos_bond = tripos_bond.sort_values(by=['bond_type', 'begin'])
        tripos_bond.index = np.arange(1, len(tripos_bond) + 1)
        tripos_bond = tripos_bond[['begin', 'end', 'bond_type']]

        bond_block = '@<TRIPOS>BOND\n' + tripos_bond.to_string(
            header=False) + '\n'

        return bond_block

    def print_mol2_file(self, filename=None):
        """
        Saves the current instance of the molecule to a mol2 file in the current directory
        """

        if filename is None:
            file_name = f'{self.name}.mol2'
        elif filename[-5:] == '.mol2':
            file_name = filename
        else:
            file_name = f'{filename}.mol2'

        block = self.molecule_block() + self.atom_block() + self.bond_block()

        with open(file_name, 'w') as file:
            file.write(block)

        return None
