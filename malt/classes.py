"""
A module for generating key properties of molecules needed to generate
molecular file formats, specifically aimed at Tripos .mol2 files.
"""

import numpy as np
import pandas as pd
from csv import reader
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

# Some user-defined variables - at the moment these are hard-coded,
# I may find a better way of doing it later
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


def open_file(filename):
    """
    opens file and returns csv_reader object
    :param filename: path of file to be opened
    :return: csv_reader object of file contents
    """
    with open(filename, 'r') as file:
        csv_reader = reader(file)
        output = list(csv_reader)

        return output


class Molecule:

    def __init__(self, smiles, charge_path=None, xyz_path=None,
                 **kwargs):
        """
        For each molecule there should be a smiles string to extract atomic
        connectivity information, a .csv file containing the atomic charges,
        and an .xyz file containing the atomic coordinates. The .xyz and
        .csv files must be indexed according to the atomic indices of the
        smiles string provided.

        :param smiles: SMILES string
        :param charge_path: Path  to charges - these
        should be .csv files
        :param xyz_path: Path to atomic coordinates -
        these must be .xyz files
        :param kwargs: any other keyword arguments - included for possible
        generalisable work later on
        """
        self.smiles = smiles
        self.charge_path = charge_path
        self.xyz_path = xyz_path
        self.index_lookup = None

        self.mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(self.mol)

        if self.charge_path is not None:
            self.charges = self.get_external_charges()
        else:
            self.charges = self.get_xyz_charges()

        self.coordinates = self.get_xyz_coordinates()

        self.set_coordinates()

    def molecule_block(self):
        """
        Calculates and returns, in the correct format, the '@<TRIPOS>MOLECULE'
        block for the instance of the molecule.
        The number of features and sets are hard-coded to 0
        """
        charge_type = 'DFT'

        molecule_block = (f'@<TRIPOS>MOLECULE\n'
                          f'{self.smiles}\n'
                          f'{self.num_atoms} {self.num_bonds} {NO_SUBSTRUCTS} 0 0\n '
                          f'{MOLECULE_TYPE}\n'
                          f'{charge_type}\n')

        return molecule_block

    def atom_block(self):
        """
        Computes and returns, in the correct format, the '@<TRIPOS>ATOM
        block for the instance of the molecule.
        """

        tripos_atom = pd.DataFrame(
            columns=['rdkit_index', 'atom_name', 'x_coord', 'y_coord',
                     'z_coord', 'sybyl', 'substruct', 'substruct_name',
                     'partial_charge', 'atom_symbol', 'atom_index_label'])

        substruct_name = self.smiles

        for atom in self.mol.GetAtoms():
            atom_idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            charge = self.charges[atom_idx]

            # Get co-ordinates for each atom
            atom_coords = self.coordinates[atom_idx]
            x, y, z = self.get_coordinates_from_list(atom_coords)

            atom_index_label = self.get_atom_index_label(symbol, atom_idx)
            atom_name = self.get_atom_name(symbol, atom_index_label)

            sybyl = self.generate_sybyl_code(atom)

            # Append atom to dataframe of atomic information
            tripos_atom = tripos_atom.append(
                {'rdkit_index': atom_idx, 'atom_name': atom_name,
                 'x_coord': x,
                 'y_coord': y, 'z_coord': z,
                 'sybyl': sybyl, 'substruct': NO_SUBSTRUCTS,
                 'substruct_name': substruct_name,
                 'partial_charge': "%.3f" % charge,
                 'atom_symbol': symbol, 'atom_index_label': atom_index_label},
                ignore_index=True)

        tripos_atom = self.sort_atom_dataframe(tripos_atom)

        # TRIPOS atom ID and RDKit index are not the same, need to generate
        # a mapping from one to t'other
        self.index_lookup = self.generate_index_mapper(tripos_atom)

        # Generate final dataframe, and return as a string
        atom_block_df = tripos_atom[['atom_name', 'x_coord',
                                     'y_coord', 'z_coord', 'sybyl',
                                     'partial_charge']]

        atom_block = '@<TRIPOS>ATOM\n' + atom_block_df.to_string(
            header=False) + '\n'

        return atom_block

    def bond_block(self):
        """
        Generates and returns the TRIPOS Bond Block for the instance of the
        molecule - ATOM Block must have been generated first for this to work
        :return: string - TRIPOS bond block for molecule
        """

        # atom block must have been created for bond block to be created
        if self.index_lookup is None:
            print('Atom Block must be generated before Bond Block')
            return

        bond_block = self.get_bond_dataframe(self.mol)

        # Convert RDKit indices to TRIPOS atom ID numbers so there is
        # consistency between ATOM block and BOND block
        sorted_block = self.sort_bond_dataframe(bond_block)

        bond_block = '@<TRIPOS>BOND\n' + sorted_block.to_string(
            header=False) + '\n'

        return bond_block

    def print_mol2_file(self, filename=None):
        """
        Saves the current instance of the molecule to a mol2 file in the
        current directory
        """

        if filename is None:
            file_name = f'{self.smiles}.mol2'
        elif filename[-5:] == '.mol2':
            file_name = filename
        else:
            file_name = f'{filename}.mol2'

        block = self.molecule_block() + self.atom_block() + self.bond_block()

        with open(file_name, 'w') as file:
            file.write(block)

        return None

    @property
    def num_atoms(self):
        return len(self.mol.GetAtoms())

    @property
    def num_bonds(self):
        return len(self.mol.GetBonds())

    def set_coordinates(self):
        """
        Assigns the coordinates extracted from the .xyz file to each atom in
        the RDKit mol object
        :return: None
        """
        conformer = self.mol.GetConformer(0)

        for atom in self.mol.GetAtoms():
            atom_idx = atom.GetIdx()
            atom_coords = self.coordinates[atom_idx]

            x, y, z = self.get_coordinates_from_list(atom_coords)

            conformer.SetAtomPosition(atom_idx, Point3D(x, y, z))

        return None

    @staticmethod
    def get_coordinates_from_list(coords_list):
        """
        extracts the x,y,z coordinates from a list generated by
        self.get_xyz_coordinates
        :param coords_list: list of coords generated by
        self.get_xyz_coordinates
        :return: x, y and z co-ordinates
        """

        x = float(coords_list[1])
        y = float(coords_list[2])
        z = float(coords_list[3])

        return x, y, z

    def get_xyz_charges(self):
        """
        If the charges are included in the title lines of the xyz file,
        this function extracts them. It is called if a path to charges is
        not provided.
        :return: list of charges, indexed by atom ID
        """
        xyz_file = self.xyz_path
        all_xyz_data = open_file(xyz_file)
        charges = all_xyz_data[1]

        # This bit is a workaround for a bug that shouldn't exist where a
        # list is printed to the xyz file including the square brackets
        if charges[0][0] == '[':
            charges[0] = charges[0][1:]
            charges[-1] = charges[-1][:-1]

        charges = [float(charge) for charge in charges]

        return charges

    def get_xyz_coordinates(self):
        """
        Extract coordinates from .xyz file to assign to atoms instantiated
        from smiles string
        :return:list of lists of atomic symbol and x, y, z coordinates.
        Lists in order of atom index
        """
        xyz_file = self.xyz_path
        all_coords = open_file(xyz_file)

        # skip the first two lines of all_coords as it doesn't contain
        # coordinate information
        coordinates = [entry[0].split() for entry in all_coords[2:]]

        return coordinates

    def get_external_charges(self):
        """
        Retrieves a list of partial charges for the molecule from the list
        of DFT VEHICLE Charges args: filepath - path to file of partial
        charges returns charges - list of atomic partial charges for the
        molecule
        """
        charge_file = self.charge_path
        all_charges = open_file(charge_file)

        charges = [float(charge) for charge in all_charges[0]]

        return charges

    @staticmethod
    def sort_elements(element_dict):
        """
        Sorts a dictionary of elements into TRIPOS format - increasing
        atomic number with H at the end
        :param element_dict: dictionary keys = elements,
        values = atomic numbers - values can be discarded
        :return: list sorted into Tripos format
        """
        element_dict = dict((sorted(element_dict.items(),
                                    key=lambda item: item[1])))
        elements = [key for key in element_dict]

        if 'H' in elements:
            elements.remove('H')
            elements.append('H')

        return elements

    def elements(self):
        """
        Returns a list of all the elements present in the input molecule,
        in increasing atomic number order with H last
        """
        element_dict = {}
        for atom in self.mol.GetAtoms():
            if atom.GetSymbol() not in element_dict:
                element_dict[atom.GetSymbol()] = atom.GetAtomicNum()

        elements = self.sort_elements(element_dict)

        return elements

    def elements_by_index(self):
        """
        Returns a dictionary with element symbols as keys, and the atom ID's
        of those elements in the molecule as values
        """
        mol_by_elements = {}
        for atom in self.mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol not in mol_by_elements:
                mol_by_elements[symbol] = []
            mol_by_elements[symbol].append(atom.GetIdx())

        return mol_by_elements

    @staticmethod
    def generate_sybyl_code(atom):
        """
        Generate the sybyl code for each atom - comprised of the symbol and
        atom type e.g. aromatic/sp3 etc.
        :param atom: instance of RDKit atom object in RDKit mol object
        :return: string of sybyl code for atom
        """
        symbol = atom.GetSymbol()

        if symbol == 'C' or symbol == 'N':
            if atom.GetIsAromatic():
                sybyl = f'{symbol}.ar'
            else:
                sybyl = f'{symbol}.{atom_types[str(atom.GetHybridization())]}'
        elif symbol == 'H':
            sybyl = 'H'
        else:
            sybyl = f'{symbol}.{atom_types[str(atom.GetHybridization())]}'

        return sybyl

    def get_atom_index_label(self, symbol, atom_index):
        """
        Generates the TRIPOS atom index label
        :param symbol: symbol of RDKit atom
        :param atom_index: index of atom in RDKit mol object
        :return: string = TRIPOS atom name
        """
        mol_by_elements = self.elements_by_index()

        return mol_by_elements[symbol].index(atom_index) + 1

    @staticmethod
    def get_atom_name(symbol, index_label):
        """
        Returns the TRIPOS atom name
        :param symbol: atom symbol
        :param index_label: TRIPOS atom index label
        :return: name of atom
        """
        return f'{symbol}{index_label}'

    def sort_atom_dataframe(self, dataframe):
        """
        Sorts TRIPOS atom dataframe by sybyl - elements first, followed by
        atom index label
        :param dataframe: pandas dataframe of TRIPOS atom information
        :return: pandas dataframe sorted by sybyl code
        """
        dataframe['atom_symbol'] = pd.Categorical(dataframe['atom_symbol'],
                                                  self.elements())

        dataframe = dataframe.sort_values(by=['atom_symbol',
                                              'atom_index_label'])
        dataframe.index = np.arange(1, len(dataframe) + 1)

        return dataframe

    @staticmethod
    def generate_index_mapper(dataframe):
        """
        TRIPOS Atom ID's and RDKit atom indices are not the same - a mapping is
        needed between the rdkit indices and the TRIPOS indices to generate
        the bond block.
        :param dataframe: pandas dataframe containing the information for
        the TRIPOS atom block
        :return: dictionary mapping rdkit atom ID's to TRIPOS atom ID's
        """
        index_mapper = pd.DataFrame(data=np.arange(1, len(dataframe) + 1),
                                    index=dataframe['rdkit_index'])

        index_lookup = index_mapper.to_dict()[0]

        return index_lookup

    @staticmethod
    def get_bond_dataframe(mol_object):
        """
        Generates a pandas dataframe containing the information necessary to
        construct the TRIPOS bond block
        :param mol_object: RDKit mol object
        :return: pandas dataframe of bond indices (beginning and end),
        and TRIPOS bond type
        """

        tripos_bond = pd.DataFrame(
            columns=['begin_atom_rdkit', 'end_atom_rdkit', 'bond_type'])

        for index, bond in enumerate(mol_object.GetBonds()):
            beginning = bond.GetBeginAtom().GetIdx()
            end = bond.GetEndAtom().GetIdx()

            # Get rdkit bond type and convert to TRIPOS bond type
            bond_type = str(bond.GetBondType())
            bond_type = bond_types[bond_type]

            tripos_bond = tripos_bond.append(
                {'begin_atom_rdkit': beginning, 'end_atom_rdkit': end,
                 'bond_type': bond_type}, ignore_index=True)

        return tripos_bond

    def sort_bond_dataframe(self, dataframe):
        """
        Converts RDKit indices in dataframe into TRIPOS atom ID numbers so
        atom ID's are consistent between atom and bond block. Then sorts by
        bond ID
        :param dataframe: pandas dataframe containing bond block information
        :return: pandas dataframe sorted by bond ID, with atom ID's replaced
        with TRIPOS ID's
        """
        # Replace beginning and end atom indices with TRIPOS atom indices
        # to match those in ATOM Block
        dataframe['begin'] = dataframe.apply(
            lambda row: self.index_lookup[row.begin_atom_rdkit], axis=1)
        dataframe['end'] = dataframe.apply(
            lambda row: self.index_lookup[row.end_atom_rdkit], axis=1)

        # Sort by bond type
        dataframe['bond_type'] = pd.Categorical(dataframe['bond_type'],
                                                list(bond_types.values()))

        tripos_bond = dataframe.sort_values(by=['bond_type', 'begin'])
        tripos_bond.index = np.arange(1, len(tripos_bond) + 1)
        tripos_bond = tripos_bond[['begin', 'end', 'bond_type']]

        return tripos_bond
