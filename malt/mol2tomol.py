import csv
import json
import pkg_resources
from rdkit import Chem


#Set path to package data
MOL2_DICT = pkg_resources.resource_filename('malt', 'Data/mol2_dictionary.json')
VEHICLE_MOL2 = pkg_resources.resource_filename('malt', 'Data/vehicle_dft.mol2')

class mol2:

    def __init__(self, mol_id, *args):
        self.mol = None
        self.charges = []
        self.mol_id = mol_id
        self.mol_block = None
        
        #Read in the position dictionary
        with open(MOL2_DICT, 'r') as file:
            self.lookup = json.load(file)

        #Read in the VEHICLE Mol2 file
        with open(VEHICLE_MOL2, 'r') as file:
            csv_reader = csv.reader(file)
            self.vehicle = list(csv_reader)

    def get_block_lists(self):
        """
        Returns the mol2 block for the molecule with ID mol_id as a list of lists - useful for isolating charges and 
        other information
        """
        beg = self.lookup[self.mol_id][0]
        end = self.lookup[self.mol_id][1]
        
        block_list = self.vehicle[beg:end]

        return block_list



    def get_block(self):
        """
        Reads in molecule block for molecule with id mol_id from full VEHICLE mol2 file
        """
        block = ''
        block_lists = self.get_block_lists()

        for list in block_lists:
            block += list[0] + '\n'

        return block

    def mol_from_mol2(self):
        """
        Returns an RDKit Mol object for an instance of the mol2 class
        """

        block = self.get_block()

        mol = Chem.MolFromMol2Block(block, sanitize=False, removeHs=False)

        return mol 

    def get_charges(self):
        """
        Returns a list of the atomic charges, ordered by atom index, for the instance of the class
        """

        charges = []
        block_lists = self.get_block_lists()

        for idx, list in enumerate(block_lists):
            if idx > 5:
                if '@<TRIPOS>BOND' in list[0]:
                    break
                else:
                    list = list[0].split()
                    charge = float(list[-1])
                    charges.append(charge)
        
        return charges
