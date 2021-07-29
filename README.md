# malt
---
Malt is a small package for converting .pdb molecular geometry files into TRIPOS .mol2 files, crucially including atomic partial charges as a fourth coordinate for each constituent atom. It is limited in its scope, and is likely to be updated  only if useful to my research. It is free to use for all that may find it useful.

## Installation
Currently this Repo is the only source of this code. To install it, clone the repository into the desired directory on your local machine using 

`git clone https://github.com/matthewtoholland/malt.git`

Then create and activate a new virtual environment. Using Conda this would be
```
conda env create -f environment.yml
conda activate malt
```

Finally install the package to your local environment using
`python setup.py install`
from the command line. 

## Usage

The process of generating a molecular block in the TRIPOS MOL2 format and saving it to a .mol2 file have been separated to allow users the flexibility to save multiple molecules to one file.

### Single Molecule

For a single molecule in a pdb file

```
from malt import methods
from malt.methods import molecule_block, write_mol2

block = molecule_block('/path/to/pdb/file.pdb')
write_mol2(block, 'file_as_mol2.mol2')
```

### Multiple Molecules

For multiple molecules the easiest way is to use a loop, adding each molecule block to a variable, which is then saved to a mol2 outside of the loop.