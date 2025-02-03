import hashlib
import sqlite3
import pandas as pd
import cclib
from rdkit import Chem
from openbabel import pybel


from rdkit import Chem


mol = pybel.readstring('smi', '[NH3]CCCC(=O)[O]')
mol = mol.write('can')

print(mol)
# Create a molecule from a SMILES string
molecule = Chem.MolFromSmiles('[NH3]CCCC(=O)[O]', sanitize=False)
print('molecule ', molecule)

# test = Chem.MolToSmiles(Chem.MolFromSmiles('[NH3]CCCC(=O)[O]'), True)
# print(test)

# Define the substructure using a SMARTS pattern
substructure = Chem.MolFromSmiles('CC')

# Find the substructure in the molecule
matches = molecule.HasSubstructMatch(substructure)

# Print the indices of the matching atoms
print("Matching atom indices:", matches)

# If you want to highlight the substructure in the molecule
# from rdkit.Chem import Draw

# # Create a highlighted image of the molecule
# highlighted_image = Draw.MolToImage(molecule, highlightAtoms=matches)
# highlighted_image.show()