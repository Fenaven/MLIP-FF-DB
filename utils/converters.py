from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from openbabel import pybel

def xyz_to_smiles(fname: str) -> str:      
    mol = next(pybel.readfile("xyz", fname)) 
    smi = mol.write(format="smi") 
    return smi.split()[0].strip()

def write_xyz_file(coords: list, out_file: str) -> None:
    with open (out_file, 'w') as file:
        file.write(f'{len(coords)}\n\n')
        for coord in coords:
            if coords.index(coord) == len(coords) - 1:
                file.write(f'{coord[0]} {coord[1]} {coord[2]} {coord[3]}')
            else:
                file.write(f'{coord[0]} {coord[1]} {coord[2]} {coord[3]}\n')