import utils
import cclib

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from openbabel import pybel

class Parser:
    def __init__(self, input_file_path: str, log_file_path: str) -> None:
        """_summary_

        Args:
            input_file_path (str): _description_
            log_file_path (str): _description_
        """
        with open(input_file_path, 'r') as file:
            self.input = file.read()
        
        with open(log_file_path, 'r') as file:
            self.log_file = file.readlines()

    def get_input_hash(self) -> str:
        return utils.get_input_hash(self.input)

    def parse_total_energy(self) -> float:
        pass

    def parse_scf_convergence(self) -> bool:
        pass

    def get_summary(self) -> list:
        prop_lst = []

        return prop_lst


class ParserORCA(Parser):
    def __init__(self, input_file_path, log_file_path):
        super().__init__(input_file_path, log_file_path)

    def __parse_total_energy(self) -> float:
        for s in self.log_file:
            if 'FINAL SINGLE POINT ENERGY' in s:
                s_split_temp = s.split(' ')
                break
        
        energy = float(s_split_temp[-1])

        return energy
    
    def __parse_scf_convergence(self) -> bool:
        for s in self.log_file:
            if 'SCF CONVERGED AFTER' in s:
                scf_conv = True
                break

            else:
                scf_conv = False

        return scf_conv
    
    def __parse_forces(self):
        count = 0
        for s in self.log_file:
            if 'CARTESIAN GRADIENT' in s:
                cartesian_index = count
                break
            count += 1

        s = self.log_file[cartesian_index + 3]
        count_test = 0

        forces_lst = []

        while s != '\n':
            forces_splitted_lst = s.split(' ')

            forces_lst_tmp = []

            for s_force in forces_splitted_lst[4:]:
                try:
                    force = float(s_force.replace('\n', ''))
                    forces_lst_tmp.append(force)
                except:
                    continue
            count_test += 1
            s = self.log_file[cartesian_index + 3 + count_test]

            forces_lst.append(forces_lst_tmp)
            
        return forces_lst

    def __parse_atomic_coordinates(self):
        
        for s in self.log_file:
            if 'CARTESIAN COORDINATES (ANGSTROEM)' in s:
                coord_start_index = self.log_file.index(s)
        
        s = self.log_file[coord_start_index + 2]
        count = 0

        coords_lst = []

        while s != '\n':
            coords_lst_tmp = []
            coordinate_splitted_lst = s.split(' ')

            while '' in coordinate_splitted_lst:
                coordinate_splitted_lst.remove('')

            coordinate_splitted_lst[-1] = coordinate_splitted_lst[-1].replace('\n', '')

            for i in range(0, 4):
                if i == 0:
                    coords_lst_tmp.append(coordinate_splitted_lst[i])
                else:
                    coords_lst_tmp.append(float(coordinate_splitted_lst[i]))
            
            coords_lst.append(coords_lst_tmp)


            # print(coordinate_splitted_lst)
            count += 1
            s = self.log_file[coord_start_index + 2 + count]

        return coords_lst

    def get_summary(self) -> dict:
        input_hash = self.get_input_hash()

        energy = self.__parse_total_energy()
        scf_conv = self.__parse_scf_convergence()

        forces = self.__parse_forces()
        coords = self.__parse_atomic_coordinates()

        data_dict = {'input_hash': input_hash, 'total_energy': energy, 'scf_convergence': scf_conv, 'forces': forces,
                     'coords': coords}

        return data_dict


class ParserCFG:
    def __init__(self, file) -> None:
        pass

    def parse_energy(self) -> float:
        pass

    def check_scf(self) -> float:
        pass


def main():
    input_file_path = 'data/input.inp'
    log_file_path = 'data/butane_cut.log'

    parser = ParserORCA(input_file_path, log_file_path)

    coords = parser.get_summary()['coords']
    # conda print(coords)

    with open('test_structure.xyz', 'w') as file:
        file.write(f'{len(coords)}\n')
        file.write('\n')
        for coord in coords:
            if coords.index(coord) == len(coords) - 1:
                file.write(coord[0] + f' {coord[1]} {coord[2]} {coord[3]}')
            else:
                file.write(coord[0] + f' {coord[1]} {coord[2]} {coord[3]}\n')


    def xyz_to_smiles(fname: str) -> str: 
     
        mol = next(pybel.readfile("xyz", fname)) 
 
        smi = mol.write(format="smi") 
 
        return smi.split()[0].strip()

    smi = xyz_to_smiles("test_structure.xyz") 
    print(smi)





if __name__ == '__main__':
    main()