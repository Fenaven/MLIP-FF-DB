from typing import Dict, List, Optional
from .base_parser import BaseParser
import numpy as np

class LammpsParser(BaseParser):
    def __init__(self, log_file_path: str, dump_file_path: str):
        super().__init__(log_file_path)
        self.dump_file_path = dump_file_path
        self.dump_data = self._read_dump_file()
        self.log_data = self._read_log_file()

    def _read_dump_file(self) -> Dict:
        """Чтение и парсинг LAMMPS dump файла"""
        data = {
            'timestep': None,
            'num_atoms': None,
            'box_bounds': None,
            'atoms': []
        }

        try:
            with open(self.dump_file_path, 'r') as f:
                lines = f.readlines()

            i = 0
            while i < len(lines):
                line = lines[i].strip()

                if line == 'ITEM: TIMESTEP':
                    i += 1
                    if i < len(lines):
                        try:
                            data['timestep'] = int(lines[i].strip())
                        except ValueError:
                            print(f"Ошибка при чтении timestep")
                
                elif line == 'ITEM: NUMBER OF ATOMS':
                    i += 1
                    if i < len(lines):
                        try:
                            data['num_atoms'] = int(lines[i].strip())
                        except ValueError:
                            print(f"Ошибка при чтении числа атомов")
                
                elif line.startswith('ITEM: BOX BOUNDS'):
                    data['box_bounds'] = []
                    i += 1
                    try:
                        for _ in range(3):  # x, y, z bounds
                            if i < len(lines):
                                values = lines[i].strip().split()
                                if len(values) >= 2:
                                    bounds = [float(values[0]), float(values[1])]
                                    data['box_bounds'].append(bounds)
                                i += 1
                    except (ValueError, IndexError) as e:
                        print(f"Ошибка при чтении границ ячейки: {e}")
                
                elif line.startswith('ITEM: ATOMS'):
                    headers = line.split()[2:]
                    i += 1
                    while i < len(lines):
                        line = lines[i].strip()
                        if not line:
                            break
                        
                        try:
                            values = line.split()
                            if len(values) == len(headers):
                                atom_dict = {}
                                for header, value in zip(headers, values):
                                    atom_dict[header] = value
                                data['atoms'].append(atom_dict)
                        except Exception as e:
                            print(f"Ошибка при чтении атома: {e}")
                        i += 1
                        
                i += 1

            if not data['atoms']:
                raise ValueError("Не найдены данные об атомах в dump файле")

        except Exception as e:
            print(f"Ошибка при чтении dump файла: {e}")
            raise

        return data

    def _read_log_file(self) -> Dict:
        """Чтение и парсинг LAMMPS log файла"""
        data = {
            'energy': 0.0,
            'forces': []
        }

        try:
            with open(self.log_file_path, 'r') as f:
                lines = f.readlines()

            energy_col = None
            headers = None
            
            for i in range(len(lines) - 1, -1, -1):
                line = lines[i].strip()
                
                if 'Step' in line and ('PotEng' in line or 'E_pair' in line):
                    headers = line.split()
                    if 'PotEng' in headers:
                        energy_col = headers.index('PotEng')
                    elif 'E_pair' in headers:
                        energy_col = headers.index('E_pair')
                    break

            if energy_col is not None:
                for j in range(i + 1, len(lines)):
                    line = lines[j].strip()
                    if line and not line.startswith('Loop'):
                        try:
                            values = line.split()
                            if len(values) > energy_col:
                                data['energy'] = float(values[energy_col])
                                break
                        except (ValueError, IndexError) as e:
                            print(f"Ошибка при чтении энергии: {e}")

        except Exception as e:
            print(f"Ошибка при чтении log файла: {e}")
            raise

        return data

    def get_data(self) -> Dict:
        """Получение всех данных из файлов"""
        try:
            coordinates = []
            elements = []
            
            for atom in self.dump_data['atoms']:
                try:
                    x = float(atom.get('x', 0))
                    y = float(atom.get('y', 0))
                    z = float(atom.get('z', 0))
                    coordinates.append([x, y, z])
                    elements.append(str(atom.get('type', '1')))
                except (ValueError, KeyError) as e:
                    print(f"Ошибка при обработке атома: {e}")
                    continue

            if not coordinates or not elements:
                raise ValueError("Не удалось получить координаты или типы атомов")

            return {
                'energy': self.log_data['energy'],
                'coordinates': coordinates,
                'elements': elements,
                'forces': None
            }

        except Exception as e:
            print(f"Ошибка при получении данных: {e}")
            raise

    def get_scf_convergence(self) -> bool:
        return True 