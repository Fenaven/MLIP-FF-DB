from typing import List, Optional
import numpy as np
from .base_parser import BaseParser

class CFGParser(BaseParser):
    """Парсер для CFG файлов"""
    
    def __init__(self, cfg_file: str):
        self.cfg_file = cfg_file
        self._data = self._read_cfg()
    
    def _read_cfg(self):
        """Читает CFG файл"""
        with open(self.cfg_file, 'r') as f:
            lines = f.readlines()
        
        data = {
            'energy': None,
            'coordinates': [],
            'elements': [],
            'forces': [],
            'cell': np.zeros((3, 3)),
            'has_forces': False
        }
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            if line.startswith('Size'):
                n_atoms = int(lines[i+1])
                i += 2
            elif line.startswith('Supercell'):
                # Читаем векторы ячейки
                for j in range(3):
                    parts = lines[i+1+j].split()
                    data['cell'][j] = [float(x) for x in parts]
                i += 4
            elif line.startswith('AtomData:'):
                # Определяем формат данных
                headers = line.split(':')[1].strip().split()
                col_indices = {
                    'type': headers.index('type'),
                    'x': headers.index('cartes_x'),
                    'y': headers.index('cartes_y'),
                    'z': headers.index('cartes_z')
                }
                if 'fx' in headers and 'fy' in headers and 'fz' in headers:
                    col_indices.update({
                        'fx': headers.index('fx'),
                        'fy': headers.index('fy'),
                        'fz': headers.index('fz')
                    })
                    data['has_forces'] = True
                
                # Читаем данные атомов
                i += 1
                while i < len(lines) and lines[i].strip():
                    parts = lines[i].split()
                    if len(parts) >= 5:  # id type x y z [fx fy fz]
                        data['elements'].append(parts[col_indices['type']])
                        coords = [float(parts[col_indices[x]]) for x in ['x', 'y', 'z']]
                        data['coordinates'].append(coords)
                        
                        if data['has_forces']:
                            forces = [float(parts[col_indices[x]]) for x in ['fx', 'fy', 'fz']]
                            data['forces'].append(forces)
                    i += 1
            elif line.startswith('Energy'):
                data['energy'] = float(lines[i+1])
                i += 2
            else:
                i += 1
        
        data['coordinates'] = np.array(data['coordinates'])
        if data['has_forces']:
            data['forces'] = np.array(data['forces'])
        
        return data
    
    def get_energy(self) -> float:
        """Получает энергию из файла"""
        return self._data.get('energy', 0.0)
    
    def get_forces(self) -> Optional[np.ndarray]:
        """Получает силы из файла"""
        return self._data['forces'] if self._data['has_forces'] else None
    
    def get_coordinates(self) -> np.ndarray:
        """Получает координаты из файла"""
        return self._data['coordinates']
    
    def get_elements(self) -> List[str]:
        """Получает список элементов"""
        return self._data['elements']
    
    def get_cell(self) -> Optional[np.ndarray]:
        """Получает параметры ячейки"""
        return self._data['cell'] 