from typing import List, Optional
import numpy as np
from .base_parser import BaseParser

class XYZParser(BaseParser):
    """Парсер для XYZ файлов"""
    
    def __init__(self, xyz_file: str):
        self.xyz_file = xyz_file
        self._data = self._read_xyz()
    
    def _read_xyz(self):
        """Читает XYZ файл"""
        with open(self.xyz_file, 'r') as f:
            lines = f.readlines()
        
        n_atoms = int(lines[0])
        comment = lines[1].strip()
        
        data = {'comment': comment, 'energy': None}
        if 'energy' in comment.lower():
            try:
                energy = float(comment.split('=')[1].strip())
                data['energy'] = energy
            except:
                pass
        
        coords = []
        elements = []
        forces = []
        has_forces = False
        
        for line in lines[2:2+n_atoms]:
            parts = line.split()
            elements.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
            if len(parts) >= 7:
                has_forces = True
                forces.append([float(x) for x in parts[4:7]])
        
        data['coordinates'] = np.array(coords)
        data['elements'] = elements
        if has_forces:
            data['forces'] = np.array(forces)
        
        return data
    
    def get_energy(self) -> float:
        """Получает энергию из комментария"""
        return self._data.get('energy', 0.0)
    
    def get_forces(self) -> Optional[np.ndarray]:
        """Получает силы из файла"""
        return self._data.get('forces')
    
    def get_coordinates(self) -> np.ndarray:
        """Получает координаты из файла"""
        return self._data['coordinates']
    
    def get_elements(self) -> List[str]:
        """Получает список элементов"""
        return self._data['elements']
    
    def get_cell(self) -> Optional[np.ndarray]:
        """не содержит параметры ячейки"""
        return None 