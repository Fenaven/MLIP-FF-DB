from typing import List, Optional
import numpy as np
from .base_parser import BaseParser

class OrcaParser(BaseParser):
    """Парсер для файлов Orca"""
    
    def __init__(self, input_file: str = None, output_file: str = None):
        self.input_file = input_file
        self.output_file = output_file
        self._output_lines = []
        self._input_lines = []
        
        if output_file:
            with open(output_file, 'r') as f:
                self._output_lines = f.readlines()
        
        if input_file:
            with open(input_file, 'r') as f:
                self._input_lines = f.readlines()
    
    def get_energy(self) -> float:
        """Получает энергию из файла"""
        for line in self._output_lines:
            if 'FINAL SINGLE POINT ENERGY' in line:
                return float(line.split()[-1])
        return 0.0
    
    def get_forces(self) -> Optional[np.ndarray]:
        """Получает силы из файла"""
        forces = []
        reading_forces = False
        
        for i, line in enumerate(self._output_lines):
            if 'CARTESIAN GRADIENT' in line:
                reading_forces = True
                continue
            
            if reading_forces:
                if line.strip() == '':
                    break
                    
                if len(line.split()) >= 7: 
                    parts = [x for x in line.split() if x]
                    forces.append([float(x) for x in parts[4:7]])
        
        return np.array(forces) if forces else None
    
    def get_coordinates(self) -> np.ndarray:
        """Получает координаты из файла"""
        coords = []
        reading_coords = False
        
        for line in self._output_lines:
            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                reading_coords = True
                continue
            
            if reading_coords:
                if line.strip() == '':
                    break
                    
                parts = [x for x in line.split() if x]
                if len(parts) == 4:
                    coords.append([float(x) for x in parts[1:4]])
        
        return np.array(coords)
    
    def get_elements(self) -> List[str]:
        """Получает список элементов"""
        elements = []
        reading_coords = False
        
        for line in self._output_lines:
            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                reading_coords = True
                continue
            
            if reading_coords:
                if line.strip() == '':
                    break
                    
                parts = [x for x in line.split() if x]
                if len(parts) == 4:
                    elements.append(parts[0])
        
        return elements
    
    def get_cell(self) -> Optional[np.ndarray]:
        """Получает параметры ячейки"""
        return None
    
    def get_scf_convergence(self) -> bool:
        """Проверяет сходимость SCF"""
        for line in self._output_lines:
            if 'SCF CONVERGED AFTER' in line:
                return True
        return False 