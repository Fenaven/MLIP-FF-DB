from abc import ABC, abstractmethod
from typing import Dict, Any, List, Optional
import numpy as np

class BaseParser(ABC):
    """Базовый класс для всех парсеров"""
    
    @abstractmethod
    def get_energy(self) -> float:
        """Получает энергию из файла"""
        pass
    
    @abstractmethod
    def get_forces(self) -> Optional[np.ndarray]:
        """Получает силы из файла"""
        pass
    
    @abstractmethod
    def get_coordinates(self) -> np.ndarray:
        """Получает координаты из файла"""
        pass
    
    @abstractmethod
    def get_elements(self) -> List[str]:
        """Получает список элементов"""
        pass
    
    @abstractmethod
    def get_cell(self) -> Optional[np.ndarray]:
        """Получает параметры ячейки"""
        pass
    
    def get_data(self) -> Dict[str, Any]:
        """Получает все данные из файла"""
        data = {
            'energy': self.get_energy(),
            'forces': self.get_forces(),
            'coordinates': self.get_coordinates(),
            'elements': self.get_elements(),
            'cell': self.get_cell()
        }
        return data 