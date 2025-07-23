import sqlite3
import pandas as pd
from typing import List, Optional, Dict, Any
import numpy as np

def get_connection():
    """Создает и возвращает соединение с базой данных"""
    return sqlite3.connect('main.db')

def upload_calculation(calc_id: int, structure_id: int, total_energy: float, input_hash: str, program_id: str, 
                      date: str, server_name: str, scf_converged: bool, input_file_path: str, log_file_path: str, 
                      atomic_file_path: str) -> None:
    """Добавление данных в таблицу расчетов"""
    connection = get_connection()
    cursor = connection.cursor()
    
    insert_prompt = """
    INSERT INTO computations (calc_id, structure_id, total_energy, input_hash, program_id, date, server_name, 
                            scf_converged, input_file, log_file, atomic_file)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    
    try:
        cursor.execute(insert_prompt, (calc_id, structure_id, total_energy, input_hash, program_id,
                                     date, server_name, scf_converged, input_file_path, 
                                     log_file_path, atomic_file_path))
        connection.commit()
    finally:
        connection.close()

def upload_atomic_properties(calc_id: int, atom_number: int, grade: float, f_x: float, f_y: float, f_z: float) -> None:
    """Добавление данных в таблицу атомных свойств"""
    connection = get_connection()
    cursor = connection.cursor()
    
    atom_id = int(f'{calc_id}' + f'{atom_number}')
    
    insert_prompt = """
    INSERT INTO atomic_properties (calc_id, atom_id, atom_number, grade, F_x, F_y, F_z)
    VALUES (?, ?, ?, ?, ?, ?, ?)
    """
    
    try:
        cursor.execute(insert_prompt, (calc_id, atom_id, atom_number, grade, f_x, f_y, f_z))
        connection.commit()
    finally:
        connection.close()

def upload_atoms(str_id: int, atom_number: int, element: str, x_coord: float, y_coord: float, z_coord: float) -> None:
    """Добавление данных в таблицу атомов"""
    connection = get_connection()
    cursor = connection.cursor()
    
    insert_prompt = """
    INSERT INTO atoms (str_id, atom_number, element, x_coord, y_coord, z_coord)
    VALUES (?, ?, ?, ?, ?, ?)
    """
    
    try:
        cursor.execute(insert_prompt, (str_id, atom_number, element, x_coord, y_coord, z_coord))
        connection.commit()
    finally:
        connection.close()

def get_table(table_type: str) -> pd.DataFrame:
    """Получение содержимого таблицы"""
    connection = get_connection()
    cursor = connection.cursor()
    
    try:
        cursor.execute(f"SELECT * FROM {table_type}")
        data = cursor.fetchall()
        column_names = [desc[0] for desc in cursor.description]
        return pd.DataFrame(data, columns=column_names)
    finally:
        connection.close()

def get_last_id(table: str, id_column: str) -> int:
    """Получение последнего ID из таблицы"""
    connection = get_connection()
    cursor = connection.cursor()
    
    try:
        cursor.execute(f'SELECT {id_column} FROM {table} ORDER BY rowid DESC LIMIT 1')
        result = cursor.fetchone()
        return result[0] if result else 0
    finally:
        connection.close()

def get_last_calcid() -> int:
    return get_last_id('computations', 'calc_id')

def get_last_strid() -> int:
    return get_last_id('computations', 'structure_id')

def get_last_molid() -> int:
    return get_last_id('molecules', 'mol_id')

def compare_structures(str_id1: int, str_id2: int) -> float:
    """
    Сравнение двух структур по RMSD
    
    Args:
        str_id1 (int): ID первой структуры
        str_id2 (int): ID второй структуры
        
    Returns:
        float: Значение RMSD
    """
    connection = get_connection()
    
    # Получаем координаты атомов для обеих структур
    query = """
    SELECT element, x_coord, y_coord, z_coord
    FROM atoms
    WHERE str_id = ?
    ORDER BY atom_number
    """
    
    try:
        struct1 = pd.read_sql_query(query, connection, params=[str_id1])
        struct2 = pd.read_sql_query(query, connection, params=[str_id2])
        
        if len(struct1) != len(struct2):
            return float('inf')
        
        if not (struct1['element'] == struct2['element']).all():
            return float('inf')
        
        # Вычисляем RMSD
        coords1 = struct1[['x_coord', 'y_coord', 'z_coord']].values
        coords2 = struct2[['x_coord', 'y_coord', 'z_coord']].values
        
        diff = coords1 - coords2
        rmsd = np.sqrt(np.mean(np.sum(diff * diff, axis=1)))
        return rmsd
    finally:
        connection.close()

def find_similar_structures(str_id: int, rmsd_threshold: float = 0.1) -> List[int]:
    """
    Поиск похожих структур по RMSD
    
    Args:
        str_id (int): ID структуры для поиска
        rmsd_threshold (float): Пороговое значение RMSD
        
    Returns:
        List[int]: Список ID похожих структур
    """
    connection = get_connection()
    cursor = connection.cursor()
    
    try:
        cursor.execute("SELECT DISTINCT str_id FROM atoms")
        all_structures = [row[0] for row in cursor.fetchall()]
        
        similar_structures = []
        for other_id in all_structures:
            if other_id != str_id:
                rmsd = compare_structures(str_id, other_id)
                if rmsd <= rmsd_threshold:
                    similar_structures.append(other_id)
        
        return similar_structures
    finally:
        connection.close() 