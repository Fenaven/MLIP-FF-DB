import streamlit as st
import sqlite3
import pandas as pd
import utils
import numpy as np
import os
from src.parsers import OrcaParser, LammpsParser, XYZParser, CFGParser
from openbabel import pybel
from rdkit import Chem

from streamlit_ketcher import st_ketcher
from src.ui.search_tab import render_search_tab


software_lst = ['Orca', 'LAMMPS', 'VASP']

software_version_dict = {'Orca': ['5.0.0', '6.0.0'],
                         'LAMMPS': ['1.1.1', '6.6.6'],
                         'VASP': ['1.1.1', '6.6.6']}

server_lst = ['Server1', 'Server2', 'Server3']

def xyz_to_smiles(fname: str) -> str: 
    """Конвертирует XYZ-файл в SMILES-код

    Args:
        fname (str): Путь к XYZ-файлу

    Returns:
        str: SMILES-код молекулы
    """
    try:
        mol = next(pybel.readfile("xyz", fname)) 
        
        smi = mol.write(format="can") 
        print(smi)
    
        return smi.split()[0].strip()
    except Exception as e:
        print(f"Ошибка при конвертации XYZ в SMILES: {e}")
        return ""


def test_insert() -> None:
    """
    Тестовая функция для добавления NoneType в таблицу расчетов
    """    
    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    insert_prompt = """
    INSERT INTO computations (calc_id, structure_id, total_energy, input_hash, program_id, date, server_name, scf_converged, input_file, log_file, atomic_file)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
     
    try: 
        cursor.execute(insert_prompt, (1, None, None, None, None, None, None))
    except:
        pass
    finally:
        connection.commit()
        connection.close()


def upload_calculation(calc_id: int, structure_id: int, total_energy: float, input_hash: str, program_id: str, 
                       date: str, server_name:str, scf_converged: bool, input_file_path: str, log_file_path: str, 
                       atomic_file_path: str) -> None:
    """Добавление данных в таблицу расчетов 

    Args:
        calc_id (int): id расчета
        structure_id (int): id структуры
        total_energy (float): общая энергия
        input_hash (str): хэш-строка инпут файла
        program_id (str): id программы для расчетов (выбирается вручную)
        date (str): дата (выбирается вручную)
        server_name (str): кластер, где производились расчеты (выбирается вручную, но возможно автоматизировать по логам через абсолютные пути к файлам)
        scf_converged (bool): свелся ли SCF
        input_file_path (str): путь к инпут файлу
        log_file_path (str): путь к лог файлу
        atomic_file_path (str): путь к atomic connectivity файлу
    """    

    insert_prompt = """
    INSERT INTO computations (calc_id, structure_id, total_energy, input_hash, program_id, date, server_name, scf_converged, input_file, log_file, atomic_file)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """

    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    cursor.execute(insert_prompt, (calc_id, structure_id, total_energy, input_hash, program_id,
                                   date, server_name, scf_converged, input_file_path, 
                                   log_file_path, atomic_file_path))
    
    connection.commit()
    connection.close()

def upload_atomic_properties(calc_id: int, atom_number: int, grade: float, f_x: float, f_y: float, f_z: float) -> None:
    """Добавление данных в таблицу атомных свойств

    Args:
        calc_id (int): id расчета
        atom_number (int): номер атома
        grade (float): грейд
        f_x (float): проекция силы на x
        f_y (float): проекция силы на y
        f_z (float): проекция силы на z
    """    

    atom_id = int(f'{calc_id}' + f'{atom_number}')

    insert_prompt = """
    INSERT INTO atomic_properties (calc_id, atom_id, atom_number, grade, f_x, f_y, f_z)
    VALUES (?, ?, ?, ?, ?, ?, ?)
    """

    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    cursor.execute(insert_prompt, (calc_id, atom_id, atom_number, grade, f_x, f_y, f_z))
    
    connection.commit()
    connection.close()


def upload_atoms(str_id: int, atom_number: int, element: str, x_coord: float, y_coord: float, z_coord: float) -> None:
    """Добавление данных в таблицу атомов

    Args:
        str_id (int): id структуры
        atom_number (int): номер атома
        element (str): элемент
        x_coord (float): координата x
        y_coord (float): координата y
        z_coord (float): координата z
    """    
    insert_prompt = """
    INSERT INTO atoms (str_id, atom_number, element, x_coord, y_coord, z_coord)
    VALUES (?, ?, ?, ?, ?, ?)
    """

    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    cursor.execute(insert_prompt, (str_id, atom_number, element, x_coord, y_coord, z_coord))

    connection.commit()
    connection.close()

def upload_structures(str_id: int, A_vector=None, B_vector=None, C_vector=None, alpha_angle=None, beta_angle=None, gamma_angle=None) -> None:
    """Добавление данных в таблицу структур

    Args:
        str_id (int): id структуры
        A_vector (_type_, optional): _description_. Defaults to None.
        B_vector (_type_, optional): _description_. Defaults to None.
        C_vector (_type_, optional): _description_. Defaults to None.
        alpha_angle (_type_, optional): _description_. Defaults to None.
        beta_angle (_type_, optional): _description_. Defaults to None.
        gamma_angle (_type_, optional): _description_. Defaults to None.
    """    

    insert_prompt = """
    INSERT INTO structures (str_id, A_vector, B_vector, C_vector, alpha_angle, beta_angle, gamma_angle)
    VALUES (?, ?, ?, ?, ?, ?, ?)
    """

    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    cursor.execute(insert_prompt, (str_id, A_vector, B_vector, C_vector, alpha_angle, beta_angle, gamma_angle))

    connection.commit()
    connection.close()


def upload_molecules(mol_id: int, coords: list):
    insert_prompt = """
    INSERT INTO molecules (mol_id, SMILES, elements, names)
    VALUES (?, ?, ?, ?)
    """

    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    elements_lst = []
    
    temp_xyz_path = f'data_uploaded/temp_structure_{mol_id}.xyz'
    
    try:
        with open(temp_xyz_path, 'w') as file:
            file.write(f'{len(coords)}\n')
            file.write('\n')
            for coord in coords:
                elements_lst.append(coord[0])
                if coords.index(coord) == len(coords) - 1:
                    file.write(coord[0] + f' {coord[1]} {coord[2]} {coord[3]}')
                else:
                    file.write(coord[0] + f' {coord[1]} {coord[2]} {coord[3]}\n')

        smi = xyz_to_smiles(temp_xyz_path)
        elements_unique = list(set(elements_lst))
        elements_str = ".".join(elements_unique)

        cursor.execute(insert_prompt, (mol_id, smi, elements_str, None))

        connection.commit()
        
    finally:
        # Удаляем временный файл
        # if os.path.exists(temp_xyz_path):
        #     os.remove(temp_xyz_path)
        connection.close()


def get_table(table_type: str) -> pd.DataFrame:
    """_summary_

    Args:
        table_type (str): _description_

    Returns:
        pd.DataFrame: _description_
    """
    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    cursor.execute(f"SELECT * FROM {table_type}")
    data = cursor.fetchall()

    column_names = [desc[0] for desc in cursor.description]
    df = pd.DataFrame(data, columns=column_names)

    connection.commit()
    connection.close()
     
    return df

def search_by_energy(min_energy: float = None, max_energy: float = None) -> pd.DataFrame:
    """Поиск расчетов по диапазону энергии

    Args:
        min_energy (float, optional): Минимальная энергия. Defaults to None.
        max_energy (float, optional): Максимальная энергия. Defaults to None.

    Returns:
        pd.DataFrame: Результаты поиска
    """
    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    query = "SELECT * FROM computations WHERE 1=1"
    params = []

    if min_energy is not None:
        query += " AND total_energy >= ?"
        params.append(min_energy)
    
    if max_energy is not None:
        query += " AND total_energy <= ?"
        params.append(max_energy)

    cursor.execute(query, params)
    data = cursor.fetchall()

    column_names = [desc[0] for desc in cursor.description]
    df = pd.DataFrame(data, columns=column_names)

    connection.close()
    return df

def search_by_elements(elements: list) -> pd.DataFrame:
    """Поиск молекул по химическим элементам

    Args:
        elements (list): Список элементов для поиска. 
            Находит все молекулы, содержащие перечисленные элементы, перечисленные в любом порядке.

    Returns:
        pd.DataFrame: Результаты поиска
    """
    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()
    
    search_elements = set(elements)
    
    cursor.execute("SELECT * FROM molecules")
    data = cursor.fetchall()
    
    filtered_data = []
    for row in data:
        # Преобразуем строку элементов из БД в множество
        mol_elements = set(row[2].split('.'))
        # Проверяем, содержит ли молекула все искомые элементы
        if search_elements.issubset(mol_elements):
            filtered_data.append(row)

    column_names = [desc[0] for desc in cursor.description]
    df = pd.DataFrame(filtered_data, columns=column_names)

    connection.close()
    return df

def search_by_substructure(substructure_smiles: str) -> pd.DataFrame:
    """Поиск молекул по подструктуре

    Args:
        substructure_smiles (str): SMILES-код подструктуры

    Returns:
        pd.DataFrame: Результаты поиска
    """
    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    cursor.execute("SELECT * FROM molecules")
    molecules = cursor.fetchall()

    results = []
    substructure_mol = Chem.MolFromSmiles(substructure_smiles, sanitize=False)

    for mol_id, smiles, elements, names in molecules:
        if '.' in smiles:
            mol_parts = smiles.split('.')
            for part in mol_parts:
                main_mol = Chem.MolFromSmiles(part, sanitize=False)
                if main_mol and main_mol.HasSubstructMatch(substructure_mol):
                    results.append((mol_id, smiles, elements, names))
                    break
        else:
            main_mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if main_mol and main_mol.HasSubstructMatch(substructure_mol):
                results.append((mol_id, smiles, elements, names))

    column_names = ['mol_id', 'SMILES', 'elements', 'names']
    df = pd.DataFrame(results, columns=column_names)

    connection.close()
    return df

def get_last_calcid() -> int:
    """Функция для определения последнего значения calc_id в таблице расчетов

    Returns:
        int: Последнее значение calc_id, которое было записано в таблицу
    """

    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    try:
        cursor.execute('SELECT calc_id FROM computations ORDER BY rowid DESC LIMIT 1')
        last_score = cursor.fetchone()[0]
    except:
        last_score = 0

    connection.commit()
    connection.close()

    return last_score

def get_last_strid() -> int:
    """Функция для определения последнего значения str_id в таблице расчетов

    Returns:
        int: Последнее значение str_id, которое было записано в таблицу
    """
    # TODO: добавить сравнение новой структуры с уже имеющимися в БД через RMSD

    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    try:
        cursor.execute('SELECT structure_id FROM computations ORDER BY rowid DESC LIMIT 1')
        last_score = cursor.fetchone()[0]
    except:
        last_score = 0

    connection.commit()
    connection.close()

    return last_score

def get_last_molid() -> int:
    """Функция для определения последнего значения mol_id в таблице молекул

    Returns:
        int: Последнее значение mol_id, которое было записано в таблицу
    """
    # TODO: добавить сравнение новой структуры с уже имеющимися в БД через RMSD

    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    try:
        cursor.execute('SELECT mol_id FROM molecules ORDER BY rowid DESC LIMIT 1')
        last_score = cursor.fetchone()[0]
    except:
        last_score = 0

    connection.commit()
    connection.close()

    return last_score

def get_main_structures(substructure: str) -> list:
    # TODO: переписать с зарядами
    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    cursor.execute('SELECT SMILES FROM molecules')
    rows = cursor.fetchall()
    column_values = [row[0] for row in rows]

    matches_lst = []

    for mol in column_values:
        substructure_mol = Chem.MolFromSmiles(substructure, sanitize=False)
        if '.' in mol:
            mol_lst_tmp = mol.split('.')
            
            matches_lst_tmp = []
            for mol_tmp in mol_lst_tmp:
                main_structure_mol = Chem.MolFromSmiles(mol_tmp, sanitize=False)

                matches = main_structure_mol.HasSubstructMatch(substructure_mol)

                matches_lst_tmp.append(matches)

            matches_lst.append(matches_lst_tmp)

        main_structure_mol = Chem.MolFromSmiles(mol, sanitize=False)
        substructure_mol = Chem.MolFromSmiles(substructure, sanitize=False)

        print(main_structure_mol)
        print(substructure_mol)

        matches = main_structure_mol.HasSubstructMatch(substructure_mol)

        matches_lst.append(matches)

    connection.close()

    return matches_lst

def search_by_calculation_params(program: str = None, include_forces: bool = False) -> pd.DataFrame:
    """
    Поиск расчётов по программе и наличию сил
    
    Args:
        program (str, optional): Название программы (Orca, LAMMPS, VASP)
        include_forces (bool, optional): Включать ли расчёты с силами
    
    Returns:
        pd.DataFrame: DataFrame с результатами поиска
    """
    connection = sqlite3.connect('main.db')
    
    program_id_map = {
        'Orca': 100,
        'LAMMPS': 200,
        'VASP': 300
    }
    
    query = """
    SELECT c.calc_id, c.program_id, c.date, c.server_name, c.scf_converged, c.total_energy,
           m.mol_id, m.SMILES as smiles
    FROM computations c
    LEFT JOIN molecules m ON c.structure_id = m.mol_id
    WHERE 1=1
    """
    
    params = []
    
    if program:
        query += " AND c.program_id = ?"
        params.append(program_id_map.get(program))
    
    if include_forces:
        query += " AND EXISTS (SELECT 1 FROM atomic_properties ap WHERE ap.calc_id = c.calc_id AND ap.f_x IS NOT NULL)"
    
    # отладочная 
    # st.write("Параметры поиска:")
    # st.write(f"- Программа: {program} (program_id: {program_id_map.get(program) if program else None})")
    # st.write(f"- Поиск с силами: {include_forces}")
    # st.write("SQL запрос:")
    # st.code(query)
    # st.write("Параметры запроса:", params)
    
    # проверим содержимое таблицы computations
    # all_computations = pd.read_sql_query("SELECT * FROM computations", connection)
    # st.write("Все записи в таблице computations:")
    # st.dataframe(all_computations)
    
    # if include_forces:
    #     # Проверим таблицу atomic_properties
    #     forces = pd.read_sql_query("SELECT * FROM atomic_properties", connection)
    #     st.write("Записи в таблице atomic_properties:")
    #     st.dataframe(forces)
    
    df = pd.read_sql_query(query, connection, params=params)
    connection.close()
    
    return df

def export_structures(calc_ids: list, format: str = 'xyz') -> str:
    """
    Экспортирует все выбранные структуры и их данные в выбранном формате
    
    Args:
        calc_ids (list): Список ID расчётов
        format (str): Формат файла (xyz, pdb, cfg)
    
    Returns:
        str: Содержимое файла в выбранном формате
    """
    connection = sqlite3.connect('main.db')
    
    # Получаем данные о расчётах
    calc_query = """
    SELECT c.*, m.SMILES
    FROM computations c
    LEFT JOIN molecules m ON c.structure_id = m.mol_id
    WHERE c.calc_id IN ({})
    """.format(','.join('?' * len(calc_ids)))
    
    calc_df = pd.read_sql_query(calc_query, connection, params=calc_ids)
    
    # Получаем координаты 
    atoms_query = """
    SELECT a.*, c.calc_id
    FROM atoms a
    JOIN computations c ON a.str_id = c.structure_id
    WHERE c.calc_id IN ({})
    ORDER BY c.calc_id, a.atom_number
    """.format(','.join('?' * len(calc_ids)))
    
    atoms_df = pd.read_sql_query(atoms_query, connection, params=calc_ids)
    
    # Получаем силы 
    forces_query = """
    SELECT ap.calc_id, ap.atom_id, ap.atom_number, ap.grade, ap.f_x, ap.f_y, ap.f_z
    FROM atomic_properties ap
    WHERE ap.calc_id IN ({})
    ORDER BY ap.calc_id, ap.atom_number
    """.format(','.join('?' * len(calc_ids)))
    
    forces_df = pd.read_sql_query(forces_query, connection, params=calc_ids)
    
    connection.close()
    
    if atoms_df.empty or calc_df.empty:
        return None
    
    # отладка 
    print("Columns in forces_df:", forces_df.columns.tolist())
    print("Sample of forces_df:")
    print(forces_df.head())
    
    content = ""
    
    if format == 'xyz':
        for calc_id in calc_ids:
            calc_data = calc_df[calc_df['calc_id'] == calc_id].iloc[0]
            atoms_data = atoms_df[atoms_df['calc_id'] == calc_id]
            forces_data = forces_df[forces_df['calc_id'] == calc_id] if not forces_df.empty else pd.DataFrame()
            
            content += f"{len(atoms_data)}\n"
            content += f"Calculation ID: {calc_id}\n"
            content += f"SMILES: {calc_data['SMILES']}\n"
            content += f"Program: {calc_data['program_id']}\n"
            content += f"Total Energy: {calc_data['total_energy']}\n"
            content += f"SCF Converged: {calc_data['scf_converged']}\n"
            content += f"Server: {calc_data['server_name']}\n"
            content += f"Date: {calc_data['date']}\n"
            content += "Coordinates and Forces:\n"
            
            for i, row in atoms_data.iterrows():
                content += f"{row['element']} {row['x_coord']:.6f} {row['y_coord']:.6f} {row['z_coord']:.6f}"
                if not forces_data.empty:
                    force = forces_data[forces_data['atom_number'] == row['atom_number']]
                    if not force.empty and not pd.isna(force['F_x'].iloc[0]):
                        content += f" {force['F_x'].iloc[0]:.6f} {force['F_y'].iloc[0]:.6f} {force['F_z'].iloc[0]:.6f}"
                content += "\n"
            content += "\n"
            
    elif format == 'pdb':
        for calc_id in calc_ids:
            calc_data = calc_df[calc_df['calc_id'] == calc_id].iloc[0]
            atoms_data = atoms_df[atoms_df['calc_id'] == calc_id]
            forces_data = forces_df[forces_df['calc_id'] == calc_id] if not forces_df.empty else pd.DataFrame()
            
            content += f"TITLE     Calculation ID: {calc_id}\n"
            content += f"REMARK    SMILES: {calc_data['SMILES']}\n"
            content += f"REMARK    Program: {calc_data['program_id']}\n"
            content += f"REMARK    Total Energy: {calc_data['total_energy']}\n"
            content += f"REMARK    SCF Converged: {calc_data['scf_converged']}\n"
            content += f"REMARK    Server: {calc_data['server_name']}\n"
            content += f"REMARK    Date: {calc_data['date']}\n"
            content += "MODEL     1\n"
            
            for i, row in atoms_data.iterrows():
                content += f"ATOM  {(i+1):5d}  {row['element']:<3s} MOL     1    {row['x_coord']:8.3f}{row['y_coord']:8.3f}{row['z_coord']:8.3f}  1.00  0.00          {row['element']:>2s}\n"
                if not forces_data.empty:
                    force = forces_data[forces_data['atom_number'] == row['atom_number']]
                    if not force.empty and not pd.isna(force['F_x'].iloc[0]):
                        content += f"REMARK    Forces for atom {i+1}: {force['F_x'].iloc[0]:.6f} {force['F_y'].iloc[0]:.6f} {force['F_z'].iloc[0]:.6f}\n"
            
            content += "ENDMDL\n"
            content += "END\n\n"
            
    elif format == 'cfg':
        for calc_id in calc_ids:
            calc_data = calc_df[calc_df['calc_id'] == calc_id].iloc[0]
            atoms_data = atoms_df[atoms_df['calc_id'] == calc_id]
            forces_data = forces_df[forces_df['calc_id'] == calc_id] if not forces_df.empty else pd.DataFrame()
            
            content += "BEGIN_CFG\n"
            content += f"Size\n{len(atoms_data)}\n"
            content += "Supercell\n"
            content += "0.0 0 0\n"
            content += "0 0.0 0\n"
            content += "0 0 0.0\n"
            
            content += f"# Calculation ID: {calc_id}\n"
            content += f"# SMILES: {calc_data['SMILES']}\n"
            content += f"# Program: {calc_data['program_id']}\n"
            content += f"# Total Energy: {calc_data['total_energy']}\n"
            content += f"# SCF Converged: {calc_data['scf_converged']}\n"
            content += f"# Server: {calc_data['server_name']}\n"
            content += f"# Date: {calc_data['date']}\n"
            
            has_forces = not forces_data.empty and not forces_data['F_x'].isna().all()
            content += "AtomData: id type cartes_x cartes_y cartes_z"
            if has_forces:
                content += " fx fy fz"
            content += "\n"
            
            for i, row in atoms_data.iterrows():
                content += f"{i+1} {row['element']} {row['x_coord']:.6f} {row['y_coord']:.6f} {row['z_coord']:.6f}"
                if has_forces:
                    force = forces_data[forces_data['atom_number'] == row['atom_number']]
                    if not force.empty and not pd.isna(force['F_x'].iloc[0]):
                        content += f" {force['F_x'].iloc[0]:.6f} {force['F_y'].iloc[0]:.6f} {force['F_z'].iloc[0]:.6f}"
                content += "\n"
            
            content += f"Energy\n{calc_data['total_energy']}\n"
            content += "END_CFG\n\n"
    
    return content

def main():
    try:
        os.mkdir('data_uploaded')
    except:
        pass
    
    tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9 = st.tabs([
        'Импорт', 'Таблица расчетов', 'Таблица программ', 'Таблица атомных свойств', 
        'Таблица атомов', 'Таблица молекул', 'Редактор', 'Поиск по молекулам', 'Поиск по параметрам'
    ])
    
    with tab1:
        st.title('Импорт log-файлов')

        col1, col2, col3 = st.columns(3)

        software_selected_option = col1.selectbox("Выберите программу", software_lst)
        version_selected_option = col2.selectbox(f"Выберите версию {software_selected_option}",
                                               software_version_dict[software_selected_option])
        
        server_selected_option = col3.selectbox("Выберите сервер", server_lst)
    
        log_file = st.file_uploader("Загрузить log-файл", type=['out', 'cfg', 'log', 'dump'])
        input_file = st.file_uploader("Загрузить input-файл", type=['inp', 'gjf'])
        atomic_connectivity_file = st.file_uploader("Загрузить atomic connectivity", type=['xyz', 'pdb', 'out', 'log'])

        if st.button('Загрузить'):
            try:
                new_index = get_last_calcid() + 1
                new_index_strid = get_last_strid() + 1
                new_index_molid = get_last_molid() + 1

                # Создаем пути к файлам
                log_file_path = f'data_uploaded/id_{new_index}.log'
                dump_file_path = f'data_uploaded/id_{new_index}.dump'
                input_file_path = f'data_uploaded/id_{new_index}.inp'
                xyz_file_path = f'data_uploaded/id_{new_index}.xyz'

                # Проверяем наличие необходимых файлов
                if software_selected_option == 'Orca':
                    if not input_file or not log_file:
                        st.error("Для Orca необходимы input и log файлы")
                        return
                elif software_selected_option == 'LAMMPS':
                    if not log_file or not atomic_connectivity_file:
                        st.error("Для LAMMPS необходимы log и dump файлы")
                        return

                # Сохраняем файлы в зависимости от выбранной программы
                if software_selected_option == 'Orca':
                    with open(input_file_path, 'wb') as f:
                        f.write(input_file.getvalue())
                    with open(log_file_path, 'wb') as f:
                        f.write(log_file.getvalue())
                    if atomic_connectivity_file:
                        with open(xyz_file_path, 'wb') as f:
                            f.write(atomic_connectivity_file.getvalue())
                elif software_selected_option == 'LAMMPS':
                    with open(log_file_path, 'wb') as f:
                        f.write(log_file.getvalue())
                    with open(dump_file_path, 'wb') as f:
                        f.write(atomic_connectivity_file.getvalue())

                # Создаем парсер в зависимости от программы
                parser = None
                try:
                    if software_selected_option == 'Orca':
                        parser = OrcaParser(input_file_path, log_file_path)
                    elif software_selected_option == 'LAMMPS':
                        parser = LammpsParser(log_file_path, dump_file_path)
                    elif software_selected_option == 'VASP':
                        st.error("Поддержка VASP пока не реализована")
                        return
                except Exception as e:
                    st.error(f"Ошибка при создании парсера: {str(e)}")
                    return

                if parser:
                    try:
                        data = parser.get_data()
                        
                        program_id_map = {
                            'Orca': 100,
                            'LAMMPS': 200,
                            'VASP': 300
                        }

                        # Определяем пути к файлам для загрузки в БД
                        input_path = input_file_path if software_selected_option == 'Orca' else None
                        atomic_path = xyz_file_path if software_selected_option == 'Orca' else dump_file_path

                        upload_calculation(
                            calc_id=new_index,
                            structure_id=new_index_strid,
                            total_energy=data.get('energy', 0.0),
                            input_hash='',  # TODO: добавить хеширование
                            program_id=program_id_map[software_selected_option],
                            date='19.03.2025',
                            server_name=server_selected_option,
                            scf_converged=True if software_selected_option != 'Orca' else parser.get_scf_convergence(),
                            input_file_path=input_path,
                            log_file_path=log_file_path,
                            atomic_file_path=atomic_path
                        )

                        coords = data.get('coordinates', [])
                        forces = data.get('forces')
                        elements = data.get('elements', [])

                        if not coords or not elements or len(coords) != len(elements):
                            st.error("Ошибка: некорректные данные о координатах или элементах")
                            return

                        for i, (element, coord) in enumerate(zip(elements, coords)):
                            upload_atoms(
                                str_id=new_index_strid,
                                atom_number=i+1,
                                element=element,
                                x_coord=coord[0],
                                y_coord=coord[1],
                                z_coord=coord[2]
                            )

                            if forces and i < len(forces):
                                force = forces[i]
                                upload_atomic_properties(
                                    calc_id=new_index,
                                    atom_number=i+1,
                                    grade=None,
                                    f_x=force[0],
                                    f_y=force[1],
                                    f_z=force[2]
                                )

                        st.success(f'Файл загружен и сохранен в БД, calc_id={new_index}')
                    except Exception as e:
                        st.error(f"Ошибка при обработке данных: {str(e)}")
            except Exception as e:
                st.error(f"Общая ошибка: {str(e)}")


    with tab2:
        st.dataframe(get_table('computations'), use_container_width=True)

    with tab3:
        st.dataframe(get_table('programs'), use_container_width=True)

    with tab4:
        st.dataframe(get_table('atomic_properties'), use_container_width=True)

    with tab5:
        st.dataframe(get_table('atoms'), use_container_width=True)

    with tab6:
        st.dataframe(get_table('molecules'), use_container_width=True)

    with tab7:
        molecule = st.text_input("Molecule", "CCO")
        smile_code = st_ketcher(molecule)
        st.markdown(f"Smile code: ``{smile_code}``")

        substructure_smiles = st.text_input('c')

        if st.button('Check'):
            st.text(get_main_structures(substructure=substructure_smiles))

    with tab8:
        st.title('Поиск по молекулам')
        
        search_type = st.selectbox(
            "Выберите тип поиска",
            ["Поиск по энергии", "Поиск по элементам", "Поиск по подструктуре"]
        )
        
        if search_type == "Поиск по энергии":
            col1, col2 = st.columns(2)
            min_energy = col1.number_input("Минимальная энергия", value=-400.0)
            max_energy = col2.number_input("Максимальная энергия", value=-200.0)
            
            if st.button("Найти"):
                results = search_by_energy(min_energy, max_energy)
                st.dataframe(results, use_container_width=True)
                
        elif search_type == "Поиск по элементам":
            elements = st.text_input("Введите элементы через запятую (например: C,H,O)")
            if st.button("Найти"):
                elements_list = [e.strip().upper() for e in elements.split(",")]
                elements_list = list(filter(None, elements_list))
                elements_list = list(dict.fromkeys(elements_list))
                
                st.write(f"Поиск молекул, содержащих только элементы: {', '.join(elements_list)}")
                
                results = search_by_elements(elements_list)
                
                if len(results) == 0:
                    st.warning("Молекулы с таким набором элементов не найдены")
                else:
                    st.success(f"Найдено молекул: {len(results)}")
                    st.dataframe(results, use_container_width=True)
                
        else:  # Поиск по подструктуре
            substructure = st.text_input("Введите SMILES-код подструктуры")
            if st.button("Найти"):
                results = search_by_substructure(substructure)
                st.dataframe(results, use_container_width=True)
    
    with tab9:
        render_search_tab()


if __name__ == '__main__':
    main()