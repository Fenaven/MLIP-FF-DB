import streamlit as st
import sqlite3
import pandas as pd
import utils
import numpy as np
import os
import parsers
from openbabel import pybel
from rdkit import Chem

from streamlit_ketcher import st_ketcher


software_lst = ['Orca', 'LAMMPS', 'VASP']

software_version_dict = {'Orca': ['5.0.0', '6.0.0'],
                         'LAMMPS': ['1.1.1', '6.6.6'],
                         'VASP': ['1.1.1', '6.6.6']}

server_lst = ['Server1', 'Server2', 'Server3']

def xyz_to_smiles(fname: str) -> str: 
     
        mol = next(pybel.readfile("xyz", fname)) 
 
        smi = mol.write(format="can") 
        print(smi)
    
        return smi.split()[0].strip()


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

# --- Далее все функции upload_ добавляют в таблицу данные

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

    with open('test_structure.xyz', 'w') as file:
        file.write(f'{len(coords)}\n')
        file.write('\n')
        for coord in coords:
            elements_lst.append(coord[0])
            if coords.index(coord) == len(coords) - 1:
                file.write(coord[0] + f' {coord[1]} {coord[2]} {coord[3]}')
            else:
                file.write(coord[0] + f' {coord[1]} {coord[2]} {coord[3]}\n')

    smi = xyz_to_smiles("test_structure.xyz") 
    elements_unique = list(set(elements_lst))
    elements_str = ".".join(elements_unique)

    cursor.execute(insert_prompt, (mol_id, smi, elements_str, None))

    connection.commit()
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


def main():
    # test_insert()

    try:
        os.mkdir('data_uploaded')
    except:
        pass

    # --- streamlit app ---
    tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs(['Импорт', 'Таблица расчетов', 'Таблица программ', 'Таблица атомных свойств', 
                                            'Таблица атомов', 'Таблица молекул', 'Редактор'])
    
    with tab1:
        st.title('Импорт log-файлов')

        col1, col2, col3 = st.columns(3)

        software_selected_option = col1.selectbox("Выберите программу", software_lst)
        version_selected_option = col2.selectbox(f"Выберите версию {software_selected_option}",
                                               software_version_dict[software_selected_option])
        
        server_selected_option = col3.selectbox("Выберите программу", server_lst)
    
        log_file = st.file_uploader("Загрузить log-файл", type=['out', 'cfg', 'log'])
        input_file = st.file_uploader("Загрузить input-файл", type=['inp', 'gjf'])
        atomic_connectivity_file = st.file_uploader("Загрузить atomic connectivity ", type=['xyz', 'pdb', 'out', 'log'])

        # --- загрузка лога в базу --- 
        if st.button('Загрузить'):
            new_index = get_last_calcid() + 1
            new_index_strid = get_last_strid() + 1
            new_index_molid = get_last_molid() + 1

            input_file_path = f'data_uploaded/id_{new_index}.inp'
            log_file_path = f'data_uploaded/id_{new_index}.log'
            atomic_connectivity_file_path = f'data_uploaded/id_{new_index}.txt'

            for file_path, file in zip([input_file_path, log_file_path, atomic_connectivity_file_path],
                                           [input_file, log_file, atomic_connectivity_file]):
                with open(file_path, 'wb') as f:
                    f.write(file.getvalue())

            if software_selected_option == 'Orca':
                parser = parsers.ParserORCA(input_file_path=input_file_path, log_file_path=log_file_path)
            else:
                pass

            parsed_data_dict = parser.get_summary()
            st.text(parsed_data_dict)
            print(new_index)

            # --- ОБНОВЛЕНИЕ ТАБЛИЦ ---

            upload_calculation(calc_id=new_index, structure_id=new_index_strid, total_energy=parsed_data_dict['total_energy'], input_hash=parsed_data_dict['input_hash'], 
                               program_id=100, date='19.03.2025', server_name=server_selected_option, scf_converged=parsed_data_dict['scf_convergence'],
                               input_file_path=input_file_path, log_file_path=log_file_path, atomic_file_path=atomic_connectivity_file_path)
            
            upload_molecules(new_index_molid, parsed_data_dict['coords'])

            for i in range(len(parsed_data_dict['forces'])):
                forces_per_atom_lst = parsed_data_dict['forces'][i]
                coords_per_atom_lst = parsed_data_dict['coords'][i]

                upload_atomic_properties(calc_id=new_index, atom_number=i+1, grade=None, f_x=forces_per_atom_lst[0],
                                         f_y=forces_per_atom_lst[1], f_z=forces_per_atom_lst[2])
                
                upload_atoms(str_id=new_index_strid, atom_number=i+1, element=coords_per_atom_lst[0], x_coord=coords_per_atom_lst[1],
                             y_coord=coords_per_atom_lst[2], z_coord=coords_per_atom_lst[3])
                
            st.text(f'Файл загружен и сохранен в БД, calc_id={new_index}')


    with tab2:
        st.dataframe(get_table('computations'))

    with tab3:
        st.dataframe(get_table('programs'))

    with tab4:
        st.dataframe(get_table('atomic_properties'))

    with tab5:
        st.dataframe(get_table('atoms'))

    with tab6:
        st.dataframe(get_table('molecules'))

    with tab7:
        molecule = st.text_input("Molecule", "CCO")
        smile_code = st_ketcher(molecule)
        st.markdown(f"Smile code: ``{smile_code}``")

        substructure_smiles = st.text_input('c')

        if st.button('Check'):
            st.text(get_main_structures(substructure=substructure_smiles))


if __name__ == '__main__':
    main()