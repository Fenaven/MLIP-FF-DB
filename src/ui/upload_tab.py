import streamlit as st
import os
from typing import Dict, Any, Optional
from src.parsers import OrcaParser, LammpsParser, XYZParser, CFGParser
from src.database.db_operations import (
    upload_calculation, upload_atomic_properties,
    upload_atoms, get_last_calcid, get_last_strid, get_last_molid
)

def save_uploaded_file(uploaded_file, directory: str, filename: str) -> Optional[str]:
    """Сохраняет загруженный файл"""
    if uploaded_file is None:
        return None
        
    try:
        os.makedirs(directory, exist_ok=True)
        file_path = os.path.join(directory, filename)
        
        with open(file_path, 'wb') as f:
            f.write(uploaded_file.getvalue())
            
        return file_path
    except Exception as e:
        st.error(f"Ошибка при сохранении файла: {e}")
        return None

def render_upload_tab():
    """Отображает вкладку загрузки файлов"""
    st.title('Импорт файлов расчётов')
    
    # Выбор программы и версии
    col1, col2, col3 = st.columns(3)
    
    software_lst = ['Orca', 'LAMMPS']
    software_version_dict = {
        'Orca': ['5.0.0', '6.0.0'],
        'LAMMPS': ['1.1.1', '6.6.6']
    }
    server_lst = ['Server1', 'Server2', 'Server3']
    
    with col1:
        program = st.selectbox("Выберите программу", software_lst)
    
    with col2:
        version = st.selectbox(
            f"Выберите версию {program}",
            software_version_dict[program]
        )
    
    with col3:
        server = st.selectbox("Выберите сервер", server_lst)
    
    # Загрузка файлов в зависимости от программы
    if program == 'Orca':
        input_file = st.file_uploader("Загрузить input файл", type=['inp'])
        output_file = st.file_uploader("Загрузить output файл", type=['out', 'log'])
        xyz_file = st.file_uploader("Загрузить XYZ файл (опционально)", type=['xyz'])
        
        if st.button('Загрузить файлы'):
            if not (input_file and output_file):
                st.error("Необходимо загрузить input и output файлы")
                return
            
            # Сохраняем файлы
            calc_id = get_last_calcid() + 1
            str_id = get_last_strid() + 1
            mol_id = get_last_molid() + 1
            
            input_path = save_uploaded_file('data_uploaded', f'id_{calc_id}.inp', input_file)
            output_path = save_uploaded_file('data_uploaded', f'id_{calc_id}.out', output_file)
            xyz_path = None
            if xyz_file:
                xyz_path = save_uploaded_file('data_uploaded', f'id_{calc_id}.xyz', xyz_file)
            
            # Парсим файлы
            parser = OrcaParser(input_path, output_path)
            data = parser.get_data()
            
            # Загружаем данные в БД
            upload_calculation(
                calc_id=calc_id,
                structure_id=str_id,
                total_energy=data['energy'],
                input_hash='',  # TODO: добавить хеширование
                program_id=100,  # Orca
                date='19.03.2025',  # TODO: добавить реальную дату
                server_name=server,
                scf_converged=True,  # TODO: добавить проверку
                input_file_path=input_path,
                log_file_path=output_path,
                atomic_file_path=xyz_path or ''
            )
            
            # Загружаем координаты и силы
            coords = data['coordinates']
            forces = data['forces']
            elements = data['elements']
            
            for i, (element, coord) in enumerate(zip(elements, coords)):
                upload_atoms(
                    str_id=str_id,
                    atom_number=i+1,
                    element=element,
                    x_coord=coord[0],
                    y_coord=coord[1],
                    z_coord=coord[2]
                )
                
                if forces is not None:
                    force = forces[i]
                    upload_atomic_properties(
                        calc_id=calc_id,
                        atom_number=i+1,
                        grade=None,
                        f_x=force[0],
                        f_y=force[1],
                        f_z=force[2]
                    )
            
            st.success(f"Файлы успешно загружены. ID расчёта: {calc_id}")
            
    elif program == 'LAMMPS':
        log_file = st.file_uploader("Загрузить log файл", type=['log'])
        dump_file = st.file_uploader("Загрузить dump файл", type=['dump'])
        cfg_file = st.file_uploader("Загрузить cfg файл (опционально)", type=['cfg'])
        
        if st.button('Загрузить файлы'):
            if not (log_file and dump_file):
                st.error("Необходимо загрузить log и dump файлы")
                return
            
            # Сохраняем файлы
            calc_id = get_last_calcid() + 1
            str_id = get_last_strid() + 1
            mol_id = get_last_molid() + 1
            
            log_path = save_uploaded_file('data_uploaded', f'id_{calc_id}.log', log_file)
            dump_path = save_uploaded_file('data_uploaded', f'id_{calc_id}.dump', dump_file)
            cfg_path = None
            if cfg_file:
                cfg_path = save_uploaded_file('data_uploaded', f'id_{calc_id}.cfg', cfg_file)
            
            # Парсим файлы
            parser = LammpsParser(dump_path, log_path)
            data = parser.get_data()
            
            # Загружаем данные в БД
            upload_calculation(
                calc_id=calc_id,
                structure_id=str_id,
                total_energy=data['energy'],
                input_hash='',  # TODO: добавить хеширование
                program_id=200,  # LAMMPS
                date='19.03.2025',  # TODO: добавить реальную дату
                server_name=server,
                scf_converged=True,
                input_file_path=log_path,
                log_file_path=dump_path,
                atomic_file_path=cfg_path or ''
            )
            
            # Загружаем координаты и силы
            coords = data['coordinates']
            forces = data['forces']
            elements = data['elements']
            
            for i, (element, coord) in enumerate(zip(elements, coords)):
                upload_atoms(
                    str_id=str_id,
                    atom_number=i+1,
                    element=element,
                    x_coord=coord[0],
                    y_coord=coord[1],
                    z_coord=coord[2]
                )
                
                if forces is not None:
                    force = forces[i]
                    upload_atomic_properties(
                        calc_id=calc_id,
                        atom_number=i+1,
                        grade=None,
                        f_x=force[0],
                        f_y=force[1],
                        f_z=force[2]
                    )
            
            st.success(f"Файлы успешно загружены. ID расчёта: {calc_id}") 