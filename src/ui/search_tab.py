import streamlit as st
import pandas as pd
from typing import List, Optional

from src.database.db_operations import get_table
from src.utils.export import export_structures

def render_search_tab():
    """Отображает вкладку поиска по параметрам расчёта"""
    
    st.header("Поиск по параметрам расчёта")
    
    # Инициализируем session_state для хранения результатов поиска
    if 'search_results' not in st.session_state:
        st.session_state.search_results = None
    
    # Используем columns для разделения интерфейса на секции
    search_col, results_col = st.columns([1, 2])
    
    with search_col:
        program = st.selectbox(
            "Выберите программу",
            options=[""] + ['Orca', 'LAMMPS', 'VASP'],
            index=0
        )
        
        data_type = st.radio(
            "Тип данных",
            options=["Только энергии", "Энергии и силы"],
            index=0
        )
        
        include_forces = data_type == "Энергии и силы"
        
        if st.button("Найти расчёты", use_container_width=True):
            if not program:
                st.warning("Пожалуйста, выберите программу для поиска")
            else:
                st.session_state.search_results = search_by_calculation_params(
                    program=program if program else None,
                    include_forces=include_forces
                )
    
    with results_col:
        if st.session_state.search_results is not None:
            if len(st.session_state.search_results) > 0:
                st.dataframe(st.session_state.search_results, use_container_width=True)
                st.success(f"Найдено расчётов: {len(st.session_state.search_results)}")
                
                # Добавляем кнопки экспорта для всех форматов
                st.subheader("Экспорт всех найденных структур")
                
                # Используем columns для размещения кнопок в ряд
                col1, col2, col3 = st.columns(3)
                
                calc_ids = st.session_state.search_results['calc_id'].tolist()
                
                with col1:
                    if st.button("Экспорт в XYZ", use_container_width=True):
                        content = export_structures(calc_ids, 'xyz')
                        if content:
                            st.download_button(
                                label="Скачать .xyz",
                                data=content,
                                file_name=f"structures.xyz",
                                mime="text/plain",
                                use_container_width=True
                            )
                
                with col2:
                    if st.button("Экспорт в PDB", use_container_width=True):
                        content = export_structures(calc_ids, 'pdb')
                        if content:
                            st.download_button(
                                label="Скачать .pdb",
                                data=content,
                                file_name=f"structures.pdb",
                                mime="text/plain",
                                use_container_width=True
                            )
                
                with col3:
                    if st.button("Экспорт в CFG", use_container_width=True):
                        content = export_structures(calc_ids, 'cfg')
                        if content:
                            st.download_button(
                                label="Скачать .cfg",
                                data=content,
                                file_name=f"structures.cfg",
                                mime="text/plain",
                                use_container_width=True
                            )
            else:
                st.info("Расчёты с указанными параметрами не найдены")

def search_by_calculation_params(program: Optional[str] = None, include_forces: bool = False) -> pd.DataFrame:
    """
    Поиск расчётов по программе и наличию сил
    
    Args:
        program (str, optional): Название программы (Orca, LAMMPS, VASP)
        include_forces (bool, optional): Включать ли расчёты с силами
    
    Returns:
        pd.DataFrame: DataFrame с результатами поиска
    """
    program_id_map = {
        'Orca': 100,
        'LAMMPS': 200,
        'VASP': 300
    }
    
    # Получаем данные из базы
    computations = get_table('computations')
    molecules = get_table('molecules')
    
    # Фильтруем по программе
    if program:
        computations = computations[computations['program_id'] == program_id_map[program]]
    
    # Фильтруем по наличию сил
    if include_forces:
        atomic_properties = get_table('atomic_properties')
        calc_ids_with_forces = atomic_properties[atomic_properties['F_x'].notna()]['calc_id'].unique()
        computations = computations[computations['calc_id'].isin(calc_ids_with_forces)]
    
    # Объединяем с данными о молекулах
    result = pd.merge(
        computations,
        molecules[['mol_id', 'SMILES']],
        left_on='structure_id',
        right_on='mol_id',
        how='left'
    )
    
    return result 