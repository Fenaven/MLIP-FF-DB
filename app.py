import streamlit as st
import sqlite3
import pandas as pd
import utils


software_lst = ['Orca', 'LAMMPS', 'VASP']

software_version_dict = {'Orca': ['5.0.0', '6.0.0'],
                         'LAMMPS': ['1.1.1', '6.6.6'],
                         'VASP': ['1.1.1', '6.6.6']}

server_lst = ['Server1', 'Server2', 'Server3']

# --- prompts ---

insert_prompt = """
    INSERT INTO computations (calc_id, structure_id, total_energy, input_hash, server_name, scf_converged, log_file)
    VALUES (?, ?, ?, ?, ?, ?, ?)
    """


def test_insert() -> None:
    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()
     
    try: 
        cursor.execute(insert_prompt, (1, None, None, None, None, None, None))
    except:
        pass
    finally:
        connection.commit()
        connection.close()


def upload_calculation(calc_id: int, structure_id: int, input_hash: str, data: dict, server_name: str, log_file_path: str) -> None:
    """_summary_

    Args:
        calc_id (int): _description_
        structure_id (int): _description_
        data (dict): _description_
        server_name (str): _description_
        log_file (str): _description_
    """

    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    cursor.execute(insert_prompt, (calc_id, structure_id, data['total_energy'], 
                                   input_hash, server_name, data['scf'], log_file_path))
    
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

def get_last_calcid():

    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    cursor.execute('SELECT calc_id FROM computations ORDER BY rowid DESC LIMIT 1')
    last_score = cursor.fetchone()[0]

    connection.commit()
    connection.close()

    return last_score

def main():
    test_insert()

    # --- streamlit app ---
    tab1, tab2, tab3 = st.tabs(['Импорт', 'Таблица расчетов', 'Таблица программ'])
    
    with tab1:
        st.title('Импорт log-файлов')

        col1, col2, col3 = st.columns(3)

        software_selected_option = col1.selectbox("Выберите программу", software_lst)
        version_selected_option = col2.selectbox(f"Выберите версию {software_selected_option}",
                                               software_version_dict[software_selected_option])
        
        server_selected_option = col3.selectbox("Выберите программу", server_lst)
    
        log_file = st.file_uploader("Загрузить log-файл", type=['out', 'cfg', 'log'])
        input_file = st.file_uploader("Загрузить input-файл", type=['inp', 'gjf'])
        atomic_connectivity_file = st.file_uploader("Загрузить atomic connectivity", type=['xyz', 'pdb'])

        # --- загрузка лога в базу --- 
        if st.button('Загрузить'):
            try:
                new_index = get_last_calcid() + 1

                log_file_path = f'uploaded_files/id_{new_index}.log'
                input_file_path = f'uploaded_files/id_{new_index}.inp'

                with open(log_file_path, 'wb') as file:
                    file.write(log_file.getvalue())

                # здесь нужно добавить парсинг с файла
                data = {'total_energy': 100, 'scf': False}
                input_hash = 'test_hash'
                upload_calculation(new_index, 1, input_hash, data, server_selected_option, log_file_path)
                
                st.text(f'Файл загружен и сохранен в БД, calc_id={new_index}')

            except:
                st.text('Нет загруженного файла')

    with tab2:
        # здесь нужно добавить возможность обновления базы? 
        st.table(get_table('computations'))

    with tab3:
        st.table(get_table('programs'))



if __name__ == '__main__':
    main()