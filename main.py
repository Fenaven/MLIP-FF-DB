import sqlite3
from app import software_lst, software_version_dict
import utils


def main():
    connection = sqlite3.connect("main.db")
    cursor = connection.cursor()

    # sql-запрос для создания таблицы расчетов
    computation_table_prompt = """
    CREATE TABLE IF NOT EXISTS computations (
    calc_id INTEGER PRIMARY KEY,
    structure_id INTEGER,
    total_energy DOUBLE,
    input_hash VARCHAR,
    server_name VARCHAR,
    scf_converged BOOLEAN,
    log_file VARCHAR,
    program_id INTEGER
    )
    """

    soft_table_prompt = """
    CREATE TABLE IF NOT EXISTS programs (
    program_id INTEGER PRIMARY KEY,
    program_name TEXT,
    program_version TEXT
    )
    """

    # --- создание таблиц ---
    cursor.execute(computation_table_prompt)
    cursor.execute(soft_table_prompt)


    # --- заполнение таблицы программ ---
    insert_program_prompt = """
    INSERT INTO programs (program_id, program_name, program_version)
    VALUES (?, ?, ?)
    """

    
    for prog in software_lst:
        version_lst = software_version_dict[prog] 

        for version in version_lst:
            prog_id = utils.get_program_id(prog, version)

            cursor.execute(insert_program_prompt, (prog_id, prog, version))


    connection.commit()
    connection.close()


if __name__ == '__main__':
    main()