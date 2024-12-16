import sqlite3
from app import software_lst, software_version_dict
import utils


def main():
    connection = sqlite3.connect("main.db")
    cursor = connection.cursor()

    # sql-запросы для создания таблиц
    computation_table_prompt = """
    CREATE TABLE IF NOT EXISTS computations (
    calc_id INTEGER PRIMARY KEY,
    structure_id INTEGER,
    total_energy DOUBLE,
    input_hash VARCHAR,
    program_id INTEGER,
    date DATETIME,
    server_name VARCHAR,
    scf_converged BOOLEAN,
    input_file VARCHAR,
    log_file VARCHAR,
    atomic_file VARCHAR
    )
    """

    programs_table_prompt = """
    CREATE TABLE IF NOT EXISTS programs (
    program_id INTEGER PRIMARY KEY,
    program_name TEXT,
    program_version TEXT
    )
    """

    atomic_properties_table_prompt = """
    CREATE TABLE IF NOT EXISTS atomic_properties (
    calc_id INTEGER,
    atom_id INTEGER,
    atom_number INTEGER,
    grade FLOAT,
    F_x FLOAT,
    F_y FLOAT,
    F_z FLOAT
    )
    """

    structures_table_prompt = """
    CREATE TABLE IF NOT EXISTS structures (
    str_id INTEGER PRIMARY KEY,
    A_vector FLOAT,
    B_vector FLOAT,
    C_vector FLOAT,
    alpha_angle FLOAT,
    beta_angle FLOAT,
    gamma_angle FLOAT
    )
    """

    atoms_table_prompt = """
    CREATE TABLE IF NOT EXISTS atoms (
    str_id INTEGER,
    atom_number INTEGER,
    element TEXT,
    x_coord FLOAT,
    y_coord FLOAT,
    z_coord FLOAT
    )
    """

    molecules_table_prompt = """
    CREATE TABLE IF NOT EXISTS molecules (
    mol_id INTEGER,
    SMILES TEXT,
    elements TEXT,
    names TEXT
    )
    """


    # --- создание таблиц ---
    cursor.execute(computation_table_prompt)
    cursor.execute(programs_table_prompt)
    cursor.execute(atomic_properties_table_prompt)
    cursor.execute(structures_table_prompt)
    cursor.execute(atoms_table_prompt)
    cursor.execute(molecules_table_prompt)

    # --- заполнение таблицы программ ---
    # TODO: уточнить какие программы использовались в расчетах, получить от исполнителей полный список

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