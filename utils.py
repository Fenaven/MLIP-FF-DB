import hashlib
import uuid

def get_input_hash(input_string: str):
    """Функция для хэширования инпут файла без геометрии

    Args:
        input_string (str): строка для хэширования, содержащая инпут файл

    Returns:
        _type_: хэш инпута
    """
    hash_object = hashlib.sha256()
    
    hash_object.update(input_string.encode('utf-8'))

    return hash_object.hexdigest()


def get_program_id(program_name: str, program_version: str) -> int:
    """_summary_

    Args:
        program_name (str): _description_
        program_version (str): _description_

    Returns:
        int: _description_
    """
    name_to_id_dict = {'Orca': '1', 'VASP': '2', 'LAMMPS': '3'}
    version_to_id_dict = {'5.0.0': '500', '6.0.0': '600', 
                          '1.1.1': '100', '6.6.6': '666'}
    
    program_id = name_to_id_dict[program_name] + version_to_id_dict[program_version]

    return program_id
