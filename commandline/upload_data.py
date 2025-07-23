#!/usr/bin/env python3
import argparse
import os
from app import upload_calculation, upload_atomic_properties, upload_atoms, upload_molecules, get_last_calcid, get_last_strid, get_last_molid
from utils.parsers import ParserORCA
import streamlit as st

def main():
    parser = argparse.ArgumentParser(description='Загрузка данных расчетов в базу данных')
    parser.add_argument('--program', choices=['Orca', 'VASP', 'LAMMPS'], required=True,
                      help='Программа, использованная для расчетов')
    parser.add_argument('--version', required=True,
                      help='Версия программы')
    parser.add_argument('--server', required=True,
                      help='Имя сервера, где производились расчеты')
    parser.add_argument('--input', required=True,
                      help='Путь к input-файлу')
    parser.add_argument('--log', required=True,
                      help='Путь к log-файлу')
    parser.add_argument('--atomic', required=True,
                      help='Путь к файлу с атомными свойствами')
    parser.add_argument('--date', required=True,
                      help='Дата проведения расчетов (формат: DD.MM.YYYY)')

    args = parser.parse_args()

    # Создаем директорию для загруженных файлов, если её нет
    try:
        os.makedirs('data_uploaded', exist_ok=True)
    except Exception as e:
        print(f"Ошибка при создании директории: {e}")
        return

    new_index = get_last_calcid() + 1
    new_index_strid = get_last_strid() + 1
    new_index_molid = get_last_molid() + 1

    input_file_path = f'data_uploaded/id_{new_index}.inp'
    log_file_path = f'data_uploaded/id_{new_index}.log'
    atomic_connectivity_file_path = f'data_uploaded/id_{new_index}.txt'

    try:
        for src, dst in [(args.input, input_file_path),
                         (args.log, log_file_path),
                         (args.atomic, atomic_connectivity_file_path)]:
            with open(src, 'rb') as fsrc, open(dst, 'wb') as fdst:
                fdst.write(fsrc.read())
    except Exception as e:
        print(f"Ошибка при копировании файлов: {e}")
        return

    if args.program == 'Orca':
        parser = ParserORCA(input_file_path=input_file_path, log_file_path=log_file_path)
    else:
        print(f"Парсер для программы {args.program} пока не реализован")
        return

    parsed_data_dict = parser.get_summary()
    print(f"Данные успешно распарсены: {parsed_data_dict}")

    try:
        upload_calculation(
            calc_id=new_index,
            structure_id=new_index_strid,
            total_energy=parsed_data_dict['total_energy'],
            input_hash=parsed_data_dict['input_hash'],
            program_id=100,  # TODO: реализовать получение program_id
            date=args.date,
            server_name=args.server,
            scf_converged=parsed_data_dict['scf_convergence'],
            input_file_path=input_file_path,
            log_file_path=log_file_path,
            atomic_file_path=atomic_connectivity_file_path
        )

        upload_molecules(new_index_molid, parsed_data_dict['coords'])

        for i in range(len(parsed_data_dict['forces'])):
            forces_per_atom_lst = parsed_data_dict['forces'][i]
            coords_per_atom_lst = parsed_data_dict['coords'][i]

            upload_atomic_properties(
                calc_id=new_index,
                atom_number=i+1,
                grade=None,
                f_x=forces_per_atom_lst[0],
                f_y=forces_per_atom_lst[1],
                f_z=forces_per_atom_lst[2]
            )

            upload_atoms(
                str_id=new_index_strid,
                atom_number=i+1,
                element=coords_per_atom_lst[0],
                x_coord=coords_per_atom_lst[1],
                y_coord=coords_per_atom_lst[2],
                z_coord=coords_per_atom_lst[3]
            )

        print(f'Данные успешно загружены в БД, calc_id={new_index}')

    except Exception as e:
        print(f"Ошибка при загрузке данных в БД: {e}")

if __name__ == '__main__':
    main() 