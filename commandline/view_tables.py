#!/usr/bin/env python3
import argparse
import sqlite3
import pandas as pd

def get_table(table_name: str) -> pd.DataFrame:
    """Получение содержимого таблицы из базы данных

    Args:
        table_name (str): Имя таблицы

    Returns:
        pd.DataFrame: Содержимое таблицы
    """
    connection = sqlite3.connect('main.db')
    cursor = connection.cursor()

    cursor.execute(f'SELECT * FROM {table_name}')
    data = cursor.fetchall()

    column_names = [desc[0] for desc in cursor.description]
    df = pd.DataFrame(data, columns=column_names)

    connection.close()
    return df

def main():
    parser = argparse.ArgumentParser(description='Просмотр таблиц базы данных')
    parser.add_argument('--table', choices=['computations', 'programs', 'atomic_properties', 
                                          'atoms', 'molecules'], required=True,
                      help='Имя таблицы для просмотра')
    parser.add_argument('--limit', type=int,
                      help='Ограничить количество выводимых строк')
    parser.add_argument('--output', choices=['console', 'csv'],
                      help='Формат вывода (по умолчанию: console)')

    args = parser.parse_args()

    try:
        df = get_table(args.table)
        
        if args.limit:
            df = df.head(args.limit)
        
        if args.output == 'csv':
            output_file = f'{args.table}.csv'
            df.to_csv(output_file, index=False)
            print(f"Данные сохранены в файл {output_file}")
        else:
            print(f"\nТаблица {args.table}:")
            print(df.to_string())

    except Exception as e:
        print(f"Ошибка при получении данных: {e}")

if __name__ == '__main__':
    main() 