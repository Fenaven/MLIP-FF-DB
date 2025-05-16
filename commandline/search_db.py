#!/usr/bin/env python3
import argparse
from app import search_by_energy, search_by_elements, search_by_substructure
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Поиск по базе данных')
    parser.add_argument('--type', choices=['energy', 'elements', 'substructure'], required=True,
                      help='Тип поиска: elements (по элементам), substructure (по подструктуре)')
    
    parser.add_argument('--elements',
                      help='Список элементов через запятую (например: C,H,O)')
    
    parser.add_argument('--smiles',
                      help='SMILES-код подструктуры')

    args = parser.parse_args()

    try:
        if args.type == 'elements':
            if not args.elements:
                print("Для поиска по элементам необходимо указать --elements")
                return
            elements_list = [e.strip().upper() for e in args.elements.split(",")]
            elements_list = list(filter(None, elements_list))
            elements_list = list(dict.fromkeys(elements_list))
            results = search_by_elements(elements_list)
            
        else:
            if not args.smiles:
                print("Для поиска по подструктуре необходимо указать --smiles")
                return
            results = search_by_substructure(args.smiles)

        if len(results) == 0:
            print("Результаты не найдены")
        else:
            print(f"Найдено результатов: {len(results)}")
            print("\nРезультаты:")
            print(results.to_string())

    except Exception as e:
        print(f"Ошибка при поиске: {e}")

if __name__ == '__main__':
    main() 