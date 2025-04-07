#!/usr/bin/env python3
import argparse
from app import search_by_energy, search_by_elements, search_by_substructure
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Поиск по базе данных')
    parser.add_argument('--type', choices=['energy', 'elements', 'substructure'], required=True,
                      help='Тип поиска: energy (по энергии), elements (по элементам), substructure (по подструктуре)')
    
    # Аргументы для поиска по энергии
    parser.add_argument('--min-energy', type=float,
                      help='Минимальная энергия для поиска')
    parser.add_argument('--max-energy', type=float,
                      help='Максимальная энергия для поиска')
    
    # Аргументы для поиска по элементам
    parser.add_argument('--elements',
                      help='Список элементов через запятую (например: C,H,O)')
    
    # Аргументы для поиска по подструктуре
    parser.add_argument('--smiles',
                      help='SMILES-код подструктуры')

    args = parser.parse_args()

    try:
        if args.type == 'energy':
            if args.min_energy is None or args.max_energy is None:
                print("Для поиска по энергии необходимо указать --min-energy и --max-energy")
                return
            results = search_by_energy(args.min_energy, args.max_energy)
            
        elif args.type == 'elements':
            if not args.elements:
                print("Для поиска по элементам необходимо указать --elements")
                return
            elements_list = [e.strip().upper() for e in args.elements.split(",")]
            elements_list = list(filter(None, elements_list))
            elements_list = list(dict.fromkeys(elements_list))
            results = search_by_elements(elements_list)
            
        else:  # substructure
            if not args.smiles:
                print("Для поиска по подструктуре необходимо указать --smiles")
                return
            results = search_by_substructure(args.smiles)

        # Выводим результаты
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