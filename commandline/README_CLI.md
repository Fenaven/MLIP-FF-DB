# Командная строка

В этом документе описаны скрипты для работы с базой данных через командную строку.

## Загрузка данных

Скрипт `upload_data.py` позволяет загружать данные расчетов в базу данных.

```bash
./upload_data.py --help
```

### Примеры использования:

```bash
# Загрузка данных расчета Orca
./upload_data.py --program Orca \
                 --version 5.0.0 \
                 --server Server1 \
                 --input data_uploaded/test.inp \
                 --log data_uploaded/test.log \
                 --atomic data_uploaded/test.log \
                 --date 19.03.2024
```

### Аргументы:
- `--program`: Программа для расчетов (Orca, VASP, LAMMPS)
- `--version`: Версия программы
- `--server`: Имя сервера, где производились расчеты
- `--input`: Путь к input-файлу
- `--log`: Путь к log-файлу
- `--atomic`: Путь к файлу atomic connectivity
- `--date`: Дата проведения расчетов (формат: DD.MM.YYYY)

## Поиск по базе данных

Скрипт `search_db.py` позволяет искать данные в базе по различным критериям.

```bash
./search_db.py --help
```

### Примеры использования:

```bash
# Поиск по элементам
./search_db.py --type elements \
               --elements C,H,O

# TODO
# Поиск по подструктуре
./search_db.py --type substructure \
               --smiles "[O]"
```

### Аргументы:
- `--type`: Тип поиска (energy, elements, substructure)
- `--min-energy`: Минимальная энергия (для поиска по энергии)
- `--max-energy`: Максимальная энергия (для поиска по энергии)
- `--elements`: Список элементов через запятую (для поиска по элементам)
- `--smiles`: SMILES-код подструктуры (для поиска по подструктуре)

## Просмотр таблиц

Скрипт `view_tables.py` позволяет просматривать содержимое таблиц базы данных.

```bash
./view_tables.py --help
```

### Примеры использования:

```bash
# Просмотр таблицы в консоли
./view_tables.py --table computations

# Просмотр с ограничением количества строк
./view_tables.py --table molecules --limit 10

# Сохранение таблицы в CSV
./view_tables.py --table atomic_properties --output csv
```

### Аргументы:
- `--table`: Имя таблицы (computations, programs, atomic_properties, atoms, molecules)
- `--limit`: Ограничение количества выводимых строк
- `--output`: Формат вывода (console, csv)

## Тестовые данные

В директории `data` находятся тестовые файлы:
- `input.inp`: Тестовый input-файл для Orca
- `butane_cut.log`: Тестовый log-файл для Orca