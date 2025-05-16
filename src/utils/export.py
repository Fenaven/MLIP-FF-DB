import pandas as pd
from typing import List, Optional
from ..database.db_operations import get_connection
import io
import zipfile

def _generate_structure_content(calc_data, atoms_data, forces_data, format, calc_id):
    content = ""
    if format == 'xyz':
        content += f"{len(atoms_data)}\n"
        content += f"Calculation ID: {calc_id}\n"
        content += f"SMILES: {calc_data['SMILES']}\n"
        content += f"Program: {calc_data['program_id']}\n"
        content += f"Total Energy: {calc_data['total_energy']}\n"
        content += f"SCF Converged: {calc_data['scf_converged']}\n"
        content += f"Server: {calc_data['server_name']}\n"
        content += f"Date: {calc_data['date']}\n"
        content += "Coordinates and Forces:\n"
        for i, row in atoms_data.iterrows():
            content += f"{row['element']} {row['x_coord']:.6f} {row['y_coord']:.6f} {row['z_coord']:.6f}"
            if not forces_data.empty:
                force = forces_data[forces_data['atom_number'] == row['atom_number']]
                if not force.empty and not pd.isna(force['F_x'].iloc[0]):
                    content += f" {force['F_x'].iloc[0]:.6f} {force['F_y'].iloc[0]:.6f} {force['F_z'].iloc[0]:.6f}"
            content += "\n"
        content += "\n"
    elif format == 'pdb':
        content += f"TITLE     Calculation ID: {calc_id}\n"
        content += f"REMARK    SMILES: {calc_data['SMILES']}\n"
        content += f"REMARK    Program: {calc_data['program_id']}\n"
        content += f"REMARK    Total Energy: {calc_data['total_energy']}\n"
        content += f"REMARK    SCF Converged: {calc_data['scf_converged']}\n"
        content += f"REMARK    Server: {calc_data['server_name']}\n"
        content += f"REMARK    Date: {calc_data['date']}\n"
        content += "MODEL     1\n"
        for i, row in atoms_data.iterrows():
            content += f"ATOM  {(i+1):5d}  {row['element']:<3s} MOL     1    {row['x_coord']:8.3f}{row['y_coord']:8.3f}{row['z_coord']:8.3f}  1.00  0.00          {row['element']:>2s}\n"
            if not forces_data.empty:
                force = forces_data[forces_data['atom_number'] == row['atom_number']]
                if not force.empty and not pd.isna(force['F_x'].iloc[0]):
                    content += f"REMARK    Forces for atom {i+1}: {force['F_x'].iloc[0]:.6f} {force['F_y'].iloc[0]:.6f} {force['F_z'].iloc[0]:.6f}\n"
        content += "ENDMDL\n"
        content += "END\n\n"
    elif format == 'cfg':
        content += "BEGIN_CFG\n"
        content += f"Size\n{len(atoms_data)}\n"
        content += "Supercell\n"
        content += "0.0 0 0\n"
        content += "0 0.0 0\n"
        content += "0 0 0.0\n"
        content += f"# Calculation ID: {calc_id}\n"
        content += f"# SMILES: {calc_data['SMILES']}\n"
        content += f"# Program: {calc_data['program_id']}\n"
        content += f"# Total Energy: {calc_data['total_energy']}\n"
        content += f"# SCF Converged: {calc_data['scf_converged']}\n"
        content += f"# Server: {calc_data['server_name']}\n"
        content += f"# Date: {calc_data['date']}\n"
        has_forces = not forces_data.empty and not forces_data['F_x'].isna().all()
        content += "AtomData: id type cartes_x cartes_y cartes_z"
        if has_forces:
            content += " fx fy fz"
        content += "\n"
        for i, row in atoms_data.iterrows():
            content += f"{i+1} {row['element']} {row['x_coord']:.6f} {row['y_coord']:.6f} {row['z_coord']:.6f}"
            if has_forces:
                force = forces_data[forces_data['atom_number'] == row['atom_number']]
                if not force.empty and not pd.isna(force['F_x'].iloc[0]):
                    content += f" {force['F_x'].iloc[0]:.6f} {force['F_y'].iloc[0]:.6f} {force['F_z'].iloc[0]:.6f}"
            content += "\n"
        content += f"Energy\n{calc_data['total_energy']}\n"
        content += "END_CFG\n\n"
    return content

def export_structures(calc_ids: List[int], format: str = 'xyz') -> bytes:
    """
    Экспортирует все выбранные структуры в отдельные файлы и упаковывает их в zip-архив
    Args:
        calc_ids (list): Список ID расчётов
        format (str): Формат файла (xyz, pdb, cfg)
    Returns:
        bytes: zip-архив с файлами
    """
    connection = get_connection()
    calc_query = """
    SELECT c.*, m.SMILES
    FROM computations c
    LEFT JOIN molecules m ON c.structure_id = m.mol_id
    WHERE c.calc_id IN ({})
    """.format(','.join('?' * len(calc_ids)))
    calc_df = pd.read_sql_query(calc_query, connection, params=calc_ids)
    atoms_query = """
    SELECT a.*, c.calc_id
    FROM atoms a
    JOIN computations c ON a.str_id = c.structure_id
    WHERE c.calc_id IN ({})
    ORDER BY c.calc_id, a.atom_number
    """.format(','.join('?' * len(calc_ids)))
    atoms_df = pd.read_sql_query(atoms_query, connection, params=calc_ids)
    forces_query = """
    SELECT ap.calc_id, ap.atom_id, ap.atom_number, ap.grade, ap.F_x, ap.F_y, ap.F_z
    FROM atomic_properties ap
    WHERE ap.calc_id IN ({})
    ORDER BY ap.calc_id, ap.atom_number
    """.format(','.join('?' * len(calc_ids)))
    forces_df = pd.read_sql_query(forces_query, connection, params=calc_ids)
    connection.close()
    if atoms_df.empty or calc_df.empty:
        return None
    mem_zip = io.BytesIO()
    with zipfile.ZipFile(mem_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for calc_id in calc_ids:
            calc_data = calc_df[calc_df['calc_id'] == calc_id].iloc[0]
            atoms_data = atoms_df[atoms_df['calc_id'] == calc_id]
            forces_data = forces_df[forces_df['calc_id'] == calc_id] if not forces_df.empty else pd.DataFrame()
            file_content = _generate_structure_content(calc_data, atoms_data, forces_data, format, calc_id)
            ext = format
            zf.writestr(f"structure_{calc_id}.{ext}", file_content)
    mem_zip.seek(0)
    return mem_zip.read() 