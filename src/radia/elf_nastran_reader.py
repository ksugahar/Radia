#!/usr/bin/env python
"""
ELF Nastran Reader - Read ELF_magic.nas files for verification

This module reads Nastran bulk data files from ELF_MAGIC and converts
them to Radia hexahedral elements.

Supported formats:
- GRID* (large field format, continuation with *)
- CHEXA (8-node hexahedron elements)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional


def parse_elf_nastran(filepath: str) -> Tuple[Dict[int, np.ndarray], List[dict]]:
    """
    Parse ELF_magic Nastran file.

    Args:
        filepath: Path to .nas file

    Returns:
        nodes: Dict mapping node ID to [x, y, z] coordinates
        elements: List of dicts with 'id', 'property_id', 'nodes' keys
    """
    nodes = {}
    elements = []

    with open(filepath, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].rstrip()

        # GRID* - Large field format (16-character fields)
        # GRID* 1234567890123456               X1234567890123456Y1234567890123456*
        # *     Z1234567890123456
        if line.startswith('GRID*'):
            # Parse node ID and X, Y from first line
            # Fields: GRID*, ID, CP, X, Y, *, (continuation)
            # Positions: 0-8, 8-24, 24-40, 40-56, 56-72, 72+
            node_id = int(line[8:24].strip())
            # CP field (24-40) is blank or coordinate system ID
            x = float(line[40:56].strip())
            y = float(line[56:72].strip())

            # Get continuation line for Z
            i += 1
            cont_line = lines[i].rstrip()
            # Continuation line starts with *, Z is in first field after *
            z = float(cont_line[8:24].strip())

            nodes[node_id] = np.array([x, y, z])

        # CHEXA - Hexahedron element
        # CHEXA  EID     PID     G1      G2      G3      G4      G5      G6      +
        # +      G7      G8
        elif line.startswith('CHEXA'):
            # Small field format (8-character fields)
            # Fields: CHEXA, EID, PID, G1, G2, G3, G4, G5, G6, +
            parts = line.split()
            elem_id = int(parts[1])
            prop_id = int(parts[2])
            # Nodes G1-G6 (first 6 of 8) - remove trailing + if present
            node_ids = []
            for j in range(3, 9):
                val = parts[j].rstrip('+')
                node_ids.append(int(val))

            # Get continuation line for G7, G8
            i += 1
            cont_line = lines[i].rstrip()
            cont_parts = cont_line.split()
            # Skip continuation marker '+', then get G7, G8
            if cont_parts[0].startswith('+'):
                node_ids.extend([int(cont_parts[1]), int(cont_parts[2])])
            else:
                node_ids.extend([int(cont_parts[0]), int(cont_parts[1])])

            elements.append({
                'id': elem_id,
                'property_id': prop_id,
                'nodes': node_ids
            })

        i += 1

    return nodes, elements


def get_hex_vertices(nodes: Dict[int, np.ndarray], element: dict) -> List[List[float]]:
    """
    Get 8 vertex coordinates for a hexahedron element.

    Args:
        nodes: Dict mapping node ID to coordinates
        element: Element dict with 'nodes' key

    Returns:
        List of 8 [x, y, z] coordinates
    """
    return [nodes[nid].tolist() for nid in element['nodes']]


def convert_to_radia(nodes: Dict[int, np.ndarray], elements: List[dict],
                     property_filter: Optional[int] = None) -> List[int]:
    """
    Convert Nastran hex elements to Radia ObjHexahedron.

    Args:
        nodes: Dict mapping node ID to coordinates
        elements: List of element dicts
        property_filter: If specified, only convert elements with this property ID

    Returns:
        List of Radia object IDs
    """
    import radia as rad

    radia_objs = []

    for elem in elements:
        if property_filter is not None and elem['property_id'] != property_filter:
            continue

        vertices = get_hex_vertices(nodes, elem)
        obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        radia_objs.append(obj)

    return radia_objs


def print_mesh_info(nodes: Dict[int, np.ndarray], elements: List[dict]):
    """Print summary of mesh contents."""
    print(f"Nodes: {len(nodes)}")
    print(f"Elements: {len(elements)}")

    # Count by property ID
    prop_counts = {}
    for elem in elements:
        pid = elem['property_id']
        prop_counts[pid] = prop_counts.get(pid, 0) + 1

    print("Elements by property ID:")
    for pid, count in sorted(prop_counts.items()):
        print(f"  Property {pid}: {count} elements")

    # Bounding box
    coords = np.array(list(nodes.values()))
    print(f"Bounding box:")
    print(f"  X: [{coords[:, 0].min():.6f}, {coords[:, 0].max():.6f}]")
    print(f"  Y: [{coords[:, 1].min():.6f}, {coords[:, 1].max():.6f}]")
    print(f"  Z: [{coords[:, 2].min():.6f}, {coords[:, 2].max():.6f}]")


def read_elf_matrix(filepath: str) -> np.ndarray:
    """
    Read ELF_magic.mat binary file.

    Args:
        filepath: Path to .mat file

    Returns:
        Interaction matrix as numpy array

    Format:
        Fortran unformatted file with N records (one per row):
        - Each record: [4-byte len][N*8 bytes data][4-byte len]
        - First record = first row of matrix
        - N records total for N x N matrix
    """
    with open(filepath, 'rb') as f:
        raw = f.read()

    # Parse Fortran unformatted records
    pos = 0
    rows = []
    while pos < len(raw) - 4:
        rec_len = np.frombuffer(raw[pos:pos+4], dtype='<i4')[0]
        if rec_len <= 0 or pos + 8 + rec_len > len(raw):
            break
        data = raw[pos+4:pos+4+rec_len]
        end_marker = np.frombuffer(raw[pos+4+rec_len:pos+8+rec_len], dtype='<i4')[0]
        if rec_len != end_marker:
            raise ValueError(f"Fortran record marker mismatch at pos {pos}")
        row = np.frombuffer(data, dtype='<f8')
        rows.append(row)
        pos += 8 + rec_len

    if len(rows) == 0:
        raise ValueError("No Fortran records found in file")

    # Verify all rows have same length
    n_dof = len(rows[0])
    if len(rows) != n_dof:
        raise ValueError(f"Expected {n_dof} rows but found {len(rows)}")

    for i, row in enumerate(rows):
        if len(row) != n_dof:
            raise ValueError(f"Row {i} has {len(row)} elements, expected {n_dof}")

    matrix = np.vstack(rows)
    return matrix


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print("Usage: python elf_nastran_reader.py <file.nas>")
        print("       python elf_nastran_reader.py <file.mat>  (read matrix)")
        sys.exit(1)

    filepath = sys.argv[1]

    if filepath.endswith('.mat'):
        # Read and display matrix info
        matrix = read_elf_matrix(filepath)
        print(f"Matrix shape: {matrix.shape}")
        print(f"Matrix range: [{matrix.min():.6f}, {matrix.max():.6f}]")
        print(f"Diagonal (first 6): {matrix.diagonal()[:6]}")
    else:
        # Read mesh
        nodes, elements = parse_elf_nastran(filepath)
        print_mesh_info(nodes, elements)
