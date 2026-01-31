#!/usr/bin/env python
"""
Convert ELF .meg mesh files to Netgen .vol format.
"""

import re
import numpy as np


def parse_meg_file(meg_path):
    """Parse ELF .meg mesh file."""
    nodes = {}
    elements = []
    scale = 1.0
    
    with open(meg_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('*'):
                continue
            
            parts = line.split()
            if not parts:
                continue
            
            keyword = parts[0]
            
            if keyword == 'MGSC':
                scale = float(parts[1])
            
            elif keyword == 'MGR1':
                # Node: MGR1 node_id flag x y z
                node_id = int(parts[1])
                x = float(parts[3])
                y = float(parts[4])
                z = float(parts[5])
                nodes[node_id] = [x * scale, y * scale, z * scale]
            
            elif keyword in ('MMB8T', 'MCL8T', 'MYO8T'):
                # Hex element: TYPE elem_id flag material_id n1 n2 n3 n4 n5 n6 n7 n8
                elem_id = int(parts[1])
                material_id = int(parts[3])
                node_ids = [int(parts[i]) for i in range(4, 12)]
                elements.append({
                    'id': elem_id,
                    'type': keyword,
                    'material': material_id,
                    'nodes': node_ids
                })
    
    return nodes, elements, scale


def write_vol_file(vol_path, nodes, elements):
    """Write Netgen .vol format mesh file."""
    
    # Convert nodes to list sorted by ID
    node_ids = sorted(nodes.keys())
    node_map = {nid: i+1 for i, nid in enumerate(node_ids)}  # 1-indexed
    
    with open(vol_path, 'w') as f:
        f.write('mesh3d\n')
        f.write('dimension\n3\n')
        f.write('geomtype\n0\n')
        
        # Points
        f.write(f'\n# Points\n')
        f.write(f'points\n{len(nodes)}\n')
        for nid in node_ids:
            x, y, z = nodes[nid]
            f.write(f'{x:.15e} {y:.15e} {z:.15e}\n')
        
        # Volume elements (hexahedra)
        hex_elements = [e for e in elements if len(e['nodes']) == 8]
        
        f.write(f'\n# Volume elements\n')
        f.write(f'volumeelements\n{len(hex_elements)}\n')
        for elem in hex_elements:
            mat = elem['material']
            # Convert node IDs and reorder for Netgen hex convention
            # ELF: 1,2,4,3,5,6,8,7 -> Netgen: 1,2,3,4,5,6,7,8
            elf_nodes = elem['nodes']
            # ELF ordering: n1,n2,n4,n3,n5,n6,n8,n7 (CHEXA-like)
            # Netgen ordering: n1,n2,n3,n4,n5,n6,n7,n8
            # Map: ELF[0,1,2,3,4,5,6,7] -> Netgen[0,1,3,2,4,5,7,6]
            netgen_order = [0, 1, 3, 2, 4, 5, 7, 6]
            reordered = [node_map[elf_nodes[i]] for i in netgen_order]
            
            f.write(f'{mat} {" ".join(map(str, reordered))}\n')
        
        # Surface elements (empty for now)
        f.write(f'\n# Surface elements\n')
        f.write(f'surfaceelements\n0\n')
        
        # Edge segments (empty)
        f.write(f'\n# Edge segments\n')
        f.write(f'edgesegmentsgi2\n0\n')
        
        # End markers
        f.write(f'\n# Materials\n')
        f.write(f'materials\n1\n')
        f.write(f'1 domain1\n')
        
        f.write(f'\nendmesh\n')


def convert_meg_to_vol(meg_path, vol_path=None):
    """Convert .meg file to .vol file."""
    if vol_path is None:
        vol_path = meg_path.replace('.meg', '.vol').replace('.MEG', '.vol')
    
    print(f'Reading: {meg_path}')
    nodes, elements, scale = parse_meg_file(meg_path)
    
    print(f'  Nodes: {len(nodes)}')
    print(f'  Elements: {len(elements)}')
    print(f'  Scale: {scale} (coordinates in meters)')
    
    print(f'Writing: {vol_path}')
    write_vol_file(vol_path, nodes, elements)
    
    print('Done.')
    return vol_path


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print('Usage: python meg_to_vol.py <input.meg> [output.vol]')
        sys.exit(1)
    
    meg_path = sys.argv[1]
    vol_path = sys.argv[2] if len(sys.argv) > 2 else None
    convert_meg_to_vol(meg_path, vol_path)
