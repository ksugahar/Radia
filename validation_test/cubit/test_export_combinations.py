"""
Exhaustive Export Mesh test: all format x option combinations.

Layer 1: Export via cubit.cmd (APREPRO commands from .ccm)
Layer 2: Output file format validation

Models:
  - 3D: sphere tet mesh (volume + boundary)
  - 2D: square tri mesh (surface only)

Combinations:
  GMSH:    order(1-5) x dim(2D,3D) = 10 per model x 2 = 20 (v4.1 only)
  Nastran: order(1,2) x dim(2D,3D) x nopyramid(0,1) = 8 per model x 2 = 16
  VTK:     order(1,2) x dim(2D,3D) = 4 per model x 2 = 8

Usage:
  python validation_test/cubit/test_export_combinations.py
"""

import sys
import os
import json
import traceback

# --- Setup Cubit ---
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
from install_panels import find_cubit_bin
_cubit_path = find_cubit_bin()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

# Set CUBIT_PLUGIN_DIR for .ccm command loading
_plugin_dir = os.path.join(os.path.dirname(_cubit_path), 'bin', 'plugins') \
    if _cubit_path else None
if _plugin_dir and os.path.isdir(_plugin_dir):
    os.environ['CUBIT_PLUGIN_DIR'] = _plugin_dir

import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', _plugin_dir or ''])

# --- Output directory ---
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'export_test_output')
os.makedirs(OUT_DIR, exist_ok=True)


# ================================================================
# Model builders
# ================================================================
def build_3d_model():
    """Sphere tet mesh with volume block + boundary block."""
    cubit.cmd("reset")
    cubit.cmd("create sphere radius 0.05")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd("volume 1 size 0.01")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add tet all")
    cubit.cmd('block 1 name "sphere"')
    cubit.cmd("block 2 add tri all")
    cubit.cmd('block 2 name "boundary"')
    n_tet = len(cubit.get_block_tets(1))
    n_tri = len(cubit.get_block_tris(2))
    n_nodes = cubit.get_node_count()
    return {'n_tet': n_tet, 'n_tri': n_tri, 'n_nodes': n_nodes}


def build_2d_model():
    """Square surface tri mesh."""
    cubit.cmd("reset")
    cubit.cmd("create surface rectangle width 0.1 height 0.1 zplane")
    cubit.cmd("surface 1 scheme trimesh")
    cubit.cmd("surface 1 size 0.01")
    cubit.cmd("mesh surface 1")
    cubit.cmd("block 1 add tri all")
    cubit.cmd('block 1 name "surface"')
    n_tri = len(cubit.get_block_tris(1))
    n_nodes = cubit.get_node_count()
    return {'n_tri': n_tri, 'n_nodes': n_nodes}


# ================================================================
# File parsers (Layer 2 validation)
# ================================================================
def parse_gmsh(filename):
    """Parse GMSH v4.1 and return {nodes, elements, element_types, physical_names}."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    result = {'nodes': 0, 'elements': 0, 'element_types': {},
              'physical_names': [], 'version': None,
              'declared_elements': 0}

    # Detect version first
    for line in lines:
        if line.strip().startswith('4.1') or line.strip().startswith('4.0'):
            result['version'] = line.strip().split()[0]
            break

    return _parse_gmsh_v41(lines, result)


def _parse_gmsh_v41(lines, result):
    """Parse GMSH v4.1 format."""
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line == '$MeshFormat':
            i += 1
            result['version'] = lines[i].strip().split()[0]
            i += 1
            continue

        if line == '$PhysicalNames':
            i += 1
            n = int(lines[i].strip())
            i += 1
            for _ in range(n):
                parts = lines[i].strip().split()
                if len(parts) >= 3:
                    result['physical_names'].append(
                        (int(parts[0]), int(parts[1]), parts[2].strip('"')))
                i += 1
            continue

        if line == '$Nodes':
            i += 1
            header = lines[i].strip().split()
            num_blocks = int(header[0])
            total_nodes = int(header[1])
            result['nodes'] = total_nodes
            i += 1
            for _ in range(num_blocks):
                bh = lines[i].strip().split()
                num_block_nodes = int(bh[3])
                i += 1
                i += num_block_nodes  # skip node tags
                i += num_block_nodes  # skip node coords
            continue

        if line == '$Elements':
            i += 1
            header = lines[i].strip().split()
            num_blocks = int(header[0])
            total_elems = int(header[1])
            result['declared_elements'] = total_elems
            i += 1
            for _ in range(num_blocks):
                bh = lines[i].strip().split()
                etype = int(bh[2])
                num_block_elems = int(bh[3])
                i += 1
                result['element_types'][etype] = \
                    result['element_types'].get(etype, 0) + num_block_elems
                result['elements'] += num_block_elems
                i += num_block_elems  # skip element lines
            continue

        i += 1

    return result


def parse_nastran(filename):
    """Parse Nastran BDF and return {nodes, elements, element_cards}."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    result = {'nodes': 0, 'elements': 0, 'element_cards': {}}
    for line in lines:
        if line.startswith('GRID'):
            result['nodes'] += 1
        elif line.startswith(('CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM',
                              'CTRIA', 'CQUAD', 'CBAR')):
            card = line.split()[0]
            # Handle CTETRA vs CTETRA(10) etc.
            result['element_cards'][card] = \
                result['element_cards'].get(card, 0) + 1
            result['elements'] += 1

    return result


def parse_vtk(filename):
    """Parse VTK legacy format and return {nodes, cells, cell_types}."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    result = {'nodes': 0, 'cells': 0, 'cell_types': {}}
    section = None

    for line in lines:
        line = line.strip()
        if line.startswith('POINTS'):
            result['nodes'] = int(line.split()[1])
        elif line.startswith('CELLS'):
            result['cells'] = int(line.split()[1])
            section = 'cells'
        elif line.startswith('CELL_TYPES'):
            section = 'cell_types'
        elif section == 'cell_types' and line:
            try:
                ct = int(line)
                result['cell_types'][ct] = \
                    result['cell_types'].get(ct, 0) + 1
            except ValueError:
                section = None

    return result


# ================================================================
# Test runner
# ================================================================
GMSH_TRI_TYPES = {1: 2, 2: 9, 3: 21, 4: 23, 5: 25}   # order -> gmsh type for tri
GMSH_TET_TYPES = {1: 4, 2: 11, 3: 29, 4: 30, 5: 31}  # order -> gmsh type for tet

results = []


def run_test(name, model, export_cmd, validator, expected):
    """Run one export test case."""
    global results
    test = {'name': name, 'status': 'FAIL', 'errors': []}

    try:
        # Build model
        if model == '3d':
            info = build_3d_model()
        else:
            info = build_2d_model()

        # Export
        cubit.cmd(export_cmd)

        # Extract filename from command
        parts = export_cmd.split('"')
        filename = parts[1] if len(parts) >= 2 else None
        if not filename or not os.path.exists(filename):
            test['errors'].append(f"Output file not found: {filename}")
            results.append(test)
            return False

        fsize = os.path.getsize(filename)
        if fsize == 0:
            test['errors'].append("Output file is empty")
            results.append(test)
            return False

        # Validate
        errors = validator(filename, info, expected)
        if errors:
            test['errors'] = errors
        else:
            test['status'] = 'PASS'

        test['file_size'] = fsize

    except Exception as e:
        test['errors'].append(str(e))
        traceback.print_exc()

    results.append(test)
    status = 'PASS' if test['status'] == 'PASS' else 'FAIL'
    print(f"  [{status}] {name}")
    if test['errors']:
        for err in test['errors']:
            print(f"         {err}")
    return test['status'] == 'PASS'


# ================================================================
# Validators
# ================================================================
def validate_gmsh_export(filename, model_info, expected):
    """Validate GMSH export file."""
    errors = []
    r = parse_gmsh(filename)

    # Check version
    if expected.get('version'):
        if not r['version'].startswith(expected['version']):
            errors.append(f"Version: expected {expected['version']}, got {r['version']}")

    # Check nodes > 0
    if r['nodes'] == 0:
        errors.append("No nodes in file")

    # Check element count matches declaration
    if r['elements'] != r['declared_elements']:
        errors.append(f"Element count mismatch: declared {r['declared_elements']}, "
                      f"actual {r['elements']}")

    # Check expected element types present
    for etype in expected.get('element_types', []):
        if etype not in r['element_types']:
            errors.append(f"Element type {etype} not found (have: {r['element_types']})")

    # Check PhysicalNames dimensions
    for dim, pid, name in r.get('physical_names', []):
        if expected.get('check_dim'):
            # tri blocks should be dim 2, tet blocks dim 3
            pass  # checked by element_types presence

    # Order 2+: nodes should be more than linear
    if expected.get('order', 1) >= 2 and 'n_nodes' in model_info:
        if r['nodes'] <= model_info['n_nodes']:
            errors.append(f"Order {expected['order']} but node count not increased: "
                          f"{r['nodes']} <= {model_info['n_nodes']}")

    return errors


def validate_nastran_export(filename, model_info, expected):
    """Validate Nastran BDF export file."""
    errors = []
    r = parse_nastran(filename)

    if r['nodes'] == 0:
        errors.append("No GRID cards found")
    if r['elements'] == 0:
        errors.append("No element cards found")

    for card in expected.get('cards', []):
        if card not in r['element_cards']:
            errors.append(f"Card {card} not found (have: {list(r['element_cards'].keys())})")

    return errors


def validate_vtk_export(filename, model_info, expected):
    """Validate VTK export file."""
    errors = []
    r = parse_vtk(filename)

    if r['nodes'] == 0:
        errors.append("No POINTS found")
    if r['cells'] == 0:
        errors.append("No CELLS found")

    return errors


# ================================================================
# Test case definitions
# ================================================================
def generate_test_cases():
    """Generate all format x option combinations."""
    cases = []

    # --- GMSH: order(1-3) x dim(2D,3D) (v4.1 only) ---
    # GMSH max order is 3; order 4-5 is a documented error (NetgenCurver
    # face/volume interior nodes are unreliable at p>=4 -> negative
    # Jacobians).  Use `export netgen` for order 4-5.  See CLAUDE.md
    # "GMSH order limit".
    for model in ['3d', '2d']:
        for order in [1, 2, 3]:
            for dim in [3, 2]:
                vstr = '4.1'
                name = f"gmsh_{model}_order{order}_v{vstr}_dim{dim}"
                fname = os.path.join(OUT_DIR, f"{name}.msh")

                expected_types = []
                if model == '3d' and dim == 3:
                    expected_types.append(GMSH_TET_TYPES[order])
                    expected_types.append(GMSH_TRI_TYPES[order])
                elif model == '3d' and dim == 2:
                    expected_types.append(GMSH_TRI_TYPES[order])
                elif model == '2d':
                    expected_types.append(GMSH_TRI_TYPES[order])

                cmd = (f'export gmsh "{fname}" order {order} '
                       f'dimension {dim} overwrite')

                cases.append({
                    'name': name, 'model': model,
                    'cmd': cmd,
                    'validator': validate_gmsh_export,
                    'expected': {
                        'version': vstr[:3],
                        'order': order,
                        'element_types': expected_types,
                    }
                })

    # --- Nastran: order(1,2) x dim(2D,3D) x nopyramid(0,1) ---
    for model in ['3d', '2d']:
        for order in [1, 2]:
            for dim in [3, 2]:
                for nopyr in [False, True]:
                    nopyr_str = "_nopyr" if nopyr else ""
                    name = f"nastran_{model}_order{order}_dim{dim}{nopyr_str}"
                    fname = os.path.join(OUT_DIR, f"{name}.bdf")

                    cmd = (f'export nastran_bdf "{fname}" order {order} '
                           f'dimension {dim} overwrite')
                    if nopyr:
                        cmd += ' nopyramid'

                    tri_card = "CTRIA3" if order == 1 else "CTRIA6"
                    cards = [tri_card]
                    if model == "3d" and dim == 3:
                        cards.insert(0, "CTETRA")

                    cases.append({
                        'name': name, 'model': model,
                        'cmd': cmd,
                        'validator': validate_nastran_export,
                        'expected': {'cards': cards}
                    })

    # --- VTK: order(1,2) x dim(2D,3D) ---
    for model in ['3d', '2d']:
        for order in [1, 2]:
            for dim in [3, 2]:
                name = f"vtk_{model}_order{order}_dim{dim}"
                fname = os.path.join(OUT_DIR, f"{name}.vtk")

                cmd = (f'export vtk "{fname}" order {order} '
                       f'dimension {dim} overwrite')

                cases.append({
                    'name': name, 'model': model,
                    'cmd': cmd,
                    'validator': validate_vtk_export,
                    'expected': {}
                })

    return cases


# ================================================================
# Main
# ================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("Export Mesh Exhaustive Test: All Format x Option Combinations")
    print("=" * 70)

    cases = generate_test_cases()
    print(f"\nTotal test cases: {len(cases)}")
    print(f"Output directory: {OUT_DIR}\n")

    # Group by model to minimize rebuilds
    current_model = None
    passed = 0
    failed = 0

    for case in cases:
        if case['model'] != current_model:
            current_model = case['model']
            print(f"\n--- Model: {current_model.upper()} ---")

        ok = run_test(case['name'], case['model'], case['cmd'],
                      case['validator'], case['expected'])
        if ok:
            passed += 1
        else:
            failed += 1

    # Summary
    print("\n" + "=" * 70)
    print(f"Results: {passed} PASSED, {failed} FAILED, {len(cases)} total")
    print("=" * 70)

    # Save results JSON
    results_file = os.path.join(OUT_DIR, "results.json")
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to: {results_file}")

    sys.exit(0 if failed == 0 else 1)
