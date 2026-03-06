#!/usr/bin/env python
"""
Compare three solvers (LU, BiCGSTAB, HACApK) for nonlinear electromagnet.

Model: C-type electromagnet with saturable iron (quarter model)
- 13 hexahedral elements with MIMA X + MIMA -Z
- Nonlinear B-H curve material
- 20000 AT racetrack coil (full model)

ELF Reference: Bz = -994.16 mT at (0,0,0)
"""

import sys
import os
import time
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
nonlinear_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(nonlinear_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, nonlinear_dir)

import radia as rad

ELF_NONLINEAR = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_1x1x1\quater"

# Import the standard coil model
from coil_model import create_racetrack_coil

scale = 0.001  # mm to m


def load_elf_geometry(path):
    """Load ELF geometry from .meg file."""
    nodes = {}
    elements = []
    with open(os.path.join(path, "ELF_magic.meg"), 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('MGR1'):
                parts = line.split()
                node_id = int(parts[1])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                nodes[node_id] = np.array([x, y, z])
            elif line.startswith('MMB8T'):
                parts = line.split()
                elem_id = int(parts[1])
                node_ids = [int(parts[i]) for i in range(4, 12)]
                elements.append((elem_id, node_ids))
    return nodes, elements


def load_bh_curve(filepath):
    """Load B-H curve from text file."""
    data = np.loadtxt(filepath, comments='#')
    return data.tolist()


def load_elf_field_at_origin(path):
    """Load ELF field at (0,0,0) from .mag file."""
    with open(os.path.join(path, "ELF_magic.mag"), 'r') as f:
        for line in f:
            if line.startswith('M3GB'):
                parts = line.split()
                if len(parts) >= 8:
                    point_id = int(parts[1])
                    if point_id == 1:
                        Bx = float(parts[5])
                        By = float(parts[6])
                        Bz = float(parts[7])
                        return np.array([Bx, By, Bz])
    return None


def build_model(nodes, elements, bh_data):
    """Build Radia model with yoke and coil."""
    rad.UtiDelAll()

    # Create nonlinear material
    mat = rad.MatSatIsoTab(bh_data)

    # Create yoke geometry (quarter)
    hex_objects = []
    for elem_id, node_ids in elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(hex_obj, mat)
        hex_objects.append(hex_obj)

    yoke = rad.ObjCnt(hex_objects)

    # Create racetrack coil (standard model from coil_model.py)
    coil = create_racetrack_coil(20000.0)

    # Combine model
    model = rad.ObjCnt([yoke, coil])

    return model, len(hex_objects)


def run_solver(model, method, precision=0.0001, max_iter=100):
    """Run solver and return results."""
    t_start = time.time()
    result = rad.Solve(model, precision, max_iter, method, image='+x-z')
    t_solve = time.time() - t_start

    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))

    return {
        'method': method,
        'max_M': result[0],
        'iterations': int(result[2]),
        'time': t_solve,
        'Bx': B[0],
        'By': B[1],
        'Bz': B[2],
        'B_mag': np.linalg.norm(B)
    }


def main():
    print("=" * 70)
    print("Solver Comparison: LU vs BiCGSTAB vs HACApK")
    print("=" * 70)
    print("Model: C-type electromagnet (quarter, nonlinear)")
    print("Coil:  Full racetrack, 20000 AT")
    print("Image: +x-z")
    print()

    # Load geometry
    nodes, elements = load_elf_geometry(ELF_NONLINEAR)
    n_elem = len(elements)
    n_dof = n_elem * 6
    print(f"Geometry: {n_elem} elements, {n_dof} DOF")

    # Load B-H curve
    bh_file = os.path.join(nonlinear_dir, "BH.txt")
    bh_data = load_bh_curve(bh_file)
    print(f"Material: {len(bh_data)} B-H points")

    # Load ELF reference
    elf_B = load_elf_field_at_origin(ELF_NONLINEAR)
    if elf_B is not None:
        print(f"ELF Bz:  {elf_B[2]*1000:.2f} mT")
    print()

    # Solver configurations
    solvers = [
        {'name': 'LU', 'method': 0, 'desc': 'Dense LU (LAPACK dgesv)'},
        {'name': 'BiCGSTAB', 'method': 1, 'desc': 'Iterative BiCGSTAB + Jacobi precond'},
        {'name': 'HACApK', 'method': 2, 'desc': 'H-matrix accelerated BiCGSTAB'},
    ]

    results = []

    for solver in solvers:
        print("-" * 70)
        print(f"Solver: {solver['name']} (Method {solver['method']})")
        print(f"  {solver['desc']}")
        print("-" * 70)

        # Build fresh model for each solver
        model, n_elem = build_model(nodes, elements, bh_data)

        # Set HACApK parameters if needed
        if solver['method'] == 2:
            rad.SetHACApKParams(1e-4, 10, 2.0)

        # Run solver
        result = run_solver(model, solver['method'])
        result['name'] = solver['name']
        results.append(result)

        print(f"  Max |M|:     {result['max_M']:.2f} A/m")
        print(f"  Iterations:  {result['iterations']}")
        print(f"  Time:        {result['time']:.3f} s")
        print(f"  Bz:          {result['Bz']*1000:.2f} mT")
        if elf_B is not None:
            error = abs(result['Bz'] - elf_B[2]) / abs(elf_B[2]) * 100
            print(f"  Error vs ELF: {error:.2f}%")
        print()

    # Summary table
    print("=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'Solver':<12} {'Method':<8} {'Iter':<6} {'Time(s)':<10} {'Bz(mT)':<12} {'Error(%)':<10}")
    print("-" * 70)

    for r in results:
        if elf_B is not None:
            error = abs(r['Bz'] - elf_B[2]) / abs(elf_B[2]) * 100
        else:
            error = 0.0
        print(f"{r['name']:<12} {r['method']:<8} {r['iterations']:<6} {r['time']:<10.3f} {r['Bz']*1000:<12.2f} {error:<10.2f}")

    print("-" * 70)
    if elf_B is not None:
        print(f"{'ELF ref':<12} {'-':<8} {'-':<6} {'-':<10} {elf_B[2]*1000:<12.2f} {0.0:<10.2f}")
    print("=" * 70)

    # Check consistency
    print("\nSOLVER CONSISTENCY CHECK:")
    bz_values = [r['Bz'] for r in results]
    bz_spread = (max(bz_values) - min(bz_values)) / abs(np.mean(bz_values)) * 100
    print(f"  Bz spread between solvers: {bz_spread:.4f}%")
    if bz_spread < 0.1:
        print("  -> EXCELLENT: All solvers agree within 0.1%")
    elif bz_spread < 1.0:
        print("  -> GOOD: All solvers agree within 1%")
    else:
        print("  -> WARNING: Significant difference between solvers")


if __name__ == '__main__':
    main()
