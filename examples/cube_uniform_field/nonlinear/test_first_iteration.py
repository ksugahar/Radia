#!/usr/bin/env python
"""
Test first iteration of LU vs HACApK to verify matrix element consistency.

This script runs each solver in a SEPARATE PROCESS to ensure clean state.
"""

import sys
import os
import subprocess
import json

MU_0 = 4 * 3.14159265358979 * 1e-7
H_EXT = 200000.0  # Same as benchmark


def run_solver_subprocess(solver_method, n_div=5, max_iter=100):
    """Run solver in a separate subprocess for clean state."""

    script = f'''
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
import radia as rad
import json

MU_0 = 4 * 3.14159265358979 * 1e-7
H_EXT = {H_EXT}

BH_DATA = [
    [0.0, 0.0], [100.0, 0.1], [200.0, 0.3], [500.0, 0.8], [1000.0, 1.2],
    [2000.0, 1.5], [5000.0, 1.7], [10000.0, 1.8], [50000.0, 2.0], [100000.0, 2.1],
]

HEX_FACES = [
    [1, 4, 3, 2], [5, 6, 7, 8], [1, 2, 6, 5], [3, 4, 8, 7], [1, 5, 8, 4], [2, 3, 7, 6]
]

rad.FldUnits('m')
rad.UtiDelAll()

n_div = {n_div}
cube_size = 1.0 / n_div
elements = []

for ix in range(n_div):
    for iy in range(n_div):
        for iz in range(n_div):
            x0 = ix * cube_size - 0.5
            y0 = iy * cube_size - 0.5
            z0 = iz * cube_size - 0.5
            vertices = [
                [x0, y0, z0], [x0 + cube_size, y0, z0],
                [x0 + cube_size, y0 + cube_size, z0], [x0, y0 + cube_size, z0],
                [x0, y0, z0 + cube_size], [x0 + cube_size, y0, z0 + cube_size],
                [x0 + cube_size, y0 + cube_size, z0 + cube_size], [x0, y0 + cube_size, z0 + cube_size]
            ]
            obj = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 0])
            elements.append(obj)

container = rad.ObjCnt(elements)
mat = rad.MatSatIsoTab(BH_DATA)
rad.MatApl(container, mat)

B_ext = MU_0 * H_EXT
ext = rad.ObjBckg([0, 0, B_ext])
grp = rad.ObjCnt([container, ext])

solver_method = {solver_method}
if solver_method == 2:
    try:
        rad.SetHACApKParams(1e-8, 10, 2.0)
    except:
        pass

result = rad.Solve(grp, 0.001, {max_iter}, solver_method)

M_data = rad.ObjM(container)
M_z_sum = 0.0
count = 0
for elem_data in M_data:
    M = elem_data[1]
    if isinstance(M, (list, tuple)) and len(M) >= 3:
        M_z_sum += M[2]
        count += 1

M_avg_z = M_z_sum / count if count > 0 else 0.0
n_iter = int(result[3]) if result[3] else 0
residual = result[0] if result[0] else 0.0

print(json.dumps({{"M_avg_z": M_avg_z, "n_iter": n_iter, "residual": residual}}))
'''

    # Write temp script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    temp_script = os.path.join(script_dir, '_temp_solver.py')
    with open(temp_script, 'w') as f:
        f.write(script)

    # Run in subprocess
    try:
        result = subprocess.run(
            ['python', temp_script],
            capture_output=True,
            text=True,
            cwd=script_dir,
            timeout=120
        )

        # Parse result from last line (JSON)
        output_lines = result.stdout.strip().split('\n')
        for line in reversed(output_lines):
            if line.startswith('{'):
                return json.loads(line)

        print(f"STDERR: {result.stderr}")
        return None
    except Exception as e:
        print(f"Error running subprocess: {e}")
        return None
    finally:
        if os.path.exists(temp_script):
            os.remove(temp_script)


def main():
    print('=' * 70)
    print('TEST: LU vs HACApK Comparison (Subprocess isolation)')
    print('=' * 70)
    print(f'H_ext = {H_EXT} A/m')

    n_div = 5
    print(f'\nMesh: {n_div}x{n_div}x{n_div} = {n_div**3} elements')

    # Run LU in subprocess
    print('\nRunning LU solver...')
    result_lu = run_solver_subprocess(0, n_div=n_div, max_iter=100)
    if result_lu:
        print(f"  LU: M_avg_z = {result_lu['M_avg_z']:.2f} A/m, iter = {result_lu['n_iter']}")

    # Run HACApK in subprocess
    print('\nRunning HACApK solver...')
    result_hacapk = run_solver_subprocess(2, n_div=n_div, max_iter=100)
    if result_hacapk:
        print(f"  HACApK: M_avg_z = {result_hacapk['M_avg_z']:.2f} A/m, iter = {result_hacapk['n_iter']}")

    if result_lu and result_hacapk:
        print('\n' + '=' * 70)
        print('COMPARISON')
        print('=' * 70)
        diff = result_hacapk['M_avg_z'] - result_lu['M_avg_z']
        pct = 100 * diff / result_lu['M_avg_z'] if result_lu['M_avg_z'] != 0 else 0
        print(f"  LU:     M_avg_z = {result_lu['M_avg_z']:.2f} A/m, iter = {result_lu['n_iter']}")
        print(f"  HACApK: M_avg_z = {result_hacapk['M_avg_z']:.2f} A/m, iter = {result_hacapk['n_iter']}")
        print(f"  Diff:   {diff:.2f} A/m ({pct:.4f}%)")


if __name__ == '__main__':
    main()
