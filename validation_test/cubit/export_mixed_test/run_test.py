"""Test Cubit mesh export: mixed elements, hex cylinder, p-convergence.

Runs journal files in Cubit batch mode and verifies output files.
Requires Coreform Cubit with the Radia plugin installed.

Usage:
    python run_test.py
"""

import os
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def find_cubit():
    """Find Cubit executable."""
    candidates = [
        r"C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.exe",
        r"C:\Program Files\Coreform Cubit 2025.12\coreform_cubit.exe",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return None


def run_jou(cubit_exe, jou_template, out_dir, expected_files):
    """Run a .jou template and verify output files. Returns (passed, failed)."""
    os.makedirs(out_dir, exist_ok=True)

    with open(jou_template, "r") as f:
        content = f.read()
    content = content.replace("{DIR}", out_dir.replace("\\", "/"))

    tmp_jou = os.path.join(out_dir, "_test.jou")
    with open(tmp_jou, "w") as f:
        f.write(content)

    result = subprocess.run(
        [cubit_exe, "-batch", "-nographics", "-nojournal", tmp_jou],
        capture_output=True, text=True, timeout=300
    )

    # Print important lines
    if result.stdout:
        for line in result.stdout.split("\n"):
            line = line.strip()
            if any(k in line.lower() for k in
                   ["error", "netgencurver", "exported", "build_high",
                    "interrupt"]):
                print(f"  {line}")

    if result.returncode != 0:
        print(f"  NOTE: Cubit exited with code {result.returncode}")

    # Verify output files
    passed = 0
    failed = 0
    for fname, min_size in expected_files:
        fpath = os.path.join(out_dir, fname)
        if not os.path.isfile(fpath):
            print(f"  FAIL: {fname} not created")
            failed += 1
        elif os.path.getsize(fpath) < min_size:
            print(f"  FAIL: {fname} too small ({os.path.getsize(fpath)} bytes)")
            failed += 1
        else:
            size = os.path.getsize(fpath)
            print(f"  OK:   {fname} ({size:,} bytes)")
            passed += 1

    try:
        os.remove(tmp_jou)
    except OSError:
        pass

    return passed, failed


def main():
    cubit_exe = find_cubit()
    if not cubit_exe:
        print("SKIP: Cubit not found")
        return 0

    out_dir = os.path.join(SCRIPT_DIR, "output")
    total_passed = 0
    total_failed = 0

    # --- Test 1: Mixed hex+tet+pyramid (order 1 and 2) ---
    print("=== Test 1: Mixed hex+tet+pyramid ===")
    p, f = run_jou(cubit_exe,
        os.path.join(SCRIPT_DIR, "mixed_hex_tet_pyramid.jou"),
        out_dir,
        [
            ("mixed_o1.msh", 1000),
            ("mixed_o1.bdf", 500),
            ("mixed_o1.vol", 500),
            ("mixed_o2.msh", 1000),
            ("mixed_o2.bdf", 500),
            ("mixed_o2.vol", 500),
        ])
    total_passed += p
    total_failed += f

    # --- Test 2: Hex cylinder (order 2, HEX20 connectivity) ---
    print("=== Test 2: Hex cylinder order 2 ===")
    p, f = run_jou(cubit_exe,
        os.path.join(SCRIPT_DIR, "hex_cylinder.jou"),
        out_dir,
        [
            ("hex_cyl_o2.msh", 1000),
            ("hex_cyl_o2.bdf", 1000),
        ])
    total_passed += p
    total_failed += f

    # Verify HEX20 has 20 nodes (not 8)
    msh_path = os.path.join(out_dir, "hex_cyl_o2.msh")
    if os.path.isfile(msh_path):
        with open(msh_path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) > 5 and parts[1] == "17":  # GMSH HEX20
                    n_nodes = len(parts) - 5
                    if n_nodes == 20:
                        print(f"  OK:   HEX20 connectivity = {n_nodes} nodes")
                        total_passed += 1
                    else:
                        print(f"  FAIL: HEX20 has {n_nodes} nodes (expected 20)")
                        total_failed += 1
                    break

    # --- Test 3: p-convergence (sphere, order 1-5 via export netgen) ---
    print("=== Test 3: p-convergence (sphere order 1-5) ===")
    p_conv_jou = os.path.join(SCRIPT_DIR, "sphere_p_convergence.jou")
    if os.path.isfile(p_conv_jou):
        expected = []
        for order in range(1, 6):
            expected.append((f"sphere_o{order}.vol", 500))
        p, f = run_jou(cubit_exe, p_conv_jou, out_dir, expected)
        total_passed += p
        total_failed += f
    else:
        print("  SKIP: sphere_p_convergence.jou not found")

    print(f"\n=== TOTAL: {total_passed} passed, {total_failed} failed ===")
    return 1 if total_failed > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
