"""Test BEM-SIBC workpiece solver with OCC-generated cylinder mesh.

Creates a cylinder workpiece via OCC (no Cubit needed), saves as .vol,
then runs calc_heating_bem.py via subprocess to verify the full pipeline.
"""

import json
import os
import subprocess
import sys
import tempfile

# Generate cylinder surface mesh via OCC
def generate_workpiece_vol(vol_path, R_wp=0.010, H_wp=0.020, maxh_factor=4):
    """Create cylindrical workpiece .vol via OCC."""
    from netgen.occ import Cylinder, Pnt, Dir, OCCGeometry
    from ngsolve import Mesh

    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    for f in wp_cyl.faces:
        f.name = "wp_surface"
    wp_cyl.name = "workpiece"

    geo = OCCGeometry(wp_cyl)
    maxh = R_wp / maxh_factor
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    with TaskManager():
        mesh.Curve(2)
        mesh.ngmesh.Save(vol_path)
        from ngsolve import BND
        from ngsolve import TaskManager
        print(f"Saved {vol_path}: {mesh.ne} vol elems, {mesh.GetNE(BND)} surf elems")
        return vol_path


def main():
    tmpdir = tempfile.mkdtemp(prefix="ih_wp_test_")
    vol_path = os.path.join(tmpdir, "workpiece.vol")

    print("=== Step 1: Generate workpiece .vol ===")
    generate_workpiece_vol(vol_path)

    print("\n=== Step 2: Run calc_heating_bem.py ===")
    calc_script = os.path.join(
        os.path.dirname(__file__), "..", "..", "..", "src", "radia", "panels",
        "calc_heating_bem.py")
    calc_script = os.path.abspath(calc_script)

    cmd = [
        sys.executable, calc_script,
        "--vol", vol_path,
        "--coil-radius", "0.030",
        "--coil-current", "1.0",
        "--gap-deg", "5",
        "--frequency", "7000",
        "--sigma", "2e6",
        "--material", "steel",
        "--half-thickness", "0.010",
    ]
    print(f"CMD: {' '.join(cmd)}")

    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    print(f"stderr:\n{proc.stderr}")

    if proc.returncode != 0:
        print(f"FAILED (returncode={proc.returncode})")
        print(f"stdout: {proc.stdout}")
        return 1

    # Parse JSON from last line of stdout
    stdout_lines = proc.stdout.strip().split("\n")
    result = json.loads(stdout_lines[-1])

    if "error" in result:
        print(f"ERROR: {result['error']}")
        return 1

    print("\n=== Results ===")
    for k, v in result.items():
        print(f"  {k}: {v}")

    # Basic sanity checks
    P = result["P_total_W"]
    H = result["H_t_rms_Am"]
    assert P > 0, f"P_total must be positive, got {P}"
    assert H > 0, f"H_t_rms must be positive, got {H}"
    assert result["ndof"] > 10, f"Too few DOFs: {result['ndof']}"
    assert result["n_iter"] > 0, f"Should take at least 1 iteration"
    assert result["Z_s_abs_Ohm"] > 0, f"Z_s must be positive"
    assert result["skin_depth_mm"] > 0, f"skin_depth must be positive"
    assert result["wp_mu_r"] > 1, f"steel mu_r should be > 1"

    # Step 3: Test with GMSH output
    print("\n=== Step 3: Test GMSH output ===")
    msh_path = os.path.join(tmpdir, "wp_bem.msh")
    cmd_msh = cmd + ["--msh-output", msh_path]
    proc2 = subprocess.run(cmd_msh, capture_output=True, text=True, timeout=300)
    assert proc2.returncode == 0, f"GMSH run failed: {proc2.stderr}"
    result2 = json.loads(proc2.stdout.strip().split("\n")[-1])
    assert os.path.exists(msh_path), "GMSH .msh file not created"
    assert os.path.getsize(msh_path) > 100, "GMSH .msh file too small"
    geo_path = os.path.splitext(msh_path)[0] + ".geo"
    assert os.path.exists(geo_path), "GMSH .geo companion not created"
    # Results should match within tolerance
    assert abs(result2["P_total_W"] - P) / P < 0.01, "P_total mismatch with GMSH"
    print(f"GMSH output OK: {os.path.getsize(msh_path)} bytes")

    # Step 4: Test copper material (non-magnetic)
    print("\n=== Step 4: Test copper material ===")
    cmd_cu = list(cmd)
    cmd_cu[cmd_cu.index("steel")] = "copper"
    cmd_cu[cmd_cu.index("2e6")] = "5.8e7"
    proc3 = subprocess.run(cmd_cu, capture_output=True, text=True, timeout=300)
    assert proc3.returncode == 0, f"Copper run failed: {proc3.stderr}"
    result3 = json.loads(proc3.stdout.strip().split("\n")[-1])
    assert result3["P_total_W"] > 0, "Copper P_total must be positive"
    assert result3["wp_mu_r"] <= 1.1, f"Copper mu_r should be ~1, got {result3['wp_mu_r']}"
    print(f"Copper: P={result3['P_total_W']:.4e} W, "
          f"mu_r={result3['wp_mu_r']:.1f}")

    print("\nAll checks passed!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
