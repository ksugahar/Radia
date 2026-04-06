"""
Accelerator magnet MSC solver for Cubit panel (ElectromagnetMSCDialog).

Pipeline:
  Cubit hex mesh (.cub5) -> Cubit API direct -> rad.ObjHexahedron -> Radia MSC Solve

User workflow:
  1. Write coil Python script (defines build_coil() -> CoilBuilder)
  2. Write .jou to create iron yoke hex mesh + named blocks
  3. Run journal in Cubit panel, click Solve

Cubit blocks (user sets in .jou):
  yoke  - volume block, hex elements (iron)

IMA symmetry (set in dialog):
  '+x-z'  = quarter model (x>0, z>0)
  '+x'    = half model (x>0)
  '-z'    = half model (z>0)
  ''      = full model (no IMA)

Coil is NOT meshed -- Biot-Savart analytical source via CoilBuilder.to_radia().

Mesh extraction uses Cubit API directly (no export_Gmesh, no .msh intermediate):
  cubit.get_block_volumes() -> cubit.get_volume_hexes/tets/wedges()
  -> cubit.get_connectivity() -> cubit.get_nodal_coordinates()
"""

import argparse
import importlib.util
import math
import os
import sys
import time

import numpy as np

# Shared utilities
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import (MU_0, setup_paths, setup_cubit,
                          get_bh_curve, progress, calc_main)


def _log(msg):
    """Write progress to stderr (panel reads these)."""
    progress("MSC", msg)


def _load_coil_script(script_path):
    """Load coil script and call build_coil() -> CoilBuilder."""
    if not script_path or not os.path.exists(script_path):
        raise FileNotFoundError(f"Coil script not found: {script_path}")

    script_dir = os.path.dirname(os.path.abspath(script_path))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    setup_paths()

    spec = importlib.util.spec_from_file_location("coil_module", script_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    if not hasattr(mod, 'build_coil'):
        raise AttributeError(
            f"Coil script must define build_coil() function: {script_path}")

    return mod.build_coil()


def _extract_elements_from_cubit(cubit, block_name="yoke", unit_scale=0.001):
    """Extract hex/tet/wedge elements from Cubit block via direct API.

    Args:
        cubit: Cubit module (initialized, model open)
        block_name: Name of the block containing yoke elements
        unit_scale: Coordinate scale (0.001 = mm->m, 1.0 = already meters)

    Returns:
        list of dicts: [{'type': 'hex'|'tet'|'wedge', 'vertices': [[x,y,z],...]}]
    """
    # Find block ID by name
    yoke_bid = None
    for bid in cubit.get_block_id_list():
        name = cubit.get_exodus_entity_name("block", bid).lower()
        if name == block_name.lower():
            yoke_bid = bid
            break

    if yoke_bid is None:
        raise ValueError(f"Block '{block_name}' not found in Cubit model")

    vids = cubit.get_block_volumes(yoke_bid)
    if not vids:
        raise ValueError(f"Block '{block_name}' has no volumes")

    elements = []
    for vid in vids:
        # Hexahedra
        for eid in cubit.get_volume_hexes(vid):
            nids = cubit.get_connectivity("hex", eid)
            verts = [[c * unit_scale for c in cubit.get_nodal_coordinates(n)]
                     for n in nids]
            elements.append({'type': 'hex', 'vertices': verts})

        # Tetrahedra
        for eid in cubit.get_volume_tets(vid):
            nids = cubit.get_connectivity("tet", eid)
            verts = [[c * unit_scale for c in cubit.get_nodal_coordinates(n)]
                     for n in nids]
            elements.append({'type': 'tet', 'vertices': verts})

        # Wedges (prisms)
        for eid in cubit.get_volume_wedges(vid):
            nids = cubit.get_connectivity("wedge", eid)
            verts = [[c * unit_scale for c in cubit.get_nodal_coordinates(n)]
                     for n in nids]
            elements.append({'type': 'wedge', 'vertices': verts})

    return elements


def solve_msc(coil_script="", cub5_file="",
              material="steel", mu_r=1000.0, bh_file=None,
              ima="", solver=0,
              max_iter=100, tol=1e-3, relax=0.0,
              unit_scale=0.001,
              msh_output=""):
    """MSC solver: Cubit hex -> Radia ObjHexahedron -> Solve.

    Args:
        coil_script: Path to Python script with build_coil()
        cub5_file: Cubit .cub5 file with hex mesh
        material: "steel" (BH curve) or "linear" (mu_r)
        mu_r: Relative permeability (linear mode only)
        bh_file: Custom BH curve file (H [A/m], B [T])
        ima: IMA string e.g. '+x-z' for quarter model ('' = no IMA)
        solver: 0=LU, 1=BiCGSTAB, 2=HACApK
        max_iter: Max nonlinear iterations
        tol: Convergence tolerance
        relax: Under-relaxation (0=full step)
        unit_scale: Cubit coordinate scale (0.001=mm->m, 1.0=meters)
        msh_output: Optional GMSH .msh output path

    Returns:
        dict with B_center, n_elem, ndof, iterations, converged, etc.
    """
    setup_paths()
    import radia as rad

    t_total_start = time.perf_counter()

    # ============================================================
    # Step 1: Load coil -> Radia objects (exact analytical)
    # ============================================================
    _log("COIL:loading script")
    coil = _load_coil_script(coil_script)
    current = coil.current
    radia_coil_objs = coil.to_radia()
    coil_container = rad.ObjCnt(radia_coil_objs)
    _log(f"COIL:{len(radia_coil_objs)} Radia objects, I={current}A")

    # ============================================================
    # Step 2: Extract elements from Cubit via direct API
    # ============================================================
    _log("MESH:loading Cubit model")

    if not cub5_file or not os.path.exists(cub5_file):
        return {"error": f"Cubit .cub5 file required: {cub5_file}"}

    cubit = setup_cubit(cub5_file)
    if cubit is None:
        return {"error": "Cubit not available"}

    try:
        elements = _extract_elements_from_cubit(
            cubit, block_name="yoke", unit_scale=unit_scale)
    except ValueError as e:
        return {"error": str(e)}

    n_hex = sum(1 for e in elements if e['type'] == 'hex')
    n_tet = sum(1 for e in elements if e['type'] == 'tet')
    n_wedge = sum(1 for e in elements if e['type'] == 'wedge')
    n_total = len(elements)
    if n_total == 0:
        return {"error": "No volume elements found in 'yoke' block"}

    _log(f"MESH:{n_hex} hex + {n_tet} tet + {n_wedge} wedge = {n_total} elem")

    # ============================================================
    # Step 3: Create Radia objects with material
    # ============================================================
    rad.UtiDelAll()
    _log("RADIA:creating objects")

    is_linear = (material == "linear")
    if is_linear:
        mat = rad.MatLin(mu_r)
        _log(f"MAT:linear mu_r={mu_r}")
    else:
        bh_data, _ = get_bh_curve(material, bh_file)
        if bh_data is None:
            return {"error": f"No BH curve for material: {material}"}
        mat = rad.MatSatIsoTab(bh_data)
        _log(f"MAT:nonlinear BH ({len(bh_data)} pts)")

    all_objs = []
    for elem in elements:
        verts = elem['vertices']
        etype = elem['type']
        if etype == 'hex':
            obj = rad.ObjHexahedron(verts, [0, 0, 0])
        elif etype == 'tet':
            obj = rad.ObjTetrahedron(verts, [0, 0, 0])
        elif etype == 'wedge':
            obj = rad.ObjWedge(verts, [0, 0, 0])
        else:
            continue
        rad.MatApl(obj, mat)
        all_objs.append(obj)

    yoke = rad.ObjCnt(all_objs)
    n_dof = n_hex * 6 + n_tet * 3 + n_wedge * 5
    _log(f"RADIA:{len(all_objs)} objects, {n_dof} DOF")

    # ============================================================
    # Step 4: Assemble model (yoke + coil) and solve
    # ============================================================
    model = rad.ObjCnt([yoke, coil_container])

    solver_names = {0: "LU", 1: "BiCGSTAB", 2: "HACApK"}
    ima_str = ima if ima else "(none)"
    _log(f"SOLVE:solver={solver_names.get(solver, solver)}, "
         f"IMA={ima_str}, tol={tol}, maxiter={max_iter}")

    if relax > 0:
        rad.SolverConfig(relax_param=relax)

    t_solve_start = time.perf_counter()
    if ima:
        result = rad.Solve(model, tol, max_iter, solver, image=ima)
    else:
        result = rad.Solve(model, tol, max_iter, solver)
    t_solve = time.perf_counter() - t_solve_start

    if relax > 0:
        rad.SolverConfig(relax_param=0.0)

    n_iter = int(result[3]) if len(result) > 3 else 0
    residual = float(result[1]) if len(result) > 1 else 0.0
    converged = residual < tol

    _log(f"SOLVE:done in {t_solve:.2f}s, {n_iter} iter, "
         f"residual={residual:.2e}")

    # ============================================================
    # Step 5: Evaluate field at origin
    # ============================================================
    B = rad.Fld(model, 'b', [0, 0, 0])
    B_mag = math.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
    _log(f"FIELD:B(0)=({B[0]*1e3:.2f}, {B[1]*1e3:.2f}, {B[2]*1e3:.2f}) mT, "
         f"|B|={B_mag*1e3:.2f} mT")

    # ============================================================
    # Step 6: Optional GMSH output (element-center B field)
    # ============================================================
    gmsh_file = ""
    if msh_output:
        try:
            _write_msc_gmsh(msh_output, elements, model, rad)
            gmsh_file = msh_output
            _log(f"GMSH:{msh_output}")
        except Exception as e:
            _log(f"GMSH:export failed: {e}")

    t_total = time.perf_counter() - t_total_start

    result_dict = {
        "B_center": [float(B[0]), float(B[1]), float(B[2])],
        "B_center_mag": float(B_mag),
        "n_elem": n_total,
        "n_hex": n_hex,
        "n_tet": n_tet,
        "n_wedge": n_wedge,
        "ndof": n_dof,
        "solver": solver_names.get(solver, str(solver)),
        "ima": ima,
        "iterations": n_iter,
        "residual": residual,
        "converged": converged,
        "material": material,
        "mu_r": mu_r if is_linear else None,
        "current": current,
        "t_solve": t_solve,
        "t_total": t_total,
        "msh_file": gmsh_file,
    }

    rad.UtiDelAll()
    return result_dict


def _write_msc_gmsh(filename, elements, model, rad):
    """Write GMSH .msh v2.2 with element-center B field as ElementData.

    Builds node/element tables directly from the Radia element list
    (no dependency on external mesh import modules).
    """
    # Collect unique nodes and build element connectivity
    node_map = {}  # (x,y,z) rounded -> node_id
    node_coords = {}  # node_id -> (x, y, z)
    elem_connectivity = []  # [(gmsh_type, [node_ids])]
    next_nid = 1

    gmsh_types = {'hex': 5, 'tet': 4, 'wedge': 6}

    for elem in elements:
        nids = []
        for v in elem['vertices']:
            key = (round(v[0], 10), round(v[1], 10), round(v[2], 10))
            if key not in node_map:
                node_map[key] = next_nid
                node_coords[next_nid] = (v[0], v[1], v[2])
                next_nid += 1
            nids.append(node_map[key])
        elem_connectivity.append((gmsh_types[elem['type']], nids))

    # Compute B at element centroids
    B_data = []
    for elem in elements:
        verts = elem['vertices']
        n = len(verts)
        cx = sum(v[0] for v in verts) / n
        cy = sum(v[1] for v in verts) / n
        cz = sum(v[2] for v in verts) / n
        B = rad.Fld(model, 'b', [cx, cy, cz])
        B_data.append(B)

    # Write GMSH v2.2
    with open(filename, 'w') as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        # Nodes
        f.write(f"$Nodes\n{len(node_coords)}\n")
        for nid in sorted(node_coords.keys()):
            x, y, z = node_coords[nid]
            f.write(f"{nid} {x} {y} {z}\n")
        f.write("$EndNodes\n")

        # Elements
        f.write(f"$Elements\n{len(elem_connectivity)}\n")
        for i, (etype, nids) in enumerate(elem_connectivity):
            eid = i + 1
            nids_str = " ".join(str(n) for n in nids)
            f.write(f"{eid} {etype} 2 1 1 {nids_str}\n")
        f.write("$EndElements\n")

        # ElementData: B vector
        f.write("$ElementData\n")
        f.write("1\n\"B [T]\"\n")
        f.write("1\n0.0\n")
        f.write("3\n0\n3\n")
        f.write(f"{len(elements)}\n")
        for i, B in enumerate(B_data):
            f.write(f"{i+1} {B[0]} {B[1]} {B[2]}\n")
        f.write("$EndElementData\n")

        # ElementData: |B| scalar
        f.write("$ElementData\n")
        f.write("1\n\"|B| [T]\"\n")
        f.write("1\n0.0\n")
        f.write("3\n0\n1\n")
        f.write(f"{len(elements)}\n")
        for i, B in enumerate(B_data):
            Bmag = math.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
            f.write(f"{i+1} {Bmag}\n")
        f.write("$EndElementData\n")


def main():
    parser = argparse.ArgumentParser(
        description="Accelerator magnet MSC solver (Radia hex)")
    parser.add_argument("--coil-script", required=True,
                        help="Python script with build_coil() -> CoilBuilder")
    parser.add_argument("--cub5", default="", help="Cubit .cub5 file")
    parser.add_argument("--material", default="steel",
                        choices=["steel", "elf_steel", "linear"],
                        help="steel (BH curve) or linear (mu_r)")
    parser.add_argument("--mu-r", type=float, default=1000,
                        help="mu_r (linear mode only)")
    parser.add_argument("--bh-file", default="",
                        help="Custom BH curve file (H B)")
    parser.add_argument("--ima", default="",
                        help="IMA string: '+x-z' (quarter), '+x' (half-x), "
                             "'-z' (half-z), '' (full)")
    parser.add_argument("--solver", type=int, default=0,
                        choices=[0, 1, 2],
                        help="0=LU, 1=BiCGSTAB, 2=HACApK")
    parser.add_argument("--max-iter", type=int, default=100,
                        help="Max nonlinear iterations")
    parser.add_argument("--tol", type=float, default=1e-3,
                        help="Convergence tolerance")
    parser.add_argument("--relax", type=float, default=0.0,
                        help="Under-relaxation (0=full step)")
    parser.add_argument("--unit-scale", type=float, default=0.001,
                        help="Cubit coordinate scale (0.001=mm->m)")
    parser.add_argument("--msh-output", default="",
                        help="GMSH .msh output path")
    parser.add_argument("--output", default="",
                        help="JSON output file")

    def run(args):
        return solve_msc(
            coil_script=args.coil_script,
            cub5_file=args.cub5,
            material=args.material,
            mu_r=args.mu_r,
            bh_file=args.bh_file if args.bh_file else None,
            ima=args.ima,
            solver=args.solver,
            max_iter=args.max_iter,
            tol=args.tol,
            relax=args.relax,
            unit_scale=args.unit_scale,
            msh_output=args.msh_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
