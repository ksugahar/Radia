"""calc_streamfunction_volume.py -- 3D VOLUME stream-function coil designer.

CLI entry point for the radia-streamfunction "Volume 3D" panel mode.  Designs
equal-current solenoid windings inside a hollow cylindrical CONDUCTOR VOLUME by
the foliated current-Clebsch volume stream function

    J = grad(lambda) x grad(mu),   mu = r  (radial foliation, FIXED)

solved as a LINEAR least-norm inverse (ACA+TSVD) for a uniform axial Bz, then
contoured into windable equal-current wires whose Biot-Savart is cross-checked
against Radia (two-codebase invariant).  See
``radia.streamfunction_volume`` and
``examples/vim/README_streamfunction_cohomology.md``.

Pure headless: NGSolve / Radia imports live inside ``run()`` only, so the panel
can introspect ``build_argparser()`` without paying the NGSolve import cost.
"""
import argparse
import os
import sys

# Allow bare imports from src/radia (streamfunction_volume, stream_function)
# when run as a subprocess.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from radia.panels.calc_common import calc_main


def build_argparser() -> argparse.ArgumentParser:
    """Argparse for the Volume-3D stream-function coil designer.

    Importable WITHOUT side effects (no NGSolve / Radia import) so the panel can
    build its widgets from this surface alone."""
    p = argparse.ArgumentParser(
        description="3D volume stream-function (foliated Clebsch) solenoid "
                    "coil designer: target axial Bz -> equal-current wires.")
    p.add_argument("--vol", required=True,
                   help="Hollow cylindrical CONDUCTOR volume mesh (.vol).")
    p.add_argument("--target-bz", type=float, default=1.0e-3,
                   help="Target uniform axial field Bz in the bore [T].")
    p.add_argument("--foliation", choices=["radial"], default="radial",
                   help="Foliation mu.  'radial' = nested cylinders mu=r.")
    p.add_argument("--n-leaves", type=int, default=3,
                   help="Number of mu-leaf cylinders (radial shells).")
    p.add_argument("--fes-order", type=int, default=2,
                   help="H1 polynomial order for the lambda field.")
    p.add_argument("--aca-eps", type=float, default=1.0e-10,
                   help="ACA+ stopping tolerance for the response factorisation.")
    p.add_argument("--rings-per-span", type=int, default=14,
                   help="Equal-current rings across the lambda range "
                        "(wire density; larger = more, finer wires).")
    p.add_argument("--n-targets", type=int, default=9,
                   help="Number of axial field-target points along the bore.")
    p.add_argument("--nphi", type=int, default=49,
                   help="Azimuthal samples for contour extraction.")
    p.add_argument("--nz", type=int, default=41,
                   help="Axial samples for contour extraction.")
    p.add_argument("--n-threads", type=int, default=4,
                   help="NGSolve TaskManager thread count.")
    p.add_argument("--msh-output", default=None,
                   help="GMSH .msh (v4.1) wire overlay to write.  Panel sets it.")
    p.add_argument("--output", default=None,
                   help="Result JSON.  calc_main writes the return dict here.")
    return p


def run(args) -> dict:
    """Design the volume coil; return the standard result dict.

    Heavy imports happen HERE, never at module top.  TaskManager is the FIRST
    ngsolve symbol imported (UnboundLocalError trap)."""
    import time
    from ngsolve import TaskManager, SetNumThreads, Mesh

    from streamfunction_volume import design_volume_coil, write_wires_msh

    SetNumThreads(int(args.n_threads))

    t0 = time.perf_counter()
    with TaskManager():
        mesh = Mesh(args.vol)
        t_mesh_s = time.perf_counter() - t0

        t1 = time.perf_counter()
        res = design_volume_coil(
            mesh,
            target_bz=args.target_bz,
            n_leaves=args.n_leaves,
            fes_order=args.fes_order,
            aca_eps=args.aca_eps,
            rings_per_span=args.rings_per_span,
            n_targets=args.n_targets,
            nphi=args.nphi,
            nz=args.nz,
        )
        t_design_s = time.perf_counter() - t1

    wires = res.pop("_wires")
    res.pop("_currents")

    gmsh_file = ""
    if args.msh_output:
        n_node, n_line = write_wires_msh(wires, args.msh_output)
        gmsh_file = args.msh_output
        res["msh_nodes"] = int(n_node)
        res["msh_lines"] = int(n_line)

    res["t_mesh_s"] = round(t_mesh_s, 3)
    res["t_design_s"] = round(t_design_s, 3)
    res["t_total_s"] = round(time.perf_counter() - t0, 3)
    res["gmsh_file"] = gmsh_file
    return res


def main():
    parser = build_argparser()
    calc_main(run, parser)


if __name__ == "__main__":
    main()
