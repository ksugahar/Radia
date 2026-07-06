"""
Analytic function evaluation on mesh: .vol + expression -> .msh v4.1

Called as subprocess from Cubit Solve -> Analytic Function:
    python calc_analytic.py --vol model.vol --expr "sin(x)*cos(y)" --name "f"

Evaluates a mathematical expression as an NGSolve CoefficientFunction
on the mesh, exports to GMSH .msh v4.1 for visualization.

IMPORTANT: NGSolve must be imported BEFORE cubit.
Outputs JSON to stdout (last line).
"""

import argparse
import json
import os
import sys

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import setup_paths, progress, calc_main


def _log(msg):
    progress("ANALYTIC", msg)


def evaluate_analytic(vol_path, expression, field_name="f",
                      fes_order=1, vector=False):
    """Evaluate expression on mesh and export to .msh v4.1."""
    setup_paths()

    from ngsolve import Mesh, x, y, z, CF, sqrt, sin, cos, tan, exp, log, TaskManager
    from ngsolve import atan, asin, acos, pi
    import numpy as np

    _log(f"Loading: {vol_path}")
    mesh = Mesh(vol_path)
    _log(f"Mesh: {mesh.ne} elements, {mesh.nv} vertices")

    # Build CoefficientFunction from expression
    # Make math functions available in eval namespace
    ns = {
        "x": x, "y": y, "z": z,
        "CF": CF, "sqrt": sqrt,
        "sin": sin, "cos": cos, "tan": tan,
        "exp": exp, "log": log,
        "atan": atan, "asin": asin, "acos": acos,
        "pi": pi, "abs": abs,
    }
    try:
        cf = eval(expression, {"__builtins__": {}}, ns)
    except Exception as e:
        return {"error": f"Invalid expression: {e}"}

    _log(f"Expression: {expression}")

    # Evaluate on GridFunction for visualization
    from ngsolve import H1, GridFunction

    if vector:
        fes = H1(mesh, order=fes_order, dim=3)
    else:
        fes = H1(mesh, order=fes_order)
    gf = GridFunction(fes, name=field_name)
    with TaskManager():
        gf.Set(cf)

        base = os.path.splitext(vol_path)[0]
        sol_path = base + "_analytic.sol"
        gf.Save(sol_path)
        _log(f"Written: {sol_path}")

        from gmsh_post_export import vol2msh
        msh_path = base + "_analytic.msh"
        vol2msh(
            msh_path, vol_path,
            [{
                "name": field_name,
                "sol": sol_path,
                "fes": "H1",
                "fes_order": fes_order,
                "fes_dim": 3 if vector else 1,
                "ncomp": 3 if vector else 1,
            }],
        )
        _log(f"Written: {msh_path}")

        return {
            "vol_path": vol_path,
            "sol_path": sol_path,
            # NOTE: key MUST be "msh_file" (not "msh_path") so notebook
            # workbenches and GMSH helpers recognize this output.
            "msh_file": msh_path,
            "field_name": field_name,
            "expression": expression,
            "n_elements": mesh.ne,
            "n_vertices": mesh.nv,
        }


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate analytic function on mesh")
    parser.add_argument("--vol", required=True,
                        help="Netgen .vol mesh file")
    parser.add_argument("--expr", required=True,
                        help="Expression (e.g. 'sin(x)*cos(y)')")
    parser.add_argument("--name", default="f",
                        help="Field name for GMSH display")
    parser.add_argument("--fes-order", type=int, default=1,
                        help="FES interpolation order")
    parser.add_argument("--vector", action="store_true",
                        help="Treat expression as 3-component vector")

    def run(args):
        return evaluate_analytic(
            args.vol, args.expr, args.name,
            args.fes_order, args.vector)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
