"""Persist the per-DOF |H_t| and Z_s of an ESIM per-panel run as
NGSolve `.sol` files alongside the per-DOF JSON.

The per-panel JSONs (sweep_data_dense/I*_f*_per_panel.json) carry the raw
arrays `esim_per_panel_H_t`, `esim_per_panel_Z_s_real`, and
`esim_per_panel_Z_s_imag`.  This script wraps them as NGSolve
GridFunctions on the workpiece mesh and writes:

    sweep_data_dense/<case>_Ht.sol    # H1 order 1, real
    sweep_data_dense/<case>_Zs.sol    # H1 order 1, complex

so the field is reloadable as a true NGSolve GridFunction --
feeds the NGSolve → GMSH (`vol2msh`/`GmshPostExport`) post pipeline
and is the persistence form the lab visualization policy expects.

Run once per case:

    python save_per_panel_sols.py I100_f50k    # 100 A / 50 kHz
    python save_per_panel_sols.py I500_f10k    # 500 A / 10 kHz (max gap)

Or omit the argument to write both committed cases.
"""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SWEEP = HERE / "sweep_data_dense"
VOL = (HERE.parent.parent / "src" / "radia" / "panels" / "samples"
       / "ih_bem_sample_p1.vol")


def save_sols_for_case(case: str) -> None:
    """case is the JSON stem, e.g. 'I100_f50k' or 'I500_f10k'."""
    from ngsolve import Mesh, H1, GridFunction, TaskManager, BND

    perp = SWEEP / f"{case}_per_panel.json"
    if not perp.exists() or not VOL.exists():
        sys.exit(f"ERROR: need {perp} and {VOL}.")
    d = json.load(open(perp))
    zr = np.asarray(d["esim_per_panel_Z_s_real"])
    zi = np.asarray(d["esim_per_panel_Z_s_imag"])
    Ht = np.asarray(d["esim_per_panel_H_t"])

    with TaskManager():
        mesh = Mesh(str(VOL))
        labs = list(mesh.GetBoundaries())
        sibc = {i for i, n in enumerate(labs) if n.lower() == "sibc"}
        bv = sorted({int(v.nr) for el in mesh.Elements(BND)
                     if int(el.index) in sibc for v in el.vertices})

        # |H_t| as a real H1 order-1 GridFunction on the full mesh
        # (zero outside the SIBC vertices).
        fes_r = H1(mesh, order=1)
        gf_Ht = GridFunction(fes_r)
        gf_Ht.vec[:] = 0.0
        for k, vnr in enumerate(bv):
            gf_Ht.vec[vnr] = float(Ht[k])
        ht_path = SWEEP / f"{case}_Ht.sol"
        gf_Ht.Save(str(ht_path))

        # Z_s as a complex H1 order-1 GridFunction.
        fes_c = H1(mesh, order=1, complex=True)
        gf_Zs = GridFunction(fes_c)
        gf_Zs.vec[:] = 0.0
        for k, vnr in enumerate(bv):
            gf_Zs.vec[vnr] = complex(float(zr[k]), float(zi[k]))
        zs_path = SWEEP / f"{case}_Zs.sol"
        gf_Zs.Save(str(zs_path))

    print(f"Saved {ht_path.name} ({ht_path.stat().st_size} B) and "
          f"{zs_path.name} ({zs_path.stat().st_size} B).")


def main():
    cases = sys.argv[1:] or ["I100_f50k", "I500_f10k"]
    for c in cases:
        save_sols_for_case(c)


if __name__ == "__main__":
    main()
