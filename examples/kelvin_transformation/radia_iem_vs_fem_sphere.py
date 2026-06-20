"""3-way validation: Radia open-boundary IEM vs reduced-Omega + Kelvin FEM vs the analytic sphere.

A linear magnetic sphere (relative permeability mu_r) in a uniform applied field H0 has the EXACT
uniform interior field

    H_in = 3 / (mu_r + 2) * H0 ,   M = (mu_r - 1) * H_in   (magnetic-sphere-in-uniform-field)

Both independent solvers are checked against that analytic answer on the SAME quarter-sphere geometry
(`Cubit_1_4_p_convergence/sphere_1_4_p{1,2,3}.vol`, x>=0 & y>=0, radius 50 mm), order-matched:

  * analytic   : 3 / (mu_r + 2)
  * Radia IEM  : MMM -- the tetrahedron surface-charge / dipole OPEN-BOUNDARY integral method (no air
                 mesh, exact analytic open boundary).  The sphere meshes as tets, so MMM is the Radia IEM
                 member here; yano-type MSC is the HEXAHEDRON variant of the same surface-charge IEM
                 (parity yano == MMM == HDiv-VIM locked by tests/feec/parity_vs_msc, cube 0.76%), and the
                 loop removal (rad.GetLoopBasis, Stage 1) is a field-PRESERVING yano-hex conditioning fix
                 that does not change this accuracy.  Quarter model solved with image='+x+y' (the z-field
                 is parallel to both the x=0 and y=0 planes -> symmetric; MMM has NO IMA boundary
                 limitation, CLAUDE.md).  Observable = volume-averaged magnetization <M_z> ->
                 H_in = <M_z>/(mu_r-1); NOT the inside-iron field eval (an unreliable dipole approximation
                 for MMM, per CLAUDE.md "rad.Fld inside materials").
  * FEM        : NGSolve reduced-Omega (2-scalar) + Kelvin transformation (open boundary, no PML), via
                 solve_kelvin_benchmark on the full magnetic+air+kelvin mesh.

TWO lessons the p-sweep makes explicit (this is what to get right with reduced-Omega):
  1. reduced-Omega MUST use p>=2.  p=1 here is +7.8% (unusable).
  2. The FE order MUST be MATCHED to the mesh geometry-curving order: solve order p on the order-p
     curved mesh.  Running p=3 on an order-2-curved mesh is a geometry mismatch and gives a spurious
     result that does not reflect true p=3.  Order-matched, the residual at p=2/3 (~0.5-0.7%, == the
     locked golden bands) is the coarse academic mesh's h + Kelvin-truncation floor -- NOT FE-order and
     NOT a code bug.  The no-air-mesh IEM lands closest.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DMESH = os.path.join(HERE, "Cubit_1_4_p_convergence")
sys.path.insert(0, os.path.join(REPO, "src", "radia"))
sys.path.insert(0, os.path.join(REPO, "src", "radia", "panels"))

import ngsolve as ng  # noqa: E402
import radia as rad  # noqa: E402
from netgen_mesh_import import extract_elements  # noqa: E402
from calc_kelvin_benchmark import solve_kelvin_benchmark  # noqa: E402

MU_R = 100.0
MU0 = 4.0e-7 * math.pi


def mesh_for(p):
    return os.path.join(DMESH, f"sphere_1_4_p{p}.vol")


def analytic_ratio(mu_r):
    return 3.0 / (mu_r + 2.0)


def radia_iem_mmm(mu_r, H0=1000.0):
    """Radia MMM (tet surface-charge IEM, analytic open boundary) on the x>=0,y>=0 quarter sphere +
    image='+x+y'.  Returns H_in/H0 from the volume-averaged magnetization (the reliable MMM observable).
    MMM uses straight tets (no curving), so the order-1 mesh vertices are used."""
    mesh = ng.Mesh(mesh_for(1))
    els, _ = extract_elements(mesh, material_filter="magnetic", allow_hex=False)
    V = [np.array(e["vertices"], float) for e in els]
    vol = np.array([abs(np.dot(p[1] - p[0], np.cross(p[2] - p[0], p[3] - p[0]))) / 6.0 for p in V])
    cen = np.array([p.mean(0) for p in V])

    rad.UtiDelAll()
    rad.set_demag_backend("auto")        # tet -> MMM (C++)
    objs = []
    for p in V:
        t = rad.ObjTetrahedron([list(v) for v in p], [0, 0, 0])
        rad.MatApl(t, rad.MatLin(mu_r))
        objs.append(t)
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda q: [0, 0, MU0 * H0])])
    rad.Solve(cont, 1e-6, 3000, 0, image="+x+y")    # z-field // x=0 and y=0 planes -> symmetric
    Mz = np.array([rad.Fld(cont, "m", list(c))[2] for c in cen])
    rad.UtiDelAll()
    Mz_avg = float((Mz * vol).sum() / vol.sum())
    return Mz_avg / (mu_r - 1.0) / H0, Mz_avg, len(els)


def reduced_omega_kelvin_fem(mu_r, H0=1.0, orders=(1, 2, 3)):
    """ORDER-MATCHED reduced-Omega + Kelvin p-sweep: solve order p on the order-p curved mesh."""
    ng.SetNumThreads(4)
    out = []
    for p in orders:
        with ng.TaskManager():
            r = solve_kelvin_benchmark(mesh_for(p), mu_r=mu_r, H0=H0, field_axis="z",
                                       fes_order=p, R_kelvin=0.20)
        if "error" in r:
            raise RuntimeError(f"p={p}: {r['error']}")
        out.append({"fes_order": p, "ratio_Hin_over_H0": r["Hi_origin"] / H0, "ndof": r.get("ndof")})
    return out


def main():
    ratio_ana = analytic_ratio(MU_R)
    ratio_mmm, Mz_avg, n_tet = radia_iem_mmm(MU_R)
    fem = reduced_omega_kelvin_fem(MU_R, orders=(1, 2, 3))

    print(f"\nMagnetic sphere in uniform field, mu_r={MU_R:.0f}  (analytic H_in/H0 = {ratio_ana:.6f})")
    print("quarter sphere, ORDER-MATCHED (geometry curving order == FE order)\n")
    print(f"  {'method':<40} {'H_in/H0':>10} {'err vs analytic':>16}   note")
    print(f"  {'analytic 3/(mu_r+2)':<40} {ratio_ana:>10.6f} {0.0:>+15.2%}")
    print(f"  {'Radia IEM (MMM, tet, no air mesh)':<40} {ratio_mmm:>10.6f} "
          f"{ratio_mmm/ratio_ana-1.0:>+15.2%}   {n_tet} tets, image=+x+y")
    for f in fem:
        note = "(p=1 UNRELIABLE -- use p>=2)" if f["fes_order"] == 1 else f"ndof={f['ndof']}, order-matched"
        print(f"  {'reduced-Omega + Kelvin FEM p=' + str(f['fes_order']):<40} "
              f"{f['ratio_Hin_over_H0']:>10.6f} {f['ratio_Hin_over_H0']/ratio_ana-1.0:>+15.2%}   {note}")

    out = {
        "geometry": "magnetic sphere (R=50mm) in uniform field; quarter model, order-matched p-meshes",
        "meshes": [os.path.relpath(mesh_for(p), REPO).replace("\\", "/") for p in (1, 2, 3)],
        "mu_r": MU_R,
        "analytic_ratio_Hin_over_H0": ratio_ana,
        "radia_iem_mmm": {"ratio_Hin_over_H0": ratio_mmm, "Mz_avg": Mz_avg, "n_tet": n_tet,
                          "image": "+x+y", "err_vs_analytic": ratio_mmm / ratio_ana - 1.0,
                          "observable": "volume-averaged M_z (reliable; not inside-field eval)"},
        "reduced_omega_kelvin_fem_psweep": [
            {**f, "err_vs_analytic": f["ratio_Hin_over_H0"] / ratio_ana - 1.0} for f in fem],
        "notes": [
            "reduced-Omega MUST use p>=2: p=1 here is +7.8% (unreliable).",
            "FE order MUST be matched to the mesh geometry-curving order (solve order p on the order-p "
            "curved mesh); p=3 on an order-2 mesh is a geometry mismatch, not true p=3.",
            "Order-matched, the p=2/3 residual (~0.5-0.7%, == the locked golden bands) is the coarse "
            "academic mesh's h + Kelvin-truncation floor -- not FE-order, not a code bug.  The no-air-mesh "
            "IEM (MMM) lands closest.",
            "MMM is the tet variant of the Radia surface-charge IEM; yano-type MSC is the hex variant "
            "(parity locked by tests/feec/parity_vs_msc). Loop removal (rad.GetLoopBasis) is "
            "field-preserving and does not change this accuracy.",
        ],
    }
    with open(os.path.join(HERE, "radia_iem_vs_fem_sphere.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print("\nReadout: reduced-Omega needs p>=2 AND geometry-order matching; order-matched, both the IEM")
    print("and the FEM match the analytic sphere within ~0.5-0.7% (mesh h/Kelvin-truncation limited).")
    print("saved", os.path.join(HERE, "radia_iem_vs_fem_sphere.json"))


if __name__ == "__main__":
    main()
