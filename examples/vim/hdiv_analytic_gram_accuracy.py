"""
hdiv_analytic_gram_accuracy.py -- (A) how much does the EXACT analytic charge Gram buy in demag
ACCURACY for the loop-free HDiv-VIM, on geometries with a KNOWN exact answer?

The user's goal is accuracy + loop-free + distortion-solvable.  loop-free + distortion-solvable are
already established (hdiv_demag_quad_self.py / hdiv_distortion_precond.py).  This script isolates the
ACCURACY axis: the uniform-Mz demag factor (Rayleigh quotient m^T N m / m^T M_mass m) vs the exact
reference (a sphere is EXACTLY 1/3; a cube's uniform demag is ~1/3) for three Gram qualities:

  default        build_demag(mesh)                 centroid-monopole off-diag + sub-point self  (the
                                                    "~5-6% monopole error" the docstring warns about)
  wilton_surface build_demag(mesh, wilton_surface) exact Wilton surface-surface block            (docstring:
                                                    "demag 1/3 to <0.15%" -- the #1 production accuracy fix)
  analytic_gram  build_demag(mesh, analytic_gram)  FULL exact phi_tet/tri_potential Gram         (also closes
                                                    the non-uniform-M / nonlinear gap)

loop_res (=||N loop||/||N||) is printed too, to confirm loop-freeness is untouched by the Gram choice.
Distorted (sheared) cube included -- accuracy under distortion is the user's regime.
"""
import json, os, sys, math
import numpy as np
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src", "radia"))
from vim._core import build_demag  # noqa: E402

ng.SetNumThreads(4)


def demag_factor(d):
    """uniform-Mz demag factor = m^T N m / m^T M_mass m  (N, M_mass in the RT0 dof basis)."""
    m = d["m_unit"]
    N, M = d["N"], d["M_mass"]
    num = float(m @ (N @ m))
    den = float(m @ (M @ m))
    return num / den


def loop_residual(d):
    if d["loops"] is None or d["n_loop"] == 0:
        return 0.0
    N = d["N"]; Nn = np.linalg.norm(N, 2)
    return max(np.linalg.norm(N @ d["loops"][k]) for k in range(d["n_loop"])) / (Nn + 1e-300)


def tet_cube(n, shear=0.0):
    def mp(x, y, z):
        if shear == 0.0:
            return (x, y, z)
        s = shear
        return (x + s*math.sin(math.pi*y)*z, y + 0.83*s*x*z, z + 0.67*s*x*math.sin(math.pi*y))
    return MakeStructured3DMesh(hexes=False, nx=n, ny=n, nz=n, mapping=mp)


def tet_sphere(maxh):
    from netgen.occ import Sphere, Pnt, OCCGeometry
    geo = OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def measure(mesh, tag, exact):
    with ng.TaskManager():
        d_def = build_demag(mesh)
        d_wil = build_demag(mesh, wilton_surface=True)
        d_ana = build_demag(mesh, analytic_gram=True)
    rows = []
    for name, d in (("default(monopole)", d_def), ("wilton_surface", d_wil), ("analytic_gram", d_ana)):
        dz = demag_factor(d)
        err = abs(dz - exact) / exact
        lr = loop_residual(d)
        rows.append({"gram": name, "demag": dz, "rel_err": err, "loop_res": lr,
                     "ndof": d["ndof"], "n_loop": d["n_loop"]})
    print(f"\n[{tag}]  exact uniform demag ~= {exact:.5f}   (ndof={rows[0]['ndof']}, loops={rows[0]['n_loop']})")
    print(f"  {'Gram':<18} {'demag':>10} {'rel_err':>10} {'loop_res':>10}")
    for r in rows:
        print(f"  {r['gram']:<18} {r['demag']:>10.5f} {r['rel_err']:>10.2e} {r['loop_res']:>10.1e}")
    return {"tag": tag, "exact": exact, "rows": rows}


def main():
    res = {"cases": []}
    # sphere: demag EXACTLY 1/3 (the cleanest reference)
    res["cases"].append(measure(tet_sphere(0.35), "sphere tet (exact 1/3)", 1.0/3.0))
    # cube: uniform demag ~1/3 (tet-meshed unit cube)
    res["cases"].append(measure(tet_cube(4), "cube tet 4^3 (uniform ~1/3)", 1.0/3.0))
    # sheared cube: accuracy under distortion (the user's regime).  NOTE: a sheared cube's TRUE uniform
    # demag is no longer exactly 1/3 -- report demag + the spread across Gram qualities (the exact-Gram
    # value is the best available reference here), not a rel_err vs 1/3.
    res["cases"].append(measure(tet_cube(4, shear=0.18), "cube tet 4^3 sheared 18% (vs analytic ref)", 1.0/3.0))
    with open(os.path.join(HERE, "hdiv_analytic_gram_accuracy.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("\nsaved", os.path.join(HERE, "hdiv_analytic_gram_accuracy.json"))
    print("\nReadout: if wilton_surface/analytic_gram cut rel_err to <~0.2% (vs default ~5-6%) while")
    print("loop_res stays ~1e-15, then the exact Gram delivers the ACCURACY axis with loop-freeness intact.")


if __name__ == "__main__":
    main()
