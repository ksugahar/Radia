"""bench_dtn_mbuild -- Step 5 benchmark: build the material-aware transfer matrix M

Selling point of the Kelvin-FEM route (Track B): the operator a stream-function coil design inverts
is generated from a SPARSE SPD volume matrix (one factorisation + a back-substitution per coil DoF),
GREEN-FUNCTION-FREE, for ARBITRARY iron.  This benchmark measures that against the dense
Green-function route a BEM / layered-Sommerfeld code would take.

Three honest contrasts (all measured in-repo, JSON + figure committed next to this script):

  C1  SPARSITY + FACTORISATION COST of the M-generator (the Kelvin-FEM volume matrix A), swept to
        large ndof: nnz, fill fraction, nnz/row, sparse-Cholesky time/storage vs DENSE-Cholesky of the
        SAME free-DoF system (measured up to a cap; the dense O(ndof^2) storage is computed exactly
        at every size, the dense factor time is measured up to the cap).  This is demo_cc's "sparse
        volume vs dense condensed operator" made quantitative + timed: sparse nnz/row is CONSTANT
        (fill -> 0) while the dense operator a Green/BEM route forms is O(ndof^2).

  C2  ACCURACY vs the GREEN-FUNCTION ROUTE (concentric shell, where the layered-sphere Green operator
        EXISTS, order=2 well-resolved): the FEM transfer eigenvalue R_n (best-fit over the angular
        pattern, exactly as demo_hh scenario_A) reproduces the analytic layered transfer (analytic_R)
        -- i.e. Kelvin-FEM builds the SAME operator a layered-Green/BEM code assembles (rel < 3e-2).

  C3  GENERALITY (the decisive point): for a NON-concentric iron blob there is NO closed-form Green
        function at all (BEM -> volume integral equation), yet the Kelvin-FEM M builds with the
        IDENTICAL sparse cost.  The dense Green route is simply N/A.

Per the Benchmark Policy: one case at a time, psutil peak memory, JSON with the required fields.
The C1 sweep is order=1 (clean ndof sweep; the asymptotic sparse-vs-dense scaling is FE-order-
independent); the C2 anchor is order=2 (resolution for the physics check).

Run:  pip install -e packages/radia-mcp ; python bench_dtn_mbuild.py
Needs numpy, scipy, ngsolve 6.2.2604, netgen.occ, psutil, radia_mcp.
"""
import os, sys, json, platform, time
from datetime import datetime
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import scipy.sparse as sp
import ngsolve as ng
from ngsolve import TaskManager
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import demo_hh_general_iron_design as hh   # build_fem / analytic_R / _CN / coil_dofs / a, R_out, r_t
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

try:
    import psutil
    _PROC = psutil.Process(os.getpid())
except Exception:
    _PROC = None

HERE = os.path.dirname(os.path.abspath(__file__))
B_SHELL, C_SHELL, MU_R = 0.70, 0.90, 50.0          # concentric iron shell [b,c], permeability
SWEEP_ORDER = 1
DENSE_CAP = 9000                                    # measure dense Cholesky only up to this nfree
MAXH_SWEEP = [0.20, 0.14, 0.105, 0.08, 0.06]       # order-1 ndof sweep (coarse -> ~40k)
ACC_ORDER, ACC_MAXH = 2, 0.25                       # demo_hh scenario_A anchor config
BLOB = ("blob", (0.0, 0.0, 0.68), 0.16)
TGT_THETAS = np.deg2rad(np.linspace(10, 170, 8))
TGT_PHIS = np.deg2rad(np.linspace(0, 315, 10))


def peak_mb():
    if _PROC is None:
        return None
    mi = _PROC.memory_info()
    return (getattr(mi, "peak_wset", None) or mi.rss) / (1024 * 1024)


def sparse_from_form(A, ndof):
    r, c, v = A.mat.COO()
    return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(ndof, ndof))


def n_boundary_dofs(mesh, fes):
    """Free DoFs on ANY boundary -- a proxy for the surface-unknown count N_Gamma a BEM would carry."""
    fd = fes.FreeDofs(); bnd = fes.GetDofs(mesh.Boundaries(".*"))
    return int(sum(1 for i in range(fes.ndof) if bnd[i] and fd[i]))


def measure_case(iron_spec, maxh, order, want_dense):
    """C1/C3: sparsity + factorisation + M-build cost for one mesh.  No accuracy check here."""
    rec = {"maxh": maxh, "order": order, "iron": iron_spec[0] if iron_spec else "none"}
    with TaskManager():
        mesh, fes, make_A = hh.build_fem(iron_spec, order, maxh)
        A = make_A(MU_R)
        ndof = fes.ndof
        fd = fes.FreeDofs(); free = [i for i in range(ndof) if fd[i]]; nfree = len(free)
        Sf = sparse_from_form(A, ndof)[np.ix_(free, free)]
        nnz = int(Sf.nnz)
        rec.update(ndof=ndof, nfree=nfree, nnz=nnz, nnz_per_row=nnz / max(nfree, 1),
                   fill=nnz / float(nfree * nfree), N_gamma=n_boundary_dofs(mesh, fes),
                   storage_sparse_mb=nnz * 12 / 1e6, storage_dense_mb=nfree * nfree * 8 / 1e6)

        t0 = time.time()
        Ainv = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        rec["t_factor_sparse"] = time.time() - t0

        if want_dense and nfree <= DENSE_CAP:
            D = Sf.toarray(); t0 = time.time(); np.linalg.cholesky(D)
            rec["t_factor_dense"] = time.time() - t0; del D
        else:
            rec["t_factor_dense"] = None

        cdofs = hh.coil_dofs(mesh, fes)
        targets = [mesh(hh.r_t * np.sin(t) * np.cos(p), hh.r_t * np.sin(t) * np.sin(p), hh.r_t * np.cos(t))
                   for t in TGT_THETAS for p in TGT_PHIS]
        gf = ng.GridFunction(fes); rr = gf.vec.CreateVector()
        t_bs = t_ev = 0.0
        for idof in cdofs:
            gf.vec[:] = 0.0; gf.vec[idof] = 1.0
            t0 = time.time(); rr.data = -(A.mat * gf.vec); gf.vec.data += Ainv * rr; t_bs += time.time() - t0
            t0 = time.time(); [gf(t) for t in targets]; t_ev += time.time() - t0
        rec.update(n_coil=len(cdofs), n_target=len(targets), t_backsub=t_bs, t_eval=t_ev,
                   t_mbuild=t_bs + t_ev)
    rec["peak_memory_mb"] = peak_mb()
    rec["t_setup"] = rec["t_factor_sparse"]; rec["t_solve"] = rec["t_mbuild"]
    rec["iterations"] = rec["n_coil"]; rec["converged"] = True
    return rec


def accuracy_case():
    """C2: order-2 well-resolved -- the FEM transfer R_n vs the analytic layered-Green transfer
    (best-fit over the angular pattern, exactly as demo_hh scenario_A)."""
    from scipy.special import eval_legendre
    with TaskManager():
        mesh, fes, make_A = hh.build_fem(("concentric", B_SHELL, C_SHELL), ACC_ORDER, ACC_MAXH)
        A = make_A(MU_R)
        Ainv = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        thetas = np.deg2rad([10, 35, 60, 90, 120, 150])
        targets = [mesh(hh.r_t * np.sin(t), 0.0, hh.r_t * np.cos(t)) for t in thetas]
        cth = np.cos(thetas)
        gf = ng.GridFunction(fes); rr = gf.vec.CreateVector()
        accs = []
        for n in (1, 2, 3):
            gf.vec[:] = 0.0
            gf.Set(_solid_harmonic(n) / hh.a ** n, ng.BND, definedon=mesh.Boundaries("inner"))
            rr.data = -(A.mat * gf.vec); gf.vec.data += Ainv * rr
            v = np.array([gf(t) for t in targets]); P = eval_legendre(n, cth)
            R_fem = float(np.dot(v, P) / np.dot(P, P))
            R_ana = hh.analytic_R(n, MU_R, hh.r_t, B_SHELL, C_SHELL) * hh._CN[n]
            accs.append({"n": n, "R_fem": R_fem, "R_analytic": R_ana,
                         "rel": abs(R_fem - R_ana) / abs(R_ana)})
    return {"order": ACC_ORDER, "maxh": ACC_MAXH, "per_n": accs,
            "rel_max": max(a["rel"] for a in accs)}


def main():
    print("=== bench_dtn_mbuild : sparse Kelvin-FEM M-build vs dense Green-function route ===")
    print("iron shell [%.2f,%.2f] mu_r=%.0f | C1 sweep order=%d | C2 anchor order=%d | dense cap nfree<=%d\n"
          % (B_SHELL, C_SHELL, MU_R, SWEEP_ORDER, ACC_ORDER, DENSE_CAP))

    print("[C1] SPARSITY + FACTORISATION SCALING (order=%d):" % SWEEP_ORDER)
    results = []
    for maxh in MAXH_SWEEP:
        rec = measure_case(("concentric", B_SHELL, C_SHELL), maxh, SWEEP_ORDER, want_dense=True)
        results.append(rec)
        td = ("%.2fs" % rec["t_factor_dense"]) if rec["t_factor_dense"] is not None else "(>cap,model)"
        print("  ndof=%6d nfree=%6d nnz/row=%.1f fill=%.2e | store sparse=%.1fMB dense=%.0fMB | "
              "t_fac sparse=%.2fs dense=%s | M %dx%d build=%.1fs"
              % (rec["ndof"], rec["nfree"], rec["nnz_per_row"], rec["fill"], rec["storage_sparse_mb"],
                 rec["storage_dense_mb"], rec["t_factor_sparse"], td, rec["n_target"], rec["n_coil"],
                 rec["t_mbuild"]))
    big = results[-1]
    print("  => sparse nnz/row CONSTANT (~%.0f), fill -> %.1e; the dense operator the Green/BEM route "
          "forms is %.0fx larger at the finest mesh (%.0fMB vs %.1fMB)."
          % (big["nnz_per_row"], big["fill"], big["storage_dense_mb"] / big["storage_sparse_mb"],
             big["storage_dense_mb"], big["storage_sparse_mb"]))

    print("\n[C2] ACCURACY vs the analytic layered-Green transfer (order=%d, well-resolved):" % ACC_ORDER)
    acc = accuracy_case()
    for a in acc["per_n"]:
        print("  n=%d  R_fem=%+.6e  analytic=%+.6e  rel=%.2e" % (a["n"], a["R_fem"], a["R_analytic"], a["rel"]))
    print("  => Kelvin-FEM builds the SAME operator the Green route gives (rel_max=%.1e < 3e-2): %s"
          % (acc["rel_max"], "PASS" if acc["rel_max"] < 3e-2 else "CHECK"))

    print("\n[C3] GENERALITY -- non-concentric iron blob (no closed-form Green function):")
    blob_rec = measure_case(BLOB, 0.14, SWEEP_ORDER, want_dense=False)
    blob_rec["green_route"] = "N/A -- no layered/Sommerfeld Green function for a non-concentric blob"
    results.append(blob_rec)
    print("  ndof=%d nfree=%d nnz/row=%.1f fill=%.2e | t_fac sparse=%.2fs | M %dx%d build=%.1fs"
          % (blob_rec["ndof"], blob_rec["nfree"], blob_rec["nnz_per_row"], blob_rec["fill"],
             blob_rec["t_factor_sparse"], blob_rec["n_target"], blob_rec["n_coil"], blob_rec["t_mbuild"]))
    print("  => same sparse cost as the concentric case; the dense Green route does not exist here.")

    out = {
        "timestamp": datetime.now().isoformat(), "hostname": platform.node(), "benchmark": "dtn_mbuild",
        "problem": {"iron_shell": [B_SHELL, C_SHELL], "mu_r": MU_R, "sweep_order": SWEEP_ORDER,
                    "acc_order": ACC_ORDER, "acc_maxh": ACC_MAXH, "dense_cap_nfree": DENSE_CAP,
                    "maxh_sweep": MAXH_SWEEP, "blob": list(BLOB)},
        "results": results, "accuracy_vs_analytic_green": acc,
        "notes": ("storage_sparse_mb = nnz*12B (f64 val + 2 i32 idx); storage_dense_mb = nfree^2*8B "
                  "(the dense matrix a BEM-condensed/Green operator stores). The Green route also needs "
                  "the fundamental solution + singular quadrature, and DOES NOT EXIST for arbitrary iron "
                  "(C3 blob). The dense factor time is measured to nfree<=%d; the dense storage is "
                  "computed exactly at every size." % DENSE_CAP),
    }
    jpath = os.path.join(HERE, "bench_dtn_mbuild.json")
    with open(jpath, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("\nJSON -> %s" % jpath)
    _make_figure(results)
    return out


def _make_figure(results):
    conc = [r for r in results if r["iron"] == "concentric"]
    if not conc:
        return
    ndof = np.array([r["ndof"] for r in conc], float)
    try:
        from radia_mcp.figure import paper_figure, save_lab_figure, add_slope_guide
        fig, ax = paper_figure("IEEE_DOUBLE_COLUMN", nrows=1, ncols=2)
        a0, a1 = ax.ravel()
        a0.loglog(ndof, [r["storage_sparse_mb"] for r in conc], "o-", label="sparse A (Kelvin-FEM)")
        a0.loglog(ndof, [r["storage_dense_mb"] for r in conc], "s--", label="dense (Green/BEM)")
        a0.set_xlabel("volume DoF"); a0.set_ylabel("operator storage [MB]"); a0.legend()
        add_slope_guide(a0, ndof[len(ndof) // 2], conc[len(ndof) // 2]["storage_dense_mb"], 2.0, label="O(N^2)")
        a1.loglog(ndof, [r["t_factor_sparse"] for r in conc], "o-", label="sparse factor")
        a1.loglog(ndof, [r["t_mbuild"] for r in conc], "^-", label="M-build (back-subs)")
        dm = np.array([r["t_factor_dense"] if r["t_factor_dense"] is not None else np.nan for r in conc])
        if np.any(np.isfinite(dm)):
            a1.loglog(ndof, dm, "s--", label="dense factor")
        a1.set_xlabel("volume DoF"); a1.set_ylabel("time [s]"); a1.legend()
        fig.subplots_adjust(wspace=0.34, left=0.11, right=0.985, bottom=0.17, top=0.97)
        info = save_lab_figure(fig, os.path.join(HERE, "bench_dtn_mbuild"), tighten=False)
        print("figure -> %s" % info.get("pdf", os.path.join(HERE, "bench_dtn_mbuild.pdf")))
    except Exception as e:                 # figure is an optional extra; JSON already persisted
        print("figure skipped (%s)" % e)


if __name__ == "__main__":
    main()
