"""Proposed surrogate tests for the historical linear volumetric F1 slab.

Guards: Breakage of the complex AMS and full-HCurl BDDC+sparsesolv_ams solve
routes that changes computed Joule power or loses finite convergent output.
Zero-output mutations of each production-native route failed the existing
true-residual contract; restoration solved normally. This does not lock the
internal shape or efficiency of a numerically equivalent preconditioner.
Assumes: surrogate test of adoption record induction_heating/20261024_ams_linear_slab (convergence level 1).
m2p1 is the lowest-order route; m2p2 is the coarsest h level of the BDDC route.
Constant real mu/sigma, exact tangent trace, peak phasors exp(+i omega t).
Historical values are from radia5.3.3/source8a988152/native4b09e365/run04.
PROPOSAL ONLY: rtol1e-9/atol0 awaits USER approval; this candidate is not cleared
for merge or numerical qualification. No analytic accuracy gate belongs here.
Runtime is metadata, not an assertion. No nonlinear BH, GPU, full-coil or thermal claim.
"""
import importlib.util
from pathlib import Path
import pytest
import cmath
import math
import time


def _compute(p):
    import numpy as np
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    import radia.sparsesolv_ngsolve as ssn
    start = time.perf_counter()
    mu, sigma, delta, H0 = 4*math.pi*1e-7, 5e6, 1e-3, 1000.0
    omega = 2/(mu*sigma*delta**2)
    a, width, k = 2*delta, delta/4, (1+1j)/delta
    Ay = mu*H0*(ng.exp(k*ng.x)-ng.exp(-k*ng.x))/(2*k*cmath.cosh(k*a))
    boundary_A = ng.CoefficientFunction((0, Ay, 0))
    mesh = MakeStructured3DMesh(hexes=False,nx=32,ny=2,nz=2,
        secondorder=False,mapping=lambda x,y,z:(4*delta*x-a,width*(y-.5),width*(z-.5)))
    if p not in (1,2): raise ValueError("admitted p1/p2 only")
    space = (ng.HCurl(mesh,order=1,nograds=True,complex=True,dirichlet=".*") if p==1
             else ng.HCurl(mesh,order=2,complex=True,dirichlet=".*"))
    u,v = space.TnT()
    form = ng.BilinearForm(space,symmetric=True)
    form += ((1/mu)*ng.InnerProduct(ng.curl(u),ng.curl(v))
             +1j*omega*sigma*ng.InnerProduct(u,v))*ng.dx(bonus_intorder=2*p+2)
    before_stats = dict(ssn.AMSCoarseStats()) if p==2 else {}
    if p == 2:
        pre = ng.Preconditioner(form,"bddc",coarsetype="sparsesolv_ams",coarseflags={"cycles":4})
    with ng.TaskManager():
        form.Assemble()
    if p == 1:
        realspace = ng.HCurl(mesh,order=1,nograds=True,dirichlet=".*")
        ur,vr = realspace.TnT()
        aux = ng.BilinearForm(realspace,symmetric=True)
        aux += ((1/mu)*ng.InnerProduct(ng.curl(ur),ng.curl(vr))
                +omega*sigma*ng.InnerProduct(ur,vr))*ng.dx
        with ng.TaskManager(): aux.Assemble()
        grad,_ = realspace.CreateGradient()
        coords = [tuple(vertex.point) for vertex in mesh.vertices]
        pre = ssn.ComplexHypreBasedAMSPreconditioner(aux.mat,grad,
            freedofs=realspace.FreeDofs(),coord_x=[x[0] for x in coords],
            coord_y=[x[1] for x in coords],coord_z=[x[2] for x in coords],
            ndof_complex=space.ndof,cycle_type=1,print_level=0)
    g = ng.GridFunction(space)
    g.Set(boundary_A,definedon=mesh.Boundaries(".*"))
    lift = g.vec.CreateVector(); lift.data = g.vec
    rhs = g.vec.CreateVector(); rhs.data = -form.mat*lift
    free = np.asarray(list(space.FreeDofs()),dtype=bool)
    solver = ssn.GMRESSolver(form.mat,pre,freedofs=space.FreeDofs(),
                             tol=1e-12,maxiter=500,restart=50)
    with ng.TaskManager():
        g.vec.data += solver*rhs
    res = rhs.CreateVector(); res.data = -form.mat*g.vec
    relative_residual = float(np.linalg.norm(res.FV().NumPy()[free])/np.linalg.norm(rhs.FV().NumPy()[free]))
    if not math.isfinite(relative_residual) or relative_residual > 1e-6:
        raise RuntimeError(f"F1 solve residual {relative_residual}; required <= 1e-6")
    if not np.array_equal(g.vec.FV().NumPy()[~free],lift.FV().NumPy()[~free]):
        raise RuntimeError("F1 Dirichlet lift changed")
    J = -1j*omega*sigma*g
    power = float(ng.Integrate(ng.InnerProduct(J,J).real,mesh,order=2*p+12))/(2*sigma)
    if not math.isfinite(power):
        raise RuntimeError(f"Nonfinite computed Joule power: {power}")
    stats = dict(ssn.AMSCoarseStats()) if p==2 else {}
    if p==2 and not (stats["builds"] > before_stats["builds"] and stats["applies"] > before_stats["applies"]):
        raise RuntimeError("BDDC coarse native route not observed")
    return {"p":p,"m":2,"native_coarse_stats":stats,"ndof":space.ndof,"power_W":power,"runtime_s":time.perf_counter()-start,
            "relative_residual":relative_residual,"iterations":solver.iterations}


@pytest.mark.parametrize("p", [1, 2])
def test_f1_historical_joule_power(p):
    value = _compute(p)
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("f1_regression_helper", root/"tools/regression_lock.py")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    baseline = helper.load_baseline(Path(__file__).parent/"fixtures"/f"ih_f1_m2_p{p}_5_3_3.json")
    conditions = {"a_over_delta": 2, "m": 2, "p": p, "phasor": "exp(+i omega t), peak",
                  "sigma_S_per_m": 5e6, "mu_H_per_m": 4*math.pi*1e-7,
                  "delta_m": 1e-3, "H0_A_per_m": 1000,
                  "route": "Complex AMS" if p==1 else "full HCurl BDDC+sparsesolv_ams",
                  "gmres_tol": 1e-12, "restart": 50, "maxiter": 500,
                  "runtime_metadata": "measured coarse case on source6d9e8fd30; historical run04 has no per-level timings"}
    helper.assert_regression({"power_W":value["power_W"]}, baseline, conditions=conditions)
