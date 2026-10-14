"""Independent Python entry-point oracle for the MATLAB sparsesolv gateway."""
import json
import hashlib
from pathlib import Path
import sys

import ngsolve as ng
import numpy as np
from netgen.occ import Box, OCCGeometry, Pnt
import radia.sparsesolv_ngsolve as ss


def generate(directory):
    # This is a small parity fixture, not a scaling benchmark. Never inherit
    # the host-wide default thread count for mesh generation and assembly.
    ng.SetNumThreads(2)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    mesh_path = directory / "mesh.vol"
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=0.6))
    mesh.ngmesh.Save(str(mesh_path))
    space = ng.HCurl(mesh, order=1, nograds=True)
    u, v = space.TnT()
    a = ng.BilinearForm(space)
    a += u*v*ng.dx
    with ng.TaskManager():
        a.Assemble()
    grad, _ = space.CreateGradient()
    coords = np.array([v.point for v in mesh.vertices])
    kwargs = dict(grad_mat=grad, freedofs=space.FreeDofs(),
                  coord_x=coords[:, 0].tolist(), coord_y=coords[:, 1].tolist(),
                  coord_z=coords[:, 2].tolist())
    ams = ss.HypreBasedAMSPreconditioner(a.mat, **kwargs)
    cs = ng.HCurl(mesh, order=1, nograds=True, complex=True)
    u, v = cs.TnT()
    ac = ng.BilinearForm(cs)
    ac += (1+0.5j)*u*v*ng.dx
    with ng.TaskManager():
        ac.Assemble()
    amsc = ss.ComplexHypreBasedAMSPreconditioner(a.mat, **kwargs)
    output = dict(mesh=str(mesh_path), python=sys.version, python_executable=sys.executable,
                  ngsolve=ng.__version__, ngsolve_threads=2,
                  mesh_sha256=hashlib.sha256(mesh_path.read_bytes()).hexdigest(),
                  sparsesolv_sha256=hashlib.sha256(Path(ss.__file__).read_bytes()).hexdigest(),
                  sparsesolv_binary=ss.__file__, cases=[])
    for name, mat, pre in [("ams", a.mat, ams), ("complex_ams", ac.mat, amsc),
                           ("ic", a.mat, ss.ICPreconditioner(a.mat)),
                           ("complex_ic", ac.mat, ss.ICPreconditioner(ac.mat))]:
        rhs = mat.CreateColVector()
        rhs.FV().NumPy()[:] = np.sin(np.arange(space.ndof)+0.2)
        if mat.is_complex:
            rhs.FV().NumPy()[:] *= 1+0.3j
        applied = rhs.CreateVector()
        with ng.TaskManager():
            applied.data = pre*rhs
        solver = ss.COCRSolver(mat, pre, tol=1e-11, maxiter=500)
        solution = rhs.CreateVector()
        with ng.TaskManager():
            solution.data = solver*rhs
        def split(vec):
            values = vec.FV().NumPy()
            return dict(real=values.real.tolist(), imag=values.imag.tolist())
        output["cases"].append(dict(name=name, rhs=split(rhs), applied=split(applied),
                                     solution=split(solution), iterations=solver.iterations))
    output["iccg"] = iccg_reference(mesh)
    path = directory / "reference.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    return path


# Option sets shared with test_sparsesolv_mex.m:testICCGPythonParity (MATLAB names)
ICCG_CASES = [
    ("default", {}, {}),
    ("fixed_unscaled", dict(Shift=1.2, AutoShift=False, DiagonalScaling=False),
     dict(shift=1.2, auto_shift=False, diagonal_scaling=False)),
    ("iteration_limit", dict(MaxIterations=4), dict(maxiter=4)),
    ("strict_stagnation", dict(DivergenceThreshold=1.0, DivergenceCount=0),
     dict(divergence_threshold=1.0, divergence_count=0)),
    ("abmc", dict(UseABMC=True), dict(use_abmc=True)),
]


def iccg_reference(mesh):
    """ICCG solves on the H1 order-2 Dirichlet Laplace system, real and complex."""
    out = []
    for complex_ in (False, True):
        space = ng.H1(mesh, order=2, dirichlet=".*", complex=complex_)
        u, v = space.TnT()
        weight = (1 + 0.5j) if complex_ else 1.0
        a = ng.BilinearForm(space)
        a += weight * ng.grad(u) * ng.grad(v) * ng.dx
        with ng.TaskManager():
            a.Assemble()
        free = np.array(list(space.FreeDofs()), dtype=bool)
        rhs = a.mat.CreateColVector()
        values = np.cos(0.7 * np.arange(space.ndof) + 0.1)
        if complex_:
            values = values * (1 - 0.4j)
        rhs.FV().NumPy()[:] = values * free
        for name, _, kwargs in ICCG_CASES:
            solver = ss.SparseSolvSolver(a.mat, method="ICCG", freedofs=space.FreeDofs(),
                                         save_residual_history=True, **kwargs)
            sol = rhs.CreateVector()
            sol[:] = 0.0
            with ng.TaskManager():
                result = solver.Solve(rhs, sol)
            x = sol.FV().NumPy()
            out.append(dict(name=name, complex=complex_, ndof=space.ndof,
                            rhs=dict(real=rhs.FV().NumPy().real.tolist(),
                                     imag=rhs.FV().NumPy().imag.tolist()),
                            solution=dict(real=x.real.tolist(), imag=x.imag.tolist()),
                            converged=bool(result.converged), iterations=result.iterations,
                            best_iteration=result.best_iteration,
                            final_residual=result.final_residual,
                            true_residual=result.true_residual,
                            actual_shift=result.actual_shift,
                            residual_history=list(result.residual_history)))
    return out


if __name__ == "__main__":
    print(generate(sys.argv[1]))
