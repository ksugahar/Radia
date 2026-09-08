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
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    mesh_path = directory / "mesh.vol"
    mesh = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=0.6))
    mesh.ngmesh.Save(str(mesh_path))
    space = ng.HCurl(mesh, order=1, nograds=True)
    u, v = space.TnT()
    a = ng.BilinearForm(space)
    a += u*v*ng.dx
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
    ac.Assemble()
    amsc = ss.ComplexHypreBasedAMSPreconditioner(a.mat, **kwargs)
    output = dict(mesh=str(mesh_path), python=sys.version, ngsolve=ng.__version__,
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
        applied.data = pre*rhs
        solver = ss.COCRSolver(mat, pre, tol=1e-11, maxiter=500)
        solution = rhs.CreateVector()
        solution.data = solver*rhs
        def split(vec):
            values = vec.FV().NumPy()
            return dict(real=values.real.tolist(), imag=values.imag.tolist())
        output["cases"].append(dict(name=name, rhs=split(rhs), applied=split(applied),
                                     solution=split(solution), iterations=solver.iterations))
    path = directory / "reference.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    return path


if __name__ == "__main__":
    print(generate(sys.argv[1]))
