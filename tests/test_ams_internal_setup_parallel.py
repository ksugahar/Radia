"""Exercise native hierarchy setup in a disposable process."""
import subprocess
import sys


def test_internal_setup_workers_and_update():
    import radia.sparsesolv_ngsolve as native
    script = r'''
import importlib.util, sys
import numpy as np
from ngsolve import Mesh, HCurl, BilinearForm, curl, dx, TaskManager, SetNumThreads
from netgen.csg import unit_cube
spec = importlib.util.spec_from_file_location('sparsesolv_ngsolve', sys.argv[1])
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)
# Bound the initial assembly too: a fresh process otherwise inherits the host's
# full core count and reserves a multi-gigabyte local heap before the loop.
SetNumThreads(1)
mesh = Mesh(unit_cube.GenerateMesh(maxh=0.12))
fes = HCurl(mesh, order=1, nograds=True)
u,v = fes.TnT()
a = BilinearForm(fes)
a += (curl(u)*curl(v)+u*v)*dx
a.Assemble()
grad,_ = fes.CreateGradient()
xyz = np.asarray(mesh.ngmesh.Coordinates())
outputs = []
for threads in (1,4):
    SetNumThreads(threads)
    pre = native.HypreBasedAMSPreconditioner(a.mat,grad,
        coord_x=xyz[:,0].tolist(),coord_y=xyz[:,1].tolist(),coord_z=xyz[:,2].tolist())
    assert (pre.setup_workers == 1 if threads == 1 else pre.setup_workers > 1)
    rhs, out = pre.CreateColVector(), pre.CreateRowVector()
    rhs.FV().NumPy()[:] = np.sin(np.arange(fes.ndof))
    for repeat in range(3):
        pre.Update(a.mat)
        with TaskManager():
            pre.Mult(rhs,out)
        result = out.FV().NumPy().copy()
        assert np.isfinite(result).all()
        if outputs:
            np.testing.assert_allclose(result,outputs[0],rtol=1e-9,atol=1e-10)
        outputs.append(result)
print('INTERNAL_SETUP_PARALLEL_OK')
'''
    result = subprocess.run([sys.executable, '-c', script, native.__file__],
        capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, (result.returncode,result.stdout,result.stderr)
    assert 'INTERNAL_SETUP_PARALLEL_OK' in result.stdout
