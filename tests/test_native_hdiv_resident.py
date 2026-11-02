"""NG-owned material staging and the shared inverse BH law on the native GPU binding.

Not FE assembly, charge Gram or complete MMM acceptance; the device solve
routes are covered by the flat-TET workflow tests.
"""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pytest


def _check_material_staging(module, ng, tmp_path, family, order):
    from netgen.csg import Pnt
    from netgen.meshing import Mesh, MeshPoint, Element3D, Element2D, FaceDescriptor
    from scipy.sparse import csc_matrix
    mesh = Mesh(dim=3)
    mesh.SetMaterial(1, "body")
    fd = mesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    if family == "tet":
        vertices = ((0,0,0),(1,0,0),(0,1,0),(0,0,1))
        faces = ((0,2,1),(0,1,3),(1,2,3),(2,0,3))
    else:
        vertices = ((0,0,0),(1,0,0),(1,1,0),(0,1,0),
                    (0,0,1),(1,0,1),(1,1,1),(0,1,1))
        faces = ((0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))
    ids = [mesh.Add(MeshPoint(Pnt(*p))) for p in vertices]
    mesh.Add(Element3D(1, ids))
    for face in faces:
        mesh.Add(Element2D(fd, [ids[i] for i in face]))
    mesh_path = tmp_path / f"material-{family}-{order}.vol"
    mesh.Save(str(mesh_path))
    with ng.TaskManager():
        owner = module._MaterialStagingContext(str(mesh_path), order)
        payload = owner.quadrature()
        assert payload["epoch"] == 1 and payload["first_element"] == 0
        E = csc_matrix((payload["evaluation_values"], payload["evaluation_row_indices"],
                        payload["evaluation_column_offsets"]),
                       shape=(payload["evaluation_rows"], payload["dofs"]))
        reference_mesh = ng.Mesh(str(mesh_path))
        fes = ng.HDiv(reference_mesh, order=order)
        u, v = fes.TnT()
        mass = ng.BilinearForm(fes)
        mass += u*v*ng.dx
        mass.Assemble()
        r, c, values = mass.mat.COO()
        expected = csc_matrix((values, (r, c)), shape=(fes.ndof, fes.ndof)).toarray()
        actual = (E.T @ E.multiply(np.tile(payload["weights"], 3)[:, None])).toarray()
        assert np.linalg.norm(actual-expected) <= 1e-12*np.linalg.norm(expected)
        owner.quadrature()["evaluation_values"][:] = 0
        assert np.array_equal(payload["evaluation_values"], owner.quadrature()["evaluation_values"])
        with pytest.raises((ValueError, RuntimeError)):
            owner.quadrature(first_element=1)
        owner.close(); owner.close()
        with pytest.raises((ValueError, RuntimeError)):
            owner.quadrature()
    assert np.all(np.isfinite(payload["evaluation_values"]))


def test_shared_hdiv_resident_native_gpu(tmp_path):
    configured = os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR")
    if not configured:
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    import ngsolve as ng
    build = Path(configured)
    from radia.gpu import _native
    module = _native("_hdiv_cuda")
    try:
        assert "matlab.engine" not in sys.modules
        # Independent analytic nonlinear law h(a)=a+a^3 with linear barrier.
        knots=np.array([0.,1.,2.]);coefficients=np.array([1.,0.,1.,0.,1.,3.,4.,2.]);base=np.array([0.,.75,6.])
        amplitudes=np.array([0.,1e-8,.2,.999999,1.,1.2,1.5,1.500001,2.,10.])
        cap=1.5;last=cap+cap**3
        bh=module._inverse_bh_samples(knots,coefficients,base,amplitudes,cap,7.,last,1.,backend="gpu")
        clipped=np.minimum(amplitudes,cap);extra=np.maximum(amplitudes-cap,0.)
        field=clipped+clipped**3+7*extra
        expected_secant=np.divide(field,amplitudes,out=np.ones_like(field),where=amplitudes>0)
        expected_tangent=np.where(amplitudes>cap,7.,1+3*clipped**2)
        expected_energy=.5*clipped**2+.25*clipped**4+last*extra+3.5*extra**2
        np.testing.assert_allclose(bh["secant"],expected_secant,rtol=2e-14,atol=1e-15)
        np.testing.assert_allclose(bh["tangent"],expected_tangent,rtol=2e-14,atol=1e-15)
        np.testing.assert_allclose(bh["energy"],expected_energy,rtol=2e-14,atol=1e-15)
        for bad in (np.array([-1.]),np.array([np.nan])):
            with pytest.raises((ValueError,RuntimeError)):
                module._inverse_bh_samples(knots,coefficients,base,bad,cap,7.,last,1.,backend="gpu")
        with pytest.raises((ValueError,RuntimeError)):
            module._inverse_bh_samples(knots,coefficients,base,amplitudes,cap,7.,last+1.,1.,backend="gpu")
        recovered=module._inverse_bh_samples(knots,coefficients,base,amplitudes,cap,7.,last,1.,backend="gpu")
        assert np.array_equal(recovered["energy"],bh["energy"])
        # NGSolve owns each type-specific rule, mapped shape and mass oracle.
        for family in ("tet", "hex"):
            for order in (1, 2):
                _check_material_staging(module, ng, tmp_path, family, order)
    finally:
        pass  # The Radia lazy runtime owns its pip library handles.
