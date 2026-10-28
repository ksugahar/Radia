"""Shared resident B^T G B component: not FE assembly or complete MMM acceptance."""
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
        # Independent analytic nonlinear material reference over NG-owned E/weights.
        law=dict(knots=np.array([0.,1.,2.]),coefficients=np.array([1.,0.,1.,0.,1.,3.,4.,2.]),
                 integral_base=np.array([0.,.75,6.]),cap=1.5,barrier_slope=7.,last_field=4.875,initial_reluctivity=1.)
        csr=E.tocsr();csr.sort_indices();rng=np.random.default_rng(27)
        x=np.ascontiguousarray(.05*rng.standard_normal(fes.ndof));direction=np.ascontiguousarray(rng.standard_normal(fes.ndof))
        args=(csr.indptr.astype(np.uint64),csr.indices.astype(np.uint64),np.ascontiguousarray(csr.data),np.ascontiguousarray(payload["weights"]))
        def evaluate(state,delta=direction):
            return module._material_bh_evaluate(*args,np.ascontiguousarray(state),np.ascontiguousarray(delta),law,backend="gpu")
        v=(E@x).reshape(3,-1).T;dv=(E@direction).reshape(3,-1).T;amplitude=np.linalg.norm(v,axis=1)
        clipped=np.minimum(amplitude,1.5);extra=np.maximum(amplitude-1.5,0.)
        sec=np.divide(clipped+clipped**3+7*extra,amplitude,out=np.ones_like(amplitude),where=amplitude>0)
        dif=np.where(amplitude>1.5,7.,1+3*clipped**2)
        unit=np.divide(v,amplitude[:,None],out=np.zeros_like(v),where=amplitude[:,None]>0)
        point_field=payload["weights"][:,None]*sec[:,None]*v
        point_action=payload["weights"][:,None]*(sec[:,None]*dv+(dif-sec)[:,None]*unit*np.sum(unit*dv,axis=1)[:,None])
        gradient=E.T@point_field.T.ravel();action=E.T@point_action.T.ravel()
        energy=np.dot(payload["weights"],.5*clipped**2+.25*clipped**4+4.875*extra+3.5*extra**2)
        material=evaluate(x)
        assert np.linalg.norm(material["gradient"]-gradient)<=1e-11*max(np.linalg.norm(gradient),1e-30)
        assert np.linalg.norm(material["tangent_action"]-action)<=1e-11*max(np.linalg.norm(action),1e-30)
        assert abs(material["energy"]-energy)<=1e-12*max(abs(energy),1e-30)
        step=1e-6;plus=evaluate(x+step*direction);minus=evaluate(x-step*direction)
        assert np.linalg.norm((plus["gradient"]-minus["gradient"])/(2*step)-action)<=1e-6*max(np.linalg.norm(action),1e-30)
        assert abs((plus["energy"]-minus["energy"])/(2*step)-gradient@direction)<=1e-6*max(abs(gradient@direction),1e-12)
        zero=evaluate(np.zeros_like(x));assert np.array_equal(zero["gradient"],np.zeros_like(x))
        assert np.linalg.norm(zero["tangent_action"]-expected@direction)<=1e-11*np.linalg.norm(expected@direction)
        assert evaluate(x,np.array([]))["tangent_action"].size==0
        resident=module._ResidentBHMaterial(*args,fes.ndof,law,backend="gpu")
        setup=resident.stats();assert setup["epoch"]==1 and setup["setup_upload_bytes"]>0
        for _ in range(5):
            reused=resident.evaluate(x,direction)
            assert np.array_equal(reused["gradient"],material["gradient"])
            assert np.array_equal(reused["tangent_action"],material["tangent_action"])
        stats=resident.stats();assert stats["setup_upload_bytes"]==setup["setup_upload_bytes"]
        assert stats["state_upload_bytes"]==5*2*fes.ndof*8 and stats["evaluations"]==5
        with pytest.raises((ValueError,RuntimeError)):
            resident.evaluate(np.full_like(x,np.nan),direction)
        with pytest.raises((ValueError,RuntimeError)):
            resident.evaluate(np.full_like(x,1e308),direction)
        assert np.array_equal(resident.evaluate(x,direction)["gradient"],material["gradient"])
        assert resident.evaluate(x,np.array([]))["tangent_action"].size==0
        frozen=resident.linearize(x);epoch=frozen["linearization_epoch"]
        np.testing.assert_allclose(resident.apply_tangent(direction,epoch),action,rtol=1e-11,atol=1e-13)
        resident.evaluate(x+.01*direction,np.array([]))
        with pytest.raises((ValueError,RuntimeError),match="epoch"):
            resident.apply_tangent(direction,epoch)
        fresh=resident.linearize(x);assert fresh["linearization_epoch"]!=epoch
        np.testing.assert_allclose(resident.apply_tangent(direction,fresh["linearization_epoch"]),action,rtol=1e-11,atol=1e-13)
        assert np.array_equal(resident.evaluate(np.zeros_like(x),direction)["gradient"],np.zeros_like(x))
        resident.close();resident.close();assert resident.stats()["closed"]
        assert resident.stats()["device_payload_bound_bytes"]==0
        with pytest.raises((ValueError,RuntimeError)):
            resident.evaluate(x,direction)
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
        leaves = [dict(row=0, column=0, rows=2, columns=2, rank=0,
                       a1=np.array([4., .2, .2, 3.]), a2=np.array([])),
                  dict(row=0, column=2, rows=2, columns=2, rank=1,
                       a1=np.array([.2, .3]), a2=np.array([.1, .4])),
                  dict(row=2, column=2, rows=2, columns=2, rank=0,
                       a1=np.array([5., -.1, -.1, 6.]), a2=np.array([]))]
        permutation = np.array([2, 0, 3, 1], dtype=np.uint64)
        scale = np.array([2., .5, 1.5, .75])
        B = np.array([[1., 0., 2.], [0., 3., 0.], [4., 0., -1.], [0., .5, 2.]])
        rr, cc = np.nonzero(B)
        offsets = np.r_[0, np.cumsum(np.bincount(rr, minlength=4))].astype(np.uint64)
        indices = cc.astype(np.uint64)
        values = B[rr, cc].copy()
        L = np.zeros((4, 4));L[:2, :2] = leaves[0]["a1"].reshape(2, 2)
        L[2:, 2:] = leaves[2]["a1"].reshape(2, 2)
        L[:2, 2:] = np.outer(leaves[1]["a2"], leaves[1]["a1"])
        L[2:, :2] = L[:2, 2:].T
        G = np.empty((4, 4));G[np.ix_(permutation, permutation)] = L
        expected = B.T @ (scale[:, None]*G*scale[None, :]) @ B
        # Rank-one PSD leaf whose dot product is 1, not the naive rounded 0.
        cancellation=np.array([1e16,1.,-1e16])
        ci=np.arange(3,dtype=np.uint64);rp=np.arange(4,dtype=np.uint64)
        precise=module._ResidentHDiv([dict(row=0,column=0,rows=3,columns=3,rank=1,
            a1=cancellation,a2=cancellation)],ci,np.ones(3),3,rp,ci,np.ones(3),backend="gpu")
        try:
            np.testing.assert_array_equal(precise.apply(np.ones(3)),cancellation)
        finally:
            precise.close()
        # Correct the rounded product: (1+2^-27)*(1-2^-27)-1 = -2^-54.
        product_factors=np.array([1.+2.**-27,-1.])
        pi=np.arange(2,dtype=np.uint64);pr=np.arange(3,dtype=np.uint64)
        accurate=module._ResidentHDiv([dict(row=0,column=0,rows=2,columns=2,rank=1,
            a1=product_factors,a2=product_factors)],pi,np.ones(2),2,pr,pi,np.ones(2),backend="gpu")
        try:
            np.testing.assert_array_equal(accurate.apply(np.array([1.-2.**-27,1.])),product_factors*(-2.**-54))
        finally:
            accurate.close()
        args = (leaves, permutation, scale, 3, offsets, indices, values)
        with pytest.raises(ValueError, match="backend"):
            module._ResidentHDiv(*args, backend="cpu")
        with pytest.raises(TypeError):
            module._ResidentHDiv(leaves, permutation.astype(np.int32), *args[2:], backend="gpu")
        obj = module._ResidentHDiv(*args, backend="gpu")
        try:
            owned = obj.stats()["owned_device_bytes"]
            for x in (np.array([.2, -.3, .7]), np.zeros(3), np.ones(3)):
                np.testing.assert_allclose(obj.apply(x), expected @ x, rtol=1e-13, atol=1e-13)
            with pytest.raises(ValueError, match="dimension"):
                obj.apply(np.ones(4))
            mass=np.array([1.,2.,3.]);rhs=np.array([.4,-.7,.2])
            mi=np.arange(4,dtype=np.uint64);mj=np.arange(3,dtype=np.uint64)
            matrix=expected+2*np.diag(mass)
            reference=np.linalg.solve(matrix,rhs)
            for preconditioner,graph in (("diagonal",False),("chebyshev",True)):
                solution=obj.solve_linear_material(mi,mj,mass,rhs,2.,1e-12,100,preconditioner,graph)
                np.testing.assert_allclose(solution["x"],reference,rtol=1e-10,atol=1e-13)
                assert np.linalg.norm(matrix@solution["x"]-rhs)/np.linalg.norm(rhs)<=1e-11
                assert solution["relative_residual"]<=1e-11
                assert solution["backend"]=="gpu"
                assert obj.stats()["owned_device_bytes"]==owned
            with pytest.raises(ValueError,match="preconditioner"):
                obj.solve_linear_material(mi,mj,mass,rhs,2.,1e-12,100,"automatic")
            with pytest.raises(ValueError,match="dimension"):
                obj.solve_linear_material(mi,mj,mass,rhs[:-1],2.,1e-12,100,"diagonal")
            with pytest.raises(RuntimeError,match=r"iterations=1.*requested=.*actual="):
                obj.solve_linear_material(mi,mj,mass,rhs,2.,1e-12,1,"diagonal")
            law=dict(knots=np.array([0.,1.,2.]),coefficients=np.array([1.,0.,1.,0.,1.,3.,4.,2.]),integral_base=np.array([0.,.75,6.]),cap=1.5,barrier_slope=7.,last_field=4.875,initial_reluctivity=1.)
            material=module._ResidentBHMaterial(mi,mj,np.ones(3),np.ones(1),3,law,backend="gpu")
            target=np.array([.2,-.3,.1]);nonlinear_rhs=expected@target+(1+target@target)*target
            try:
                for pre in ("reference-diagonal","reference-chebyshev-graph"):
                    solution=obj.solve_nonlinear_material(material,mi,mj,np.ones(3),nonlinear_rhs,np.zeros(3),1e-11,5e-14,30,1000,pre)
                    independent=expected@solution["x"]+(1+solution["x"]@solution["x"])*solution["x"]-nonlinear_rhs
                    assert np.linalg.norm(independent)/np.linalg.norm(nonlinear_rhs)<=1e-11
                    np.testing.assert_allclose(solution["x"],target,rtol=1e-10,atol=1e-12)
                with pytest.raises(RuntimeError,match="Newton"):
                    obj.solve_nonlinear_material(material,mi,mj,np.ones(3),nonlinear_rhs,np.zeros(3),1e-11,5e-14,1,1000,"reference-diagonal")
                recovered=obj.solve_nonlinear_material(material,mi,mj,np.ones(3),nonlinear_rhs,np.zeros(3),1e-11,5e-14,30,1000,"reference-diagonal")
                np.testing.assert_allclose(recovered["x"],target,rtol=1e-10,atol=1e-12)
            finally:material.close()
            for _ in range(100):
                np.testing.assert_allclose(obj.apply(np.ones(3)), expected @ np.ones(3), rtol=1e-13, atol=1e-13)
                assert obj.stats()["owned_device_bytes"] == owned
        finally:
            obj.close()
        obj.close()
        assert obj.stats()["closed"] and obj.stats()["owned_device_bytes"] == 0
        with pytest.raises(RuntimeError, match="closed"):
            obj.apply(np.ones(3))
    finally:
        pass  # The Radia lazy runtime owns its pip library handles.
