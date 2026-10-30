"""NG-owned pyramid material evaluation; never a mixed charge-Gram certificate.

Expected failures on pinned NGSolve 6.2.2608: connected tet-pyramid-hex orders
1/2 fail the unchanged 1e-11 constant-projection criterion (1.59e-11/1.78e-11).
They remain hard failures so earlier mass/checker defects cannot be hidden by
an overly broad xfail. GPU action is not executed after that failure. Rule-path
parity is not geometric exactness; the NG pyramid rule loses ~1e-10 volume.
"""
import json
import subprocess
import sys
import numpy as np
import pytest
from scipy.sparse import coo_matrix, csc_matrix, diags
from scipy.sparse.linalg import norm as sparse_norm


ng = pytest.importorskip("ngsolve")
gpu = pytest.importorskip("radia._hdiv_cuda", reason="Optional isolated native GPU binding")
from test_hdiv_connected_meshes import connected_mesh
from validation_test.hmatrix_cuda.hdiv_hdf5 import write_hdiv_h5


@pytest.mark.parametrize("family", ["pyramid", "tet-pyramid-hex"])
@pytest.mark.parametrize("order", [1, 2])
def test_ng_pyramid_mapped_material_mass_and_gpu_action(tmp_path, record_property, family, order):
    mesh = connected_mesh(family)
    path = tmp_path / "material.vol"
    mesh.ngmesh.Save(str(path))
    contract = tmp_path / "labels.json"
    contract.write_text(json.dumps(dict(schema="radia.vol-label-contract.v1",
        required=dict(materials=["body"],boundaries=["exterior"]),
        allowed=dict(materials=["body"],boundaries=["exterior"]))))
    checked = subprocess.run([sys.executable,"-m","cubit_mesh_export.check",str(path),
        "--strict-labels","--contract",str(contract),"--report-json",str(tmp_path/"check-vol.json"),"--format","json"],
        capture_output=True,text=True)
    assert checked.returncode == 0, checked.stderr + checked.stdout
    assert json.loads(checked.stdout)["passed"]
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        u, v = fes.TnT()
        mass = ng.BilinearForm(fes, symmetric=True)
        mass += ng.InnerProduct(u, v)*ng.dx(bonus_intorder=12)
        mass.Assemble()
        rows, cols, values = mass.mat.COO()
        M = coo_matrix((values, (rows, cols)), shape=(fes.ndof, fes.ndof)).tocsr()
        source = ng.LinearForm(fes)
        source += ng.CF((.3, -.7, 1.1))*v*ng.dx(bonus_intorder=12)
        source.Assemble()
        projected = ng.GridFunction(fes, novisual=True)
        projected.vec.data = mass.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")*source.vec
        x = np.ascontiguousarray(projected.vec.FV().NumPy())
        staging = gpu._MaterialStagingContext(str(path), order)
        try:
            payload = staging.quadrature()
            repeated = staging.quadrature()
            for key in ("evaluation_values","evaluation_row_indices","evaluation_column_offsets","weights"):
                np.testing.assert_array_equal(payload[key],repeated[key])
        finally:
            staging.close()
        with pytest.raises(ValueError,match="closed"):
            staging.quadrature()
        record_property("dofs",fes.ndof)
        record_property("points",len(payload["weights"]))
        E = csc_matrix((payload["evaluation_values"], payload["evaluation_row_indices"],
            payload["evaluation_column_offsets"]), shape=(payload["evaluation_rows"], fes.ndof)).tocsr()
        E.sort_indices()
        weights = np.ascontiguousarray(payload["weights"])
        assert np.all(weights > 0) and np.isfinite(weights).all()
        relative_mass = float(sparse_norm(E.T@diags(np.tile(weights, 3))@E-M)/sparse_norm(M))
        record_property("relative_mass_error",relative_mass)
        assert relative_mass <= 1e-11
        assert np.max(np.abs((E@x).reshape(3,-1).T-[.3,-.7,1.1])) <= 1e-11
        # Analytic H(M)=2*M: exact energy/gradient/tangent from independent mass.
        law = dict(knots=np.array([0.,100.]),coefficients=np.array([0.,0.,2.,0.]),
                   integral_base=np.array([0.,10000.]),cap=90.,barrier_slope=1000.,
                   last_field=180.,initial_reluctivity=2.)
        direction = np.ascontiguousarray(np.random.default_rng(20261020).standard_normal(fes.ndof)*.01)
        material = gpu._ResidentBHMaterial(E.indptr.astype(np.uint64),E.indices.astype(np.uint64),
            np.ascontiguousarray(E.data),weights,fes.ndof,law,backend="gpu")
        try:
            before = material.stats()
            evaluated = material.evaluate(x, direction)
            gradient_error = np.linalg.norm(evaluated["gradient"]-2*(M@x))/np.linalg.norm(2*(M@x))
            tangent_error = np.linalg.norm(evaluated["tangent_action"]-2*(M@direction))/np.linalg.norm(2*(M@direction))
            energy_error = abs(evaluated["energy"]-float(x@(M@x)))/abs(float(x@(M@x)))
            record_property("relative_gradient_error",gradient_error)
            record_property("relative_tangent_error",tangent_error)
            record_property("relative_energy_error",energy_error)
            assert max(gradient_error,tangent_error,energy_error) <= 1e-11
            assert material.stats()["setup_upload_bytes"] == before["setup_upload_bytes"]
        finally:
            material.close()
        assert material.stats()["closed"]
        # Validation-only exchange for the independent standalone MEX publisher.
        # This reference mass projection is not a GPU-user setup prerequisite.
        fields = dict(E=E.tocsr(),weights=weights,x=x,direction=direction,
            expected_gradient=2*(M@x),expected_tangent=2*(M@direction),
            expected_energy=np.array([float(x@(M@x))]),law_knots=law["knots"],
            law_coefficients=law["coefficients"].reshape(-1,4),law_integral_base=law["integral_base"])
        metadata=dict(units="SI",labels=["body","exterior"],family=family,order=order,
            law={key:value for key,value in law.items() if key not in ("knots","coefficients","integral_base")})
        write_hdiv_h5(tmp_path/"material-input.h5",fields,metadata)
