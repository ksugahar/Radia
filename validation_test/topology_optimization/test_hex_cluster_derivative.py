"""Required heavier validation of native low-rank HEX derivative contractions."""
import os
import numpy as np


def test_native_hex_cluster_leaf_contraction_matches_analytic_dense_derivative():
    import ngsolve as ng
    ng.SetNumThreads(int(os.environ.get("RADIA_TEST_NGSOLVE_THREADS", "1")))
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram
    # The support-radius admissibility of 2026-09-05 widened the near family,
    # and 2x2x1 at leafsize 8 stopped admitting anything at all: every one of
    # its 187 leaves came back dense, so the leaf contraction this test exists
    # to check was only ever compared on dense blocks.  Measured 2026-09-20:
    # leafsize 8 admits nothing up to eta 6 on 3x3x1 and nothing at all up to
    # 4x4x2; 3x3x1 at leafsize 4 and eta 3 gives 24 low-rank leaves of rank 8
    # with real low-rank work, so the assertion guards a mesh that can satisfy
    # it.
    mesh=MakeStructured3DMesh(hexes=True,nx=3,ny=3,nz=1,
        mapping=lambda x,y,z:(x+.031*y*z,y+.019*x*z,z+.027*x*y))
    fes=ng.HDiv(mesh,order=1)
    with ng.TaskManager():
        cb=_charge_basis_hex(fes,cob_quad=3)
        _,gram,_=build_charge_gram(fes,eps=1e-8,leafsize=4,eta=3.0)
    assert gram.stats()["n_lowrank"]>0
    cells=np.asarray(cb["cell_nodes"]).reshape(-1,27,3)
    faces=np.asarray(cb["face_nodes"]).reshape(-1,9,3)
    def velocity(x):
        weight=np.maximum(0.,1.-2.*np.abs(x[...,0]-.5))
        return weight[...,None]*np.array([.021,-.017,.013])
    vc,vf=velocity(cells),velocity(faces)
    dense=np.asarray(gram.hex_charge_gram_directional_derivative(vc,vf))
    left=np.linspace(-.4,.7,dense.shape[0]);right=np.cos(np.arange(dense.shape[0]))
    observed=gram.directional_derivative_contractions("hex",vc[None],vf[None],left,right)[0]
    left_many=np.ascontiguousarray(np.stack((left,2*left-right)))
    observed_many=np.asarray(gram.directional_derivative_contractions_many(
        "hex",np.ascontiguousarray(np.stack((vc,2*vc))),
        np.ascontiguousarray(np.stack((vf,2*vf))),left_many,right))
    np.testing.assert_allclose(observed,left@dense@right,rtol=2e-6,atol=2e-9)
    np.testing.assert_allclose(observed_many,
        [[left@dense@right,2*left@dense@right],
         [(2*left-right)@dense@right,2*(2*left-right)@dense@right]],
        rtol=2e-6,atol=2e-9)

