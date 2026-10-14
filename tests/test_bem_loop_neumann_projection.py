"""Weak Neumann trace: exact constant reproduction on a planar P1 patch."""
import numpy as np
from radia.bem_loop_extension import _project_ring_neumann


def test_constant_normal_field_is_reproduced():
    pts = np.array([[0.,0.,0.], [2.,0.,0.], [2.,1.,0.], [0.,1.,0.]])
    tris = np.array([[0,1,2], [0,2,3]])
    areas = np.ones(2)
    normals = np.tile([0.,0.,1.], (2,1))
    mass = np.zeros((4,4))
    for t, a in zip(tris, areas):
        mass[np.ix_(t,t)] += a*(np.ones((3,3))+np.eye(3))/12
    q = _project_ring_neumann(pts, tris, areas, normals,
        lambda x: np.tile([2.,-3.,7.], (len(x),1)), np.linalg.inv(mass))
    np.testing.assert_allclose(q, -7., rtol=0, atol=2e-14)
