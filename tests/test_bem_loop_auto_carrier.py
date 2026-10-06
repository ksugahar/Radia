"""Automatic carrier on a stepped tube whose surface mean is outside material."""
import numpy as np
import pytest
from radia.bem_loop_extension import _auto_material_ring


def stepped_surface():
    profile=np.array([[.005,-.012],[.030,-.012],[.030,.002],[.012,.002],[.012,.050],[.005,.050]])
    count=64
    pts=np.array([[r*np.cos(t),r*np.sin(t),z] for t in np.arange(count)*2*np.pi/count for r,z in profile])
    tris=[]
    for i in range(count):
        for j in range(len(profile)):
            a=i*6+j;b=((i+1)%count)*6+j;c=((i+1)%count)*6+(j+1)%6;d=i*6+(j+1)%6
            tris.extend([(a,b,c),(a,c,d)])
    return pts,np.array(tris)


def test_stepped_tube_carrier_is_inside_without_cad_or_manual_parameters():
    pts,tris=stepped_surface()
    mean_r=np.linalg.norm(pts[:,:2],axis=1).mean()
    assert pts[:,2].mean()>.002 and mean_r>.012  # legacy mean is outside
    (r,z),ring=_auto_material_ring(pts,tris)
    assert -.012<z<.002 and .005<r<.030
    assert ring.shape==(256,3)
    assert np.allclose(np.linalg.norm(ring[:,:2],axis=1),r)
    assert np.allclose(ring[:,2],z)


def test_off_axis_hole_fails_instead_of_inventing_a_carrier():
    pts,tris=stepped_surface();pts[:,0]+=.1
    with pytest.raises(ValueError,match='through-hole around z'):
        _auto_material_ring(pts,tris)
