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


def revolved(profile, count=64):
    profile = np.asarray(profile)
    m = len(profile)
    pts = np.array([[r*np.cos(t), r*np.sin(t), z]
                    for t in np.arange(count)*2*np.pi/count for r,z in profile])
    tris = []
    for i in range(count):
        for j in range(m):
            a, b = i*m+j, ((i+1)%count)*m+j
            c, d = ((i+1)%count)*m+(j+1)%m, i*m+(j+1)%m
            tris.extend([(a,b,c), (a,c,d)])
    return pts, np.array(tris)


def test_step_plane_tangent_rays_are_retried_and_reported():
    # Both shelf heights coincide exactly with the original uniform z grid.
    pts, tris = revolved([[.005,-.016],[.030,-.016],[.030,-.012],
                         [.026,-.012],[.026,-.010],[.022,-.010],
                         [.022,.016],[.005,.016]])
    anchor, ring, diag = _auto_material_ring(pts,tris,return_diagnostics=True)
    rays = diag['ray_attempts']
    for z in (-.012,-.010):
        assert any(abs(r['z']-z)<1e-14 and r['reason']=='ambiguous_tangent_edge_or_vertex' for r in rays)
    assert diag['candidate_count'] > 0
    assert diag['candidates'][-1]['reason'] == 'accepted'
    assert diag['minimum_wall_distance'] > diag['clearance_tolerance']
    assert not np.isclose(anchor[1],[-.012,-.010]).any()
    # The rejected ray levels really do admit a carrier: odd/tangent
    # crossing counts there are not evidence that the material is absent.
    from radia.bem_loop_extension import _segment_surface_distance, _surface_winding
    theta=np.arange(32)*2*np.pi/32
    for z in (-.012,-.010):
        at_step=np.c_[.020*np.cos(theta),.020*np.sin(theta),np.full(32,z)]
        assert abs(_surface_winding(at_step[0],pts[tris]))==pytest.approx(1.)
        assert min(_segment_surface_distance(a,b,pts[tris])
                   for a,b in zip(at_step,np.roll(at_step,-1,axis=0)))>1e-3


@pytest.mark.parametrize('elliptic', [False, True])
def test_thin_and_non_axisymmetric_wall(elliptic):
    # A polygonal thin annulus; ellipse still has a common circular interior.
    inner, outer = ((.010,.016) if elliptic else (.010,.0101))
    pts,tris = revolved([[inner,-.01],[outer,-.01],[outer,.01],[inner,.01]],128)
    if elliptic:
        pts[:,0] *= 1.25
    (r,z),ring,diag = _auto_material_ring(pts,tris,return_diagnostics=True)
    assert diag['minimum_wall_distance'] > 0
    assert inner < r < outer
    if elliptic:
        assert r > 1.25*inner


@pytest.mark.parametrize('rotate', [False, True])
def test_segment_clearance_detects_crossing_between_sample_points(rotate):
    from radia.bem_loop_extension import _segment_surface_distance
    tri=np.array([[[.2,-1.,-1.],[.2,1.,-1.],[.2,0.,1.]]])
    rotation = np.linalg.qr(np.array([[1.,2.,3.],[4.,2.,1.],[2.,5.,1.]]))[0] if rotate else np.eye(3)
    distance = lambda a,b: _segment_surface_distance(np.array(a)@rotation,np.array(b)@rotation,tri@rotation)
    # Endpoints and midpoint all miss a wall at one fifth of the edge.
    assert distance([0.,0.,0.],[1.,0.,0.])==pytest.approx(0.,abs=1e-14)
    assert distance([0.,0.,0.],[.1,0.,0.])==pytest.approx(.1)
    # Coplanar edge contact and a coplanar segment disjoint from the triangle.
    assert distance([.2,-2.,-1.],[.2,2.,-1.])==pytest.approx(0.,abs=1e-14)
    assert distance([.2,-2.,-2.],[.2,2.,-2.])==pytest.approx(1.)


def test_failure_contains_candidate_reasons_and_clearance_diagnostics():
    import json
    pts,tris=stepped_surface();pts[:,0]+=.1
    with pytest.raises(ValueError) as caught:
        _auto_material_ring(pts,tris)
    diag=json.loads(str(caught.value).split('Diagnostics: ')[1])
    assert 'candidate_count' in diag and 'minimum_wall_distance' in diag
    assert all('reason' in c and 'minimum_wall_distance' in c for c in diag['candidates'])
    assert len({r['angle'] for r in diag['ray_attempts']}) > 1


def test_containment_certificate_requires_a_closed_oriented_surface():
    pts,tris=stepped_surface()
    with pytest.raises(ValueError,match='closed consistently oriented'):
        _auto_material_ring(pts,tris[:-1])
    tris[0]=tris[0,::-1]
    with pytest.raises(ValueError,match='closed consistently oriented'):
        _auto_material_ring(pts,tris)
