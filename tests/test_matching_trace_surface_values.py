"""Nonzero P1 surface values survive exact matching-trace condensation."""
import numpy as np
import pytest
import ngsolve as ng
from netgen.occ import Box, Glue, OCCGeometry, Pnt, Z
from radia.kelvin_solver import (
    solve_magnetostatic_matching_trace_total_reduced_omega as condensed,
    solve_magnetostatic_mixed_total_reduced_omega_kelvin as multiplier,
)


@pytest.fixture(scope='module')
def case():
    total=Box(Pnt(-1,-1,-1),Pnt(1,1,0));total.mat('total')
    total.faces.name='outer_total';total.faces.Max(Z).name='interface'
    reduced=Box(Pnt(-1,-1,0),Pnt(1,1,1));reduced.mat('reduced')
    reduced.faces.name='outer_reduced';reduced.faces.Min(Z).name='interface'
    mesh=ng.Mesh(OCCGeometry(Glue([total,reduced])).GenerateMesh(maxh=.65))
    space=ng.H1(mesh,order=1)
    fields=[]
    for cf in [-ng.x-.4*ng.y-ng.z,-ng.x-.4*ng.y-.01*ng.z+.05*(ng.x**2-ng.y**2),.05*(ng.x**2-ng.y**2)]:
        g=ng.GridFunction(space)
        for vertex in mesh.vertices:
            d=space.GetDofNrs(ng.NodeId(ng.VERTEX,vertex.nr))[0]
            g.vec[d]=float(cf(mesh(*vertex.point)))
        fields.append(g)
    trace,bt,br=fields
    kwargs=dict(mu_r_by_material={'total':100.,'reduced':1.},reduced_materials=('reduced',),total_materials=('total',),interface_boundary='interface',order=1,
                surface_dirichlet={'total':{'outer_total':bt},'reduced':{'outer_reduced':br}})
    return mesh,trace,kwargs


@pytest.mark.parametrize('solver,pre',[('direct','local'),('cg','local'),('cg','h1amg'),('cg','bddc_h1amg')])
def test_nonzero_surface_fields_match_original_saddle(case,solver,pre):
    mesh,trace,kwargs=case
    with ng.TaskManager():
        reference=multiplier(mesh,ng.CF((1,.4,1)),trace,1.,(0,0,0),kelvin_mats=(),**kwargs)
        result=condensed(mesh,ng.CF((1,.4,1)),trace,solver=solver,cg_preconditioner=pre,cg_tolerance=1e-13,**kwargs)
    assert result['linear_residual_relative']<1e-10
    for e in mesh.Elements(ng.VOL):
        xyz=np.mean([mesh[v].point for v in e.vertices],axis=0)
        p=mesh(*xyz)
        assert np.linalg.norm(np.array(result['H_cf'](p))-reference['H_cf'](p))<1e-9
    assert result['cg_preconditioner']==(pre if solver=='cg' else None)


def test_surface_junction_conflict_fails(case):
    mesh,trace,kwargs=case
    bad=dict(kwargs);bad['surface_dirichlet']={'total':{'outer_total':0},'reduced':{'outer_reduced':0}}
    with ng.TaskManager(),pytest.raises(ValueError,match='junction'):
        condensed(mesh,ng.CF((1,.4,1)),trace,**bad)


def test_invalid_cg_selection_fails(case):
    mesh,trace,kwargs=case
    with pytest.raises(ValueError,match='cg_preconditioner'):
        condensed(mesh,ng.CF((1,.4,1)),trace,solver='cg',cg_preconditioner='unknown',**kwargs)


def test_zero_volume_load_with_nonzero_boundary_is_accepted(case):
    mesh,trace,kwargs=case
    zero=ng.GridFunction(trace.space)
    boundary=ng.GridFunction(trace.space);boundary.Set(-ng.x-.4*ng.y)
    data=dict(kwargs,surface_dirichlet={'total':{'outer_total':boundary},'reduced':{'outer_reduced':boundary}})
    with ng.TaskManager():
        result=condensed(mesh,ng.CF((0,0,0)),zero,solver='cg',cg_tolerance=1e-13,**data)
    residual=result['linear_residual']['free_dofs']
    assert residual['rhs_l2']==0.
    assert residual['effective_load_l2']>0.
    assert residual['relative']<1e-10
    assert np.linalg.norm(np.array(result['H_cf'](mesh(.1,.2,.3)))-[1,.4,0])<1e-10


def test_iteration_exhaustion_never_returns_a_solution(case):
    from radia.kelvin_solver import LinearSolveNotConverged
    mesh,trace,kwargs=case
    with ng.TaskManager(),pytest.raises(LinearSolveNotConverged):
        condensed(mesh,ng.CF((1,.4,1)),trace,solver='cg',cg_max_iterations=1,cg_tolerance=1e-13,**kwargs)


@pytest.mark.parametrize('order',[2,3])
@pytest.mark.parametrize('solver,pre',[('direct','local'),('cg','local'),('cg','h1amg'),('cg','bddc_h1amg')])
def test_high_order_nonzero_boundary_and_source_modes(case,order,solver,pre):
    """Representable harmonic fields expose dropped edge/face trace modes."""
    mesh,_,_=case
    f=-ng.x-.4*ng.y+.02*(ng.x**2-ng.y**2)
    source=-ng.x-.4*ng.y-ng.z+.011*(ng.x**2-ng.y**2)
    if order==3:
        f+=.03*(ng.x**3-3*ng.x*ng.y**2)
        source+=.013*(ng.x**3-3*ng.x*ng.y**2)
    normal=1+.05*ng.x+.02*ng.y
    pt=f-ng.z*normal/100
    pa=f-ng.z*normal
    hs=-ng.CF(tuple(source.Diff(v) for v in (ng.x,ng.y,ng.z)))
    exact=mesh.MaterialCF({name:-ng.CF(tuple(cf.Diff(v) for v in (ng.x,ng.y,ng.z)))
                          for name,cf in [('total',pt),('reduced',pa)]})
    space=ng.H1(mesh,order=order)
    trace,bt,br=[ng.GridFunction(space) for _ in range(3)]
    with ng.TaskManager():
        for gf,cf in [(trace,source),(bt,pt),(br,pa-source)]:
            gf.Set(cf,bonus_intorder=6)
        result=condensed(mesh,hs,trace,order=order,solver=solver,cg_preconditioner=pre,
            cg_tolerance=1e-13,bonus_intorder=6,mu_r_by_material={'total':100.,'reduced':1.},
            reduced_materials=('reduced',),total_materials=('total',),interface_boundary='interface',
            surface_dirichlet={'total':{'outer_total':bt},'reduced':{'outer_reduced':br}})
        err=ng.Integrate(ng.InnerProduct(result['H_cf']-exact,result['H_cf']-exact),mesh,order=8)
    assert result['fes'].globalorder==order
    assert result['linear_residual_relative']<1e-10
    assert err<1e-17


def test_p3_reduced_load_maps_all_face_modes(case):
    from radia.kelvin_solver import _map_matching_h1_load_to_global
    mesh,_,_=case
    region=mesh.Materials('reduced')
    local=ng.H1(mesh,order=3,definedon=region)
    global_space=ng.H1(mesh,order=3)
    loads=[]
    with ng.TaskManager():
        for space in (local,global_space):
            load=ng.LinearForm(space)
            load+=(1+ng.x*ng.y)*space.TestFunction()*ng.dx(definedon=region)
            load.Assemble();loads.append(load.vec.FV().NumPy().copy())
    actual=_map_matching_h1_load_to_global(mesh,local,global_space,loads[0],order=3)
    assert np.linalg.norm(actual-loads[1])<1e-12
