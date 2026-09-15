"""Coincident hex interfaces: preserve boundary flags and complete face maps."""
import itertools
import pytest
import ngsolve as ng
from netgen.meshing import Mesh, MeshPoint, Pnt, Element3D, Element2D, FaceDescriptor, IdentificationType
from radia.periodic_h1 import periodic_h1_single_interface as _periodic_h1_single_interface


def split_mesh(identify=True, scale=1., shape='box'):
    native = Mesh(dim=3)
    corners = ((0,0,0),(1,0,0),(1,1,0),(0,1,0),
               (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    faces = ((0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))
    sides = []
    for side in range(2):
        native.SetMaterial(side+1, 'left' if side == 0 else 'right')
        def position(q):
            x,y,z = side+q[0]/2, q[1]/2, q[2]/2
            if shape == 'skew':
                x,y = x+.3*y+.2*z, y+.1*z
            elif shape == 'rotated':
                x,y = -y,x
            return Pnt(scale*x,scale*y,scale*z)
        points = {q:native.Add(MeshPoint(position(q)))
                  for q in itertools.product(range(3),repeat=3)}
        sides.append(points)
        fd_outer = native.Add(FaceDescriptor(surfnr=side*2+1,domin=side+1,bc=1))
        fd_inner = native.Add(FaceDescriptor(surfnr=side*2+2,domin=side+1,bc=2))
        for cell in itertools.product(range(2),repeat=3):
            keys = [tuple(a+b for a,b in zip(cell,c)) for c in corners]
            native.Add(Element3D(side+1,[points[q] for q in keys]))
            for face in faces:
                nodes = [keys[i] for i in face]
                if not any(len({q[k] for q in nodes}) == 1 and nodes[0][k] in (0,2) for k in range(3)):
                    continue
                interface = all(q[0] == (2 if side == 0 else 0) for q in nodes)
                native.Add(Element2D(fd_inner if interface else fd_outer,[points[q] for q in nodes]))
    native.SetBCName(0,'outer'); native.SetBCName(1,'interface')
    if identify:
        for y,z in itertools.product(range(3), repeat=2):
            native.AddPointIdentification(sides[0][2,y,z].nr,sides[1][0,y,z].nr,1,IdentificationType.PERIODIC)
    return ng.Mesh(native)


@pytest.mark.parametrize('order',[1,2,3])
def test_uniform_field_and_face_dofs(order):
    mesh = split_mesh()
    base = ng.H1(mesh,order=order,dirichlet='outer')
    periodic = _periodic_h1_single_interface(base)
    vertices = {b:a for (a,b),_ in mesh.GetPeriodicNodePairs(ng.VERTEX)}
    by_vertices = {tuple(sorted(v.nr for v in f.vertices)):f for f in mesh.faces}
    for face in mesh.faces:
        if not all(v.nr in vertices for v in face.vertices):
            continue
        master = by_vertices[tuple(sorted(vertices[v.nr] for v in face.vertices))]
        assert periodic.GetDofNrs(face) == periodic.GetDofNrs(master)
    fes = ng.Compress(periodic)
    u,v = fes.TnT()
    a = ng.BilinearForm(fes,symmetric=True)
    a += ng.InnerProduct(ng.grad(u),ng.grad(v))*ng.dx
    f = ng.LinearForm(fes)
    source = ng.CoefficientFunction((.2,-.3,1.))
    f += ng.InnerProduct(source,ng.grad(v))*ng.dx
    a.Assemble(); f.Assemble()
    g = ng.GridFunction(fes)
    g.vec.data = a.mat.Inverse(fes.FreeDofs(),inverse='pardiso')*f.vec
    assert ng.Integrate(ng.InnerProduct(ng.grad(g),ng.grad(g)),mesh) < 1e-20


def test_missing_interface_rejected():
    with pytest.raises(Exception,match='Exactly one'):
        _periodic_h1_single_interface(ng.H1(split_mesh(False),order=2))


def test_other_space_rejected():
    with pytest.raises(Exception,match='scalar H1'):
        _periodic_h1_single_interface(ng.HCurl(split_mesh(),order=1))


@pytest.mark.parametrize('order',[2,3])
@pytest.mark.parametrize('shape',['box','skew','rotated'])
@pytest.mark.parametrize('scale',[.001,1.,10.])
def test_nonuniform_harmonic_field(order, shape, scale, record_property):
    mesh = split_mesh(scale=scale,shape=shape)
    fes = ng.Compress(_periodic_h1_single_interface(ng.H1(mesh,order=order,dirichlet='outer')))
    u,v = fes.TnT()
    exact = (ng.x**2-ng.y**2+.3*ng.x*ng.y)/scale**2 + .2*ng.z/scale
    gradient = ng.CoefficientFunction(((2*ng.x+.3*ng.y)/scale**2,
                                      (-2*ng.y+.3*ng.x)/scale**2,.2/scale))
    a = ng.BilinearForm(fes,symmetric=True)
    a += ng.InnerProduct(ng.grad(u),ng.grad(v))*ng.dx(bonus_intorder=4)
    a.Assemble()
    g = ng.GridFunction(fes)
    g.Set(exact,definedon=mesh.Boundaries('outer'))
    g.vec.data -= a.mat.Inverse(fes.FreeDofs(),inverse='pardiso')*(a.mat*g.vec)
    diff = ng.grad(g)-gradient
    error = (ng.Integrate(ng.InnerProduct(diff,diff),mesh)/ng.Integrate(ng.InnerProduct(gradient,gradient),mesh))**.5
    record_property('relative_field_error',float(error))
    assert error < 1e-11
