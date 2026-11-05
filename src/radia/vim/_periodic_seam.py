"""Orientation contract of a rotational HDiv seam tied by NGSolve ``Periodic``.

NGSolve orients every face, and the face basis on it, by the global numbers of
the face's vertices.  ``Periodic(HDiv)`` copies the DOFs of a slave seam face
from its master face without a sign or mode map.  The copy is the rotation's
trace map only when each slave face orders its vertices like the image of its
master face.  On other numberings the tie targets the wrong face or the wrong
face-mode orientation, and the reduced space silently loses the exact periodic
field.  ``check_periodic_seam_orientation`` rejects such meshes.
``renumber_periodic_seam`` returns a copy whose seam numbering satisfies the
contract on every face type.
"""

import numpy as np


def periodic_seam_vertex_map(ngmesh, boundary_vertices, boundaries):
    """Master -> slave vertex map (1-based) of the PERIODIC identifications between the two seam boundaries."""
    first, second = boundaries
    vertex_map = {}
    for endpoint_a, endpoint_b, *_rest in ngmesh.GetIdentifications():
        endpoint_a = int(getattr(endpoint_a, "nr", endpoint_a))
        endpoint_b = int(getattr(endpoint_b, "nr", endpoint_b))
        if endpoint_a in boundary_vertices[first] and endpoint_b in boundary_vertices[second]:
            master, slave = endpoint_a, endpoint_b
        elif endpoint_b in boundary_vertices[first] and endpoint_a in boundary_vertices[second]:
            master, slave = endpoint_b, endpoint_a
        else:
            continue
        if vertex_map.get(master, slave) != slave:
            raise ValueError("vim.Solve: seam vertex %d is identified with two vertices" % master)
        vertex_map[master] = slave
    return vertex_map


def _seam_faces(mesh, boundaries):
    import ngsolve as ng
    available = tuple(mesh.GetBoundaries())
    faces = {name: [] for name in boundaries}
    for element in mesh.Elements(ng.BND):
        name = available[element.index]
        if name in faces:
            faces[name].append(tuple(int(vertex.nr) + 1 for vertex in element.vertices))
    return faces


def periodic_seam_misoriented_faces(mesh, boundaries):
    """(misoriented, total) seam face pairs; a pair is consistent when its vertex numbers keep their order under the map."""
    first, second = boundaries
    faces = _seam_faces(mesh, boundaries)
    vertex_map = periodic_seam_vertex_map(
        mesh.ngmesh, {name: {v for face in faces[name] for v in face} for name in boundaries}, boundaries)
    slave_faces = {frozenset(face) for face in faces[second]}
    misoriented = 0
    for face in faces[first]:
        try:
            image = [vertex_map[v] for v in face]
        except KeyError as missing:
            raise ValueError("vim.Solve: seam vertex %s has no PERIODIC partner" % missing) from None
        if frozenset(image) not in slave_faces:
            raise ValueError("vim.Solve: the image of seam face %s is not a face of %r" % (face, second))
        order_master, order_image = np.argsort(face), np.argsort(image)
        if not np.array_equal(order_master, order_image):
            misoriented += 1
    if len(faces[first]) != len(faces[second]):
        raise ValueError("vim.Solve: the seam boundaries %r and %r have %d and %d faces"
                         % (first, second, len(faces[first]), len(faces[second])))
    return misoriented, len(faces[first])


def check_periodic_seam_orientation(mesh, boundaries):
    """Fail loudly unless every seam face pair is orientation-consistent for ``Periodic(HDiv)``."""
    misoriented, total = periodic_seam_misoriented_faces(mesh, boundaries)
    if misoriented:
        raise ValueError(
            "vim.Solve: %d of %d periodic seam face pairs order their vertex numbers differently from their image; "
            "NGSolve Periodic(HDiv) would tie those faces to the wrong face or face-mode orientation and silently "
            "change the reduced space. Renumber the mesh with radia.vim.renumber_periodic_seam(mesh, %r) first."
            % (misoriented, total, tuple(boundaries)))


def renumber_periodic_seam(mesh, boundaries=("periodic_min", "periodic_max")):
    """A copy of a linear netgen mesh whose slave seam vertices follow their masters' order.

    Every non-slave vertex keeps its number.  The slave vertices keep their set of
    numbers, reassigned in the order of their masters, so the master -> slave map is
    monotone and every seam face pair is orientation-consistent.  Elements keep their local
    vertex order, so geometry and element orientation do not change.  Curved,
    second-order or geometry-bound meshes are rejected.
    Returns an ``ngsolve.Mesh``.
    """
    import ngsolve as ng
    from netgen.meshing import Element1D, Element2D, Element3D, FaceDescriptor, Mesh, MeshPoint, Pnt
    boundaries = tuple(str(name) for name in boundaries)
    if not isinstance(mesh, ng.Mesh):
        mesh = ng.Mesh(mesh)
    source = mesh.ngmesh
    if source.dim != 3:
        raise ValueError("renumber_periodic_seam: needs a 3D mesh")
    if source.GetCurveOrder() > 1:
        raise ValueError("renumber_periodic_seam: curved meshes are not supported; renumber before curving")
    if len(source.Elements0D()):
        raise ValueError("renumber_periodic_seam: 0D elements are not supported")
    faces = _seam_faces(mesh, boundaries)
    vertex_map = periodic_seam_vertex_map(
        source, {name: {v for face in faces[name] for v in face} for name in boundaries}, boundaries)
    if not vertex_map:
        raise ValueError("renumber_periodic_seam: no PERIODIC identification joins %r" % (boundaries,))
    slaves = set(vertex_map.values())
    if slaves & set(vertex_map):
        raise ValueError("renumber_periodic_seam: a seam vertex is both master and slave")
    count = len(source.Points())
    new = {v: v for v in range(1, count + 1)}
    for number, master in zip(sorted(slaves), sorted(vertex_map)):
        new[vertex_map[master]] = number
    old_of = sorted(new, key=new.__getitem__)

    def renumbered(element):
        if len(element.points) != len(element.vertices):
            raise ValueError("renumber_periodic_seam: second-order elements are not supported")
        return [new[int(p.nr)] for p in element.vertices]

    target = Mesh(dim=3)
    coordinates = source.Coordinates()
    for old in old_of:
        target.Add(MeshPoint(Pnt(*map(float, coordinates[old - 1]))))
    for index in range(1, source.GetNDomains() + 1):
        target.SetMaterial(index, source.GetMaterial(index))
    for index in range(1, source.GetNFaceDescriptors() + 1):
        fd = source.FaceDescriptor(index)
        target.Add(FaceDescriptor(surfnr=fd.surfnr, domin=fd.domin, domout=fd.domout, bc=fd.bc))
        target.SetBCName(index - 1, source.GetBCName(index - 1))
    for index in range(source.GetNCD2Names()):
        target.SetCD2Name(index + 1, source.GetCD2Name(index))
    for element in source.Elements3D():
        target.Add(Element3D(element.index, renumbered(element)))
    for element in source.Elements2D():
        target.Add(Element2D(element.index, renumbered(element)))
    for element in source.Elements1D():
        target.Add(Element1D(renumbered(element), index=element.index, edgenr=element.edgenr))
    for endpoint_a, endpoint_b, identnr in source.GetIdentifications():
        a, b = int(getattr(endpoint_a, "nr", endpoint_a)), int(getattr(endpoint_b, "nr", endpoint_b))
        target.AddPointIdentification(new[a], new[b], identnr=int(identnr), type=source.GetIdentificationType(int(identnr)))
    result = ng.Mesh(target)
    if periodic_seam_misoriented_faces(result, boundaries)[0]:
        raise AssertionError("renumber_periodic_seam: the renumbered seam is still misoriented")
    return result
