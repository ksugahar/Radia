"""Test helpers for programmatic meshes under NGSolve 6.2.2606."""

from netgen.meshing import EdgeDescriptor


def add_programmatic_edge_descriptors(mesh):
    """Populate descriptors omitted by NGSolve's structured-mesh helpers."""
    ngmesh = mesh.ngmesh
    if ngmesh.EdgeDescriptors():
        return

    if ngmesh.dim == 2:
        for index, name in enumerate(mesh.GetBoundaries(), start=1):
            descriptor = EdgeDescriptor()
            descriptor.edgenr = index
            descriptor.surfnr = (1, -1)
            descriptor.domin = 1
            descriptor.domout = 0
            descriptor.name = name
            assert ngmesh.Add(descriptor) == index
        return

    if ngmesh.dim == 3:
        for index, face in enumerate(ngmesh.FaceDescriptors(), start=1):
            descriptor = EdgeDescriptor()
            descriptor.edgenr = index
            descriptor.surfnr = (face.surfnr, -1)
            descriptor.domin = face.domin
            descriptor.domout = face.domout
            descriptor.name = face.bcname or f"boundary_{index}"
            assert ngmesh.Add(descriptor) == index
        return

    raise ValueError(f"unsupported programmatic mesh dimension: {ngmesh.dim}")


def curve_mesh(mesh, order):
    """Apply Curve after satisfying the 6.2.2606 descriptor contract."""
    add_programmatic_edge_descriptors(mesh)
    mesh.Curve(order)
    return mesh
