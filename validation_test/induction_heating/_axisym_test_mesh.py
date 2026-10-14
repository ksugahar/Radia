"""In-memory algebra fixture, not a production Cubit/VOL export route."""


def make_axisymmetric_mesh():
    import ngsolve
    from netgen.meshing import (Mesh, MeshPoint, Pnt, FaceDescriptor,
                               EdgeDescriptor, Element1D, Element2D)

    radius, height, count = .025, .025, 9
    mesh = Mesh(dim=2)
    points = {(i,j): mesh.Add(MeshPoint(Pnt(radius*i/count,
              -height/2+height*j/count,0)))
              for j in range(count+1) for i in range(count+1)}
    face = mesh.Add(FaceDescriptor(bc=1, domin=1, surfnr=1))
    mesh.SetMaterial(1,'workpiece')
    for j in range(count):
        for i in range(count):
            mesh.Add(Element2D(face,[points[i,j],points[i+1,j],
                                     points[i+1,j+1],points[i,j+1]]))
    boundaries = [
        ('bot',[(i,0,i+1,0) for i in range(count)]),
        ('outer',[(count,j,count,j+1) for j in range(count)]),
        ('top',[(i,count,i+1,count) for i in range(count)]),
        ('axis',[(0,j,0,j+1) for j in range(count)]),
    ]
    for index,(name,edges) in enumerate(boundaries,1):
        mesh.SetBCName(index-1,name)
        descriptor = EdgeDescriptor()
        descriptor.edgenr=index; descriptor.surfnr=(index,-1)
        descriptor.domin=1; descriptor.domout=0; descriptor.name=name
        mesh.Add(descriptor)
        for i,j,k,l in edges:
            mesh.Add(Element1D([points[i,j],points[k,l]],index=index))
    return ngsolve.Mesh(mesh)
