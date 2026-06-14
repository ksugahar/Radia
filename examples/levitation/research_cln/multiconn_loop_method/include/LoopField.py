from ngsolve import *
import sys
sys.path.append('..\include')
from MatrixSolver import MatrixSolver as solver 

def surface_mesh_from_boundary(boundary ):
    m=boundary
    bfaces=set()
    bedges=set()
    bnodes=set()
    for f in m.Elements():
        bfaces.add(f)
        nv=len(f.vertices)
        for i in range(nv):
            e=f.edges[i]
            bedges.add(e)
            n=f.vertices[i]
            bnodes.add(n) 
    return (bfaces,bedges, bnodes)

def surface_genus(smesh, p):
    nf=len(smesh[0])
    print("Surface face count= ", nf)  
    ne=len(smesh[1])
    nv=len(smesh[2])
    print("Surface edge count= ", ne)
    print("Surface node count= ", nv)
    nx=nv-ne+nf
    g=int(p-nx/2)
    return g

def random_edge(smesh):
    import random
    edges=smesh[1]
    ne=len(smesh[1])
    number = random.randint(1, ne)
    return list(edges)[number-1]  

def LoopFields(model, **kwargs):
    feOrder=1
    mesh=model.GetMesh()
    default_values = {"connected": 1}
    default_values.update(kwargs)
    connected=default_values["connected"]
    
    boundary=mesh.Boundaries(model.conductor_boundary)
    smesh=surface_mesh_from_boundary(boundary)
    g=surface_genus(smesh, connected)
    print("genus= ", g)

    fes = HCurl(mesh, order=feOrder, nograds=True, definedon=model.total_air_region, dirichlet=model.total_boundary)
    u,v = fes.TnT()
    fesPhi = H1(mesh, order=feOrder, definedon=model.total_air_region, dirichlet=model.total_boundary)
    phi,psi= fesPhi.TnT()     
    loops=[]
    for k in range(g):
        gfu = GridFunction(fes)
        id=random_edge(smesh)
        edge_dofs = fes.GetDofNrs(id)[0]   
            
        gfu.vec[:] = 0
        gfu.vec[edge_dofs] = 1
        print("edge DOF = ", edge_dofs)

        fes.FreeDofs().__setitem__(edge_dofs,False)

        a = BilinearForm(fes)
        a += curl(u)*curl(v)*dx(model.total_air_region)
        f=LinearForm(fes)
        with TaskManager():
            a.Assemble()
            f.Assemble()
        fr=-a.mat*gfu.vec

        gfu=solver.iccg_solve(fes, gfu, a, fr.Evaluate(), tol=1.e-16, max_iter=200, accel_factor=0)
        #Draw(gfu, mesh)
        
        gfPhi = GridFunction(fesPhi)
        a = BilinearForm(fesPhi)
        a += grad(phi)*grad(psi)*dx
        f=LinearForm(fesPhi)
        f += grad(psi)*gfu*dx
        with TaskManager():
            a.Assemble()
            f.Assemble()

        gfPhi=solver.iccg_solve(fesPhi, gfPhi, a, f.vec.FV(), tol=1.e-16, max_iter=200, accel_factor=1.0)  
        gfw=gfu-grad(gfPhi)

        gft=gfw
        for kd in range(len(loops)):
            #prod=Integrate(gfw*loops[kd]*dx, mesh, definedon=model.total_air_region)
            prod=Integrate(gfw*loops[kd]*dx, mesh)
            print("k=", k, "  kd=", kd, "   prod=", prod)
            gft=gft-prod*loops[kd]
            #for i in range(len(gft.vec)):
            #    gft.vec[i] -=prod*loops[kd].vec[i]
      
        norm2=Integrate(gft*gft*dx, mesh)
        print("k=", k,  "   norm2=", norm2)
        print
        norm=sqrt(norm2)
        gft=gft/norm
        #for i in range(len(gft.vec)):
        #    gft.vec[i]/=norm

        #Draw(gft, mesh)
        loops.append(gft)       
    return loops

def loopFieldCouplings(loopFields, s, model, fesTOmega):
    mesh=model.mesh
    mu=model.Mu
    sig=model.Sigma
    cond=model.conductive_region
    (T,omega),(W,psi) = fesTOmega.TnT()
    fv=[]
    fafv=[]
    gfs=[]
    g=len(loopFields)
    for n in range(g):
        loopField=loopFields[n]

        gfTOmega=GridFunction(fesTOmega)
        gfT, gfOmega=gfTOmega.components
        gfT.Set(loopField, BND, mesh.Boundaries(model.conductor_boundary))
        gfs.append(gfT)
    
        f=LinearForm(fesTOmega)
        f += 1./(s*sig)*curl(gfT)*curl(W)*dx(cond)
        f += mu*gfT*(W+grad(psi))*dx(cond)
        f += mu*loopField*grad(psi)*dx(model.total_air_region)

        with TaskManager():
            f.Assemble()
        fv.append(f)
  
        tmp=[]
        for n2 in range(n+1):
            gfT2=gfs[n2]
            faf =Integrate(1./(s*sig)*curl(gfT)*curl(gfT2)*dx(cond), mesh) 
            faf +=Integrate(mu*gfT*gfT2*dx(cond), mesh)
            faf +=Integrate(mu*loopField*loopFields[n2]*dx(model.total_air_region), mesh)   
            tmp.append(faf)
            if n !=n2: fafv[n2].append(faf)
        fafv.append(tmp)
    return fv, fafv