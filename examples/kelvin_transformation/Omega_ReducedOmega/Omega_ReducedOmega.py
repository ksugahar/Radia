import math
from ngsolve import *
from ngsolve.webgui import Draw
import numpy as np
import sys
sys.path.append(r'..\include')
from MatrixSolver import MatrixSolver as solver 
sys.path.append(r'..\bin\Release') 
from EMPY_Field import *
from Static_Method import Static_Method

class Omega_ReducedOmega_Method(Static_Method):
    def __init__(self,  model, **kwargs):  
        super().__init__(model,  **kwargs)
        
    def Calc(self, **kwargs):
        default_values = {"feOrder":1,
                          "boundaryCD":"Bn0", 
                          "Kelvin": "off"
                         }
        default_values.update(kwargs)
        self.feOrder=default_values["feOrder"]
        boundaryCD=default_values["boundaryCD"]
        Kelvin=default_values["Kelvin"]
        self.Kelvin=Kelvin

        feOrder=self.feOrder
        model=self.model
        mesh=model.mesh
        self.mesh=mesh
        
        import time
        start_time = time.perf_counter()


        total_region=model.total_region
        reduced_region=model.reduced_region
        total_boundary=model.total_boundary
        reduced_boundary=model.reduced_boundary
        Bn0_boundary=model.Bn0_boundary
        Ht0_boundary=model.Ht0_boundary
        self.Mu=model.Mu
        Mu=self.Mu

        self.total_region=total_region
        self.reduced_region=reduced_region

        #field=UNIF(0,0,1,0)
        coil=model.coil.field
        Ov=Ofield(coil)
        Bv=Bfield(coil)
        mu0=4.e-7*math.pi
        Hv=Bv/mu0
        zero=(0,0,0)
        Bs_dic = {"iron":zero, "A_domain":zero, "Omega_domain":Bv, "Kelvin":zero, 'default':zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in mesh.GetMaterials()])
        Os_dic = {"iron":0, "A_domain":0, "Omega_domain":Ov, "Kelvin":0, 'default':zero}
        Os = CoefficientFunction([Os_dic[mat] for mat in mesh.GetMaterials()])

        if Kelvin=="off":
            if boundaryCD=="Ht0" :
                fes=H1(mesh, order=feOrder, dirichlet=reduced_boundary+"|"+Ht0_boundary)
            else:
                fes=H1(mesh, order=feOrder, dirichlet=Ht0_boundary)
            
        else:
            fes=H1(mesh, order=feOrder, dirichlet=Ht0_boundary)
            fes=Periodic(fes)

        #dofLimit=1.0e6
        if fes.ndof > self.dofLimit:
            print(" fespace Dof >  dofLimit, DOF=", fes.ndof)
            return 0
            
        self.fes=fes
        omega,psi = fes.TnT() 
        gfOmega = GridFunction(fes)
        a= BilinearForm(fes)
        a +=Mu*(grad(omega)*grad(psi))*dx(total_region)
        a +=Mu*(grad(omega)*grad(psi))*dx(reduced_region) 

        if Kelvin=="on":
            rs=self.model.rKelvin
            xs=2*rs
            r=sqrt((x-xs)*(x-xs)+y*y+z*z)
            fac=rs*rs/r/r
            a +=Mu*fac*(grad(omega)*grad(psi))*dx("Kelvin") 
            #a +=Mu*(grad(omega)*grad(psi))*dx("Kelvin") 
        with TaskManager():
            a.Assemble()
        normal = specialcf.normal(mesh.dim)

        #surfaceOmega=HtoOmega(mesh, total_boundary, feOrder, Hv)
        #surfaceOmega=y
        # Calculate Dirichlet condition terms
        gfOmega.Set(Ov, BND, mesh.Boundaries(total_boundary))
        #gfOmega.Set(surfaceOmega, BND, mesh.Boundaries(total_boundary))

        f = LinearForm(fes)
        f +=Mu*grad(gfOmega)*grad(psi)*dx(reduced_region)
        with TaskManager():
            f.Assemble() 
        #remove components of the Dirichlet boundary
        fcut = np.array(f.vec.FV())[fes.FreeDofs()]
        np.array(f.vec.FV(), copy=False)[fes.FreeDofs()] = fcut

        # Add Neumann condition terms
        f += (normal*Bv)*psi*ds(total_boundary)
        with TaskManager():
            f.Assemble()
        
        gfOmega=self.Solve(fes, a, f)
        """
        gfOmega = GridFunction(fes)   #Clear gfOmega
        gfOmega=solver.iccg_solve(fes, gfOmega, a, f.vec.FV(), tol=1.e-8, max_iter=1000, accel_factor=0, divfac=100, diviter=10,
                     scaling=True, complex=False,logplot=True)
        """
        fesOt=H1(mesh, order=feOrder, definedon=total_region)
        fesOr=H1(mesh, order=feOrder, definedon=reduced_region+"|Kelvin")
        Ot=GridFunction(fesOt)
        Orr=GridFunction(fesOr)
        Oxr=GridFunction(fesOr)

        Ot.Set(gfOmega,VOL, definedon=total_region)
        Orr.Set(gfOmega,VOL, definedon=reduced_region+"|Kelvin")
        """
        if Kelvin=="on":
            #Orr.Set(gfOmega,VOL, definedon=reduced_region+"|"+"Kelvin")
            print(reduced_region+"|"+"Kelvin")
            #Orr.Set(1, definedon=reduced_region+"|"+"Kelvin")
            Orr.Set(1, definedon="Kelvin")
            Draw(Orr,mesh)
        """
        #Oxr.Set(surfaceOmega, BND, mesh.Boundaries(total_boundary))
        Oxr.Set(Ov, BND, mesh.Boundaries(total_boundary))

        Bt=grad(Ot)*Mu
        Or=Orr-Oxr
        #print("*** Omega ***")
        #Draw(Or, mesh)
        
        Br=(grad(Orr)-grad(Oxr))*mu0
        BField=Bt+Br+Bs
        self.BField=BField
        self.JField=0

        print("feOrder=", feOrder,"  ", "ndof=",fes.ndof,"  ")
        #self.CalcResult(model, BField)

        #Draw (Ot, mesh, order=3)  
        #Draw (Or, mesh, order=3) 
        
        self.Br=Br
        self.Bt=Bt
        #Draw ((Ot+Or+Os)*Mu, mesh, order=3, min=0., max=1.0, deformation=False)  
        """   
        mip = mesh(0,0,0)
        print("center magnetic field = ", BField(mip),"  ")
        Wm=Integrate(BField*BField/Mu*dx("iron"), mesh)
        print("magnetic energy=", Wm,"  ")


        Draw ((Ot+Or+Os)*mu, mesh, order=3, min=0., max=1.0, deformation=False)  
        print("**** B field ****")
        Draw (BField, mesh, order=feOrder, min=0, max=5, deformation=False)
        """
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print(f"経過時間: {elapsed_time:.4f} 秒  ")
        return 1
        
    def CalcError(self):
        mesh=self.mesh
        feOrder=self.feOrder
        total_region=self.total_region
        reduced_region=self.reduced_region
        Mu=self.Mu
        
        #gfu=self.Ot
        fesBd = HDiv(mesh, order=feOrder-1, definedon=total_region)
        Bd=GridFunction(fesBd, "flux")
        B = self.Bt
        Bd.Set(B)
        err=1/Mu*(B-Bd)*(B-Bd)
        eta = Integrate(err, mesh, VOL, element_wise=True)

        #gfur=self.Or
        fsBdr = HDiv(mesh, order=feOrder-1, definedon=reduced_region)
        Bdr=GridFunction(fsBdr, "flux")
        Br = self.Br
        Bdr.Set(Br)
        err_r=1/Mu*(Br-Bdr)*(Br-Bdr)
        eta_r = Integrate(err_r, mesh, VOL, element_wise=True) 
 
        if self.Kelvin=="on":
            fkBdr = HDiv(mesh, order=feOrder-1, definedon="Kelvin")
            Bdr=GridFunction(fkBdr, "flux")
            
            rs=self.model.rKelvin
            xs=2*rs
            r=sqrt((x-xs)*(x-xs)+y*y+z*z)
            fac=rs*rs/r/r
            Br=Br*fac
            Bdr.Set(Br)
            err_rk=1/(Mu*fac)*(Br-Bdr)*(Br-Bdr)
            eta_rk =Integrate(err_rk, mesh, VOL, element_wise=True) 
            print(" maxerr = ", max(eta), max(eta_r), max(eta_rk) )
            eta_r = eta_r+eta_rk
 
        error=eta+eta_r
        maxerr = max(error)
        print ("ndof=", self.fes.ndof, " maxerr = ", maxerr)
        return maxerr, error



        