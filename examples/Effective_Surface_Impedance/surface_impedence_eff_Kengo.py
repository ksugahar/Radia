from Compumag2025.SIBC.stray_meshes import *

from ngsolve import *
from mylibcem import *
import numpy as np
import pickle
import json
import time

linear=True
# linear=False

# Richardson = True
Richardson = False
FixedPoint = True
# FixedPoint = False

Z_FP = (0.7+0.6j)*1.0e-3*CF(1) # Ohm

err_rel = 1e-2
N_it_nl_max = 100  # max. nonlin. iterations 

muAir=4*pi*1e-7
muIron=1/420
sigmaIron=2e6#*1e-6
sigmaAir=1

BiotSavart = True
# BiotSavart = False
J0 = 1e6 # A/m^2
Z1 = 0.4 # m
Z2 = 0.45 # m
Ri = 0.3 # m
Ro = 0.4 # m


dir = 1000
if BiotSavart:
	dir = 0

order = 1


f = 50
omega = 2*pi*f

delta = sqrt(2/(sigmaIron*muIron*omega))

H_KL_ref=[0,42,53,62,70,79,88,100,113,132,157,193,255,376,677,1624,1e9]
B_KL_ref=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1e9*muAir]

# Z = 0.6139e-3
Z = 0.4325e-3*(1+1j)
Pfac = 0.000210
Qfac = 0.000210

Pcurve = BSpline (2, [0]+[0,1e12], [0,Pfac*1e12])
Qcurve = BSpline (2, [0]+[0,1e12], [0,Qfac*1e12])


time_start = time.time()


if not linear:
	# filename = "eff_values_nonlinear.pkl"
	# filename = "eff_values_nonlinear_complex.pkl"
	filename = "eff_values_nonlinear_complex_big.pkl"
	
	loaded = pickle.load(open( filename, "rb" ))
	H0_amp_i = loaded["H0_amp_i"]
	p_eddy_i = loaded["p_eddy_i"]
	p_hyst_i = loaded["p_hyst_i"]
	q_i = loaded["q_i"]
	P_eddy_i = loaded["P_eddy_i"]
	P_hyst_i = loaded["P_hyst_i"]
	Q_i = loaded["Q_i"]
	ti = loaded["ti"]
	Z_i = loaded["Z_i"]
	
	Z_i[0] = Z_i[1]
	
	Hvec = H0_amp_i + [1e6]
	Zvec_real = [Z_i[i].real*Hvec[i] for i in range(len(Z_i))] + [Z_i[-1].real]
	Zvec_imag = [Z_i[i].imag*Hvec[i] for i in range(len(Z_i))] + [Z_i[-1].imag]
	Pvec = P_eddy_i +[Pfac*1e12]
	Qvec = Q_i +[Qfac*1e12]
	
	Pcurve = BSpline (2, [0]+[Hvec[i]**2 for i in range(len(Hvec))], Pvec)
	Qcurve = BSpline (2, [0]+[Hvec[i]**2 for i in range(len(Hvec))], Qvec)
	
	
	Zreal = [Z_i[i].real for i in range(len(Z_i))]
	Zimag = [Z_i[i].imag for i in range(len(Z_i))]
	Zrealcurve = BSpline (2, [0]+Hvec, Zreal)
	Zimagcurve = BSpline (2, [0]+Hvec, Zimag)

	
if linear:
	
	Hvec = [0,1e6]
	Zvec_real = [0,0.4325e-3*1e6]
	Zvec_imag = [0,0.4325e-3*1e6]

	Pvec = [0,Pfac*1e12]
	Qvec = [0,Qfac*1e12]
	
	Pcurve = BSpline (2, [0]+[Hvec[i]**2 for i in range(len(Hvec))], Pvec)
	Qcurve = BSpline (2, [0]+[Hvec[i]**2 for i in range(len(Hvec))], Qvec)
	



ngsglobals.msg_level = 5

with TaskManager():
	
	farBND = 1
	
	a = 0.5
	b = 0.5
	c = 0.1
	
	a2 = 0.4
	b2 = 0.4
	
	SICube("SICube",farBND,a,b,c,delta)
	mesh = Mesh("SICube.vol")
	
	Draw(CF([1,2]),mesh,'test')
	

	
	pR,wR=np.polynomial.legendre.leggauss(5)
	nPhi = 10
	pZ,wZ=np.polynomial.legendre.leggauss(3)
	
	BStemp = J0*BiotSavartCylinder(1,[0,0,Z1],[0,0,Z2],Ri,Ro,[pR,pZ],[wR,wZ],nPhi)
	BStemp += J0*BiotSavartCylinder(1,[0,0,-Z2],[0,0,-Z1],Ri,Ro,[pR,pZ],[wR,wZ],nPhi)
	
	BSorder = 2
	VBS = HCurl(mesh,order=BSorder)
	BS = GridFunction(VBS)
	
	BS.Set(BStemp)
	BS.Save("BSShieldingHole.sol")
	BS.Load("BSShieldingHole.sol")
	

	if not BiotSavart:
		BS = 0*1000*CF((0,0,1))

	muvals={mat:muAir for mat in mesh.GetMaterials()}
	muvals["iron"]=muIron
	mu = CoefficientFunction([muvals[mat] for mat in mesh.GetMaterials()])

	sigmavals={mat:sigmaAir for mat in mesh.GetMaterials()}
	sigmavals["iron"]=sigmaIron
	sigma = CoefficientFunction([sigmavals[mat] for mat in mesh.GetMaterials()])
	
	VSpace = H1(mesh, order=order, dirichlet="top|right|back|bottom", definedon = "air|hole", complex=True)
	print("ndof:",sum(VSpace.FreeDofs()))
	
	sol = GridFunction(VSpace)
	sol.Set(dir*z,BND)
	
	H = - grad(sol) + BS
	
	if not linear:
		Z = Zrealcurve(H.Norm()) + 1j*Zimagcurve(H.Norm())
	
	uPhi=VSpace.TrialFunction()
	vPhi=VSpace.TestFunction()
	
	a = BilinearForm(VSpace, symmetric = True)
	a += SymbolicBFI(1j*omega*mu*grad(uPhi)*grad(vPhi))
	
	f = LinearForm(VSpace)
	f += SymbolicLFI(1j*omega*mu*BS*grad(vPhi))
	
	eltype = QUAD
	
	if linear:
		a += SymbolicBFI(Z*uPhi.Trace().Deriv()*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
		f += SymbolicLFI(Z*BS*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
		
	else:
		intrule = IntegrationRule(eltype,2*(order+2))

		if Richardson:
			a += SymbolicBFI(Z*uPhi.Trace().Deriv()*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
			f += SymbolicLFI(Z*BS*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
		if FixedPoint:
			a += SymbolicBFI(Z_FP*uPhi.Trace().Deriv()*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
			f += SymbolicLFI((Z_FP-Z)*grad(sol).Trace()*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
			f += SymbolicLFI(Z*BS*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
		
	
	c = Preconditioner(a, type="direct")
	
	testvals={mat:0 for mat in mesh.GetMaterials()}
	testvals["iron"]=1
	test = CoefficientFunction([testvals[mat] for mat in mesh.GetMaterials()])
	

	if linear:
	
		a.Assemble()
		f.Assemble()
		
		solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=2, needsassembling=False)

	else:

		it=0
		solold = GridFunction(VSpace)

		if Richardson:

			while True:
				it += 1
				print ("Iteration",it)
				
				
				solold.vec.data=sol.vec

				a.Assemble()
				f.Assemble()
				
				solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=2, needsassembling=False)
				
				err = sum([abs(sol.vec[i]-solold.vec[i]) for i in range(len(sol.vec))])/sum([abs(sol.vec[i]) for i in range(len(sol.vec))])
				if it == 1:
					err0 = err
				print("error:",err)
				print("error/error0:",err/err0)
				if 100*err/err0 < err_rel:
					break
					
				if it == N_it_nl_max:
					print("too many iterations")
					break
				
		if FixedPoint:

			a.Assemble()

			while True:
				it += 1
				print ("Iteration",it)
				
				
				solold.vec.data=sol.vec

				f.Assemble()
				
				solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=2, needsassembling=False)
				
				err = sum([abs(sol.vec[i]-solold.vec[i]) for i in range(len(sol.vec))])/sum([abs(sol.vec[i]) for i in range(len(sol.vec))])
				if it == 1:
					err0 = err
				print("error:",err)
				print("error/error0:",err/err0)
				if 100*err/err0 < err_rel:
					break
					
				if it == N_it_nl_max:
					print("too many iterations")
					break
				
	B = mu*H

	Draw(H*CF(1),mesh,'H')
	Draw(B,mesh,'B')
	Bred = -mu*grad(sol)
	Draw(Bred,mesh,'Bred')
	
	print("Htest:",Integrate(H,mesh,BND,definedon=mesh.Boundaries("sTop")))
	
	Losses1=Integrate(Pcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sTop"))
	Losses2=Integrate(Pcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sRight"))
	Losses3=Integrate(Pcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sBack"))
	print("Losses:")
	print("oben:",Losses1)
	print("rechts:",Losses2)
	print("hinten:",Losses3)
	Losses=Losses1+Losses2+Losses3 #+Losses4+Losses5
	
	ReactivePower1=Integrate(Qcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sTop"))
	ReactivePower2=Integrate(Qcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sRight"))
	ReactivePower3=Integrate(Qcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sBack"))
	print("Reactive power:")
	print("oben:",ReactivePower1)
	print("rechts:",ReactivePower2)
	print("hinten:",ReactivePower3)
	ReactivePower=ReactivePower1+ReactivePower2+ReactivePower3

	print("")
	print("Losses:",Losses)
	print("Reactive power:",ReactivePower)

	print('')
	time_end = time.time()
	print('Computation Time = ',time_end-time_start)
	print('')

	
	