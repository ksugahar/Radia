from Compumag2025.SIBC.stray_meshes import *

from ngsolve import *
from mylibcem import *
import numpy as np
import pickle
import json
import time

# import matplotlib.pyplot as plt
# import sys
# sys.argv=["argv"]

# sys.argv=["temp"]
# font = {'size'   : 20}
# matplotlib.rc('font', **font)

linear=True
# linear=False

# Richardson = True
Richardson = False
FixedPoint = True
# FixedPoint = False

# Z_FP = (1+1j)*0.6e-3*CF(1) # Ohm
Z_FP = (0.7+0.6j)*1.0e-3*CF(1) # Ohm
# Z_FP = 0.9e-3 # Ohm

err_rel = 1e-2
N_it_nl_max = 100  # max. nonlin. iterations 

muAir=4*pi*1e-7
# muIron=muAir*1000
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
	
	# plt.plot(H0_amp_i,[Z_i[i].real for i in range(len(Z_i))],"b")
	# plt.plot(H0_amp_i,[Z_i[i].imag for i in range(len(Z_i))],"b--")
	
	# filename = "myZ.p"
	# with open(filename, "r") as fp:
		# mydata = json.load(fp)
		
		# H0_amp_i = mydata["Hvec"]
		# Z_i_r = mydata["Zvecreal"]
		# Z_i_i = mydata["Zvecimag"]
		# Z_i = [Z_i_r[i] + 1j*Z_i_i[i] for i in range(len(Z_i_r))]
		# P_eddy_i = mydata["Pvec"]
		# Q_i = mydata["Qvec"]
	
	# plt.plot(H0_amp_i,[Z_i[i].real for i in range(len(Z_i))],"r")
	# plt.plot(H0_amp_i,[Z_i[i].imag for i in range(len(Z_i))],"r--")
	# plt.show()
	# print(sfsff)
	
	
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
	# Zreal = Zimag
	# Zreal = [0.0005254077432022083 for i in range(len(Z_i))]
	# Zimag = [0.00040473429476440637 for i in range(len(Z_i))]
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
	# SICubeHole("SICube",farBND,a,b,c,delta,a2,b2)
	mesh = Mesh("SICube.vol")
	
	Draw(CF([1,2]),mesh,'test')
	

	
	pR,wR=np.polynomial.legendre.leggauss(5)
	nPhi = 10
	pZ,wZ=np.polynomial.legendre.leggauss(3)
	
	BStemp = J0*BiotSavartCylinder(1,[0,0,Z1],[0,0,Z2],Ri,Ro,[pR,pZ],[wR,wZ],nPhi)
	BStemp += J0*BiotSavartCylinder(1,[0,0,-Z2],[0,0,-Z1],Ri,Ro,[pR,pZ],[wR,wZ],nPhi)
	
	BSorder = 2
	# VBS = H1(mesh,order=BSorder)
	# BS=GridFunction(VBS**3)
	VBS = HCurl(mesh,order=BSorder)
	BS = GridFunction(VBS)
	
	BS.Set(BStemp)
	BS.Save("BSShieldingHole.sol")
	BS.Load("BSShieldingHole.sol")
	
	
	# Draw(BStemp,mesh,'BS')
	
	# print(sfsfsf)
	
	if not BiotSavart:
		BS = 0*1000*CF((0,0,1))
	
	# BStemp = 1000*CF((0,0,1))
	# BS.Set(BStemp)

	
	muvals={mat:muAir for mat in mesh.GetMaterials()}
	muvals["iron"]=muIron
	mu = CoefficientFunction([muvals[mat] for mat in mesh.GetMaterials()])

	sigmavals={mat:sigmaAir for mat in mesh.GetMaterials()}
	sigmavals["iron"]=sigmaIron
	sigma = CoefficientFunction([sigmavals[mat] for mat in mesh.GetMaterials()])
	
	VSpace = H1(mesh, order=order, dirichlet="top|right|back|bottom", definedon = "air|hole", complex=True)
	# VSpace = H1(mesh, order=order, dirichlet="top|bottom", definedon = "air|hole", complex=True)
	# VSpace = H1(mesh, order=order, dirichlet="top|bottom", definedon = "air", complex=True)
	print("ndof:",sum(VSpace.FreeDofs()))
	
	# print(sfsfsf)
	
	sol = GridFunction(VSpace)
	sol.Set(dir*z,BND)
	
	H = - grad(sol) + BS
	
	if not linear:
		Z = Zrealcurve(H.Norm()) + 1j*Zimagcurve(H.Norm())
	
	# Ztop = Zrealcurve(CF((grad(sol)[0],grad(sol)[1],0)).Norm()) + 1j*Zrealcurve(CF((grad(sol)[0],grad(sol)[1],0)).Norm())
	# Zright = Zrealcurve(CF((0,grad(sol)[1],grad(sol)[2])).Norm()) + 1j*Zrealcurve(CF((0,grad(sol)[1],grad(sol)[2])).Norm())
	# Zback = Zrealcurve(CF((grad(sol)[0],0,grad(sol)[2])).Norm()) + 1j*Zrealcurve(CF((grad(sol)[0],0,grad(sol)[2])).Norm())
	
	
	uPhi=VSpace.TrialFunction()
	vPhi=VSpace.TestFunction()
	
	a = BilinearForm(VSpace, symmetric = True)
	a += SymbolicBFI(1j*omega*mu*grad(uPhi)*grad(vPhi))
	# a += SymbolicBFI(1e-6*(uPhi*vPhi+uT*vT))
	
	f = LinearForm(VSpace)
	f += SymbolicLFI(1j*omega*mu*BS*grad(vPhi))
	
	eltype = QUAD
	
	if linear:
		a += SymbolicBFI(Z*uPhi.Trace().Deriv()*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
		f += SymbolicLFI(Z*BS*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
		
	else:
		intrule = IntegrationRule(eltype,2*(order+2))
		
		# munonlinTop_real=MuNonLinBiro3D(grad(sol).Norm(),Hvec,Zvec_real,0,0,differential=False,field="H")
		# munonlinTop_imag=MuNonLinBiro3D(grad(sol).Norm(),Hvec,Zvec_imag,0,0,differential=False,field="H")
		# munonlinSide_real=MuNonLinBiro3D((-grad(sol)+BS).Norm(),Hvec,Zvec_real,0,0,differential=False,field="H")
		# munonlinSide_imag=MuNonLinBiro3D((-grad(sol)+BS).Norm(),Hvec,Zvec_imag,0,0,differential=False,field="H")
		
		
		# munonlinTop_real.SetUp(mesh,intrule)
		# munonlinTop_real.SetActiveElements(mesh,["sTop"])
		# munonlinTop_imag.SetUp(mesh,intrule)
		# munonlinTop_imag.SetActiveElements(mesh,["sTop"])
		# munonlinSide_real.SetUp(mesh,intrule)
		# munonlinSide_real.SetActiveElements(mesh,["sRight","sBack","sLeft","sFront"])
		# munonlinSide_imag.SetUp(mesh,intrule)
		# munonlinSide_imag.SetActiveElements(mesh,["sRight","sBack","sLeft","sFront"])
		
		# ZTop = munonlinTop_real + 1j*munonlinTop_imag
		# ZSide = munonlinSide_real +1j*munonlinSide_imag
		
		# a += SymbolicBFI(Zright*(uPhi.Trace().Deriv()[1]*vPhi.Trace().Deriv()[1]+uPhi.Trace().Deriv()[2]*vPhi.Trace().Deriv()[2]),definedon=mesh.Boundaries("sRight"))
		# a += SymbolicBFI(Zback*(uPhi.Trace().Deriv()[0]*vPhi.Trace().Deriv()[0]+uPhi.Trace().Deriv()[2]*vPhi.Trace().Deriv()[2]),definedon=mesh.Boundaries("sBack"))
		# a += SymbolicBFI(Ztop*(uPhi.Trace().Deriv()[0]*vPhi.Trace().Deriv()[0]+uPhi.Trace().Deriv()[1]*vPhi.Trace().Deriv()[1]),definedon=mesh.Boundaries("sTop"))
		
		# a += SymbolicBFI(ZTop*uPhi.Trace().Deriv()*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sTop")).SetIntegrationRule(eltype,intrule)
		# a += SymbolicBFI(ZSide*uPhi.Trace().Deriv()*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack")).SetIntegrationRule(eltype,intrule)
		# f += SymbolicLFI(ZSide*BS*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack")).SetIntegrationRule(eltype,intrule)
	
		if Richardson:
			a += SymbolicBFI(Z*uPhi.Trace().Deriv()*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
			f += SymbolicLFI(Z*BS*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
		if FixedPoint:
			a += SymbolicBFI(Z_FP*uPhi.Trace().Deriv()*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
			f += SymbolicLFI((Z_FP-Z)*grad(sol).Trace()*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
			f += SymbolicLFI(Z*BS*vPhi.Trace().Deriv(),definedon=mesh.Boundaries("sRight|sBack|sTop|sLeft|sFront"))
		
	
	c = Preconditioner(a, type="direct")
	# c = Preconditioner(a, type="bddc")
	
	testvals={mat:0 for mat in mesh.GetMaterials()}
	testvals["iron"]=1
	test = CoefficientFunction([testvals[mat] for mat in mesh.GetMaterials()])
	

	if linear:
	
		a.Assemble()
		f.Assemble()
		
		solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=2, needsassembling=False)
		# solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=1000, needsassembling=False)
		
		# muCF = mu

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
				
				# sol.vec.data = 0.5*(sol.vec+solold.vec)
				
				err = sum([abs(sol.vec[i]-solold.vec[i]) for i in range(len(sol.vec))])/sum([abs(sol.vec[i]) for i in range(len(sol.vec))])
				if it == 1:
					err0 = err
				print("error:",err)
				print("error/error0:",err/err0)
				if 100*err/err0 < err_rel:
					break
					
				# if it>1 and munonlinTop_real.HasConverged(1,0.1) and munonlinTop_imag.HasConverged(1,0.1) and munonlinSide_real.HasConverged(1,0.1) and munonlinSide_imag.HasConverged(1,0.1):
					# break
					
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
				
				# sol.vec.data = 0.5*(sol.vec+solold.vec)
				
				err = sum([abs(sol.vec[i]-solold.vec[i]) for i in range(len(sol.vec))])/sum([abs(sol.vec[i]) for i in range(len(sol.vec))])
				if it == 1:
					err0 = err
				print("error:",err)
				print("error/error0:",err/err0)
				if 100*err/err0 < err_rel:
					break
					
				# if it>1 and munonlinTop_real.HasConverged(1,0.1) and munonlinTop_imag.HasConverged(1,0.1) and munonlinSide_real.HasConverged(1,0.1) and munonlinSide_imag.HasConverged(1,0.1):
					# break
					
				if it == N_it_nl_max:
					print("too many iterations")
					break
				
		# muCF=test*munonlin.GetGlobalFunction() + (CoefficientFunction(1)-test)*muAir
	
	# sol.Save("sol_N7_O2_ref_z.sol")
	B = mu*H
	
	# Draw(BS,mesh,'BS')
	Draw(H*CF(1),mesh,'H')
	Draw(B,mesh,'B')
	Bred = -mu*grad(sol)
	Draw(Bred,mesh,'Bred')
	
	print("Htest:",Integrate(H,mesh,BND,definedon=mesh.Boundaries("sTop")))
	
	# lossesActiveCoef = 0.5*1/sigma*J*Conj(J)
	# lossesReactiveCoef = omega/2*muCF*H*Conj(H)
	
	# print("side")
	# temp = Integrate(1,mesh,definedon=mesh.Boundaries("sRight|sBack"))
	# print(Integrate(H/temp,mesh,definedon=mesh.Boundaries("sRight|sBack")))
	# print("top")
	# temp = Integrate(1,mesh,definedon=mesh.Boundaries("sTop"))
	# print(Integrate(grad(sol)/temp,mesh,definedon=mesh.Boundaries("sTop")))
	
	# S1Ht=Integrate(0.5*Z*H.Norm()**2,mesh,definedon=mesh.Boundaries("sTop"))
	Losses1=Integrate(Pcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sTop"))
	Losses2=Integrate(Pcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sRight"))
	Losses3=Integrate(Pcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sBack"))
	# Losses4=Integrate(Pcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sLeft"))
	# Losses5=Integrate(Pcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sFront"))
	print("Losses:")
	# print("S1Ht:",S1Ht)
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
	ReactivePower=ReactivePower1+ReactivePower2+ReactivePower3 #+ReactivePower4+ReactivePower5

	# ReactivePower1=Integrate(Qcurve((-grad(sol)+BS).Norm()**2),mesh,definedon=mesh.Boundaries("sTop"))
	# ReactivePower2=Integrate(Qcurve((-grad(sol)+BS).Norm()**2),mesh,definedon=mesh.Boundaries("sRight|sBack|sLeft|sFront"))
	# ReactivePower=ReactivePower1+ReactivePower2
	# losses=ReactivePower-1j*lossesActive
	
	print("")
	print("Losses:",Losses)
	print("Reactive power:",ReactivePower)

	print('')
	time_end = time.time()
	print('Computation Time = ',time_end-time_start)
	print('')

	
	