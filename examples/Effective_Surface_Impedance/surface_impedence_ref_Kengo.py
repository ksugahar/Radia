from Compumag2025.SIBC.stray_meshes import *

from ngsolve import *
from mylibcem import *
import numpy as np
import pickle
import time

# sys.argv=["temp"]
# font = {'size'   : 20}
# matplotlib.rc('font', **font)

linear=True
# linear=False

# PreisBiro = True
PreisBiro = False
FixedPoint = True
# FixedPoint = False
# Newton = True
Newton = False

if Newton:
	print(Newton_funktioniert_im_Komplexen_nicht)

mu_FP = 1.0/300.0 # !!!
# mu_FP = 1.0/420.0 # !!!
# mu_FP = 1.0/800.0 # !!!

maxits_Newton = 5
tol_Newton = 1e-3

err_nl_max = 1.0 # max. error nonlinear material
err_nl_av = 0.1 # average error nonlinear material

N_it_nl_max = 100  # max. nonlin. iterations 


# loadsol=True
loadsol=False

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


f = 50
omega = 2*pi*f

delta = sqrt(2/(sigmaIron*muIron*omega))

H_KL_ref=[0,42,53,62,70,79,88,100,113,132,157,193,255,376,677,1624,1e9]
B_KL_ref=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1e9*muAir]

# H_KL_ref=[0,1e9]
# B_KL_ref=[0,1e9*muIron]


time_start = time.time()


ngsglobals.msg_level = 5

with TaskManager():

	order = 1
	
	farBND = 1
	
	a = 0.5
	b = 0.5
	c = 0.1
	
	a2 = 0.4
	b2 = 0.4
	
	SICube("SICube",farBND,a,b,c,delta)
	# SICubeHollow("SICube",farBND,a,b,c,delta)
	# SICubeHole("SICube",farBND,a,b,c,delta,a2,b2)
	mesh = Mesh("SICube.vol")
	
	Draw(CF([1,2,3]),mesh,'test')
	# print(sfsfsf)
	
	
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
	
	# hp = IfPos(x-(a-15*delta),0,1)*IfPos(y-(b-15*delta),0,1)*IfPos(z-(c-15*delta),0,1)
	# Draw(hp,mesh,'hollow')
	
	
	muvals={mat:muAir for mat in mesh.GetMaterials()}
	muvals["iron"]=muIron
	mu = CoefficientFunction([muvals[mat] for mat in mesh.GetMaterials()])

	sigmavals={mat:sigmaAir for mat in mesh.GetMaterials()}
	sigmavals["iron"]=sigmaIron
	sigma = CoefficientFunction([sigmavals[mat] for mat in mesh.GetMaterials()])
	
	VCurl = HCurl(mesh, order=order, dirichlet="sRight|sBack|sTop|hollow", definedon = "iron|hole", complex=True, nograds=True)
	# VH1 = H1(mesh, order=order+1, dirichlet="top|right|back|bottom|left|front", complex=True)
	VH1 = H1(mesh, order=order+1, dirichlet="top|right|back|bottom|sBottom|hollow", definedon = "iron|air|hole", complex=True)
	VSpace = FESpace([VCurl,VH1])
	print("ndof:",sum(VSpace.FreeDofs()))
	
	
	hollowvals={bnd:1 for bnd in mesh.GetBoundaries()}
	hollowvals["hollow"]=0
	hollowcut = CoefficientFunction([hollowvals[bnd] for bnd in mesh.GetBoundaries()])
	
	sol = GridFunction(VSpace)
	sol_old = GridFunction(VSpace)
	sol.components[1].Set(dir*hollowcut*z,BND)
	
	H = sol.components[0] - grad(sol.components[1]) + BS
	
	uT,uPhi=VSpace.TrialFunction()
	vT,vPhi=VSpace.TestFunction()
	
	a = BilinearForm(VSpace, symmetric = True)
	a += SymbolicBFI(1/sigma*curl(uT)*curl(vT),definedon=mesh.Materials("iron|hole"))
	a += SymbolicBFI(1j*omega*mu*grad(uPhi)*grad(vPhi),definedon=~mesh.Materials("iron|hole"))
	a += SymbolicBFI(1e-6*(uPhi*vPhi+uT*vT)) # Für Tensor erforderlich lt. Schöbinger
	
	f = LinearForm(VSpace)
	f += SymbolicLFI(1j*omega*mu*BS*grad(vPhi), definedon=~mesh.Materials("iron|hole"))
	
	eltype = HEX
	
	if linear:
	
		a += SymbolicBFI(1j*omega*mu*(uT - grad(uPhi))*(vT - grad(vPhi)),definedon=mesh.Materials("iron|hole"))
		f += SymbolicLFI(-1j*omega*mu*BS*(vT-grad(vPhi)), definedon=mesh.Materials("iron|hole"))
	
	else:

		intrule = IntegrationRule(eltype,2*(order+2))
		
		# if not Newton:
		if PreisBiro:
			munonlin=MuNonLinBiro3D(H.Norm(),H_KL_ref,B_KL_ref,113,0.8,differential=False,field="H")
		if FixedPoint:
			munonlin=MuNonLinBiro3D(H.Norm(),H_KL_ref,B_KL_ref,0,0,differential=False,field="H")
		if not Newton:	
			munonlin.SetUp(mesh,intrule)
			munonlin.SetActiveElements(mesh,["iron"])

		if Newton:
			SPLorder = 2
			# munonlin = (BSpline(SPLorder, [0 for i in range(SPLorder-1)]+list(Bcurve), list(Hcurve))(sqrt(B*B+1e-6)))/(sqrt(B*B+1e-6))
			munonlinNM = (BSpline(SPLorder, [0 for i in range(SPLorder-1)]+list(H_KL_ref), list(B_KL_ref))(sqrt(H*H+1e-6)))/(sqrt(H*H+1e-6))
		
		if PreisBiro:
			a += SymbolicBFI(1j*omega*munonlin*(uT - grad(uPhi))*(vT - grad(vPhi)),definedon=mesh.Materials("iron")).SetIntegrationRule(eltype,intrule)
			f += SymbolicLFI(-1j*omega*munonlin*BS*(vT-grad(vPhi)), definedon=mesh.Materials("iron")).SetIntegrationRule(eltype,intrule)
	
		if FixedPoint:
			a += SymbolicBFI(1j*omega*mu_FP*(uT - grad(uPhi))*(vT - grad(vPhi)),definedon=mesh.Materials("iron"))
			f += SymbolicLFI(1j*omega*(mu_FP-munonlin)*H*(vT - grad(vPhi)),definedon=mesh.Materials("iron")).SetIntegrationRule(eltype,intrule)
			f += SymbolicLFI(-1j*omega*mu_FP*BS*(vT-grad(vPhi)), definedon=mesh.Materials("iron"))

		if Newton:
			a += SymbolicBFI(1j*omega*munonlinNM*(uT - grad(uPhi))*(vT - grad(vPhi)),definedon=mesh.Materials("iron")).SetIntegrationRule(eltype,intrule)
			a += SymbolicBFI(1j*omega*munonlinNM*BS*(vT-grad(vPhi)), definedon=mesh.Materials("iron")).SetIntegrationRule(eltype,intrule)
	
	c = Preconditioner(a, type="direct")
	# c = Preconditioner(a, type="bddc")
	
	testvals={mat:0 for mat in mesh.GetMaterials()}
	testvals["iron"]=1
	test = CoefficientFunction([testvals[mat] for mat in mesh.GetMaterials()])
	
	if loadsol:
		sol.Load("SI_nonline_o1.sol")
		if linear:
			muCF = mu
		else:
			munonlin.Update()
			muCF=test*munonlin.GetGlobalFunction() + (CoefficientFunction(1)-test)*muAir
	else:
	
		if linear:

			a.Assemble()
			f.Assemble()
			
			solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=2, needsassembling=False)
			# solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=1000, needsassembling=False)
			
			muCF = mu

		else:

			if PreisBiro:

				it=0
				while True:
					it += 1
					print ("Iteration",it)
					
					if it > 1:
						munonlin.Update()

					a.Assemble()
					f.Assemble()
					
					solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=2, needsassembling=False)
						
					if it > 1 and munonlin.HasConverged(err_nl_max,err_nl_av,True):
						break
						
					if it == N_it_nl_max:
						print("too many iterations")
						break
					
			if FixedPoint:
			# if False:

				it=0
				a.Assemble()
				while True:
					it += 1
					print ("Iteration",it)
					
					if it > 1:
						munonlin.Update()

					f.Assemble()
					
					sol_old.vec.data = sol.vec
					solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=2, needsassembling=False)
						
					sol.vec.data = 0.5*(sol.vec+sol_old.vec)

					if it > 1 and munonlin.HasConverged(err_nl_max,err_nl_av,True):
						break
						
					if it == N_it_nl_max:
						print("too many iterations")
						break
					

			# if FixedPoint:
			if False:

				it=0
				a.Assemble()
				inva = a.mat.Inverse(freedofs=VSpace.FreeDofs(coupling=True))
				while True:
					it += 1
					print ("Iteration",it)
					
					if it > 1:
						munonlin.Update()

					f.Assemble()
					
					sol_old.vec.data = sol.vec
					# solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=2, needsassembling=False)
					sol.vec.data = inva * f.vec
					
					sol.vec.data = 0.5*(sol.vec+sol_old.vec)

					if it > 1 and munonlin.HasConverged(err_nl_max,err_nl_av,True):
						break
						
					if it == N_it_nl_max:
						print("too many iterations")
						break
					

			if Newton:

				def SimpleNewtonSolve(gfu,a,tol=tol_Newton,maxits=maxits_Newton):
					# print('sum sol', sum(gfu.vec))
					res = gfu.vec.CreateVector()
					du = gfu.vec.CreateVector()
					fes = gfu.space
					for it in range(maxits):
						print ('')
						print ("Iteration {:3}  ".format(it),end="")
						a.Apply(gfu.vec, res)
						a.AssembleLinearization(gfu.vec)
						du.data = a.mat.Inverse(fes.FreeDofs()) * res
						# gfu.vec.data -= 0.5*du
						gfu.vec.data -= 0.75*du
						# gfu.vec.data -= 0.99*du
						# print('sum du', sum(du))
						# print('sum res', sum(res))
						#stopping criteria
						stopcritval = sqrt(abs(InnerProduct(du,res)))
						print ("<A u",it,", A u",it,">_{-1}^0.5 = ", stopcritval)
						if stopcritval < tol:
							break
				
				# gfu = GridFunction(VSpace)
				# Draw(gfu,mesh,"u")
				SimpleNewtonSolve(sol,a)

			if not Newton:
				muCF=test*munonlin.GetGlobalFunction() + (CoefficientFunction(1)-test)*muAir
			if Newton:
				muCF=test*munonlinNM + (CoefficientFunction(1)-test)*muAir
	
	# sol.Save("SI_nonline_o1.sol")
	B = muCF*H
	J = curl(sol.components[0])
	
	Draw(BS,mesh,'BS')
	Draw(H*CF(1),mesh,'H')
	Draw(B,mesh,'B')
	Draw(J,mesh,'J')
	
	airvals={mat:0 for mat in mesh.GetMaterials()}
	airvals["air"]=1
	airdom = CoefficientFunction([airvals[mat] for mat in mesh.GetMaterials()])
	
	Hred = - grad(sol.components[1])#sol.components[0] - grad(sol.components[1])
	Bred = muAir*Hred*airdom
	Draw(Bred,mesh,'Bred')

	# --------- alt:
	
	print("Htest:",Integrate(H,mesh,BND,definedon=mesh.Boundaries("sTop")))
	
	lossesActiveCoef = 0.5*1/sigma*J*Conj(J)
	lossesReactiveCoef = omega/2*muCF*H*Conj(H)
	
	lossesActive=Integrate(test*lossesActiveCoef,mesh)
	lossesReactive=Integrate(test*lossesReactiveCoef,mesh)
	losses=lossesReactive-1j*lossesActive
	
	print("side")
	temp = Integrate(1,mesh,definedon=mesh.Boundaries("sRight|sBack"))
	print(Integrate(H/temp,mesh,definedon=mesh.Boundaries("sRight|sBack")))
	print("top")
	temp = Integrate(1,mesh,definedon=mesh.Boundaries("sTop"))
	print(Integrate(grad(sol.components[1])/temp,mesh,definedon=mesh.Boundaries("sTop")))
	
	# print("Losses:",lossesActive)
	# print("Energy:",lossesReactive)
	temp = [lossesActive,lossesReactive]
	
	# --------- neu:
	
	losses_specific = 0.5*1/sigma*J*Conj(J)
	reacpower_specific = omega/2*muCF*H*Conj(H)

	Losses=Integrate(losses_specific,mesh,definedon=mesh.Materials("iron"))
	# Losses2=Integrate(losses_specific,mesh,definedon=mesh.Boundaries("sRight"))
	# Losses3=Integrate(losses_specific,mesh,definedon=mesh.Boundaries("sBack"))
	print()
	print("- - - - - - - -")
	print()
	# print("Losses:")
	# print("oben:",Losses1)
	# print("rechts:",Losses2)
	# print("hinten:",Losses3)
	# Losses=Losses1+Losses2+Losses3
	print("Losses:",Losses)

	ReactivePower=Integrate(reacpower_specific,mesh,definedon=mesh.Materials("iron"))
	# ReactivePower2=Integrate(losses_specific,mesh,definedon=mesh.Boundaries("sRight"))
	# ReactivePower3=Integrate(losses_specific,mesh,definedon=mesh.Boundaries("sBack"))
	print()
	# print("Reactive Power:")
	# print("oben:",ReactivePower1)
	# print("rechts:",ReactivePower2)
	# print("hinten:",ReactivePower3)
	# ReactivePower=ReactivePower1+ReactivePower2+ReactivePower3
	print("ReactivePower:",ReactivePower)

	
	Z = 0.4325e-3*(1+1j)
	Pfac = 0.000210
	Qfac = 0.000210
	
	Hvec = [0,1e6]
	Zvec_real = [0,0.4325e-3*1e6]
	Zvec_imag = [0,0.4325e-3*1e6]

	Pvec = [0,Pfac*1e12]
	Qvec = [0,Qfac*1e12]
	
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
	
	Hvec = H0_amp_i# + [1e6]
	Zvec_real = [Z_i[i].real*Hvec[i] for i in range(len(Z_i))]# + [Z_i[-1].real]
	Zvec_imag = [Z_i[i].imag*Hvec[i] for i in range(len(Z_i))]# + [Z_i[-1].imag]
	Pvec = P_eddy_i# +[Pfac*1e12]
	Qvec = Q_i# +[Qfac*1e12]
	
	Pcurve = BSpline (2, [0]+[Hvec[i]**2 for i in range(len(Hvec))], Pvec)
	Qcurve = BSpline (2, [0]+[Hvec[i]**2 for i in range(len(Hvec))], Qvec)
	
	if linear:
		Hvec = [0,1e6]
		Zvec_real = [0,0.4325e-3*1e6]
		Zvec_imag = [0,0.4325e-3*1e6]

		Pvec = [0,Pfac*1e12]
		Qvec = [0,Qfac*1e12]
		
		Pcurve = BSpline (2, [0]+[Hvec[i]**2 for i in range(len(Hvec))], Pvec)
		Qcurve = BSpline (2, [0]+[Hvec[i]**2 for i in range(len(Hvec))], Qvec)
	
	Losses1=Integrate(Pcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sTop"))
	Losses2=Integrate(Pcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sRight"))
	Losses3=Integrate(Pcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sBack"))
	print()
	print("- - - Comparison with effectiv curves for P' and Q' - - -")
	print()
	print("Losses:")
	print("oben:",Losses1)
	print("rechts:",Losses2)
	print("hinten:",Losses3)
	Losses=Losses1+Losses2+Losses3
	print("Losses:",Losses)
	
	ReactivePower1=Integrate(Qcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sTop"))
	ReactivePower2=Integrate(Qcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sRight"))
	ReactivePower3=Integrate(Qcurve(H.Norm()**2),mesh,definedon=mesh.Boundaries("sBack"))
	print()
	print("Reactive power:")
	print("oben:",ReactivePower1)
	print("rechts:",ReactivePower2)
	print("hinten:",ReactivePower3)
	ReactivePower=ReactivePower1+ReactivePower2+ReactivePower3
	print("ReactivePower:",ReactivePower)

	# losses=ReactivePower-1j*lossesActive
		# print("Losses:",Losses)
	# print("Energy:",ReactivePower)
	

	if False:
		print()
		print("- alt - -")
		print()

		print()
		print("Losses:",temp[0])
		print("Energy:",temp[1])

	print('')
	time_end = time.time()
	print('Computation Time = ',time_end-time_start)
	print('')
	