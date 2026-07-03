//-------------------------------------------------------------------------
// rad_moment2d.cpp -- 2D planar collocation MMMM (see rad_moment2d.h).
//
// Faithful 2D port of the 3D moment kernel (rad_interaction.cpp).  Validated
// against an independent numpy PoC: disk demag 1/2, ellipse 2:1 -> 1/3,2/3,
// linear chi-sweep chi/(1+chi/2).
//-------------------------------------------------------------------------
#include "rad_moment2d.h"
#include "rad_parallel.h"

#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>

#ifdef HAVE_LAPACK
#include "mkl_lapack.h"
#endif

namespace rad_moment2d {

static const double TWO_PI = 6.283185307179586476925286766559;

//--- field H at P (Px,Py) from a UNIT uniform line charge on a segment whose
//    frame is P0 + local coords (tangent that[2], outward normal nn[2], length L).
static inline void SegField(const double* P0, const double* that, const double* nn,
                            double L, double Px, double Py, double& Hx, double& Hy)
{
	double dx = Px - P0[0], dy = Py - P0[1];
	double x = dx * that[0] + dy * that[1];
	double y = dx * nn[0] + dy * nn[1];
	double r1 = std::sqrt(x * x + y * y);
	double r2 = std::sqrt((x - L) * (x - L) + y * y);
	double Et = std::log(r1 / r2) / TWO_PI;
	// E_n = (1/2pi)[atan((L-x)/y) + atan(x/y)] -- principal atan (NOT atan2): the
	// eval point is an INTERIOR centroid (y<0 vs the outward normal); atan2 would
	// jump by +-pi there.  y != 0 always (centroid never lies on an edge).
	double En = (std::atan((L - x) / y) + std::atan(x / y)) / TWO_PI;
	Hx = Et * that[0] + En * nn[0];
	Hy = Et * that[1] + En * nn[1];
}

//--- field Hessian G[2][2] (dH_i/dx_j, symmetric traceless) at P from the same
//    unit segment charge.  Local: Guu=(1/2pi)[u/r1^2-(u-L)/r2^2],
//    Guv=(1/2pi) v (1/r1^2-1/r2^2), Gvv=-Guu; global G = Q Gloc Q^T, Q=[that|nn].
static inline void SegGrad(const double* P0, const double* that, const double* nn,
                           double L, double Px, double Py, double G[2][2])
{
	double dx = Px - P0[0], dy = Py - P0[1];
	double x = dx * that[0] + dy * that[1];
	double y = dx * nn[0] + dy * nn[1];
	double r1s = x * x + y * y;
	double r2s = (x - L) * (x - L) + y * y;
	double Guu = (x / r1s - (x - L) / r2s) / TWO_PI;
	double Guv = y * (1.0 / r1s - 1.0 / r2s) / TWO_PI;
	double Gvv = -Guu;
	// Q = [[that0, nn0],[that1, nn1]] ; QG = Q*Gloc ; G = QG*Q^T
	double QG00 = that[0] * Guu + nn[0] * Guv;
	double QG01 = that[0] * Guv + nn[0] * Gvv;
	double QG10 = that[1] * Guu + nn[1] * Guv;
	double QG11 = that[1] * Guv + nn[1] * Gvv;
	G[0][0] = QG00 * that[0] + QG01 * nn[0];
	G[0][1] = QG00 * that[1] + QG01 * nn[1];
	G[1][0] = QG10 * that[0] + QG11 * nn[0];
	G[1][1] = QG10 * that[1] + QG11 * nn[1];
}

//--- (nE-3) residual quad eigenmodes: unit R^nE vectors orthogonal to the
//    monopole+dipole functionals {L, L*dx, L*dy}.  Gram-Schmidt on the standard
//    basis complement (the 2D twin of momentResidualEigenmodes).
static void ResidualModes(int nE, const double* Lf, const double d[][2],
                          std::vector<std::vector<double> >& modes)
{
	// functional vectors g0=L, g1=L*dx, g2=L*dy (columns in R^nE)
	std::vector<std::vector<double> > used;               // orthonormal basis of span{g0,g1,g2}
	double gcols[3][8];
	for(int f = 0; f < nE; f++){ gcols[0][f] = Lf[f]; gcols[1][f] = Lf[f] * d[f][0]; gcols[2][f] = Lf[f] * d[f][1]; }
	for(int c = 0; c < 3; c++){
		std::vector<double> v(gcols[c], gcols[c] + nE);
		for(size_t u = 0; u < used.size(); u++){
			double dot = 0; for(int f = 0; f < nE; f++) dot += v[f] * used[u][f];
			for(int f = 0; f < nE; f++) v[f] -= dot * used[u][f];
		}
		double nrm = 0; for(int f = 0; f < nE; f++) nrm += v[f] * v[f]; nrm = std::sqrt(nrm);
		if(nrm > 1e-10){ for(int f = 0; f < nE; f++) v[f] /= nrm; used.push_back(v); }
	}
	// residual modes = standard basis vectors orthogonalized against used + prior modes
	for(int e = 0; e < nE && (int)modes.size() < nE - 3; e++){
		std::vector<double> v(nE, 0.0); v[e] = 1.0;
		for(size_t u = 0; u < used.size(); u++){
			double dot = 0; for(int f = 0; f < nE; f++) dot += v[f] * used[u][f];
			for(int f = 0; f < nE; f++) v[f] -= dot * used[u][f];
		}
		for(size_t u = 0; u < modes.size(); u++){
			double dot = 0; for(int f = 0; f < nE; f++) dot += v[f] * modes[u][f];
			for(int f = 0; f < nE; f++) v[f] -= dot * modes[u][f];
		}
		double nrm = 0; for(int f = 0; f < nE; f++) nrm += v[f] * v[f]; nrm = std::sqrt(nrm);
		if(nrm > 1e-10){ for(int f = 0; f < nE; f++) v[f] /= nrm; modes.push_back(v); }
	}
}

//--- dense LU solve A x = b (A row-major n*n, b length n overwritten with x).
static int DenseSolve(std::vector<double>& A, std::vector<double>& b, int n)
{
#ifdef HAVE_LAPACK
	std::vector<double> Acol(n * n);
	for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) Acol[j * n + i] = A[i * n + j];
	std::vector<int> ipiv(n); int nrhs = 1, info = 0;
	{
		ngcore::SuspendTaskManager stm;
		radia::MKLThreadGuard mkl_guard(radia::GetNumThreads());
		dgesv_(&n, &nrhs, Acol.data(), &n, ipiv.data(), b.data(), &n, &info);
	}
	return (info == 0) ? 0 : -1;
#else
	// Gaussian elimination with partial pivoting (fallback)
	for(int k = 0; k < n - 1; k++){
		int piv = k; double best = std::fabs(A[k * n + k]);
		for(int i = k + 1; i < n; i++){ double v = std::fabs(A[i * n + k]); if(v > best){ best = v; piv = i; } }
		if(best < 1e-300) return -1;
		if(piv != k){ for(int j = 0; j < n; j++) std::swap(A[k * n + j], A[piv * n + j]); std::swap(b[k], b[piv]); }
		for(int i = k + 1; i < n; i++){
			double m = A[i * n + k] / A[k * n + k];
			for(int j = k; j < n; j++) A[i * n + j] -= m * A[k * n + j];
			b[i] -= m * b[k];
		}
	}
	for(int i = n - 1; i >= 0; i--){
		double s = b[i]; for(int j = i + 1; j < n; j++) s -= A[i * n + j] * b[j];
		b[i] = s / A[i * n + i];
	}
	return 0;
#endif
}

int SolveLinear(int nElem, const int* voff, const double* vxy,
                const double* chi, const double* Hext, double* Mout)
{
	if(nElem <= 0) return -1;

	// ---- per-element geometry (CCW-forced) + flat global-edge table ----
	std::vector<double> cx(nElem), cy(nElem), area(nElem);
	std::vector<int> nEd(nElem), dofOff(nElem + 1, 0);
	// per-element vertex coords (CCW), stored packed for the 2nd-moment tensor
	std::vector<std::vector<double> > VX(nElem), VY(nElem);
	// global edge geometry
	std::vector<double> eP0x, eP0y, eTx, eTy, eNx, eNy, eL, eMx, eMy;

	for(int k = 0; k < nElem; k++){
		int nv = voff[k + 1] - voff[k];
		if(nv < 3 || nv > 8){ return -2; }   // tri/quad (fixed local caps sized for <=8 edges)
		std::vector<double> X(nv), Y(nv);
		for(int i = 0; i < nv; i++){ X[i] = vxy[2 * (voff[k] + i)]; Y[i] = vxy[2 * (voff[k] + i) + 1]; }
		// signed 2*area (shoelace); reverse to CCW if negative
		double A2 = 0; for(int i = 0; i < nv; i++){ int j = (i + 1) % nv; A2 += X[i] * Y[j] - X[j] * Y[i]; }
		if(A2 < 0){ std::reverse(X.begin(), X.end()); std::reverse(Y.begin(), Y.end()); A2 = -A2; }
		double ar = 0.5 * A2;
		// area centroid
		double ccx = 0, ccy = 0;
		for(int i = 0; i < nv; i++){ int j = (i + 1) % nv; double cr = X[i] * Y[j] - X[j] * Y[i]; ccx += (X[i] + X[j]) * cr; ccy += (Y[i] + Y[j]) * cr; }
		ccx /= (3.0 * A2); ccy /= (3.0 * A2);
		cx[k] = ccx; cy[k] = ccy; area[k] = ar; nEd[k] = nv; dofOff[k + 1] = dofOff[k] + nv;
		VX[k] = X; VY[k] = Y;
		for(int i = 0; i < nv; i++){
			int j = (i + 1) % nv;
			double tx = X[j] - X[i], ty = Y[j] - Y[i];
			double L = std::sqrt(tx * tx + ty * ty);
			double thx = tx / L, thy = ty / L;
			double mx = 0.5 * (X[i] + X[j]), my = 0.5 * (Y[i] + Y[j]);
			double nx = thy, ny = -thx;                       // CCW outward candidate
			if(nx * (mx - ccx) + ny * (my - ccy) < 0){ nx = -nx; ny = -ny; }
			eP0x.push_back(X[i]); eP0y.push_back(Y[i]); eTx.push_back(thx); eTy.push_back(thy);
			eNx.push_back(nx); eNy.push_back(ny); eL.push_back(L); eMx.push_back(mx); eMy.push_back(my);
		}
	}
	int nDOF = dofOff[nElem];

	// ---- assemble dense system A (row-major) + rhs ----
	std::vector<double> A((size_t)nDOF * nDOF, 0.0), rhs(nDOF, 0.0);
	std::vector<double> Fx(nDOF), Fy(nDOF);                     // field at ck from each edge
	std::vector<double> G00(nDOF), G01(nDOF), G11(nDOF);        // Hessian (sym) at ck from each edge
	int row = 0;
	for(int k = 0; k < nElem; k++){
		double ck[2] = { cx[k], cy[k] }; double Ak = area[k]; int nE = nEd[k]; int o = dofOff[k]; double chik = chi[k];
		// field + Hessian at ck from every global edge
		for(int g = 0; g < nDOF; g++){
			double P0[2] = { eP0x[g], eP0y[g] }, th[2] = { eTx[g], eTy[g] }, nn[2] = { eNx[g], eNy[g] };
			double hx, hy; SegField(P0, th, nn, eL[g], ck[0], ck[1], hx, hy); Fx[g] = hx; Fy[g] = hy;
			double GG[2][2]; SegGrad(P0, th, nn, eL[g], ck[0], ck[1], GG);
			G00[g] = GG[0][0]; G01[g] = GG[0][1]; G11[g] = GG[1][1];
		}
		// ---- monopole row: sum_f sigma_f L_f = 0 ----
		for(int f = 0; f < nE; f++) A[(size_t)row * nDOF + (o + f)] = eL[o + f];
		row++;
		// ---- dipole rows (comp = x,y): (1/A) sum L_f (mid_f-c) - chi sum_g F = chi Hext ----
		for(int comp = 0; comp < 2; comp++){
			for(int f = 0; f < nE; f++){
				double dcomp = (comp == 0 ? eMx[o + f] - ck[0] : eMy[o + f] - ck[1]);
				A[(size_t)row * nDOF + (o + f)] += eL[o + f] * dcomp;
			}
			const std::vector<double>& Fc = (comp == 0 ? Fx : Fy);
			for(int g = 0; g < nDOF; g++) A[(size_t)row * nDOF + g] += -chik * Ak * Fc[g];
			rhs[row] = chik * Ak * Hext[2 * k + comp];
			row++;
		}
		// ---- quad rows (nE-3) ----
		if(nE > 3){
			// per-edge traceless 2nd moment mtil (store mxx, mxy; myy=-mxx)
			std::vector<double> mxx(nE), mxy(nE); double dloc[8][2];
			for(int f = 0; f < nE; f++){
				double dxm = eMx[o + f] - ck[0], dym = eMy[o + f] - ck[1];
				double tx = eTx[o + f], ty = eTy[o + f], L = eL[o + f];
				double Mxx = L * (dxm * dxm + (L * L / 12.0) * tx * tx);
				double Myy = L * (dym * dym + (L * L / 12.0) * ty * ty);
				double Mxy = L * (dxm * dym + (L * L / 12.0) * tx * ty);
				mxx[f] = 0.5 * (Mxx - Myy); mxy[f] = Mxy;
				dloc[f][0] = dxm; dloc[f][1] = dym;
			}
			// cell geometric 2nd-moment tensor I (about centroid) via triangle fan
			double Ixx = 0, Ixy = 0, Iyy = 0;
			for(int i = 0; i < nE; i++){
				int j = (i + 1) % nE;
				double p0x = VX[k][i] - ck[0], p0y = VY[k][i] - ck[1];
				double p1x = VX[k][j] - ck[0], p1y = VY[k][j] - ck[1];
				double cr = p0x * p1y - p0y * p1x, aTri = 0.5 * cr;
				Ixx += (aTri / 12.0) * ((p0x + p1x) * (p0x + p1x) + p0x * p0x + p1x * p1x);
				Iyy += (aTri / 12.0) * ((p0y + p1y) * (p0y + p1y) + p0y * p0y + p1y * p1y);
				Ixy += (aTri / 12.0) * ((p0x + p1x) * (p0y + p1y) + p0x * p0y + p1x * p1y);
			}
			std::vector<std::vector<double> > modes; ResidualModes(nE, &eL[o], dloc, modes);
			for(size_t mm = 0; mm < modes.size(); mm++){
				const std::vector<double>& phi = modes[mm];
				double Txx = 0, Txy = 0;
				for(int f = 0; f < nE; f++){ Txx += phi[f] * mxx[f]; Txy += phi[f] * mxy[f]; }
				// (Txx already traceless-xx; Tyy=-Txx)
				// local geometry coeff on own edges: Tq:mtil_f = 2 Txx mxx_f + 2 Txy mxy_f
				for(int f = 0; f < nE; f++)
					A[(size_t)row * nDOF + (o + f)] += 2.0 * Txx * mxx[f] + 2.0 * Txy * mxy[f];
				// field-gradient coupling: -chi sum_g Tq:(G_g I + I G_g)
				for(int g = 0; g < nDOF; g++){
					double g00 = G00[g], g01 = G01[g], g11 = G11[g];   // Hessian (sym, traceless: g11=-g00)
					// M2 = G I + I G  (both sym) -> sym
					double M2xx = 2.0 * (g00 * Ixx + g01 * Ixy);
					double M2yy = 2.0 * (g01 * Ixy + g11 * Iyy);
					double M2xy = g00 * Ixy + g01 * Iyy + g01 * Ixx + g11 * Ixy;
					double contr = Txx * (M2xx - M2yy) + 2.0 * Txy * M2xy;
					A[(size_t)row * nDOF + g] += -chik * contr;
				}
				rhs[row] = 0.0;
				row++;
			}
		}
	}

	// ---- solve ----
	std::vector<double> sol = rhs;
	if(DenseSolve(A, sol, nDOF) != 0) return -3;

	// ---- recover per-element M = (1/A) sum_f sigma_f L_f (mid_f - c) ----
	for(int k = 0; k < nElem; k++){
		int nE = nEd[k], o = dofOff[k]; double ck[2] = { cx[k], cy[k] }; double mx = 0, my = 0;
		for(int f = 0; f < nE; f++){
			double s = sol[o + f];
			mx += s * eL[o + f] * (eMx[o + f] - ck[0]);
			my += s * eL[o + f] * (eMy[o + f] - ck[1]);
		}
		Mout[2 * k] = mx / area[k]; Mout[2 * k + 1] = my / area[k];
	}
	return 0;
}

} // namespace rad_moment2d
