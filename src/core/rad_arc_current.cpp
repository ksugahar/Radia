/*-------------------------------------------------------------------------
*
* File name:      radarccu.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 rectangular cross-section arc with azimuthal current
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
*
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
* References for analytical formulas:
*   [1] Full circular coil B-field: See rad_elliptic_integral.cpp
*       - J.C. Simpson et al., NASA Tech Report, 2001
*       - J.C. Maxwell, "A Treatise on Electricity and Magnetism", 1873
*
*   [2] Arc coil B-field (analytical azimuthal integration):
*       - A. Kameari, "Calculation of Transient 3D Eddy Current using
*         Edge-Elements", IEEE Trans. Magn., Vol. 26, No. 2, 1990
*       - T. Nakata, N. Takahashi, K. Fujiwara, "Summary of Results for
*         TEAM Problem 7", COMPEL, Vol. 9, No. 2, 1990
*
*   [3] Rectangular cross-section integration (1/r integral over rectangle):
*       - A. Kameari, "Transient Eddy Current Analysis on the Two-Dimensional
*         Finite Element Method with Edge Elements", IEEJ Trans., Vol. 110-D,
*         No. 5, 1990 (in Japanese, known as "Kameari formula")
*       - D.R. Wilton, S.M. Rao, A.W. Glisson, "Potential Integrals for
*         Uniform and Linear Source Distributions on Polygonal and Polyhedral
*         Domains", IEEE Trans. AP, Vol. 32, No. 3, March 1984
*
*   [4] Adaptive quadrature for near-field azimuthal integration:
*       - R. Piessens, E. de Doncker-Kapenga, C.W. Ueberhuber, D.K. Kahaner,
*         "QUADPACK: A Subroutine Package for Automatic Integration",
*         Springer-Verlag, 1983 - Gauss-Kronrod 7-15 adaptive quadrature
*
-------------------------------------------------------------------------*/

#include "rad_application.h"
#include "rad_arc_current.h"
#include "rad_subdivided_arc_current.h"
#include "rad_elliptic_integral.h"

#include <math.h>
#include <sstream>

//-------------------------------------------------------------------------
// Analytical 1/r integral over a line segment (Wilton formula)
// Reference: Wilton et al., IEEE Trans. AP, 1984
//
// Computes: integral_{line AB} 1/|r-r'| dl'
// where r = (0, 0, z) and the line goes from A to B in the xy-plane
//-------------------------------------------------------------------------
static double InverseRIntegralLine(const double A[2], const double B[2], double z)
{
	const double EPS = 1.0e-15;

	if (z < 0.0) z = -z;  // Use |z| for symmetry

	double ax = B[0] - A[0];
	double ay = B[1] - A[1];
	double AB = sqrt(ax*ax + ay*ay);

	if (AB < EPS) return 0.0;

	// Project A and B onto the line direction
	double dA = (A[0]*ax + A[1]*ay) / AB;  // Projection of A
	double dB = (B[0]*ax + B[1]*ay) / AB;  // Projection of B
	double d = (A[1]*ax - A[0]*ay) / AB;   // Perpendicular distance to line

	double rA = sqrt(A[0]*A[0] + A[1]*A[1] + z*z);  // Distance to A
	double rB = sqrt(B[0]*B[0] + B[1]*B[1] + z*z);  // Distance to B

	// Logarithmic term: d * log((rB+dB)/(rA+dA)) with numerical safeguards
	double log_term = 0.0;
	double ff;
	if (dA > 0.0) {
		double num = rB + dB;
		double den = rA + dA;
		if (fabs(den) > EPS) {
			log_term = log(fabs(num / den));
		}
	} else if ((ff = rB - dB) <= 0.0) {
		if (dA * dB <= 0.0) {
			log_term = 0.0;
		} else if (fabs(dB) > EPS) {
			log_term = log(fabs(dA / dB));
		}
	} else if (fabs(ff) > EPS) {
		log_term = log(fabs((rA - dA) / ff));
	}

	// Arctan term: z * [atan2(d,dA) - atan2(d,dB) - atan2(rA*d, dA*z) + atan2(rB*d, dB*z)]
	double ang = 0.0;
	if (fabs(d) > EPS) {
		ang = atan2(d, dA) - atan2(d, dB) - atan2(rA*d, dA*z) + atan2(rB*d, dB*z);
	}

	return -d * log_term + z * ang;
}

//-------------------------------------------------------------------------
// Analytical 1/r integral over a rectangle (Kameari/Wilton formula)
// Reference: Kameari (1990), Wilton et al. (1984)
//
// Computes: integral_{rectangle} 1/|r-r'| dS'
// where r = (x, y, z) and rectangle has half-sizes (s1, s2) centered at origin
//-------------------------------------------------------------------------
static double InverseRIntegralRectangle(double s1, double s2, double x, double y, double z)
{
	double A[2], B[2];
	double sum = 0.0;

	// Four edges of the rectangle, traversed in order
	// Edge 1: (s1, -s2) -> (s1, s2)
	A[0] = s1 - x;  A[1] = -s2 - y;
	B[0] = s1 - x;  B[1] = s2 - y;
	sum += InverseRIntegralLine(A, B, z);

	// Edge 2: (s1, s2) -> (-s1, s2)
	A[0] = B[0];  A[1] = B[1];
	B[0] = -s1 - x;  B[1] = s2 - y;
	sum += InverseRIntegralLine(A, B, z);

	// Edge 3: (-s1, s2) -> (-s1, -s2)
	A[0] = B[0];  A[1] = B[1];
	B[0] = -s1 - x;  B[1] = -s2 - y;
	sum += InverseRIntegralLine(A, B, z);

	// Edge 4: (-s1, -s2) -> (s1, -s2)
	A[0] = B[0];  A[1] = B[1];
	B[0] = s1 - x;  B[1] = -s2 - y;
	sum += InverseRIntegralLine(A, B, z);

	return sum;
}

//-------------------------------------------------------------------------
// Analytical method using elliptic integrals for full circular coils
// and Kameari/Wilton analytical integration for arc coils
//-------------------------------------------------------------------------

void radTArcCur::B_compElliptic(radTField* FieldPtr)
{
	const double Pi = 3.141592653589793238;
	const double TwoPi = 2.0 * Pi;
	// Biot-Savart constant: B = (mu_0/4*pi) * integral(I dl x r / r^3)
	// Radia now uses SI units (meters) internally, matching ELF.
	// ConstForJ = mu_0/(4*pi) = 1e-7 H/m
	const double ConstForJ = 1.0e-7;

	TVector3d P_mi_CenPo = FieldPtr->P - CircleCentrPoint;

	const double SmallPositive = 1.E-10;
	double r2_xy = P_mi_CenPo.x*P_mi_CenPo.x + P_mi_CenPo.y*P_mi_CenPo.y;
	double r = sqrt(r2_xy + SmallPositive);  // SmallPositive avoids division by zero in B-field
	double r_exact = sqrt(r2_xy);  // Exact r for solid angle (no offset)
	double phi_obs = ((P_mi_CenPo.y < 0)? (TwoPi - acos(P_mi_CenPo.x/r)) : (acos(P_mi_CenPo.x/r)));
	double z = P_mi_CenPo.z;

	// Check if this is a full circular coil or an arc
	double delta_phi = Phi_max - Phi_min;
	bool is_full_circle = (fabs(delta_phi - TwoPi) < 1.0e-6);

	double IntForAx = 0.0, IntForAy = 0.0;
	double IntForBx = 0.0, IntForBy = 0.0, IntForBz = 0.0;
	double IntForPhi = 0.0;  // Magnetic scalar potential [A]

	if (is_full_circle) {
		// Full circular coil: use direct elliptic integral formulas
		// Integrate over the rectangular cross-section

		// Integration over cross-section using Gaussian quadrature (4x4)
		// 4-point Gauss-Legendre quadrature gives excellent accuracy for thick coils
		// Reference: Abramowitz & Stegun, Table 25.4
		static const int GAUSS_ORDER = 4;
		static const double gp[] = {
			-0.8611363115940526,
			-0.3399810435848563,
			 0.3399810435848563,
			 0.8611363115940526
		};
		static const double gw[] = {
			0.3478548451374538,
			0.6521451548625461,
			0.6521451548625461,
			0.3478548451374538
		};

		double r_mid = 0.5 * (R_max + R_min);
		double r_half = 0.5 * (R_max - R_min);
		double z_half = 0.5 * Height;

		for (int ir = 0; ir < GAUSS_ORDER; ++ir) {
			double r_coil = r_mid + r_half * gp[ir];
			double w_r = gw[ir] * r_half;

			for (int iz = 0; iz < GAUSS_ORDER; ++iz) {
				double z_coil = z_half * gp[iz];  // z relative to coil center
				double w_z = gw[iz] * z_half;

				double w_total = w_r * w_z;

				// Current density J_azim is in A/m^2 (SI units)
				// Integration weights w_r and w_z are in meters
				// dI = J_azim [A/m^2] * w_total [m^2] = J_azim * w_total [A]
				double dI = J_azim * w_total;

				// Compute B-field from this circular loop using elliptic integrals
				double dBR = 0.0, dBZ = 0.0;
				RadElliptic::CircularLoopBField(r, z, r_coil, z_coil, dI, dBR, dBZ);

				// BR and BZ are in cylindrical coordinates
				// Convert to Cartesian at observation point
				double cos_phi = cos(phi_obs);
				double sin_phi = sin(phi_obs);

				if (FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_) {
					IntForBx += dBR * cos_phi;
					IntForBy += dBR * sin_phi;
					IntForBz += dBZ;
				}

				// Vector potential
				if (FieldPtr->FieldKey.A_) {
					double dAphi = RadElliptic::CircularLoopAPhi(r, z, r_coil, z_coil, dI);
					// A = A_phi * (-sin(phi), cos(phi), 0) in Cartesian
					IntForAx += dAphi * (-sin_phi);
					IntForAy += dAphi * cos_phi;
				}

				// Magnetic scalar potential: Phi = I * Omega / (4*pi)
				// Omega is the solid angle subtended by the loop
				if (FieldPtr->FieldKey.Phi_) {
					double dOmega = RadElliptic::CircularLoopSolidAngle(r_exact, z, r_coil, z_coil);
					IntForPhi += dI * dOmega / (4.0 * Pi);
				}
			}
		}
	} else {
		// Arc coil: use Biot-Savart with Gauss quadrature over cross-section
		// This replaces the analytical Kameari formula which had accuracy issues
		// for thick coils with large aspect ratios.

		double r_mid = 0.5 * (R_max + R_min);
		double r_half = 0.5 * (R_max - R_min);  // Half radial width (a)
		double z_half = 0.5 * Height;           // Half axial width (b)

		// B-field calculation using Biot-Savart with Gauss quadrature
		if (FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_) {
			// 4-point Gauss-Legendre quadrature for cross-section integration
			static const int GAUSS_ORDER = 4;
			static const double gp[] = {
				-0.8611363115940526,
				-0.3399810435848563,
				 0.3399810435848563,
				 0.8611363115940526
			};
			static const double gw[] = {
				0.3478548451374538,
				0.6521451548625461,
				0.6521451548625461,
				0.3478548451374538
			};

			// Number of phi segments for arc integration
			// Use at least 4 segments per 90 degrees for accuracy
			int n_phi = NumberOfSectors;
			if (n_phi < 4) n_phi = 4;
			// Scale by arc angle relative to full circle
			n_phi = std::max(4, (int)(n_phi * delta_phi / TwoPi + 0.5));

			double dphi = delta_phi / n_phi;

			// Biot-Savart constant: mu_0 / (4 * pi)
			const double mu0_over_4pi = 1.0e-7;

			// Integrate over cross-section using Gauss quadrature
			for (int ir = 0; ir < GAUSS_ORDER; ++ir) {
				double r_coil = r_mid + r_half * gp[ir];
				double w_r = gw[ir] * r_half;

				for (int iz = 0; iz < GAUSS_ORDER; ++iz) {
					double z_coil = z_half * gp[iz];
					double w_z = gw[iz] * z_half;

					// Current for this cross-section element
					double dI = J_azim * w_r * w_z;

					// Integrate over arc angle using midpoint rule
					for (int iphi = 0; iphi < n_phi; ++iphi) {
						double phi_coil = Phi_min + (iphi + 0.5) * dphi;

						// Position of current element (in coil-centered coordinates)
						double x_src = r_coil * cos(phi_coil);
						double y_src = r_coil * sin(phi_coil);
						double z_src = z_coil;

						// Current direction (tangent to arc)
						// dl = r_coil * dphi * (-sin(phi), cos(phi), 0)
						double dlx = r_coil * dphi * (-sin(phi_coil));
						double dly = r_coil * dphi * cos(phi_coil);
						double dlz = 0.0;

						// Vector from source to observation point
						double Rx = P_mi_CenPo.x - x_src;
						double Ry = P_mi_CenPo.y - y_src;
						double Rz = P_mi_CenPo.z - z_src;

						double R_mag_sq = Rx*Rx + Ry*Ry + Rz*Rz;
						double R_mag = sqrt(R_mag_sq);

						if (R_mag > 1.0e-15) {
							double R_mag_cubed = R_mag_sq * R_mag;

							// dB = (mu_0 / 4pi) * dI * (dl x R) / |R|^3
							// Cross product: dl x R
							double crossX = dly * Rz - dlz * Ry;
							double crossY = dlz * Rx - dlx * Rz;
							double crossZ = dlx * Ry - dly * Rx;

							double factor = mu0_over_4pi * dI / R_mag_cubed;

							IntForBx += factor * crossX;
							IntForBy += factor * crossY;
							IntForBz += factor * crossZ;
						}
					}
				}
			}
		}

		// Vector potential and scalar potential still use numerical integration
		if (FieldPtr->FieldKey.A_ || FieldPtr->FieldKey.Phi_) {
			// Split arc into segments and integrate using analytical 1/r integral for cross-section

			int n_phi = NumberOfSectors;
			if (n_phi < 4) n_phi = 4;

			double dphi = delta_phi / n_phi;

			for (int iphi = 0; iphi < n_phi; ++iphi) {
				double phi_coil = Phi_min + (iphi + 0.5) * dphi;
				double w_phi = dphi;

				double cos_phi_coil = cos(phi_coil);
				double sin_phi_coil = sin(phi_coil);

				// Distance from observation point to coil center at this phi
				double x_coil_center = r_mid * cos_phi_coil;
				double y_coil_center = r_mid * sin_phi_coil;

				double dx = P_mi_CenPo.x - x_coil_center;
				double dy = P_mi_CenPo.y - y_coil_center;

				// Local coordinates relative to coil cross-section center
				double x_local = dx * cos_phi_coil + dy * sin_phi_coil;  // Radial component
				double y_local = -dx * sin_phi_coil + dy * cos_phi_coil; // Tangential component

				// Use analytical 1/r integral over rectangle for vector potential
				if (FieldPtr->FieldKey.A_) {
					double I_rect = InverseRIntegralRectangle(r_half, z_half, x_local, y_local, z);
					double dl = r_mid * w_phi;
					double factor_A = ConstForJ * J_azim * dl * I_rect;

					IntForAx += factor_A * (-sin_phi_coil);
					IntForAy += factor_A * cos_phi_coil;
				}

				// Scalar potential (numerical integration)
				if (FieldPtr->FieldKey.Phi_) {
					static const double gp[] = {-0.5773502691896257, 0.5773502691896257};
					static const double gw[] = {1.0, 1.0};

					for (int ir = 0; ir < 2; ++ir) {
						double r_coil = r_mid + r_half * gp[ir];
						double w_r = gw[ir] * r_half;

						for (int iz = 0; iz < 2; ++iz) {
							double z_coil = z_half * gp[iz];
							double w_z = gw[iz] * z_half;

							double dI = J_azim * w_r * w_z;
							double dl = r_coil * w_phi;

							double rx = P_mi_CenPo.x - r_coil * cos_phi_coil;
							double ry = P_mi_CenPo.y - r_coil * sin_phi_coil;
							double rz = z - z_coil;
							double dist2 = rx*rx + ry*ry + rz*rz;
							double dist = sqrt(dist2 + SmallPositive);

							double jx = -sin_phi_coil;
							double jy = cos_phi_coil;
							double cross_z = jx * ry - jy * rx;

							double r_perp = sqrt(rx*rx + ry*ry + SmallPositive);
							double dOmega = cross_z * dl / (dist * dist * r_perp);
							IntForPhi += dI * dOmega / (4.0 * Pi);
						}
					}
				}
			}
		}
	}

	// Apply results
	// Note: CircularLoopBField outputs B in Tesla.
	// However, Radia stores B internally as "equivalent magnetization" (B/mu_0 in A/m)
	// because OutFieldCompRes multiplies by mu_0 when returning B values.
	// So we divide by mu_0 here to store B/mu_0, then OutFieldCompRes will multiply
	// by mu_0 to give the correct B value in Tesla.
	const double InvMu0 = 1.0 / (4.0 * 3.141592653589793238 * 1.0e-7);  // 1/mu_0 = 795774.715 (A/m)/T

	if (FieldPtr->FieldKey.A_) {
		TVector3d BufA(IntForAx, IntForAy, 0.0);
		FieldPtr->A += BufA;
	}
	if (FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_) {
		// Convert from Tesla to A/m (internal format)
		TVector3d BufB(IntForBx * InvMu0, IntForBy * InvMu0, IntForBz * InvMu0);
		FieldPtr->B += BufB;
		FieldPtr->H += BufB;  // For coils in air, H = B/mu_0 = same as internal representation
	}
	if (FieldPtr->FieldKey.Phi_) {
		// Magnetic scalar potential in Amperes
		// Only computed for full circular coils; arcs return 0
		FieldPtr->Phi += IntForPhi;
	}
}

//-------------------------------------------------------------------------

void radTArcCur::Dump(std::ostream& o, int ShortSign)
{
	radTg3d::Dump(o);
	DumpPureObjInfo(o, ShortSign);
	if(ShortSign==1) return;

	DumpTransApplied(o);

	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

void radTArcCur::DumpPureObjInfo(std::ostream& o, int ShortSign)
{
	o << "Current carrying: ";
	o << "ArcCur";

	if(ShortSign==1) return;

	o << endl;
	o << "   {x,y,z}= {" << CircleCentrPoint.x << ','
						 << CircleCentrPoint.y << ','
						 << CircleCentrPoint.z << "}" << endl
	  << "   {rmin,rmax}= {" << R_min << ',' << R_max << "}" << endl
	  << "   {phimin,phimax}= {" << Phi_min << ',' << Phi_max << "}" << endl
	  << "   h= " << Height << endl
	  << "   nseg= " << NumberOfSectors << endl
	  << "   j= " << J_azim << endl
	  << "   Field computation mode: " << (BasedOnPrecLevel? "\"auto\"" : "\"man\"");
}

//-------------------------------------------------------------------------

void radTArcCur::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	//Dumping objects that may be used by this object
	vector<pair<int, int> > vTrfKeys;
	DumpBin_g3d_TreatTrfs(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vTrfKeys);

	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)Type_g();
	oStr << (char)Type_g3d();
	oStr << (char)0;
	oStr << (char)0;
	oStr << (char)0;

	//Members of radTg3d
	DumpBin_g3d(oStr, vTrfKeys);

	//Members of radTArcCur
	//TVector3d CircleCentrPoint;
	oStr << CircleCentrPoint;

	//double R_min, R_max;
	oStr << R_min << R_max;

	//double Phi_min, Phi_max;
	oStr << Phi_min << Phi_max;

	//double Height;
	oStr << Height;

	//double J_azim;
	oStr << J_azim;

	//int NumberOfSectors;
	oStr << NumberOfSectors;

	//short BasedOnPrecLevel;
	oStr << BasedOnPrecLevel;

	//short InternalFacesAfterCut;
	oStr << InternalFacesAfterCut;

	//char J_IsNotZero;
	oStr << J_IsNotZero;
}

//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

int radTArcCur::SubdivideItself(double* SubdivArray, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	char SubdivideCoils = pSubdivOptions->SubdivideCoils;
	char PutNewStuffIntoGenCont = pSubdivOptions->PutNewStuffIntoGenCont;

	if(!SubdivideCoils) return 1;
	radTSend Send;
	if((pSubdivOptions->SubdivisionFrame != 0) && (!g3dListOfTransform.empty())) 
	{
		Send.ErrorMessage("Radia::Error108"); return 0;
	}

	const double ZeroTol = 1.E-10;

	double kr = SubdivArray[0], kPhi = SubdivArray[2], kz = SubdivArray[4];
	double qr = SubdivArray[1], qPhi = SubdivArray[3], qz = SubdivArray[5];

	double DelPhi = Phi_max - Phi_min;
	double Del_r = R_max - R_min;

	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		double DelPhiL = DelPhi*0.666666666667*(R_max*R_max + R_max*R_min + R_min*R_min)/(R_max + R_min);
		kPhi = (kPhi < DelPhiL)? Round(DelPhiL/kPhi) : 1.;

		kr = (kr < Del_r)? Round(Del_r/kr) : 1.;
		kz = (kz < Height)? Round(Height/kz) : 1.;
	}

	if((fabs(kPhi-1.)<ZeroTol) && (fabs(kr-1.)<ZeroTol) && (fabs(kz-1.)<ZeroTol)) return 1;

	radTGroup* GroupInPlaceOfThisPtr = new radTSubdividedArcCur(this);
	radThg NewHandle(GroupInPlaceOfThisPtr);

	const double AbsZeroTol = 5.E-12;
	double q0Phi = (fabs(kPhi-1.)>AbsZeroTol)? pow(qPhi, 1./(kPhi-1.)) : qPhi;
	double q0r = (fabs(kr-1.)>AbsZeroTol)? pow(qr, 1./(kr-1.)) : qr;
	double q0z = (fabs(kz-1.)>AbsZeroTol)? pow(qz, 1./(kz-1.)) : qz;
	double BufPhi = qPhi*q0Phi - 1., BufR = qr*q0r - 1., BufZ = qz*q0z - 1.;

	double a1Phi = (fabs(BufPhi) > AbsZeroTol)? DelPhi*(q0Phi - 1.)/BufPhi : DelPhi/kPhi;
	double a1r = (fabs(BufR) > AbsZeroTol)? Del_r*(q0r - 1.)/BufR : Del_r/kr;
	double a1z = (fabs(BufZ) > AbsZeroTol)? Height*(q0z - 1.)/BufZ : Height/kz;

	TVector3d InitNewDims(a1Phi, a1r, a1z);
	TVector3d NewDims = InitNewDims;

	short NewFacesState[6], ParentFacesState[6];
	ListFacesInternalAfterCut(ParentFacesState);

	int kPhiInt = int(kPhi), krInt = int(kr), kzInt = int(kz);
	int kPhi_mi_1 = kPhiInt-1, kr_mi_1 = krInt-1, kz_mi_1 = kzInt-1;

	TVector3d InitNewCircleCenPoi = TVector3d(CircleCentrPoint.x, CircleCentrPoint.y, CircleCentrPoint.z - 0.5*(Height - a1z));
	TVector3d NewCircleCenPoi = InitNewCircleCenPoi;

	double NewAngles[2], NewRadii[2], NewHeight = a1z;

	double &NewStartAngle = *NewAngles, &NewFinAngle = *(NewAngles+1), &NewStartRad = *NewRadii, &NewFinRad = *(NewRadii+1);
	NewStartAngle = Phi_min;
	NewFinAngle = Phi_min + a1Phi;
	double SmallDelPhi = a1Phi;

	NewStartRad = R_min;
	NewFinRad = R_min + a1r;
	double SmallDel_r = a1r;

	int NewStuffCounter = 0;
	for(int iPhi=0; iPhi<kPhiInt; iPhi++)
	{
		NewFacesState[0] = NewFacesState[1] = 1;
		if(iPhi==0) NewFacesState[0] = ParentFacesState[0];
		if(iPhi==kPhi_mi_1) NewFacesState[1] = ParentFacesState[1];

		float FloatNumOfSect = (float)((SmallDelPhi/DelPhi)*NumberOfSectors);
		int IntNumOfSect = int(FloatNumOfSect);
		int NewNumberOfSectors = ((FloatNumOfSect - IntNumOfSect) < 0.5)? IntNumOfSect : IntNumOfSect + 1;

		if(NewNumberOfSectors < 1) NewNumberOfSectors = 1;

		for(int ir=0; ir<krInt; ir++)
		{
			NewFacesState[2] = NewFacesState[3] = 1;
			if(ir==0) NewFacesState[2] = ParentFacesState[2];
			if(ir==kr_mi_1) NewFacesState[3] = ParentFacesState[3];

			for(int iz=0; iz<kzInt; iz++)
			{
				NewFacesState[4] = NewFacesState[5] = 1;
				if(iz==0) NewFacesState[4] = ParentFacesState[4];
				if(iz==kz_mi_1) NewFacesState[5] = ParentFacesState[5];

				radTArcCur* ArcCurPtr = new radTArcCur(NewCircleCenPoi, NewRadii, NewAngles, NewHeight, J_azim, NewNumberOfSectors, BasedOnPrecLevel);
				if(ArcCurPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

				radThg hg(ArcCurPtr);
				if(PutNewStuffIntoGenCont) GroupInPlaceOfThisPtr->AddElement(radPtr->AddElementToContainer(hg), hg);
				else GroupInPlaceOfThisPtr->AddElement(++NewStuffCounter, hg);

				ArcCurPtr->SetFacesInternalAfterCut(NewFacesState);

				NewCircleCenPoi.z += 0.5*NewHeight;
				NewHeight *= q0z;
				NewCircleCenPoi.z += 0.5*NewHeight;
			}
			NewStartRad = NewFinRad;
			SmallDel_r *= q0r;
			NewFinRad += SmallDel_r;

			NewCircleCenPoi.z = InitNewCircleCenPoi.z;
			NewHeight = a1z;
		}
		NewStartAngle = NewFinAngle;
		SmallDelPhi *= q0Phi;
		NewFinAngle += SmallDelPhi;

		NewStartRad = R_min;
		NewFinRad = R_min + a1r;
		SmallDel_r = a1r;
	}
	In_hg = NewHandle;
	return 1;
}

//-------------------------------------------------------------------------

void radTArcCur::Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl)
{// Attention: this assumes no more than one transformation with mult. no more than 1 !!!
	if(pFieldKey->M_) return;
	else
	{
		TVector3d CP = CentrPoint;
		radTrans* pTrans = (g3dListOfTransform.empty())? 0 : (radTrans*)((*(g3dListOfTransform.begin())).Handler_g.rep);

		radTrans TotTrans;
		if(pTrans != 0)
		{
			if(pBaseTrans != 0) 
			{
				TrProduct(pBaseTrans, pTrans, TotTrans);
				pTrans = &TotTrans;
			}
		}
		else
		{
			if(pBaseTrans != 0) pTrans = pBaseTrans;
		}

		if(pTrans != 0) CP = pTrans->TrPoint(CP);
		radTPairOfVect3d Pair(CP);
		if(pFieldKey->J_)
		{
			double Phic = 0.5*(Phi_min + Phi_max);
			TVector3d J(-J_azim*sin(Phic), J_azim*cos(Phic), 0.);
			if(pTrans != 0) J = pTrans->TrVectField(J);
			Pair.V2 = J;
		}
		else
		{
			radTCompCriterium CompCriterium;
			TVector3d ZeroVect(0.,0.,0.);
			radTField Field(*pFieldKey, CompCriterium, CP, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
			g3dSrcPtr->B_genComp(&Field);
			Pair.V2 = (pFieldKey->M_)? Field.M : ((pFieldKey->B_)? Field.B : ((pFieldKey->H_)? Field.H : ((pFieldKey->A_)? Field.A : ZeroVect)));
		}
		pVectPairOfVect3d->push_back(Pair);
	}
}

//-------------------------------------------------------------------------

void radTArcCur::VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique) 
{//adds vertices of sectors 
	
	double HalfHeight = 0.5*Height;
	double dPhi = 0.;
	if(NumberOfSectors > 0) dPhi = (Phi_max - Phi_min)/NumberOfSectors;
	double CurPhi = Phi_min;
	
	for(int i=0; i<=NumberOfSectors; i++)
	{
		double CosPhi = cos(CurPhi), SinPhi = sin(CurPhi);
		double RelInnerX = R_min*CosPhi, RelOuterX = R_max*CosPhi;
		double RelInnerY = R_min*SinPhi, RelOuterY = R_max*SinPhi;

		double x = CircleCentrPoint.x + RelInnerX;
		double y = CircleCentrPoint.y + RelInnerY;
		double z = CircleCentrPoint.z - HalfHeight;
		TVector3d LowerInner(x, y, z);

		TVector3d UpperInner = LowerInner;
		UpperInner.z += Height;

		x = CircleCentrPoint.x + RelOuterX;
		y = CircleCentrPoint.y + RelOuterY;
		TVector3d LowerOuter(x, y, z);

		TVector3d UpperOuter = LowerOuter;
		UpperOuter.z += Height;

		OutVect.push_back(LowerInner);
		OutVect.push_back(UpperInner);
		OutVect.push_back(LowerOuter);
		OutVect.push_back(UpperOuter);

		CurPhi += dPhi;
	}
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
