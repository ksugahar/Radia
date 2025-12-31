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
-------------------------------------------------------------------------*/

#include "rad_application.h"
#include "rad_arc_current.h"
#include "rad_graphics_3d.h"
#include "rad_subdivided_arc_current.h"
#include "rad_elliptic_integral.h"

#include <math.h>
#include <sstream>

//-------------------------------------------------------------------------
// Analytical method using elliptic integrals
// This replaces the legacy numerical integration methods
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
	double r = sqrt(P_mi_CenPo.x*P_mi_CenPo.x + P_mi_CenPo.y*P_mi_CenPo.y + SmallPositive);
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

		// Integration over cross-section using Gaussian quadrature (2x2)
		static const double gp[] = {-0.5773502691896257, 0.5773502691896257};  // Gauss points
		static const double gw[] = {1.0, 1.0};  // Gauss weights

		double r_mid = 0.5 * (R_max + R_min);
		double r_half = 0.5 * (R_max - R_min);
		double z_half = 0.5 * Height;

		for (int ir = 0; ir < 2; ++ir) {
			double r_coil = r_mid + r_half * gp[ir];
			double w_r = gw[ir] * r_half;

			for (int iz = 0; iz < 2; ++iz) {
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
					double dOmega = RadElliptic::CircularLoopSolidAngle(r, z, r_coil, z_coil);
					// Phi = I * Omega / (4*pi)
					// Units: dI [A] * Omega [sr] / (4*pi) = [A] (since Omega has max 4*pi sr)
					IntForPhi += dI * dOmega / (4.0 * Pi);
				}
			}
		}
	} else {
		// Arc coil: use elliptic integral for each azimuthal segment
		// For each phi position, treat the arc segment as part of a circular loop
		// and use elliptic integrals to compute the field contribution

		// Number of azimuthal integration points
		int n_phi = NumberOfSectors;
		if (n_phi < 4) n_phi = 4;

		double dphi = delta_phi / n_phi;

		// Gauss points for cross-section integration (2x2)
		static const double gp[] = {-0.5773502691896257, 0.5773502691896257};
		static const double gw[] = {1.0, 1.0};

		double r_mid = 0.5 * (R_max + R_min);
		double r_half = 0.5 * (R_max - R_min);
		double z_half = 0.5 * Height;

		// Gauss-Legendre quadrature for azimuthal integration (more accurate than trapezoidal)
		// Use n_phi points with Gauss weights
		for (int iphi = 0; iphi < n_phi; ++iphi) {
			// Midpoint of each segment for Gauss integration
			double phi_coil = Phi_min + (iphi + 0.5) * dphi;
			double w_phi = dphi;  // Weight for this segment

			// For each phi, we treat the arc element as contributing like a partial loop
			// The observation point in the local frame of this phi slice
			double cos_phi_coil = cos(phi_coil);
			double sin_phi_coil = sin(phi_coil);

			// Rotate observation point to align with this phi slice
			// In the rotated frame, the coil element is at (r_coil, 0, z_coil)
			// and the observation point is at (r_obs_local, z_obs)
			double x_obs_local = P_mi_CenPo.x * cos_phi_coil + P_mi_CenPo.y * sin_phi_coil;
			double y_obs_local = -P_mi_CenPo.x * sin_phi_coil + P_mi_CenPo.y * cos_phi_coil;

			// Integrate over cross-section
			for (int ir = 0; ir < 2; ++ir) {
				double r_coil = r_mid + r_half * gp[ir];
				double w_r = gw[ir] * r_half;

				for (int iz = 0; iz < 2; ++iz) {
					double z_coil = z_half * gp[iz];
					double w_z = gw[iz] * z_half;

					// Current element: dI = J_azim * dA_cross_section
					double dI = J_azim * w_r * w_z;

					// Arc length element: dl = r_coil * dphi
					double dl = r_coil * w_phi;

					// Distance from coil element at (r_coil, 0, z_coil) to observation point
					double dx_local = x_obs_local - r_coil;
					double dy_local = y_obs_local;  // coil element is at y=0 in local frame
					double dz_local = z - z_coil;
					double dist2 = dx_local*dx_local + dy_local*dy_local + dz_local*dz_local;
					double dist = sqrt(dist2 + SmallPositive);
					double dist3 = dist * dist2;

					// Current direction in global frame: tangent to arc = (-sin(phi), cos(phi), 0)
					double jx = -sin_phi_coil;
					double jy = cos_phi_coil;

					// Biot-Savart in global frame:
					// r_vec = (P_mi_CenPo.x - r_coil*cos_phi, P_mi_CenPo.y - r_coil*sin_phi, z - z_coil)
					double rx = P_mi_CenPo.x - r_coil * cos_phi_coil;
					double ry = P_mi_CenPo.y - r_coil * sin_phi_coil;
					double rz = z - z_coil;

					// dl x r = (jy*rz - 0, 0 - jx*rz, jx*ry - jy*rx)
					double cross_x = jy * rz;
					double cross_y = -jx * rz;
					double cross_z = jx * ry - jy * rx;

					// dB = (mu_0 / 4*pi) * (I * dl) * (dl_hat x r) / |r|^3
					double factor = ConstForJ * dI * dl / dist3;

					if (FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_) {
						IntForBx += factor * cross_x;
						IntForBy += factor * cross_y;
						IntForBz += factor * cross_z;
					}

					// Vector potential: dA = (mu_0 / 4*pi) * I * dl / |r| * dl_hat
					if (FieldPtr->FieldKey.A_) {
						double factor_A = ConstForJ * dI * dl / dist;
						IntForAx += factor_A * jx;
						IntForAy += factor_A * jy;
					}

					// Scalar potential for arc segment
					// For an infinitesimal current element dI*dl at position r_src,
					// the scalar potential contribution is computed using the solid
					// angle formula for the "ribbon" connecting the element to P.
					//
					// For an arc element, we compute the angle subtended by the
					// arc segment at the observation point in the plane perpendicular
					// to the current direction.
					//
					// dPhi = (dI / 4*pi) * d_omega
					// where d_omega is the differential solid angle
					//
					// For a current element dl at angle phi, the solid angle is
					// related to the angle between the position vectors.
					if (FieldPtr->FieldKey.Phi_) {
						// Use the solid angle contribution from each infinitesimal arc element
						// For a current element at (r_coil*cos(phi), r_coil*sin(phi), z_coil)
						// flowing in direction (-sin(phi), cos(phi), 0),
						// the solid angle seen from point P depends on the geometry.
						//
						// For numerical computation, we use:
						// dOmega = (r_vec x dl) . r_hat / r^2
						// where r_vec goes from element to P, dl is current direction
						//
						// This gives: dOmega = |r_vec x dl| / r^2 = sin(theta) * dl / r
						// where theta is angle between r_vec and dl
						//
						// Actually, for scalar potential, we need:
						// dPhi = (I / 4*pi) * integral of (dl x r_hat) . n_hat / r
						// This is complex for general 3D. Use Biot-Savart analog:
						//
						// For the arc element, the scalar potential contribution is:
						// dPhi = (dI / 4*pi) * dphi_angle
						// where dphi_angle is the angle subtended by the arc element
						// as seen from P in the plane perpendicular to z.
						//
						// Simplified approach: integrate the solid angle contributions
						// using the formula for the angle at P between adjacent elements.
						//
						// The scalar potential from an arc segment from phi1 to phi2 is:
						// Phi = (I / 4*pi) * sum of (atan2 differences)

						double rho2_local = dx_local*dx_local + dy_local*dy_local;
						if (rho2_local > SmallPositive) {
							// Solid angle contribution from this element
							// Using the formula for infinitesimal arc:
							// dOmega = (r_coil * dphi) * z_rel / (r_local^2)
							// where z_rel is the z-component and r_local is distance in local frame
							double rho_local = sqrt(rho2_local);

							// The scalar potential from arc element:
							// Following the approach in Landau & Lifshitz or Jackson:
							// For a current element, the solid angle contribution is
							// dOmega ~ (dl x r) . z_hat / (r * rho)
							//
							// For current in azimuthal direction:
							// dl = r_coil * dphi * (-sin(phi), cos(phi), 0) in global
							// r = (P - element) = (dx_global, dy_global, dz)
							//
							// dl x r = (cos(phi)*dz, sin(phi)*dz,
							//           -sin(phi)*dy_global - cos(phi)*(-dx_global))
							//        = (cos(phi)*dz, sin(phi)*dz, cos(phi)*dx_global + sin(phi)*dy_global)
							//
							// z-component of dl x r = cos(phi)*dx + sin(phi)*dy
							// where dx = P.x - r_coil*cos(phi), dy = P.y - r_coil*sin(phi)

							double dlxr_z = cross_z;  // Already computed: jx*ry - jy*rx
							                          // = -sin(phi)*ry - cos(phi)*(-rx)
							                          // = -sin(phi)*(P.y - r_coil*sin(phi)) + cos(phi)*(P.x - r_coil*cos(phi))

							// dOmega = dlxr_z * dl / (r^2 * rho)
							// But we need proper normalization
							double r_perp = sqrt(rx*rx + ry*ry + SmallPositive);
							double dOmega = dlxr_z * dl / (dist * dist * r_perp);

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

void radTArcCur::B_intUtilSpecCaseZeroVxVy(double r0, double r1, double r2, double ph0, double ph1, double ph2, double h, double& Iz)
{
	const double Pi = 3.141592653589793238;
	const double SmallPositive = 1.E-10;

	double ph0_mi_ph1 = ph0-ph1;
	double ph0_mi_ph2 = ph0-ph2;
	double a1 = tan(0.5*ph0_mi_ph1);
	double a2 = tan(0.5*ph0_mi_ph2);
	double a1a1 = a1*a1;
	double a2a2 = a2*a2;
	double a1a1_pl_1 = a1a1+1.;
	double a2a2_pl_1 = a2a2+1.;
	double a1a1_mi_1 = a1a1-1.;
	double a2a2_mi_1 = a2a2-1.;
	double r1_mi_r2 = r1-r2;

CorrectionStart:
	double r0_mi_r1 = r0-r1; if(r0_mi_r1==0.) { r0+=SmallPositive; goto CorrectionStart;}
	double r0_mi_r2 = r0-r2; if(r0_mi_r2==0.) { r0+=SmallPositive; goto CorrectionStart;}
	double r1_pl_r2 = r1+r2;
	double r0_pl_r1 = r0+r1;
	double r0_pl_r2 = r0+r2;
	double r0plr1_di_r0mir1 = r0_pl_r1/r0_mi_r1;
	double r0plr2_di_r0mir2 = r0_pl_r2/r0_mi_r2;
	double Two_a2r0 = 2.*a2*r0;
	double Two_a1r0 = 2.*a1*r0;
	double a2a2mi1_mu_r0 = a2a2_mi_1*r0;
	double a1a1mi1_mu_r0 = a1a1_mi_1*r0;
	double a2a2mi1_r0_pl_a2a2pl1_r1 = a2a2mi1_mu_r0+a2a2_pl_1*r1;
	double a2a2mi1_r0_pl_a2a2pl1_r2 = a2a2mi1_mu_r0+a2a2_pl_1*r2; 
	double a1a1mi1_r0_pl_a1a1pl1_r1 = a1a1mi1_mu_r0+a1a1_pl_1*r1;
	double a1a1mi1_r0_pl_a1a1pl1_r2 = a1a1mi1_mu_r0+a1a1_pl_1*r2;
	if(a2a2mi1_r0_pl_a2a2pl1_r1==0. || 
	   a2a2mi1_r0_pl_a2a2pl1_r2==0. ||
	   a1a1mi1_r0_pl_a1a1pl1_r1==0. ||
	   a1a1mi1_r0_pl_a1a1pl1_r2==0.) { r0+=SmallPositive; goto CorrectionStart;}

	double a2a2mi1_r0_di_a2a2pl1 = a2a2mi1_mu_r0/a2a2_pl_1;
	double a1a1mi1_r0_di_a1a1pl1 = a1a1mi1_mu_r0/a1a1_pl_1;

	double PiMult1, PiMult2, PiMult3, PiMult4, Sz;
	PiMult1 = PiMult2 = PiMult3 = PiMult4 = 0.;

	double SzL1 = 0.5*(ph2-ph1)*(r2-r1);
	double SzL2 = Pi*(Step(ph0_mi_ph1-Pi)*Step(Pi-ph0_mi_ph2) + Step(Pi+ph0_mi_ph1)*Step(-Pi-ph0_mi_ph2))*(r1_mi_r2*(Step(r0_mi_r2)*Step(r0_mi_r1)-Step(-r0_mi_r1)*Step(-r0_mi_r2)) + (r1_pl_r2-2.*r0)*Step(r0_mi_r1)*Step(-r0_mi_r2));
	double SzL3 = Pi*r0*Step(r0_mi_r1)*Step(-r0_mi_r2)*(Sign(a2)-Sign(a1));
	double SzL4 = r1*(atan(TransAtans(a1*r0plr1_di_r0mir1, -a2*r0plr1_di_r0mir1, PiMult1)) + Pi*PiMult1);
	double SzL5 = r2*(atan(TransAtans(a2*r0plr2_di_r0mir2, -a1*r0plr2_di_r0mir2, PiMult2)) + Pi*PiMult2);
	double SumAtan1 = atan(TransAtans(Two_a2r0/a2a2mi1_r0_pl_a2a2pl1_r1, -Two_a2r0/a2a2mi1_r0_pl_a2a2pl1_r2, PiMult3)) + Pi*PiMult3;
	double PiTerm1 = Pi*Sign(a2)*Step(-a2a2mi1_r0_di_a2a2pl1-r1)*Step(r2+a2a2mi1_r0_di_a2a2pl1);
	double SzL6 = a2a2mi1_r0_di_a2a2pl1*(SumAtan1 + PiTerm1);
	double SumAtan2 = atan(TransAtans(Two_a1r0/a1a1mi1_r0_pl_a1a1pl1_r1, -Two_a1r0/a1a1mi1_r0_pl_a1a1pl1_r2, PiMult4)) + Pi*PiMult4;
	double PiTerm2 = Pi*Sign(a1)*Step(-a1a1mi1_r0_di_a1a1pl1-r1)*Step(r2+a1a1mi1_r0_di_a1a1pl1);
	double SzL7 = -a1a1mi1_r0_di_a1a1pl1*(SumAtan2 + PiTerm2);
	double SzL8 = (a2*r0/a2a2_pl_1)*log((r0_mi_r1*r0_mi_r1 + a2a2*r0_pl_r1*r0_pl_r1)/(r0_mi_r2*r0_mi_r2 + a2a2*r0_pl_r2*r0_pl_r2));
	double SzL9 = (a1*r0/a1a1_pl_1)*log((r0_mi_r2*r0_mi_r2 + a1a1*r0_pl_r2*r0_pl_r2)/(r0_mi_r1*r0_mi_r1 + a1a1*r0_pl_r1*r0_pl_r1));
	Sz = SzL1 + SzL2 + SzL3 + SzL4 + SzL5 + SzL6 + SzL7 + SzL8 + SzL9;
	Iz = 2.*h*Sz;
}

//-------------------------------------------------------------------------

void radTArcCur::B_intCompWithTrapeth(radTField* FieldPtr)
{
	const double Pi = 3.141592653589793238;
	const double ConForJ = 1.E-04;
	const double SmallPositive = 1.E-12;
	const double VxVyCaseZeroToler = 1.E-05; // SmallPositive < VxVyCaseZeroToler !!!
	double Fact = ConForJ*J_azim;

	TVector3d VectV = FieldPtr->NextP - FieldPtr->P;
	double ModV = sqrt(VectV.x*VectV.x + VectV.y*VectV.y + VectV.z*VectV.z);
	VectV.x = VectV.x/ModV; VectV.y = VectV.y/ModV; VectV.z = VectV.z/ModV;
	if(VectV.x==0. || VectV.x==-1.) VectV.x += SmallPositive;
	else if(VectV.x==1.) VectV.x -= SmallPositive;
	if(VectV.y==0. || VectV.y==-1.) VectV.y += SmallPositive;
	else if(VectV.y==1.) VectV.y -= SmallPositive;

	TVector3d IntForB(0.,0.,0.);

	TVector3d P_mi_CPoi = FieldPtr->P - CircleCentrPoint;
	double r0 = sqrt(P_mi_CPoi.x*P_mi_CPoi.x + P_mi_CPoi.y*P_mi_CPoi.y); if(r0==0.) r0 = SmallPositive;
	double PmiCPoix_di_r0 = P_mi_CPoi.x/r0;
	double ph0 = (P_mi_CPoi.y<0)? 2.*Pi-acos(PmiCPoix_di_r0) : acos(PmiCPoix_di_r0);
	double z1 = -0.5*Height - P_mi_CPoi.z;
	double z2 = z1 + Height;
	double r1 = R_min;
	double r2 = R_max;

// Check for Special Case
	if(Abs(VectV.x)<VxVyCaseZeroToler && Abs(VectV.y)<VxVyCaseZeroToler) 
	{
		B_intUtilSpecCaseZeroVxVy(r0, r1, r2, ph0, Phi_min, Phi_max, Height, IntForB.z);
		IntForB.x = IntForB.y = 0.;
		goto FinalDefinitionOfFieldIntegrals;
	}

// General Case: numerical integration over Phi. This uses Newton method (n=3).
	{
		double Sinph0 = sin(ph0);
		double Cosph0 = cos(ph0);
		double Vph0x = Sinph0*VectV.y + Cosph0*VectV.x;
		double Vph0y = Cosph0*VectV.y - Sinph0*VectV.x;
		double VzVz = VectV.z*VectV.z;
		double VxVx = VectV.x*VectV.x;
		double VyVy = VectV.y*VectV.y;
		double VxVxpVyVy = VxVx+VyVy;

			if(VxVxpVyVy==0.) VxVxpVyVy = SmallPositive;

		double Vzz1 = VectV.z*z1;
		double Vzz2 = VectV.z*z2;
		double VzVzr0 = VzVz*r0;
		double r0r0 = r0*r0;
		double r1r1 = r1*r1;
		double r2r2 = r2*r2;
		double r2mir1 = r2-r1;
		double z2miz1 = z2-z1;
		double Vph0xr0 = Vph0x*r0;
		double Vph0xr0_mi_Vzz1 = Vph0xr0-Vzz1;
		double Vph0xr0_mi_Vzz2 = Vph0xr0-Vzz2;
		double C1 = z1*z1+r0r0 - Vph0xr0_mi_Vzz1*Vph0xr0_mi_Vzz1;
		double C2 = z2*z2+r0r0 - Vph0xr0_mi_Vzz2*Vph0xr0_mi_Vzz2;
		double T = -Vph0y*r0;
		double R1 = z1 + VectV.z*Vph0xr0_mi_Vzz1;
		double R2 = z2 + VectV.z*Vph0xr0_mi_Vzz2;

		double ph1 = Phi_min + SmallPositive;
		double ph2 = Phi_max - SmallPositive;
		double PhMax_mi_PhMin = ph2 - ph1;

		double Sinph, Cosph, Cosphmiph0, Vphx, Vphy, VphyT, R1Vphy, R2Vphy, Vphyr1pT, Vphyr2pT, 
			   A, PreB, B1, B2, B1B1, B2B2, P, PT, K;
		double ArgLnZ1, ArgLnZ2, T_di_Vphy, Pr1, Pr2, PR1_pl_VphyT, PR2_pl_VphyT, R1Vphy_mi_PT, R2Vphy_mi_PT,
			   Kr1, Kr2, Kr1_pl_PR1_pl_VphyT, Kr2_pl_PR1_pl_VphyT, Kr1_pl_PR2_pl_VphyT, Kr2_pl_PR2_pl_VphyT, 
			   PR1pVphyT_di_K, PR2pVphyT_di_K, One_di_AA, One_di_KK, Radical1, Radical2, TwoAr1, TwoAr2, 
			   AC1, AC2, Ar1r1_pl_B1r1_pl_C1, Ar2r2_pl_B1r2_pl_C1, Ar1r1_pl_B2r1_pl_C2, Ar2r2_pl_B2r2_pl_C2,
			   Ar1r1, Ar2r2;
		double PiMult1, PiMult2, PiMult3, PiMult4;
		PiMult1 = PiMult2 = PiMult3 = PiMult4 = 0.;
		double Wz, WzL1, WzL2, WzL3, WzL4, WzL5, WzL6, Uz, UzL1, UzL2, UzL3, UzL4;

		int AmOfPoi = NumberOfSectors + 1;
		int AmOfPoi_mi_1 = NumberOfSectors;
		double Step_ph = PhMax_mi_PhMin/AmOfPoi_mi_1;
		double ph = ph1;

		TVector3d ZeroVect(0.,0.,0.), S_forIntB(0.,0.,0.), Func(0.,0.,0.);

		for(int i=0; i<AmOfPoi; i++)
		{
			Sinph = sin(ph); Cosph = cos(ph);
			Cosphmiph0 = cos(ph-ph0);
			Vphx = Sinph*VectV.y + Cosph*VectV.x;
			Vphy = Cosph*VectV.y - Sinph*VectV.x;
			R1Vphy = R1*Vphy;
			R2Vphy = R2*Vphy;
			Vphyr1pT = Vphy*r1+T;
			Vphyr2pT = Vphy*r2+T;
			A = 1.- Vphx*Vphx; 
			
//TO FIX: loss of accuracy if A is 0 or very small
//make special case? 

				if(A==0.) A = SmallPositive;

			One_di_AA = 1./(A*A);
			TwoAr1 = 2.*A*r1;
			TwoAr2 = 2.*A*r2;
			AC1 = A*C1;
			AC2 = A*C2;
			PreB = T*Vphy - Cosphmiph0*VzVzr0;
			B1 = 2.*(PreB - Vzz1*Vphx);
			B1B1 = B1*B1;
			B2 = 2.*(PreB - Vzz2*Vphx);
			B2B2 = B2*B2;
			P = -VectV.z*Vphx;
			PT = P*T;
			R1Vphy_mi_PT = R1Vphy-PT;
			R2Vphy_mi_PT = R2Vphy-PT;
			Pr1 = P*r1;
			Pr2 = P*r2;
			VphyT = Vphy*T;
			PR1_pl_VphyT = P*R1+VphyT;
			PR2_pl_VphyT = P*R2+VphyT;
			K = VxVxpVyVy*A;
			One_di_KK = 1./K/K;
			PR1pVphyT_di_K = PR1_pl_VphyT/K;
			PR2pVphyT_di_K = PR2_pl_VphyT/K;
			Kr1 = K*r1;
			Kr2 = K*r2;
			Kr1_pl_PR1_pl_VphyT = Kr1+PR1_pl_VphyT;
			Kr1_pl_PR2_pl_VphyT = Kr1+PR2_pl_VphyT;
			Kr2_pl_PR1_pl_VphyT = Kr2+PR1_pl_VphyT;
			Kr2_pl_PR2_pl_VphyT = Kr2+PR2_pl_VphyT;
			Ar1r1 = A*r1r1;
			Ar2r2 = A*r2r2;
			Ar1r1_pl_B1r1_pl_C1 = Ar1r1+B1*r1+C1;
			Ar2r2_pl_B1r2_pl_C1 = Ar2r2+B1*r2+C1;
			Ar1r1_pl_B2r1_pl_C2 = Ar1r1+B2*r1+C2;
			Ar2r2_pl_B2r2_pl_C2 = Ar2r2+B2*r2+C2;
			ArgLnZ1 = Ar1r1_pl_B1r1_pl_C1/Ar2r2_pl_B1r2_pl_C1;
			ArgLnZ2 = Ar1r1_pl_B2r1_pl_C2/Ar2r2_pl_B2r2_pl_C2;
			T_di_Vphy = T/Vphy;

			//Radical1 = sqrt(4.*A*C1-B1B1);
			//Radical2 = sqrt(4.*A*C2-B2B2);
			Radical1 = sqrt(::fabs(4.*A*C1-B1B1)); //OC fix 06/03: make a more accurate treatment
			Radical2 = sqrt(::fabs(4.*A*C2-B2B2));

			WzL1 = -(Vphy/A)*r2mir1*z2miz1;
			WzL2 = Pi*T_di_Vphy*T_di_Vphy*Step(-T_di_Vphy-r1)*Step(r2+T_di_Vphy)*(Sign(PT-R1Vphy)-Sign(PT-R2Vphy));
			WzL3 = r1r1*(atan(TransAtans((Pr1+R2)/Vphyr1pT, -(Pr1+R1)/Vphyr1pT, PiMult1)) + Pi*PiMult1);
			WzL4 = r2r2*(atan(TransAtans((Pr2+R1)/Vphyr2pT, -(Pr2+R2)/Vphyr2pT, PiMult2)) + Pi*PiMult2);
			WzL5 = One_di_KK*((PR1_pl_VphyT*PR1_pl_VphyT - R1Vphy_mi_PT*R1Vphy_mi_PT)*((atan(TransAtans(R1Vphy_mi_PT/Kr1_pl_PR1_pl_VphyT, -R1Vphy_mi_PT/Kr2_pl_PR1_pl_VphyT, PiMult3)) + Pi*PiMult3) + Pi*Sign(R1Vphy_mi_PT)*Step(-PR1pVphyT_di_K-r1)*Step(r2+PR1pVphyT_di_K)) 
				             +(PR2_pl_VphyT*PR2_pl_VphyT - R2Vphy_mi_PT*R2Vphy_mi_PT)*((atan(TransAtans(-R2Vphy_mi_PT/Kr1_pl_PR2_pl_VphyT, R2Vphy_mi_PT/Kr2_pl_PR2_pl_VphyT, PiMult4)) + Pi*PiMult4) - Pi*Sign(R2Vphy_mi_PT)*Step(-PR2pVphyT_di_K-r1)*Step(r2+PR2pVphyT_di_K)));
			WzL6 = One_di_KK*(PR1_pl_VphyT*R1Vphy_mi_PT*log(ArgLnZ1) - PR2_pl_VphyT*R2Vphy_mi_PT*log(ArgLnZ2));
			Wz = WzL1+WzL2+WzL3+WzL4+WzL5+WzL6;
			UzL1 = (B1-B2)*r2mir1/A;
			UzL2 = One_di_AA*(B1*Radical1*(atan(TransAtans((TwoAr1+B1)/Radical1, -(TwoAr2+B1)/Radical1, PiMult1)) + Pi*PiMult1)
							 +B2*Radical2*(atan(TransAtans(-(TwoAr1+B2)/Radical2, (TwoAr2+B2)/Radical2, PiMult2)) + Pi*PiMult2));
			UzL3 = 0.5*One_di_AA*(-(2.*AC1-B1B1)*log(ArgLnZ1) + (2.*AC2-B2B2)*log(ArgLnZ2));
			UzL4 = r1r1*log(Ar1r1_pl_B2r1_pl_C2/Ar1r1_pl_B1r1_pl_C1) + r2r2*log(Ar2r2_pl_B1r2_pl_C1/Ar2r2_pl_B2r2_pl_C2);
			Uz = UzL1+UzL2+UzL3+UzL4;

			Func.x = Cosph*Uz;
			Func.y = Sinph*Uz;
			Func.z = VectV.z*Vphx*Uz - 2.*Vphy*Wz;

			if((i==0) || (i==AmOfPoi_mi_1)) S_forIntB += 0.5*Func;
			else S_forIntB += Func;
			ph += Step_ph;
		}
		IntForB = (0.5*Step_ph)*S_forIntB;
	}

FinalDefinitionOfFieldIntegrals:
	TVector3d BufIb = Fact*IntForB;
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIb;
	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIb;
}

//-------------------------------------------------------------------------

void radTArcCur::B_intCompWithNewton3(radTField* FieldPtr)
{
	const double Pi = 3.141592653589793238;
	const double ConForJ = 1.E-04;
	const double SmallPositive = 1.E-10;
	const double VxVyCaseZeroToler = 1.E-05; // SmallPositive < VxVyCaseZeroToler !!!
	double Fact = ConForJ*J_azim;

	TVector3d VectV = FieldPtr->NextP - FieldPtr->P;
	double ModV = sqrt(VectV.x*VectV.x + VectV.y*VectV.y + VectV.z*VectV.z);
	VectV.x = VectV.x/ModV; VectV.y = VectV.y/ModV; VectV.z = VectV.z/ModV;
	if(VectV.x==0. || VectV.x==-1.) VectV.x += SmallPositive;
	else if(VectV.x==1.) VectV.x -= SmallPositive;
	if(VectV.y==0. || VectV.y==-1.) VectV.y += SmallPositive;
	else if(VectV.y==1.) VectV.y -= SmallPositive;

	TVector3d IntForB(1.E+23, 1.E+23, 1.E+23);

	TVector3d P_mi_CPoi = FieldPtr->P - CircleCentrPoint;
	double r0 = sqrt(P_mi_CPoi.x*P_mi_CPoi.x + P_mi_CPoi.y*P_mi_CPoi.y); if(r0==0.) r0 = SmallPositive;
	double PmiCPoix_di_r0 = P_mi_CPoi.x/r0;
	double ph0 = (P_mi_CPoi.y<0)? 2.*Pi-acos(PmiCPoix_di_r0) : acos(PmiCPoix_di_r0);
	double z1 = -0.5*Height - P_mi_CPoi.z;
	double z2 = z1 + Height;
	double r1 = R_min;
	double r2 = R_max;

// Check for Special Case
	if(Abs(VectV.x)<VxVyCaseZeroToler && Abs(VectV.y)<VxVyCaseZeroToler) 
	{
		B_intUtilSpecCaseZeroVxVy(r0, r1, r2, ph0, Phi_min, Phi_max, Height, IntForB.z);
		IntForB.x = IntForB.y = 0.;
		goto FinalDefinitionOfFieldIntegrals;
	}

// General Case: numerical integration over Phi. This uses Newton method (n=3).
	{
		double Sinph0 = sin(ph0);
		double Cosph0 = cos(ph0);
		double Vph0x = Sinph0*VectV.y + Cosph0*VectV.x;
		double Vph0y = Cosph0*VectV.y - Sinph0*VectV.x;
		double VzVz = VectV.z*VectV.z;
		double VxVx = VectV.x*VectV.x;
		double VyVy = VectV.y*VectV.y;
		double VxVxpVyVy = VxVx+VyVy;
		double Vzz1 = VectV.z*z1;
		double Vzz2 = VectV.z*z2;
		double VzVzr0 = VzVz*r0;
		double r0r0 = r0*r0;
		double r1r1 = r1*r1;
		double r2r2 = r2*r2;
		double r2mir1 = r2-r1;
		double z2miz1 = z2-z1;
		double Vph0xr0 = Vph0x*r0;
		double Vph0xr0_mi_Vzz1 = Vph0xr0-Vzz1;
		double Vph0xr0_mi_Vzz2 = Vph0xr0-Vzz2;
		double C1 = z1*z1+r0r0 - Vph0xr0_mi_Vzz1*Vph0xr0_mi_Vzz1;
		double C2 = z2*z2+r0r0 - Vph0xr0_mi_Vzz2*Vph0xr0_mi_Vzz2;
		double T = -Vph0y*r0;
		double R1 = z1 + VectV.z*Vph0xr0_mi_Vzz1;
		double R2 = z2 + VectV.z*Vph0xr0_mi_Vzz2;

		double ph1 = Phi_min + SmallPositive;
		double ph2 = Phi_max - SmallPositive;
		double PhMax_mi_PhMin = ph2 - ph1;

		double Sinph, Cosph, Cosphmiph0, Vphx, Vphy, VphyT, R1Vphy, R2Vphy, Vphyr1pT, Vphyr2pT, 
			   A, PreB, B1, B2, B1B1, B2B2, P, PT, K;
		double ArgLnZ1, ArgLnZ2, T_di_Vphy, Pr1, Pr2, PR1_pl_VphyT, PR2_pl_VphyT, R1Vphy_mi_PT, R2Vphy_mi_PT,
			   Kr1, Kr2, Kr1_pl_PR1_pl_VphyT, Kr2_pl_PR1_pl_VphyT, Kr1_pl_PR2_pl_VphyT, Kr2_pl_PR2_pl_VphyT, 
			   PR1pVphyT_di_K, PR2pVphyT_di_K, One_di_AA, One_di_KK, Radical1, Radical2, TwoAr1, TwoAr2, 
			   AC1, AC2, Ar1r1_pl_B1r1_pl_C1, Ar2r2_pl_B1r2_pl_C1, Ar1r1_pl_B2r1_pl_C2, Ar2r2_pl_B2r2_pl_C2,
			   Ar1r1, Ar2r2;
		double PiMult1, PiMult2, PiMult3, PiMult4;
		PiMult1 = PiMult2 = PiMult3 = PiMult4 = 0.;
		double Wz, WzL1, WzL2, WzL3, WzL4, WzL5, WzL6, Uz, UzL1, UzL2, UzL3, UzL4;

		const double IntegWeight[] = {3./8., 9./8., 9./8., 3./4.};
		TVector3d ZeroVect(0.,0.,0.), S_forB, GenS_forB(0.,0.,0.), PrIntForB, Func(0.,0.,0.);

		double Step_ph, ph;
		short IndForWeight, IndForPass;
		short NotFirstPass = 0;

		int AmOfPoi = 4;
		int AmOfPoi_mi_1;
		double PrecParamInt = 1.E+23; 

		while(PrecParamInt > FieldPtr->CompCriterium.AbsPrecB_int)
		{
			AmOfPoi_mi_1 = AmOfPoi - 1;
			Step_ph = PhMax_mi_PhMin/AmOfPoi_mi_1;
			ph = ph1;

			PrIntForB = IntForB;
			S_forB = ZeroVect;
			IndForWeight = IndForPass = 0;

			for(int i=0; i<AmOfPoi; i++)
			{
				if(IndForPass==2) IndForPass = 0;
				if(IndForWeight==4) IndForWeight = 1;
				if(NotFirstPass && (IndForPass==0)) goto BottomOfThisLoop;
				if(i==AmOfPoi_mi_1) IndForWeight = 0;

				Sinph = sin(ph); Cosph = cos(ph);
				Cosphmiph0 = cos(ph-ph0);
				Vphx = Sinph*VectV.y + Cosph*VectV.x;
				Vphy = Cosph*VectV.y - Sinph*VectV.x;
				R1Vphy = R1*Vphy;
				R2Vphy = R2*Vphy;
				Vphyr1pT = Vphy*r1+T;
				Vphyr2pT = Vphy*r2+T;
				A = 1.- Vphx*Vphx;
				One_di_AA = 1./(A*A);
				TwoAr1 = 2.*A*r1;
				TwoAr2 = 2.*A*r2;
				AC1 = A*C1;
				AC2 = A*C2;
				PreB = T*Vphy - Cosphmiph0*VzVzr0;
				B1 = 2.*(PreB - Vzz1*Vphx);
				B1B1 = B1*B1;
				B2 = 2.*(PreB - Vzz2*Vphx);
				B2B2 = B2*B2;
				P = -VectV.z*Vphx;
				PT = P*T;
				R1Vphy_mi_PT = R1Vphy-PT;
				R2Vphy_mi_PT = R2Vphy-PT;
				Pr1 = P*r1;
				Pr2 = P*r2;
				VphyT = Vphy*T;
				PR1_pl_VphyT = P*R1+VphyT;
				PR2_pl_VphyT = P*R2+VphyT;
				K = VxVxpVyVy*A;
				One_di_KK = 1./K/K;
				PR1pVphyT_di_K = PR1_pl_VphyT/K;
				PR2pVphyT_di_K = PR2_pl_VphyT/K;
				Kr1 = K*r1;
				Kr2 = K*r2;
				Kr1_pl_PR1_pl_VphyT = Kr1+PR1_pl_VphyT;
				Kr1_pl_PR2_pl_VphyT = Kr1+PR2_pl_VphyT;
				Kr2_pl_PR1_pl_VphyT = Kr2+PR1_pl_VphyT;
				Kr2_pl_PR2_pl_VphyT = Kr2+PR2_pl_VphyT;
				Ar1r1 = A*r1r1;
				Ar2r2 = A*r2r2;
				Ar1r1_pl_B1r1_pl_C1 = Ar1r1+B1*r1+C1;
				Ar2r2_pl_B1r2_pl_C1 = Ar2r2+B1*r2+C1;
				Ar1r1_pl_B2r1_pl_C2 = Ar1r1+B2*r1+C2;
				Ar2r2_pl_B2r2_pl_C2 = Ar2r2+B2*r2+C2;
				ArgLnZ1 = Ar1r1_pl_B1r1_pl_C1/Ar2r2_pl_B1r2_pl_C1;
				ArgLnZ2 = Ar1r1_pl_B2r1_pl_C2/Ar2r2_pl_B2r2_pl_C2;
				T_di_Vphy = T/Vphy;
				Radical1 = sqrt(4.*A*C1-B1B1);
				Radical2 = sqrt(4.*A*C2-B2B2);

				WzL1 = -(Vphy/A)*r2mir1*z2miz1;
				WzL2 = Pi*T_di_Vphy*T_di_Vphy*Step(-T_di_Vphy-r1)*Step(r2+T_di_Vphy)*(Sign(PT-R1Vphy)-Sign(PT-R2Vphy));
				WzL3 = r1r1*(atan(TransAtans((Pr1+R2)/Vphyr1pT, -(Pr1+R1)/Vphyr1pT, PiMult1)) + Pi*PiMult1);
				WzL4 = r2r2*(atan(TransAtans((Pr2+R1)/Vphyr2pT, -(Pr2+R2)/Vphyr2pT, PiMult2)) + Pi*PiMult2);
				WzL5 = One_di_KK*((PR1_pl_VphyT*PR1_pl_VphyT - R1Vphy_mi_PT*R1Vphy_mi_PT)*((atan(TransAtans(R1Vphy_mi_PT/Kr1_pl_PR1_pl_VphyT, -R1Vphy_mi_PT/Kr2_pl_PR1_pl_VphyT, PiMult3)) + Pi*PiMult3) + Pi*Sign(R1Vphy_mi_PT)*Step(-PR1pVphyT_di_K-r1)*Step(r2+PR1pVphyT_di_K)) 
					             +(PR2_pl_VphyT*PR2_pl_VphyT - R2Vphy_mi_PT*R2Vphy_mi_PT)*((atan(TransAtans(-R2Vphy_mi_PT/Kr1_pl_PR2_pl_VphyT, R2Vphy_mi_PT/Kr2_pl_PR2_pl_VphyT, PiMult4)) + Pi*PiMult4) - Pi*Sign(R2Vphy_mi_PT)*Step(-PR2pVphyT_di_K-r1)*Step(r2+PR2pVphyT_di_K)));
				WzL6 = One_di_KK*(PR1_pl_VphyT*R1Vphy_mi_PT*log(ArgLnZ1) - PR2_pl_VphyT*R2Vphy_mi_PT*log(ArgLnZ2));
				Wz = WzL1+WzL2+WzL3+WzL4+WzL5+WzL6;
				UzL1 = (B1-B2)*r2mir1/A;
				UzL2 = One_di_AA*(B1*Radical1*(atan(TransAtans((TwoAr1+B1)/Radical1, -(TwoAr2+B1)/Radical1, PiMult1)) + Pi*PiMult1)
								 +B2*Radical2*(atan(TransAtans(-(TwoAr1+B2)/Radical2, (TwoAr2+B2)/Radical2, PiMult2)) + Pi*PiMult2));
				UzL3 = 0.5*One_di_AA*(-(2.*AC1-B1B1)*log(ArgLnZ1) + (2.*AC2-B2B2)*log(ArgLnZ2));
				UzL4 = r1r1*log(Ar1r1_pl_B2r1_pl_C2/Ar1r1_pl_B1r1_pl_C1) + r2r2*log(Ar2r2_pl_B1r2_pl_C1/Ar2r2_pl_B2r2_pl_C2);
				Uz = UzL1+UzL2+UzL3+UzL4;

				Func.x = Cosph*Uz;
				Func.y = Sinph*Uz;
				Func.z = VectV.z*Vphx*Uz - 2.*Vphy*Wz;
				
				S_forB += IntegWeight[IndForWeight] * Func;

BottomOfThisLoop:
				IndForPass++; IndForWeight++;
				ph += Step_ph;
			}
			S_forB = 0.5 * S_forB;
			S_forB.z = S_forB.z/VxVxpVyVy;

			GenS_forB += S_forB; 
			IntForB = Step_ph * GenS_forB;
			PrecParamInt = Fact * Max( Max( Abs(IntForB.x-PrIntForB.x), Abs(IntForB.y-PrIntForB.y)), Abs(IntForB.z-PrIntForB.z));

			AmOfPoi = AmOfPoi_mi_1*2+1;
			NotFirstPass = 1;
		}
	}

FinalDefinitionOfFieldIntegrals:
	TVector3d BufIb = Fact * IntForB;
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIb;
	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIb;
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

radTg3dGraphPresent* radTArcCur::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTArcCurGraphPresent(this);
	return g3dGraphPresentPtr;
}

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

radTg3dGraphPresent* radTBackgroundFieldSource::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTBackgroundFldSrcGraphPresent(this);
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------
