//-------------------------------------------------------------------------
// rad_planar_charges.cpp -- shared 2D planar field + Maxwell torque (see .h).
//-------------------------------------------------------------------------
#include "rad_planar_charges.h"
#include "rad_parallel.h"

#include <vector>
#include <cmath>

namespace rad_planar_charges {

static const double TWO_PI = 6.283185307179586476925286766559;
static const double MU0 = 1.2566370614359172953850573533118e-6;   // 4 pi x 1e-7

void Field(int nq, const double* Xq, const double* Q,
           int nP, const double* P, double* Hout)
{
	// parallel over observation points (each is an independent O(nq) reduction)
	ngcore::ParallelFor(ngcore::IntRange(nP), [&](int i){
		double px = P[2 * i], py = P[2 * i + 1];
		double hx = 0.0, hy = 0.0;
		for(int a = 0; a < nq; a++){
			double dx = px - Xq[2 * a], dy = py - Xq[2 * a + 1];
			double r2 = dx * dx + dy * dy;
			if(r2 <= 0.0) continue;                 // observation point ON a charge -> skip
			double f = Q[a] / r2;
			hx += f * dx; hy += f * dy;
		}
		Hout[2 * i] = hx / TWO_PI;
		Hout[2 * i + 1] = hy / TWO_PI;
	});
}

double MaxwellTorqueCircle(int nq, const double* Xq, const double* Q,
                           double Rc, double cx, double cy, int n,
                           double hextx, double hexty)
{
	if(n < 8) n = 8;
	std::vector<double> P(2 * n), H(2 * n);
	for(int i = 0; i < n; i++){
		double phi = TWO_PI * i / n;
		P[2 * i] = cx + Rc * std::cos(phi);
		P[2 * i + 1] = cy + Rc * std::sin(phi);
	}
	Field(nq, Xq, Q, n, P.data(), H.data());
	double acc = 0.0;
	for(int i = 0; i < n; i++){
		double phi = TWO_PI * i / n, c = std::cos(phi), s = std::sin(phi);
		double Hx = H[2 * i] + hextx, Hy = H[2 * i + 1] + hexty;  // total = body + uniform applied
		double Hr = Hx * c + Hy * s;      // radial
		double Hp = -Hx * s + Hy * c;     // azimuthal
		acc += Hr * Hp;
	}
	return MU0 * Rc * Rc * (TWO_PI / n) * acc;
}

} // namespace rad_planar_charges
