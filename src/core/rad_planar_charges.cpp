//-------------------------------------------------------------------------
// rad_planar_charges.cpp -- shared 2D planar field + Maxwell torque (see .h).
//-------------------------------------------------------------------------
#include "rad_planar_charges.h"
#include "rad_parallel.h"

#include <vector>
#include <cmath>
#include <stdexcept>

namespace rad_planar_charges {

static const double TWO_PI = 6.283185307179586476925286766559;
static const double MU0 = 1.2566370614359172953850573533118e-6;   // 4 pi x 1e-7

void Field(int nq, const double* Xq, const double* Q,
           int nP, const double* P, double* Hout)
{
	// parallel over observation points (each is an independent O(nq) reduction).
	// Self-wrap so a BARE call (no caller `with TaskManager()`) still runs parallel; RegionTaskManager
	// reuses the caller's pool when one is already active (nested = no-op).  This kernel is
	// allocation-free -> scales cleanly (matches the 3D rad.Fld self-wrap in
	// radTApplication::ComputeFieldBatch).
	ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
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

void FieldAz(int nq, const double* Xq, const double* Q,
             int nP, const double* P, double* Azout)
{
	ngcore::RegionTaskManager rtm(radia::GetMaxThreads());   // self-wrap -> bare call parallel (see Field)
	ngcore::ParallelFor(ngcore::IntRange(nP), [&](int i){
		double px = P[2 * i], py = P[2 * i + 1];
		double az = 0.0;
		for(int a = 0; a < nq; a++)
			az += Q[a] * std::atan2(py - Xq[2 * a + 1], px - Xq[2 * a]);
		Azout[i] = MU0 / TWO_PI * az;
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

void MaxwellForceCircle(int nq, const double* Xq, const double* Q,
                        double Rc, double cx, double cy, int n,
                        double hextx, double hexty, double* Fout)
{
	if(n < 8) n = 8;
	std::vector<double> P(2 * n), H(2 * n);
	for(int i = 0; i < n; i++){
		double phi = TWO_PI * i / n;
		P[2 * i] = cx + Rc * std::cos(phi);
		P[2 * i + 1] = cy + Rc * std::sin(phi);
	}
	Field(nq, Xq, Q, n, P.data(), H.data());
	double fx = 0.0, fy = 0.0;
	for(int i = 0; i < n; i++){
		double phi = TWO_PI * i / n, c = std::cos(phi), s = std::sin(phi);
		double Hx = H[2 * i] + hextx, Hy = H[2 * i + 1] + hexty;
		double Hr = Hx * c + Hy * s;             // radial = H . n
		double H2 = Hx * Hx + Hy * Hy;
		fx += Hr * Hx - 0.5 * H2 * c;            // T . n , x-component
		fy += Hr * Hy - 0.5 * H2 * s;
	}
	Fout[0] = MU0 * Rc * (TWO_PI / n) * fx;
	Fout[1] = MU0 * Rc * (TWO_PI / n) * fy;
}

PlanarFieldEvaluator::PlanarFieldEvaluator(
    std::vector<double> positions, std::vector<double> strengths,
    std::vector<int> image_masks, std::vector<double> image_signs)
{
	if(positions.size() != 2*strengths.size())
		throw std::invalid_argument("PlanarFieldEvaluator: position/strength size mismatch");
	if(image_masks.size() != image_signs.size())
		throw std::invalid_argument("PlanarFieldEvaluator: image mask/sign size mismatch");
	m_baseSourceCount = strengths.size();
	m_imageCount = image_masks.size();
	m_positions.reserve(positions.size()*(image_masks.size()+1));
	m_strengths.reserve(strengths.size()*(image_masks.size()+1));
	m_positions.insert(m_positions.end(), positions.begin(), positions.end());
	m_strengths.insert(m_strengths.end(), strengths.begin(), strengths.end());
	for(std::size_t image = 0; image < image_masks.size(); ++image){
		const int mask = image_masks[image];
		if(mask < 1 || mask > 3)
			throw std::invalid_argument("PlanarFieldEvaluator: 2D image mask must be in [1,3]");
		if(!std::isfinite(image_signs[image]))
			throw std::invalid_argument("PlanarFieldEvaluator: image sign must be finite");
		for(std::size_t source = 0; source < strengths.size(); ++source){
			double x = positions[2*source], y = positions[2*source+1];
			if(mask & 1) x = -x;
			if(mask & 2) y = -y;
			m_positions.push_back(x);
			m_positions.push_back(y);
			m_strengths.push_back(image_signs[image]*strengths[source]);
		}
	}
}

void PlanarFieldEvaluator::EvaluateField(
    const double* points, std::size_t count, double* output) const
{
	Field(static_cast<int>(m_strengths.size()), m_positions.data(), m_strengths.data(),
	      static_cast<int>(count), points, output);
}

void PlanarFieldEvaluator::EvaluateAz(
    const double* points, std::size_t count, double* output) const
{
	FieldAz(static_cast<int>(m_strengths.size()), m_positions.data(), m_strengths.data(),
	        static_cast<int>(count), points, output);
}

} // namespace rad_planar_charges
