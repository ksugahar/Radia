#ifndef dipole_h
#define dipole_h
#include "exafmm_t.h"
#include "fmm_scale_invariant.h"
#include "geometry.h"
#include "intrinsics.h"
#include "timer.h"

namespace exafmm_t {
  /**
   * @brief A derived FMM class for magnetic dipole kernel.
   *
   * Magnetic dipole potential:
   *   phi_m(r) = (1/4*pi) * (m . r) / |r|^3
   *
   * H-field (gradient):
   *   H = -grad(phi_m) = (1/4*pi) * [3*(m.r)*r/r^5 - m/r^3]
   *
   * Source value: 3 components (mx, my, mz) per source
   * Target value: 4 components (phi, Hx, Hy, Hz) per target
   */
  class DipoleFmm : public FmmScaleInvariant<real_t> {
  public:
    DipoleFmm() {}

    DipoleFmm(int p_, int ncrit_, std::string filename_=std::string()) :
      FmmScaleInvariant<real_t>(p_, ncrit_, filename_)
    {
      if (this->filename.empty()) {
        this->filename = std::string("dipole_") + (std::is_same<real_t, float>::value ? "f" : "d")
                       + std::string("_p") + std::to_string(p) + std::string(".dat");
      }
    }

    /**
     * @brief Compute potentials at targets induced by dipole sources directly.
     *
     * phi_m = (1/4*pi) * (m . r) / r^3
     *
     * @param src_coord Vector of source coordinates [x1,y1,z1, x2,y2,z2, ...]
     * @param src_value Vector of dipole moments [mx1,my1,mz1, mx2,my2,mz2, ...]
     * @param trg_coord Vector of target coordinates
     * @param trg_value Vector of potentials at targets (scalar, 1 per target)
     */
    void potential_P2P(RealVec& src_coord, RealVec& src_value, RealVec& trg_coord, RealVec& trg_value) {
      int nsrcs = src_coord.size() / 3;
      int ntrgs = trg_coord.size() / 3;

      for (int t = 0; t < ntrgs; t++) {
        real_t potential = 0;
        real_t tx = trg_coord[3*t+0];
        real_t ty = trg_coord[3*t+1];
        real_t tz = trg_coord[3*t+2];

        for (int s = 0; s < nsrcs; ++s) {
          // Vector from source to target
          real_t rx = tx - src_coord[3*s+0];
          real_t ry = ty - src_coord[3*s+1];
          real_t rz = tz - src_coord[3*s+2];

          real_t r2 = rx*rx + ry*ry + rz*rz;
          if (r2 != 0) {
            real_t r = std::sqrt(r2);
            real_t r3 = r2 * r;

            // Dipole moment
            real_t mx = src_value[3*s+0];
            real_t my = src_value[3*s+1];
            real_t mz = src_value[3*s+2];

            // m . r
            real_t m_dot_r = mx*rx + my*ry + mz*rz;

            potential += m_dot_r / r3;
          }
        }
        trg_value[t] += potential / (4*PI);
      }
      add_flop((long long)ntrgs*(long long)nsrcs*20);
    }

    /**
     * @brief Compute potentials and H-field at targets induced by dipole sources.
     *
     * phi_m = (1/4*pi) * (m . r) / r^3
     * H = -grad(phi_m) = (1/4*pi) * [3*(m.r)*r/r^5 - m/r^3]
     *
     * @param src_coord Vector of source coordinates [x1,y1,z1, ...]
     * @param src_value Vector of dipole moments [mx1,my1,mz1, ...]
     * @param trg_coord Vector of target coordinates
     * @param trg_value Vector of [phi, Hx, Hy, Hz] at targets (4 per target)
     */
    void gradient_P2P(RealVec& src_coord, RealVec& src_value, RealVec& trg_coord, RealVec& trg_value) {
      int nsrcs = src_coord.size() / 3;
      int ntrgs = trg_coord.size() / 3;

      for (int t = 0; t < ntrgs; t++) {
        real_t potential = 0;
        real_t Hx = 0, Hy = 0, Hz = 0;

        real_t tx = trg_coord[3*t+0];
        real_t ty = trg_coord[3*t+1];
        real_t tz = trg_coord[3*t+2];

        for (int s = 0; s < nsrcs; ++s) {
          // Vector from source to target
          real_t rx = tx - src_coord[3*s+0];
          real_t ry = ty - src_coord[3*s+1];
          real_t rz = tz - src_coord[3*s+2];

          real_t r2 = rx*rx + ry*ry + rz*rz;
          if (r2 != 0) {
            real_t r = std::sqrt(r2);
            real_t r3 = r2 * r;
            real_t r5 = r3 * r2;

            // Dipole moment
            real_t mx = src_value[3*s+0];
            real_t my = src_value[3*s+1];
            real_t mz = src_value[3*s+2];

            // m . r
            real_t m_dot_r = mx*rx + my*ry + mz*rz;

            // Scalar potential: phi = (m.r) / r^3
            potential += m_dot_r / r3;

            // H-field: H = [3*(m.r)*r/r^5 - m/r^3]
            real_t coeff = 3.0 * m_dot_r / r5;
            Hx += coeff * rx - mx / r3;
            Hy += coeff * ry - my / r3;
            Hz += coeff * rz - mz / r3;
          }
        }

        // Apply 1/(4*pi) factor
        real_t inv_4pi = 1.0 / (4*PI);
        trg_value[4*t+0] += potential * inv_4pi;
        trg_value[4*t+1] += Hx * inv_4pi;
        trg_value[4*t+2] += Hy * inv_4pi;
        trg_value[4*t+3] += Hz * inv_4pi;
      }
      add_flop((long long)ntrgs*(long long)nsrcs*35);
    }
  };
}  // end namespace exafmm_t
#endif
