/* rad_evrs_tmethod.h -- EVRS/T-method algebra kernel.
 *
 * This is the C++ production seed for the high-order HCurl eddy-current path:
 *
 *   phi --G--> T --C--> J --D--> rho,     J = curl T.
 *
 * NGSolve owns the finite-element spaces and matrices.  Radia owns the reduced
 * VIM/CLN algebra once the discrete de Rham maps are supplied.
 */
#ifndef RAD_EVRS_TMETHOD_H
#define RAD_EVRS_TMETHOD_H

#include <vector>

namespace radia {
namespace evrs {

struct TMethodAlgebraResult {
    int n_current = 0;
    int n_t = 0;
    int n_phi = 0;
    int n_evrs = 0;
    int n_ports = 0;
    int n_rho = 0;

    // Row-major matrices.
    std::vector<double> current_evrs;      // C Q, n_current x n_evrs
    std::vector<double> resistance_t;      // C^T M_R C, n_t x n_t
    std::vector<double> inductance_t;      // C^T M_L C, n_t x n_t
    std::vector<double> resistance_evrs;   // Q^T R_T Q, n_evrs x n_evrs
    std::vector<double> inductance_evrs;   // Q^T L_T Q, n_evrs x n_evrs
    std::vector<double> port_t;            // C^T P, n_t x n_ports
    std::vector<double> port_evrs;         // Q^T C^T P, n_evrs x n_ports

    double div_curl_norm = 0.0;            // ||D C||_F
    double div_evrs_norm = 0.0;            // ||D C Q||_F
    double resistance_gauge_norm = 0.0;    // ||R_T G||_F
    double inductance_gauge_norm = 0.0;    // ||L_T G||_F
    double port_gauge_norm = 0.0;          // ||G^T C^T P||_F
    double resistance_symmetry_norm = 0.0; // ||R_T - R_T^T||_F
    double inductance_symmetry_norm = 0.0; // ||L_T - L_T^T||_F
    double evrs_resistance_symmetry_norm = 0.0;
    double evrs_inductance_symmetry_norm = 0.0;
    double evrs_resistance_galerkin_residual = 0.0; // ||Q^T R_T Q - (CQ)^T M_R (CQ)||_F
    double evrs_inductance_galerkin_residual = 0.0; // same for M_L
};

TMethodAlgebraResult BuildTMethodAlgebra(
    const std::vector<double>& curl_map, int n_current, int n_t,
    const std::vector<double>& div_map, int n_rho, int div_cols,
    const std::vector<double>& grad_map, int grad_rows, int n_phi,
    const std::vector<double>& evrs_map, int evrs_rows, int n_evrs,
    const std::vector<double>& resistance_current, int resistance_rows, int resistance_cols,
    const std::vector<double>& inductance_current, int inductance_rows, int inductance_cols,
    const std::vector<double>& port_current, int port_rows, int n_ports);

} // namespace evrs
} // namespace radia

#endif // RAD_EVRS_TMETHOD_H
