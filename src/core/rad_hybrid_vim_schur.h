/* rad_hybrid_vim_schur.h -- dense mixed Galerkin / SIBC helper kernels.
 *
 * These are small production kernels for the IGTE mixed reduction algebra:
 *
 *   S = K_kk - K_ke K_ee^{-1} K_ek
 *
 * The matrices are dense, row-major, and complex-valued because the reduced
 * VIM/CLN system is normally evaluated at s = j omega.
 */
#ifndef RAD_HYBRID_VIM_SCHUR_H
#define RAD_HYBRID_VIM_SCHUR_H

#include <complex>
#include <vector>

namespace radia {
namespace hybrid_vim {

using Complex = std::complex<double>;

std::vector<Complex> DenseSolve(
    const std::vector<Complex>& matrix, int matrix_rows, int matrix_cols,
    const std::vector<Complex>& rhs, int rhs_rows, int rhs_cols);

std::vector<Complex> DenseSchurComplement(
    const std::vector<Complex>& keep_keep, int n_keep, int kk_cols,
    const std::vector<Complex>& keep_eliminate, int ke_rows, int n_eliminate,
    const std::vector<Complex>& eliminate_keep, int ek_rows, int ek_cols,
    const std::vector<Complex>& eliminate_eliminate, int ee_rows, int ee_cols);

Complex SkinImpedance(Complex s, double sigma, double mu);

Complex SIBCAdmittanceTail(Complex s, double surface_measure, double sigma, double mu);

Complex SIBCSchurTerminationImpedance(Complex s, double k_sibc, double d);

Complex SIBCSchurTerminationAdmittance(Complex s, double k_sibc, double d);

} // namespace hybrid_vim
} // namespace radia

#endif // RAD_HYBRID_VIM_SCHUR_H
