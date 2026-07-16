#include "rad_hybrid_vim_schur.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace radia {
namespace hybrid_vim {
namespace {

void CheckSize(const std::vector<Complex>& a, int rows, int cols, const char* name) {
    if (rows < 0 || cols < 0)
        throw std::runtime_error(std::string(name) + ": negative dimension");
    const auto expected = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
    if (a.size() != expected)
        throw std::runtime_error(
            std::string(name) + ": expected " + std::to_string(rows) + " x " +
            std::to_string(cols) + " entries, got " + std::to_string(a.size()));
}

std::vector<Complex> MatMul(
        const std::vector<Complex>& a, int ar, int ac,
        const std::vector<Complex>& b, int br, int bc,
        const char* name) {
    if (ac != br)
        throw std::runtime_error(std::string(name) + ": inner dimension mismatch");
    std::vector<Complex> out(static_cast<std::size_t>(ar) * static_cast<std::size_t>(bc), Complex(0.0, 0.0));
    for (int i = 0; i < ar; ++i) {
        for (int k = 0; k < ac; ++k) {
            const Complex aik = a[static_cast<std::size_t>(i) * ac + k];
            if (aik == Complex(0.0, 0.0)) continue;
            for (int j = 0; j < bc; ++j)
                out[static_cast<std::size_t>(i) * bc + j] +=
                    aik * b[static_cast<std::size_t>(k) * bc + j];
        }
    }
    return out;
}

std::vector<Complex> Subtract(const std::vector<Complex>& a, const std::vector<Complex>& b, const char* name) {
    if (a.size() != b.size())
        throw std::runtime_error(std::string(name) + ": size mismatch");
    std::vector<Complex> out(a.size(), Complex(0.0, 0.0));
    for (std::size_t i = 0; i < a.size(); ++i)
        out[i] = a[i] - b[i];
    return out;
}

std::vector<Complex> SolveDenseSquare(
        std::vector<Complex> a,
        std::vector<Complex> b,
        int n,
        int nrhs) {
    CheckSize(a, n, n, "eliminate_eliminate");
    CheckSize(b, n, nrhs, "schur rhs");
    const double eps = 1.0e-30;

    for (int k = 0; k < n; ++k) {
        int pivot = k;
        double best = std::abs(a[static_cast<std::size_t>(k) * n + k]);
        for (int i = k + 1; i < n; ++i) {
            const double cand = std::abs(a[static_cast<std::size_t>(i) * n + k]);
            if (cand > best) {
                best = cand;
                pivot = i;
            }
        }
        if (best <= eps)
            throw std::runtime_error("eliminate_eliminate is singular to working precision");
        if (pivot != k) {
            for (int j = 0; j < n; ++j)
                std::swap(a[static_cast<std::size_t>(k) * n + j],
                          a[static_cast<std::size_t>(pivot) * n + j]);
            for (int j = 0; j < nrhs; ++j)
                std::swap(b[static_cast<std::size_t>(k) * nrhs + j],
                          b[static_cast<std::size_t>(pivot) * nrhs + j]);
        }

        for (int i = k + 1; i < n; ++i) {
            const Complex factor = a[static_cast<std::size_t>(i) * n + k] /
                                   a[static_cast<std::size_t>(k) * n + k];
            a[static_cast<std::size_t>(i) * n + k] = Complex(0.0, 0.0);
            for (int j = k + 1; j < n; ++j)
                a[static_cast<std::size_t>(i) * n + j] -=
                    factor * a[static_cast<std::size_t>(k) * n + j];
            for (int j = 0; j < nrhs; ++j)
                b[static_cast<std::size_t>(i) * nrhs + j] -=
                    factor * b[static_cast<std::size_t>(k) * nrhs + j];
        }
    }

    std::vector<Complex> x(static_cast<std::size_t>(n) * nrhs, Complex(0.0, 0.0));
    for (int col = 0; col < nrhs; ++col) {
        for (int i = n - 1; i >= 0; --i) {
            Complex sum = b[static_cast<std::size_t>(i) * nrhs + col];
            for (int j = i + 1; j < n; ++j)
                sum -= a[static_cast<std::size_t>(i) * n + j] *
                       x[static_cast<std::size_t>(j) * nrhs + col];
            x[static_cast<std::size_t>(i) * nrhs + col] =
                sum / a[static_cast<std::size_t>(i) * n + i];
        }
    }
    return x;
}

void CheckPositive(double value, const char* name) {
    if (!(std::isfinite(value) && value > 0.0))
        throw std::runtime_error(std::string(name) + " must be positive");
}

} // namespace

std::vector<Complex> DenseSolve(
    const std::vector<Complex>& matrix, int matrix_rows, int matrix_cols,
    const std::vector<Complex>& rhs, int rhs_rows, int rhs_cols) {

    CheckSize(matrix, matrix_rows, matrix_cols, "matrix");
    CheckSize(rhs, rhs_rows, rhs_cols, "rhs");
    if (matrix_rows != matrix_cols)
        throw std::runtime_error("matrix must be square");
    if (matrix_rows <= 0)
        throw std::runtime_error("matrix must not be empty");
    if (rhs_rows != matrix_rows)
        throw std::runtime_error("rhs rows must match matrix rows");
    if (rhs_cols <= 0)
        throw std::runtime_error("rhs must contain at least one column");
    return SolveDenseSquare(matrix, rhs, matrix_rows, rhs_cols);
}

std::vector<Complex> DenseSchurComplement(
    const std::vector<Complex>& keep_keep, int n_keep, int kk_cols,
    const std::vector<Complex>& keep_eliminate, int ke_rows, int n_eliminate,
    const std::vector<Complex>& eliminate_keep, int ek_rows, int ek_cols,
    const std::vector<Complex>& eliminate_eliminate, int ee_rows, int ee_cols) {

    CheckSize(keep_keep, n_keep, kk_cols, "keep_keep");
    CheckSize(keep_eliminate, ke_rows, n_eliminate, "keep_eliminate");
    CheckSize(eliminate_keep, ek_rows, ek_cols, "eliminate_keep");
    CheckSize(eliminate_eliminate, ee_rows, ee_cols, "eliminate_eliminate");
    if (kk_cols != n_keep)
        throw std::runtime_error("keep_keep must be square");
    if (ke_rows != n_keep)
        throw std::runtime_error("keep_eliminate rows must match keep_keep rows");
    if (ek_rows != n_eliminate)
        throw std::runtime_error("eliminate_keep rows must match eliminate dimension");
    if (ek_cols != n_keep)
        throw std::runtime_error("eliminate_keep columns must match keep dimension");
    if (ee_rows != n_eliminate || ee_cols != n_eliminate)
        throw std::runtime_error("eliminate_eliminate must be square in eliminate space");

    const auto solved = SolveDenseSquare(eliminate_eliminate, eliminate_keep, n_eliminate, n_keep);
    const auto correction = MatMul(keep_eliminate, n_keep, n_eliminate,
                                   solved, n_eliminate, n_keep,
                                   "K_ke solve(K_ee, K_ek)");
    return Subtract(keep_keep, correction, "Schur complement");
}

Complex SkinImpedance(Complex s, double sigma, double mu) {
    CheckPositive(sigma, "sigma");
    CheckPositive(mu, "mu");
    return std::sqrt(mu * s / sigma);
}

Complex SIBCAdmittanceTail(Complex s, double surface_measure, double sigma, double mu) {
    CheckPositive(surface_measure, "surface_measure");
    CheckPositive(sigma, "sigma");
    CheckPositive(mu, "mu");
    return surface_measure * std::sqrt(sigma / (mu * s));
}

Complex SIBCSchurTerminationImpedance(Complex s, double k_sibc, double d) {
    CheckPositive(k_sibc, "k_sibc");
    if (!(std::isfinite(d) && d >= 0.0))
        throw std::runtime_error("d must be non-negative");
    return (s + d) / (k_sibc * std::sqrt(s));
}

Complex SIBCSchurTerminationAdmittance(Complex s, double k_sibc, double d) {
    return Complex(1.0, 0.0) / SIBCSchurTerminationImpedance(s, k_sibc, d);
}

} // namespace hybrid_vim
} // namespace radia
