#include "rad_lie_map_kernel.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace rad_lie {
namespace {

constexpr int kDim = 6;
constexpr int kMaxDegree = 5;
constexpr int kR = 36;       // 6^2
constexpr int kT = 216;      // 6^3
constexpr int kU = 1296;     // 6^4
constexpr int kV = 7776;     // 6^5

// ---------------------------------------------------------------------------
// Degree-5 monomial basis in six variables (mirror of Python
// _degree5_monomial_basis: keys ordered by (total degree, lexicographic),
// product pairs enumerated in the same order so accumulation order -- and
// therefore roundoff -- matches np.bincount exactly).
// ---------------------------------------------------------------------------

struct MonomialBasis {
    std::vector<std::array<int, kDim>> keys;
    std::vector<int> degrees;
    std::vector<int> pair_left;
    std::vector<int> pair_right;
    std::vector<int> pair_out;
    int size = 0;

    int IndexOf(const std::array<int, kDim>& key) const {
        const auto found = std::lower_bound(
            keys.begin(), keys.end(), key,
            [](const std::array<int, kDim>& a, const std::array<int, kDim>& b) {
                const int da = a[0]+a[1]+a[2]+a[3]+a[4]+a[5];
                const int db = b[0]+b[1]+b[2]+b[3]+b[4]+b[5];
                if (da != db) return da < db;
                return a < b;
            });
        return static_cast<int>(found - keys.begin());
    }
};

const MonomialBasis& Basis() {
    static MonomialBasis basis;
    static std::once_flag flag;
    std::call_once(flag, []() {
        for (int total = 0; total <= kMaxDegree; ++total) {
            std::array<int, kDim> key{};
            // Enumerate lexicographically within the fixed total degree.
            for (key[0] = 0; key[0] <= total; ++key[0])
            for (key[1] = 0; key[1] <= total-key[0]; ++key[1])
            for (key[2] = 0; key[2] <= total-key[0]-key[1]; ++key[2])
            for (key[3] = 0; key[3] <= total-key[0]-key[1]-key[2]; ++key[3])
            for (key[4] = 0; key[4] <= total-key[0]-key[1]-key[2]-key[3];
                 ++key[4]) {
                key[5] = total-key[0]-key[1]-key[2]-key[3]-key[4];
                basis.keys.push_back(key);
                basis.degrees.push_back(total);
            }
        }
        basis.size = static_cast<int>(basis.keys.size());
        for (int left = 0; left < basis.size; ++left) {
            for (int right = 0; right < basis.size; ++right) {
                if (basis.degrees[left] + basis.degrees[right] > kMaxDegree)
                    continue;
                std::array<int, kDim> sum;
                for (int axis = 0; axis < kDim; ++axis)
                    sum[axis] = basis.keys[left][axis]
                        + basis.keys[right][axis];
                basis.pair_left.push_back(left);
                basis.pair_right.push_back(right);
                basis.pair_out.push_back(basis.IndexOf(sum));
            }
        }
    });
    return basis;
}

// Truncated polynomial on the shared basis (value path, no tangents).
struct Polynomial {
    std::vector<double> values;
    explicit Polynomial(int size) : values(size, 0.0) {}
};

Polynomial PolynomialConstant(double value) {
    Polynomial result(Basis().size);
    result.values[0] = value;
    return result;
}

void AddInPlace(Polynomial& left, const Polynomial& right) {
    for (std::size_t index = 0; index < left.values.size(); ++index)
        left.values[index] += right.values[index];
}

Polynomial Scale(const Polynomial& polynomial, double factor) {
    Polynomial result(Basis().size);
    for (std::size_t index = 0; index < result.values.size(); ++index)
        result.values[index] = factor * polynomial.values[index];
    return result;
}

Polynomial Multiply(const Polynomial& left, const Polynomial& right) {
    const MonomialBasis& basis = Basis();
    Polynomial result(basis.size);
    const std::size_t pair_count = basis.pair_left.size();
    for (std::size_t pair = 0; pair < pair_count; ++pair) {
        result.values[basis.pair_out[pair]] +=
            left.values[basis.pair_left[pair]]
            * right.values[basis.pair_right[pair]];
    }
    return result;
}

Polynomial Sqrt(const Polynomial& polynomial) {
    const double constant = polynomial.values[0];
    if (!(constant > 0.0))
        throw std::invalid_argument(
            "Hamiltonian square-root expansion point must be positive");
    Polynomial remainder = polynomial;
    remainder.values[0] -= constant;
    Polynomial normalized = Scale(remainder, 1.0 / constant);
    Polynomial result = PolynomialConstant(std::sqrt(constant));
    Polynomial power = PolynomialConstant(1.0);
    double coefficient = 1.0;
    for (int order = 1; order <= kMaxDegree; ++order) {
        power = Multiply(power, normalized);
        coefficient *= (0.5 - (order - 1)) / order;
        AddInPlace(result, Scale(power, std::sqrt(constant) * coefficient));
    }
    return result;
}

// ---------------------------------------------------------------------------
// Hamiltonian jet: H2..H5 symmetric tensors -> Poisson-applied A/F2/F3/F4
// plus the constant and linear terms (mirror of
// _canonical_vector_potential_hamiltonian_jet, value path).
// ---------------------------------------------------------------------------

struct Jet {
    double A[kR];
    double F2[kT];
    double F3[kU];
    double F4[kV];
    double linear[kDim];
    double constant;
};

long long Factorial(int value) {
    long long result = 1;
    for (int factor = 2; factor <= value; ++factor) result *= factor;
    return result;
}

// Scatter one monomial into the rank-`degree` symmetric tensor: every
// DISTINCT ordering of the coordinate multiset receives
// coefficient * prod(power!), exactly like _add_hamiltonian_monomial.
void ScatterMonomial(double* tensor, int degree,
                     const std::array<int, kDim>& key, double coefficient) {
    double derivative = coefficient;
    std::vector<int> slots;
    slots.reserve(degree);
    for (int coordinate = 0; coordinate < kDim; ++coordinate) {
        derivative *= static_cast<double>(Factorial(key[coordinate]));
        for (int repeat = 0; repeat < key[coordinate]; ++repeat)
            slots.push_back(coordinate);
    }
    std::sort(slots.begin(), slots.end());
    do {
        std::size_t flat = 0;
        for (int position = 0; position < degree; ++position)
            flat = flat * kDim + slots[position];
        tensor[flat] += derivative;
    } while (std::next_permutation(slots.begin(), slots.end()));
}

// out[i, rest] = sum_a P[i,a] X[a, rest]  (apply a 6x6 matrix on axis 0).
void ApplyFirstAxis(const double* matrix, const double* source, double* out,
                    int rest) {
    for (int i = 0; i < kDim; ++i) {
        double* row = out + static_cast<std::size_t>(i) * rest;
        std::memset(row, 0, sizeof(double) * rest);
        for (int a = 0; a < kDim; ++a) {
            const double factor = matrix[i*kDim + a];
            if (factor == 0.0) continue;
            const double* source_row = source + static_cast<std::size_t>(a) * rest;
            for (int r = 0; r < rest; ++r) row[r] += factor * source_row[r];
        }
    }
}

Jet BuildJet(const double* Ay, const double* As, int order_count,
             double magnetic_rigidity, double curvature,
             double curvature_sign, double reference_beta,
             bool longitudinal_covariant, const double* poisson) {
    const MonomialBasis& basis = Basis();
    const int degree_limit = order_count - 1;
    if (degree_limit < 1 || degree_limit > kMaxDegree)
        throw std::invalid_argument(
            "Ay/As coefficients need matching square shape (d+1,d+1), 1<=d<=5");
    if (!(std::isfinite(magnetic_rigidity)) || magnetic_rigidity == 0.0
        || !std::isfinite(curvature) || curvature_sign == 0.0
        || !(reference_beta > 0.0) || reference_beta > 1.0)
        throw std::invalid_argument(
            "rigidity, curvature, sign, or beta is invalid");
    double gauge_scale = 1.0;
    for (int index = 0; index < order_count * order_count; ++index) {
        gauge_scale = std::max(gauge_scale, std::fabs(Ay[index]));
        gauge_scale = std::max(gauge_scale, std::fabs(As[index]));
    }
    if (std::max(std::fabs(Ay[0]), std::fabs(As[0]))
            > 64.0 * std::numeric_limits<double>::epsilon() * gauge_scale)
        throw std::invalid_argument(
            "design-orbit gauge requires Ay[0,0]=As[0,0]=0");

    const double normalization = curvature_sign / magnetic_rigidity;

    auto variable = [&basis](int coordinate) {
        Polynomial result(basis.size);
        std::array<int, kDim> key{};
        key[coordinate] = 1;
        result.values[basis.IndexOf(key)] = 1.0;
        return result;
    };
    const Polynomial x_polynomial = variable(0);
    const Polynomial px_polynomial = variable(1);
    const Polynomial py_polynomial = variable(3);
    const Polynomial delta_polynomial = variable(5);

    auto input_polynomial = [&](const double* values) {
        Polynomial result(basis.size);
        for (int x_power = 0; x_power <= degree_limit; ++x_power) {
            for (int y_power = 0; y_power + x_power <= degree_limit;
                 ++y_power) {
                if (x_power == 0 && y_power == 0) continue;
                std::array<int, kDim> key{};
                key[0] = x_power;
                key[2] = y_power;
                result.values[basis.IndexOf(key)] +=
                    normalization * values[x_power*order_count + y_power];
            }
        }
        return result;
    };
    Polynomial Ay_polynomial = input_polynomial(Ay);
    Polynomial As_polynomial = input_polynomial(As);

    Polynomial metric = PolynomialConstant(1.0);
    AddInPlace(metric, Scale(x_polynomial, curvature));
    if (!longitudinal_covariant)
        As_polynomial = Multiply(metric, As_polynomial);

    Polynomial one_plus_delta = PolynomialConstant(1.0);
    AddInPlace(one_plus_delta, delta_polynomial);
    Polynomial mechanical_y = py_polynomial;
    AddInPlace(mechanical_y, Scale(Ay_polynomial, -1.0));
    Polynomial radicand = Multiply(one_plus_delta, one_plus_delta);
    AddInPlace(radicand, Scale(Multiply(px_polynomial, px_polynomial), -1.0));
    AddInPlace(radicand, Scale(Multiply(mechanical_y, mechanical_y), -1.0));
    Polynomial root = Sqrt(radicand);
    const double reference_mass_square =
        1.0 / (reference_beta * reference_beta) - 1.0;
    Polynomial reference_radicand = Multiply(one_plus_delta, one_plus_delta);
    reference_radicand.values[0] += reference_mass_square;
    Polynomial reference = Scale(Sqrt(reference_radicand),
                                 1.0 / reference_beta);
    Polynomial hamiltonian = Scale(Multiply(metric, root), -1.0);
    AddInPlace(hamiltonian, Scale(As_polynomial, -1.0));
    AddInPlace(hamiltonian, reference);

    Jet jet{};
    jet.constant = hamiltonian.values[0];
    double H2[kR] = {};
    double H3[kT] = {};
    double H4[kU] = {};
    double H5[kV] = {};
    for (int flat = 0; flat < basis.size; ++flat) {
        const double coefficient = hamiltonian.values[flat];
        if (coefficient == 0.0) continue;
        const int degree = basis.degrees[flat];
        const std::array<int, kDim>& key = basis.keys[flat];
        if (degree == 1) {
            for (int coordinate = 0; coordinate < kDim; ++coordinate)
                if (key[coordinate] == 1) jet.linear[coordinate] += coefficient;
        } else if (degree == 2) {
            ScatterMonomial(H2, 2, key, coefficient);
        } else if (degree == 3) {
            ScatterMonomial(H3, 3, key, coefficient);
        } else if (degree == 4) {
            ScatterMonomial(H4, 4, key, coefficient);
        } else if (degree == 5) {
            ScatterMonomial(H5, 5, key, coefficient);
        }
    }
    ApplyFirstAxis(poisson, H2, jet.A, kDim);
    ApplyFirstAxis(poisson, H3, jet.F2, kR);
    ApplyFirstAxis(poisson, H4, jet.F3, kT);
    ApplyFirstAxis(poisson, H5, jet.F4, kU);
    return jet;
}

// ---------------------------------------------------------------------------
// Fixed-shape contraction helpers of the factorial map ODE (value path).
// Each partitioned einsum family is computed once as a core product and
// scatter-added at its output-index permutations.
// ---------------------------------------------------------------------------

// out[i,j,k] += sum_ab F2[i,a,b] R[a,j] S[b,k]
void AddQuadratic(const double* F2, const double* R, const double* S,
                  double* out) {
    for (int i = 0; i < kDim; ++i) {
        double G[kR];  // G[a,k] = sum_b F2[i,a,b] S[b,k]
        for (int a = 0; a < kDim; ++a)
            for (int k = 0; k < kDim; ++k) {
                double sum = 0.0;
                for (int b = 0; b < kDim; ++b)
                    sum += F2[(i*kDim + a)*kDim + b] * S[b*kDim + k];
                G[a*kDim + k] = sum;
            }
        for (int j = 0; j < kDim; ++j)
            for (int k = 0; k < kDim; ++k) {
                double sum = 0.0;
                for (int a = 0; a < kDim; ++a)
                    sum += R[a*kDim + j] * G[a*kDim + k];
                out[(i*kDim + j)*kDim + k] += sum;
            }
    }
}

// cross2: out[i,j,k,l] += sum over the three partitions of
// F2[i,a,b] R[a,.] T[b,..]  ("iab,aj,bkl" + "iab,ak,bjl" + "iab,al,bjk").
void AddCrossSecond(const double* F2, const double* R, const double* T,
                    double* out) {
    for (int i = 0; i < kDim; ++i) {
        double core[kDim][kR];  // core[j][kl] = sum_ab F2[iab] R[aj] T[b,kl]
        double G[kDim][kR];     // G[a][kl] = sum_b F2[iab] T[b,kl]
        for (int a = 0; a < kDim; ++a)
            for (int kl = 0; kl < kR; ++kl) {
                double sum = 0.0;
                for (int b = 0; b < kDim; ++b)
                    sum += F2[(i*kDim + a)*kDim + b] * T[b*kR + kl];
                G[a][kl] = sum;
            }
        for (int j = 0; j < kDim; ++j)
            for (int kl = 0; kl < kR; ++kl) {
                double sum = 0.0;
                for (int a = 0; a < kDim; ++a)
                    sum += R[a*kDim + j] * G[a][kl];
                core[j][kl] = sum;
            }
        for (int j = 0; j < kDim; ++j)
            for (int k = 0; k < kDim; ++k)
                for (int l = 0; l < kDim; ++l)
                    out[((i*kDim + j)*kDim + k)*kDim + l] +=
                        core[j][k*kDim + l]      // (j)(kl)
                        + core[k][j*kDim + l]    // (k)(jl)
                        + core[l][j*kDim + k];   // (l)(jk)
    }
}

// cubic: out[i,j,k,l] += sum_abc F3[i,a,b,c] R[a,j] R[b,k] R[c,l]
void AddTransformCubic(const double* F3, const double* R, double* out) {
    // Successive axis transforms keep every intermediate small.
    static thread_local std::vector<double> stage1(kU), stage2(kU);
    // stage1[i,j,b,c] = sum_a F3[i,a,b,c] R[a,j]
    for (int i = 0; i < kDim; ++i)
        for (int j = 0; j < kDim; ++j)
            for (int bc = 0; bc < kR; ++bc) {
                double sum = 0.0;
                for (int a = 0; a < kDim; ++a)
                    sum += F3[(i*kDim + a)*kR + bc] * R[a*kDim + j];
                stage1[(i*kDim + j)*kR + bc] = sum;
            }
    // stage2[i,j,k,c] = sum_b stage1[i,j,b,c] R[b,k]
    for (int ij = 0; ij < kR; ++ij)
        for (int k = 0; k < kDim; ++k)
            for (int c = 0; c < kDim; ++c) {
                double sum = 0.0;
                for (int b = 0; b < kDim; ++b)
                    sum += stage1[(ij*kDim + b)*kDim + c] * R[b*kDim + k];
                stage2[(ij*kDim + k)*kDim + c] = sum;
            }
    for (int ijk = 0; ijk < kT; ++ijk)
        for (int l = 0; l < kDim; ++l) {
            double sum = 0.0;
            for (int c = 0; c < kDim; ++c)
                sum += stage2[ijk*kDim + c] * R[c*kDim + l];
            out[ijk*kDim + l] += sum;
        }
}

// cross42: out[i,j,k,l,m] += the four partitions of
// F2[i,a,b] R[a,.] U[b,...]  ("iab,aj,bklm" + k + l + m variants).
void AddCrossFourthSecond(const double* F2, const double* R, const double* U,
                          double* out) {
    for (int i = 0; i < kDim; ++i) {
        double G[kDim][kT];    // G[a][klm] = sum_b F2[iab] U[b,klm]
        double core[kDim][kT]; // core[j][klm] = sum_a R[aj] G[a][klm]
        for (int a = 0; a < kDim; ++a)
            for (int klm = 0; klm < kT; ++klm) {
                double sum = 0.0;
                for (int b = 0; b < kDim; ++b)
                    sum += F2[(i*kDim + a)*kDim + b] * U[b*kT + klm];
                G[a][klm] = sum;
            }
        for (int j = 0; j < kDim; ++j)
            for (int klm = 0; klm < kT; ++klm) {
                double sum = 0.0;
                for (int a = 0; a < kDim; ++a)
                    sum += R[a*kDim + j] * G[a][klm];
                core[j][klm] = sum;
            }
        for (int j = 0; j < kDim; ++j)
            for (int k = 0; k < kDim; ++k)
                for (int l = 0; l < kDim; ++l)
                    for (int m = 0; m < kDim; ++m)
                        out[(((i*kDim + j)*kDim + k)*kDim + l)*kDim + m] +=
                            core[j][(k*kDim + l)*kDim + m]
                            + core[k][(j*kDim + l)*kDim + m]
                            + core[l][(j*kDim + k)*kDim + m]
                            + core[m][(j*kDim + k)*kDim + l];
    }
}

// pair42: out[i,j,k,l,m] += the three pairings of F2[i,a,b] T[a,..] T[b,..]
// ("iab,ajk,blm" + "iab,ajl,bkm" + "iab,ajm,bkl").
void AddPairFourthSecond(const double* F2, const double* T, double* out) {
    for (int i = 0; i < kDim; ++i) {
        double G[kDim][kR];     // G[a][lm] = sum_b F2[iab] T[b,lm]
        double core[kR][kR];    // core[jk][lm] = sum_a T[a,jk] G[a][lm]
        for (int a = 0; a < kDim; ++a)
            for (int lm = 0; lm < kR; ++lm) {
                double sum = 0.0;
                for (int b = 0; b < kDim; ++b)
                    sum += F2[(i*kDim + a)*kDim + b] * T[b*kR + lm];
                G[a][lm] = sum;
            }
        for (int jk = 0; jk < kR; ++jk)
            for (int lm = 0; lm < kR; ++lm) {
                double sum = 0.0;
                for (int a = 0; a < kDim; ++a)
                    sum += T[a*kR + jk] * G[a][lm];
                core[jk][lm] = sum;
            }
        for (int j = 0; j < kDim; ++j)
            for (int k = 0; k < kDim; ++k)
                for (int l = 0; l < kDim; ++l)
                    for (int m = 0; m < kDim; ++m)
                        out[(((i*kDim + j)*kDim + k)*kDim + l)*kDim + m] +=
                            core[j*kDim + k][l*kDim + m]
                            + core[j*kDim + l][k*kDim + m]
                            + core[j*kDim + m][k*kDim + l];
    }
}

// cross43: out[i,j,k,l,m] += the six set partitions of
// F3[i,a,b,c] T[a,..] R[b,.] R[c,.]  (matching _cross_fourth_third).
void AddCrossFourthThird(const double* F3, const double* R, const double* T,
                         double* out) {
    for (int i = 0; i < kDim; ++i) {
        double E[kDim][kR];     // E[a][lm] = sum_bc F3[iabc] R[bl] R[cm]
        double core[kR][kR];    // core[jk][lm] = sum_a T[a,jk] E[a][lm]
        for (int a = 0; a < kDim; ++a) {
            double partial[kR];  // partial[b*6+m] = sum_c F3[iabc] R[cm]
            for (int b = 0; b < kDim; ++b)
                for (int m = 0; m < kDim; ++m) {
                    double sum = 0.0;
                    for (int c = 0; c < kDim; ++c)
                        sum += F3[((i*kDim + a)*kDim + b)*kDim + c]
                            * R[c*kDim + m];
                    partial[b*kDim + m] = sum;
                }
            for (int l = 0; l < kDim; ++l)
                for (int m = 0; m < kDim; ++m) {
                    double sum = 0.0;
                    for (int b = 0; b < kDim; ++b)
                        sum += R[b*kDim + l] * partial[b*kDim + m];
                    E[a][l*kDim + m] = sum;
                }
        }
        for (int jk = 0; jk < kR; ++jk)
            for (int lm = 0; lm < kR; ++lm) {
                double sum = 0.0;
                for (int a = 0; a < kDim; ++a)
                    sum += T[a*kR + jk] * E[a][lm];
                core[jk][lm] = sum;
            }
        for (int j = 0; j < kDim; ++j)
            for (int k = 0; k < kDim; ++k)
                for (int l = 0; l < kDim; ++l)
                    for (int m = 0; m < kDim; ++m)
                        out[(((i*kDim + j)*kDim + k)*kDim + l)*kDim + m] +=
                            core[j*kDim + k][l*kDim + m]     // (jk)(l)(m)
                            + core[j*kDim + l][k*kDim + m]   // (jl)(k)(m)
                            + core[j*kDim + m][k*kDim + l]   // (jm)(k)(l)
                            + core[k*kDim + l][j*kDim + m]   // (kl)(j)(m)
                            + core[k*kDim + m][j*kDim + l]   // (km)(j)(l)
                            + core[l*kDim + m][j*kDim + k];  // (lm)(j)(k)
    }
}

// quartic: out[i,j,k,l,m] += sum_abcd F4[i,a,b,c,d] R[aj] R[bk] R[cl] R[dm]
void AddTransformQuartic(const double* F4, const double* R, double* out) {
    static thread_local std::vector<double> stage1(kV), stage2(kV);
    // stage1[i,j,bcd] = sum_a F4[i,a,bcd] R[aj]
    for (int i = 0; i < kDim; ++i)
        for (int j = 0; j < kDim; ++j)
            for (int bcd = 0; bcd < kT; ++bcd) {
                double sum = 0.0;
                for (int a = 0; a < kDim; ++a)
                    sum += F4[(i*kDim + a)*kT + bcd] * R[a*kDim + j];
                stage1[(i*kDim + j)*kT + bcd] = sum;
            }
    // stage2[ij,k,cd] = sum_b stage1[ij,b,cd] R[bk]
    for (int ij = 0; ij < kR; ++ij)
        for (int k = 0; k < kDim; ++k)
            for (int cd = 0; cd < kR; ++cd) {
                double sum = 0.0;
                for (int b = 0; b < kDim; ++b)
                    sum += stage1[(ij*kDim + b)*kR + cd] * R[b*kDim + k];
                stage2[(ij*kDim + k)*kR + cd] = sum;
            }
    // stage1 reused: stage1[ijk,l,d] = sum_c stage2[ijk,c,d] R[cl]
    for (int ijk = 0; ijk < kT; ++ijk)
        for (int l = 0; l < kDim; ++l)
            for (int d = 0; d < kDim; ++d) {
                double sum = 0.0;
                for (int c = 0; c < kDim; ++c)
                    sum += stage2[(ijk*kDim + c)*kDim + d] * R[c*kDim + l];
                stage1[(ijk*kDim + l)*kDim + d] = sum;
            }
    for (int ijkl = 0; ijkl < kU; ++ijkl)
        for (int m = 0; m < kDim; ++m) {
            double sum = 0.0;
            for (int d = 0; d < kDim; ++d)
                sum += stage1[ijkl*kDim + d] * R[d*kDim + m];
            out[ijkl*kDim + m] += sum;
        }
}

// ---------------------------------------------------------------------------
// Factorial map ODE right-hand side, RK4 stage step, and composition.
// ---------------------------------------------------------------------------

struct MapState {
    std::vector<double> R, T, U, V;
    MapState() : R(kR, 0.0), T(kT, 0.0), U(kU, 0.0), V(kV, 0.0) {}
    static MapState Identity() {
        MapState state;
        for (int i = 0; i < kDim; ++i) state.R[i*kDim + i] = 1.0;
        return state;
    }
};

void EvaluateRhs(const Jet& jet, const MapState& state, MapState& rate) {
    ApplyFirstAxis(jet.A, state.R.data(), rate.R.data(), kDim);
    ApplyFirstAxis(jet.A, state.T.data(), rate.T.data(), kR);
    AddQuadratic(jet.F2, state.R.data(), state.R.data(), rate.T.data());
    ApplyFirstAxis(jet.A, state.U.data(), rate.U.data(), kT);
    AddCrossSecond(jet.F2, state.R.data(), state.T.data(), rate.U.data());
    AddTransformCubic(jet.F3, state.R.data(), rate.U.data());
    ApplyFirstAxis(jet.A, state.V.data(), rate.V.data(), kU);
    AddCrossFourthSecond(jet.F2, state.R.data(), state.U.data(),
                         rate.V.data());
    AddPairFourthSecond(jet.F2, state.T.data(), rate.V.data());
    AddCrossFourthThird(jet.F3, state.R.data(), state.T.data(),
                        rate.V.data());
    AddTransformQuartic(jet.F4, state.R.data(), rate.V.data());
}

void Shifted(const MapState& state, const MapState& rate, double scale,
             MapState& out) {
    for (int index = 0; index < kR; ++index)
        out.R[index] = state.R[index] + scale * rate.R[index];
    for (int index = 0; index < kT; ++index)
        out.T[index] = state.T[index] + scale * rate.T[index];
    for (int index = 0; index < kU; ++index)
        out.U[index] = state.U[index] + scale * rate.U[index];
    for (int index = 0; index < kV; ++index)
        out.V[index] = state.V[index] + scale * rate.V[index];
}

// One nonautonomous RK4 step from the identity with stage jets.
MapState Rk4StepStages(const Jet& start, const Jet& middle, const Jet& end,
                       double length) {
    const MapState identity = MapState::Identity();
    MapState k1, k2, k3, k4, shifted;
    EvaluateRhs(start, identity, k1);
    Shifted(identity, k1, 0.5 * length, shifted);
    EvaluateRhs(middle, shifted, k2);
    Shifted(identity, k2, 0.5 * length, shifted);
    EvaluateRhs(middle, shifted, k3);
    Shifted(identity, k3, length, shifted);
    EvaluateRhs(end, shifted, k4);
    MapState result;
    const double sixth = length / 6.0;
    for (int index = 0; index < kR; ++index)
        result.R[index] = identity.R[index] + sixth * (k1.R[index]
            + 2.0*k2.R[index] + 2.0*k3.R[index] + k4.R[index]);
    for (int index = 0; index < kT; ++index)
        result.T[index] = sixth * (k1.T[index]
            + 2.0*k2.T[index] + 2.0*k3.T[index] + k4.T[index]);
    for (int index = 0; index < kU; ++index)
        result.U[index] = sixth * (k1.U[index]
            + 2.0*k2.U[index] + 2.0*k3.U[index] + k4.U[index]);
    for (int index = 0; index < kV; ++index)
        result.V[index] = sixth * (k1.V[index]
            + 2.0*k2.V[index] + 2.0*k3.V[index] + k4.V[index]);
    return result;
}

// Composition outer(inner) of factorial fourth-order maps (value path).
MapState Compose(const MapState& outer, const MapState& inner) {
    MapState result;
    // R = Ro Ri
    ApplyFirstAxis(outer.R.data(), inner.R.data(), result.R.data(), kDim);
    // T = Ro Ti + To[Ri, Ri]
    ApplyFirstAxis(outer.R.data(), inner.T.data(), result.T.data(), kR);
    AddQuadratic(outer.T.data(), inner.R.data(), inner.R.data(),
                 result.T.data());
    // U = Ro Ui + cross2(To, Ri, Ti) + cubic(Uo, Ri)
    ApplyFirstAxis(outer.R.data(), inner.U.data(), result.U.data(), kT);
    AddCrossSecond(outer.T.data(), inner.R.data(), inner.T.data(),
                   result.U.data());
    AddTransformCubic(outer.U.data(), inner.R.data(), result.U.data());
    // V = Ro Vi + cross42(To, Ri, Ui) + pair42(To, Ti) + cross43(Uo, Ri, Ti)
    //     + quartic(Vo, Ri)
    ApplyFirstAxis(outer.R.data(), inner.V.data(), result.V.data(), kU);
    AddCrossFourthSecond(outer.T.data(), inner.R.data(), inner.U.data(),
                         result.V.data());
    AddPairFourthSecond(outer.T.data(), inner.T.data(), result.V.data());
    AddCrossFourthThird(outer.U.data(), inner.R.data(), inner.T.data(),
                        result.V.data());
    AddTransformQuartic(outer.V.data(), inner.R.data(), result.V.data());
    return result;
}

}  // namespace

// ---------------------------------------------------------------------------
// Public entry: the per-segment stage-jet integration loop.
// ---------------------------------------------------------------------------

void LieMapTensorsFromSpolyArrays(
    const double* Ay,
    const double* As,
    std::size_t n_segments,
    std::size_t s_order_count,
    std::size_t transverse_order_count,
    const double* lengths,
    const double* curvature,
    std::size_t curvature_columns,
    double magnetic_rigidity,
    double curvature_sign,
    double reference_beta,
    int longitudinal_covariant,
    const double* poisson,
    double maximum_step_m,
    long long maximum_steps,
    double reference_orbit_tolerance,
    double* R_out,
    double* T_out,
    double* U_out,
    double* V_out,
    double* linear_out,
    double* worst_linear_out) {
    if (!Ay || !As || !lengths || !curvature || !poisson
        || !R_out || !T_out || !U_out || !V_out || !linear_out
        || !worst_linear_out)
        throw std::invalid_argument("lie map kernel: null array pointer");
    if (n_segments == 0 || s_order_count == 0 || curvature_columns == 0)
        throw std::invalid_argument("lie map kernel: empty inputs");
    if (transverse_order_count < 2 || transverse_order_count > 6)
        throw std::invalid_argument(
            "Ay/As coefficients need matching square shape (d+1,d+1), 1<=d<=5");
    if (!(maximum_step_m > 0.0) || maximum_steps < 1
        || !(reference_orbit_tolerance > 0.0))
        throw std::invalid_argument(
            "A-map LIE integration/factorization limits are invalid");

    const std::size_t transverse_block =
        transverse_order_count * transverse_order_count;
    std::vector<double> Ay_stage(transverse_block);
    std::vector<double> As_stage(transverse_block);

    auto stage_arrays = [&](std::size_t segment, double zeta) {
        std::fill(Ay_stage.begin(), Ay_stage.end(), 0.0);
        std::fill(As_stage.begin(), As_stage.end(), 0.0);
        double power = 1.0;
        for (std::size_t order = 0; order < s_order_count; ++order) {
            const double* Ay_block = Ay
                + (segment*s_order_count + order) * transverse_block;
            const double* As_block = As
                + (segment*s_order_count + order) * transverse_block;
            for (std::size_t index = 0; index < transverse_block; ++index) {
                Ay_stage[index] += power * Ay_block[index];
                As_stage[index] += power * As_block[index];
            }
            power *= zeta;
        }
    };
    auto stage_curvature = [&](std::size_t segment, double zeta) {
        double value = 0.0;
        double power = 1.0;
        for (std::size_t column = 0; column < curvature_columns; ++column) {
            value += power * curvature[segment*curvature_columns + column];
            power *= zeta;
        }
        return value;
    };
    auto stage_jet = [&](std::size_t segment, double zeta) {
        stage_arrays(segment, zeta);
        return BuildJet(Ay_stage.data(), As_stage.data(),
                        static_cast<int>(transverse_order_count),
                        magnetic_rigidity, stage_curvature(segment, zeta),
                        curvature_sign, reference_beta,
                        longitudinal_covariant != 0, poisson);
    };
    auto gate = [&](std::size_t segment, const Jet& jet) {
        double maximum_linear = 0.0;
        for (int coordinate = 0; coordinate < kDim; ++coordinate)
            maximum_linear = std::max(maximum_linear,
                                      std::fabs(jet.linear[coordinate]));
        if (maximum_linear > reference_orbit_tolerance) {
            std::ostringstream message;
            message.precision(6);
            message << std::scientific
                    << "reference orbit is not a Hamiltonian fixed "
                    << "trajectory at segment " << segment << ": max |H1|="
                    << maximum_linear << " exceeds "
                    << reference_orbit_tolerance
                    << "; H1[x,px,y,py,zeta,delta]=[";
            for (int coordinate = 0; coordinate < kDim; ++coordinate)
                message << (coordinate ? "," : "") << jet.linear[coordinate];
            message << "]";
            throw std::invalid_argument(message.str());
        }
        return maximum_linear;
    };

    MapState accumulated = MapState::Identity();
    long long total_steps = 0;
    double worst_linear_global = 0.0;
    for (std::size_t segment = 0; segment < n_segments; ++segment) {
        Jet start = stage_jet(segment, -1.0);
        double worst_linear = gate(segment, start);
        for (int coordinate = 0; coordinate < kDim; ++coordinate)
            linear_out[segment*kDim + coordinate] = start.linear[coordinate];
        const double length = lengths[segment];
        if (!(length > 0.0))
            throw std::invalid_argument(
                "lie map kernel: segment lengths must be positive");
        const long long step_count = static_cast<long long>(
            std::ceil(length / maximum_step_m));
        total_steps += step_count;
        if (total_steps > maximum_steps)
            throw std::invalid_argument(
                "A-map LIE integration exceeds maximum_steps");
        const double step_length = length / static_cast<double>(step_count);
        MapState local = MapState::Identity();
        for (long long sub_step = 0; sub_step < step_count; ++sub_step) {
            const double zeta_mid = -1.0
                + 2.0 * (static_cast<double>(sub_step) + 0.5)
                / static_cast<double>(step_count);
            const double zeta_end = std::min(
                -1.0 + 2.0 * (static_cast<double>(sub_step) + 1.0)
                / static_cast<double>(step_count), 1.0);
            Jet middle = stage_jet(segment, zeta_mid);
            Jet end = stage_jet(segment, zeta_end);
            for (const Jet* stage : {&middle, &end}) {
                const double stage_linear = gate(segment, *stage);
                if (stage_linear > worst_linear) {
                    worst_linear = stage_linear;
                    for (int coordinate = 0; coordinate < kDim; ++coordinate)
                        linear_out[segment*kDim + coordinate] =
                            stage->linear[coordinate];
                }
            }
            const MapState step = Rk4StepStages(start, middle, end,
                                                step_length);
            local = Compose(step, local);
            start = end;
        }
        accumulated = Compose(local, accumulated);
        worst_linear_global = std::max(worst_linear_global, worst_linear);
    }
    std::memcpy(R_out, accumulated.R.data(), sizeof(double) * kR);
    std::memcpy(T_out, accumulated.T.data(), sizeof(double) * kT);
    std::memcpy(U_out, accumulated.U.data(), sizeof(double) * kU);
    std::memcpy(V_out, accumulated.V.data(), sizeof(double) * kV);
    *worst_linear_out = worst_linear_global;
}

}  // namespace rad_lie
