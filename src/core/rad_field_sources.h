// Numeric source records shared by the existing CPU evaluator and CUDA.
// Coefficients retain the CPU monomial order and origin convention.
#pragma once
#include <cstddef>
#include <cstdint>
#include <type_traits>
namespace rad_hdiv {
struct FieldTetSource {
    double v[4][3]{};
    double coefficient[84]{};
    int degree = 3;
    double origin[3]{};
};

struct FieldTriSource {
    double v[3][3]{};
    double sigma0 = 0.0;
    double slope[3]{};
    double hessian[3][3]{};
    double origin[3]{};
};

struct FieldCubicTriSource {
    double v[3][3]{};
    double coefficient[35]{};
    int degree = 3;
    double origin[3]{};
};

struct FieldCurvedTetSource {
    double nodes[10][3]{};
    double coefficient[4]{}; // 1, xi, eta, zeta
};

struct FieldCurvedTriSource {
    double nodes[6][3]{};
    double coefficient[6]{}; // 1, eta, eta^2, xi, xi*eta, xi^2
};

enum class FieldSourceKind : std::uint8_t {
    Tet, Triangle, CubicTriangle, CurvedTet, CurvedTriangle, Point
};
struct FieldSourceRef { FieldSourceKind kind; std::size_t index; };
static_assert(std::is_trivially_copyable<FieldTetSource>::value, "CUDA source must be trivially copyable");
static_assert(std::is_trivially_copyable<FieldTriSource>::value, "CUDA source must be trivially copyable");
static_assert(std::is_trivially_copyable<FieldCubicTriSource>::value, "CUDA source must be trivially copyable");
static_assert(std::is_trivially_copyable<FieldCurvedTetSource>::value, "CUDA source must be trivially copyable");
static_assert(std::is_trivially_copyable<FieldCurvedTriSource>::value, "CUDA source must be trivially copyable");
} // namespace rad_hdiv
