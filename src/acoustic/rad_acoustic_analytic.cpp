#include "rad_acoustic_analytic.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace radia::acoustics {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr int kMaximumTerms = 512;

bool Finite(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

void CheckPositive(double value, const char* name) {
    if (!std::isfinite(value) || value <= 0.0)
        throw std::invalid_argument(std::string(name) + " must be positive and finite");
}

int CheckTerms(int terms, bool allow_default) {
    if ((allow_default && terms == -1) || terms == 0)
        return terms;
    if (terms < 0 || terms > kMaximumTerms) {
        if (allow_default)
            throw std::invalid_argument(
                "terms must be -1 for automatic selection or lie between 0 and 512");
        throw std::invalid_argument("terms must lie between 0 and 512");
    }
    return terms;
}

Complex PowI(int order) {
    switch (order & 3) {
        case 0: return Complex(1.0, 0.0);
        case 1: return Complex(0.0, 1.0);
        case 2: return Complex(-1.0, 0.0);
        default: return Complex(0.0, -1.0);
    }
}

Complex SphericalJSeries(int order, Complex z) {
    if (z == Complex(0.0, 0.0))
        return order == 0 ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
    Complex numerator(1.0, 0.0);
    for (int i = 0; i < order; ++i)
        numerator *= z;
    long double denominator = 1.0L;
    for (int odd = 1; odd <= 2 * order + 1; odd += 2)
        denominator *= static_cast<long double>(odd);
    Complex term = numerator / static_cast<double>(denominator);
    Complex sum = term;
    const Complex z2 = z * z;
    for (int index = 0; index < 256; ++index) {
        term *= -z2 /
            (2.0 * static_cast<double>(index + 1) *
             static_cast<double>(2 * order + 2 * index + 3));
        sum += term;
        if (std::abs(term) <= 4.0 * std::numeric_limits<double>::epsilon() *
                                  std::max(1.0, std::abs(sum)))
            break;
    }
    return sum;
}

std::vector<Complex> SphericalJSequence(Complex z, int maximum_order) {
    if (!Finite(z) || maximum_order < 0 || maximum_order > kMaximumTerms + 1)
        throw std::invalid_argument("invalid spherical-Bessel request");
    std::vector<Complex> result(static_cast<std::size_t>(maximum_order + 1));
    if (std::abs(z) < 0.25) {
        for (int order = 0; order <= maximum_order; ++order)
            result[static_cast<std::size_t>(order)] = SphericalJSeries(order, z);
        return result;
    }

    const int start = maximum_order + 50 + static_cast<int>(std::ceil(std::abs(z)));
    std::vector<Complex> work(static_cast<std::size_t>(start + 2));
    work[static_cast<std::size_t>(start)] = Complex(1.0, 0.0);
    for (int order = start; order > 0; --order) {
        work[static_cast<std::size_t>(order - 1)] =
            (static_cast<double>(2 * order + 1) / z) *
                work[static_cast<std::size_t>(order)] -
            work[static_cast<std::size_t>(order + 1)];
        if (std::abs(work[static_cast<std::size_t>(order - 1)]) > 1.0e150) {
            for (int index = order - 1; index <= start + 1; ++index)
                work[static_cast<std::size_t>(index)] *= 1.0e-150;
        }
    }
    const Complex exact0 = std::sin(z) / z;
    const Complex exact1 = std::sin(z) / (z * z) - std::cos(z) / z;
    const double denominator = std::norm(work[0]) + std::norm(work[1]);
    if (!(denominator > 0.0) || !std::isfinite(denominator))
        throw std::runtime_error("spherical-Bessel normalization failed");
    const Complex scale =
        (std::conj(work[0]) * exact0 + std::conj(work[1]) * exact1) /
        denominator;
    for (int order = 0; order <= maximum_order; ++order)
        result[static_cast<std::size_t>(order)] =
            work[static_cast<std::size_t>(order)] * scale;
    return result;
}

std::vector<Complex> SphericalHankelSequence(Complex z, int maximum_order) {
    if (!Finite(z) || z == Complex(0.0, 0.0) || maximum_order < 0 ||
        maximum_order > kMaximumTerms + 1)
        throw std::invalid_argument("invalid spherical-Hankel request");
    std::vector<Complex> result(static_cast<std::size_t>(maximum_order + 1));
    const Complex exponential = std::exp(Complex(0.0, 1.0) * z);
    result[0] = Complex(0.0, -1.0) * exponential / z;
    if (maximum_order == 0)
        return result;
    result[1] = -exponential * (z + Complex(0.0, 1.0)) / (z * z);
    for (int order = 1; order < maximum_order; ++order)
        result[static_cast<std::size_t>(order + 1)] =
            (static_cast<double>(2 * order + 1) / z) *
                result[static_cast<std::size_t>(order)] -
            result[static_cast<std::size_t>(order - 1)];
    return result;
}

Complex FirstDerivative(const std::vector<Complex>& values, int order, Complex z) {
    if (order == 0)
        return -values[1];
    return values[static_cast<std::size_t>(order - 1)] -
           (static_cast<double>(order + 1) / z) *
               values[static_cast<std::size_t>(order)];
}

Complex SecondDerivative(const std::vector<Complex>& values, int order, Complex z) {
    const Complex first = FirstDerivative(values, order, z);
    return -(2.0 / z) * first -
           (1.0 - static_cast<double>(order * (order + 1)) / (z * z)) *
               values[static_cast<std::size_t>(order)];
}

std::vector<double> LegendreSequence(double x, int maximum_order) {
    x = std::clamp(x, -1.0, 1.0);
    std::vector<double> result(static_cast<std::size_t>(maximum_order + 1));
    result[0] = 1.0;
    if (maximum_order == 0)
        return result;
    result[1] = x;
    for (int order = 1; order < maximum_order; ++order)
        result[static_cast<std::size_t>(order + 1)] =
            ((2.0 * order + 1.0) * x * result[static_cast<std::size_t>(order)] -
             order * result[static_cast<std::size_t>(order - 1)]) /
            (order + 1.0);
    return result;
}

struct PreparedPoints {
    std::size_t count = 0;
    std::vector<double> radius;
    std::vector<double> safe_radius;
    std::vector<double> cosine;
    double maximum_radius = 0.0;
};

PreparedPoints PreparePoints(const std::vector<double>& points, double sphere_radius,
                             bool allow_interior) {
    if (points.empty() || points.size() % 3 != 0)
        throw std::invalid_argument("points must be a nonempty N-by-3 array");
    PreparedPoints prepared;
    prepared.count = points.size() / 3;
    prepared.radius.resize(prepared.count);
    prepared.safe_radius.resize(prepared.count);
    prepared.cosine.resize(prepared.count);
    for (std::size_t index = 0; index < prepared.count; ++index) {
        const double x = points[3 * index];
        const double y = points[3 * index + 1];
        const double z = points[3 * index + 2];
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
            throw std::invalid_argument("points must contain finite coordinates");
        const double r = std::sqrt(x * x + y * y + z * z);
        if (!allow_interior && r < sphere_radius * (1.0 - 1.0e-9))
            throw std::invalid_argument(
                "evaluation points must lie on or outside the sphere r >= R");
        prepared.radius[index] = r;
        prepared.safe_radius[index] = std::max(r, 1.0e-30);
        prepared.cosine[index] = z / prepared.safe_radius[index];
        prepared.maximum_radius = std::max(prepared.maximum_radius, r);
    }
    return prepared;
}

std::vector<std::vector<double>> PointLegendre(
    const PreparedPoints& points, int terms) {
    std::vector<std::vector<double>> result(points.count);
    for (std::size_t index = 0; index < points.count; ++index)
        result[index] = LegendreSequence(points.cosine[index], terms);
    return result;
}

std::vector<Complex> IncidentField(
    double wavenumber, const std::vector<double>& points) {
    const std::size_t count = points.size() / 3;
    std::vector<Complex> result(count);
    for (std::size_t index = 0; index < count; ++index)
        result[index] = std::exp(
            Complex(0.0, wavenumber * points[3 * index + 2]));
    return result;
}

double TailMagnitude(const std::vector<Complex>& values) {
    double result = 0.0;
    for (Complex value : values)
        result = std::max(result, std::abs(value));
    return result;
}

template <std::size_t Size>
std::array<Complex, Size> SolveSmall(
    std::array<std::array<Complex, Size>, Size> matrix,
    std::array<Complex, Size> rhs) {
    std::array<double, Size> column_scale{};
    for (std::size_t column = 0; column < Size; ++column) {
        for (std::size_t row = 0; row < Size; ++row)
            column_scale[column] = std::max(
                column_scale[column], std::abs(matrix[row][column]));
        if (!(column_scale[column] > 0.0) || !std::isfinite(column_scale[column]))
            throw std::runtime_error("acoustic modal system is singular");
        for (std::size_t row = 0; row < Size; ++row)
            matrix[row][column] /= column_scale[column];
    }
    for (std::size_t row = 0; row < Size; ++row) {
        double row_scale = 0.0;
        for (std::size_t column = 0; column < Size; ++column)
            row_scale = std::max(row_scale, std::abs(matrix[row][column]));
        if (!(row_scale > 0.0) || !std::isfinite(row_scale))
            throw std::runtime_error("acoustic modal system is singular");
        for (std::size_t column = 0; column < Size; ++column)
            matrix[row][column] /= row_scale;
        rhs[row] /= row_scale;
    }
    for (std::size_t pivot = 0; pivot < Size; ++pivot) {
        std::size_t selected = pivot;
        for (std::size_t row = pivot + 1; row < Size; ++row)
            if (std::abs(matrix[row][pivot]) > std::abs(matrix[selected][pivot]))
                selected = row;
        double remaining_scale = 0.0;
        for (std::size_t row = pivot; row < Size; ++row)
            for (std::size_t column = pivot; column < Size; ++column)
                remaining_scale = std::max(
                    remaining_scale, std::abs(matrix[row][column]));
        if (std::abs(matrix[selected][pivot]) <=
            64.0 * std::numeric_limits<double>::epsilon() * remaining_scale)
            throw std::runtime_error("acoustic modal system is singular");
        if (selected != pivot) {
            std::swap(matrix[selected], matrix[pivot]);
            std::swap(rhs[selected], rhs[pivot]);
        }
        for (std::size_t row = pivot + 1; row < Size; ++row) {
            const Complex factor = matrix[row][pivot] / matrix[pivot][pivot];
            matrix[row][pivot] = Complex(0.0, 0.0);
            for (std::size_t column = pivot + 1; column < Size; ++column)
                matrix[row][column] -= factor * matrix[pivot][column];
            rhs[row] -= factor * rhs[pivot];
        }
    }
    std::array<Complex, Size> solution{};
    for (std::size_t reverse = 0; reverse < Size; ++reverse) {
        const std::size_t row = Size - reverse - 1;
        Complex value = rhs[row];
        for (std::size_t column = row + 1; column < Size; ++column)
            value -= matrix[row][column] * solution[column];
        solution[row] = value / matrix[row][row];
    }
    for (std::size_t column = 0; column < Size; ++column)
        solution[column] /= column_scale[column];
    return solution;
}

Complex ElasticCoefficient(
    int order, double k, double radius, double longitudinal_speed,
    double shear_speed, double density_ratio,
    const std::vector<Complex>& j_fluid,
    const std::vector<Complex>& h_fluid,
    const std::vector<Complex>& j_longitudinal,
    const std::vector<Complex>& j_transverse) {
    const double omega = k;
    const double k_longitudinal = omega / longitudinal_speed;
    const double x = k * radius;
    const double xl = k_longitudinal * radius;
    const double mu = density_ratio * shear_speed * shear_speed;
    const double lambda = density_ratio *
        (longitudinal_speed * longitudinal_speed -
         2.0 * shear_speed * shear_speed);
    const double fluid_factor = k / (omega * omega);
    const Complex jx = j_fluid[static_cast<std::size_t>(order)];
    const Complex hx = h_fluid[static_cast<std::size_t>(order)];
    const Complex jxl = j_longitudinal[static_cast<std::size_t>(order)];
    const Complex djx = FirstDerivative(j_fluid, order, x);
    const Complex dhx = FirstDerivative(h_fluid, order, x);
    const Complex djxl = FirstDerivative(j_longitudinal, order, xl);

    if (shear_speed == 0.0) {
        std::array<std::array<Complex, 2>, 2> matrix{{
            {{fluid_factor * dhx, -k_longitudinal * djxl}},
            {{hx, -lambda * k_longitudinal * k_longitudinal * jxl}},
        }};
        std::array<Complex, 2> rhs{{-fluid_factor * djx, -jx}};
        return SolveSmall<2>(matrix, rhs)[0];
    }

    const double k_transverse = omega / shear_speed;
    const double xt = k_transverse * radius;
    const Complex jxt = j_transverse[static_cast<std::size_t>(order)];
    const Complex djxt = FirstDerivative(j_transverse, order, xt);
    const Complex ddjxl = SecondDerivative(j_longitudinal, order, xl);
    const Complex ddjxt = SecondDerivative(j_transverse, order, xt);
    const double angular = static_cast<double>(order * (order + 1));
    const Complex ur_a = k_longitudinal * djxl;
    const Complex ur_b = angular / radius * jxt;
    const Complex dur_a = k_longitudinal * k_longitudinal * ddjxl;
    const Complex dur_b = angular *
        (k_transverse * djxt / radius - jxt / (radius * radius));
    const Complex srr_a =
        -lambda * k_longitudinal * k_longitudinal * jxl + 2.0 * mu * dur_a;
    const Complex srr_b = 2.0 * mu * dur_b;
    const Complex va_a = jxl / radius;
    const Complex va_b = jxt / radius + k_transverse * djxt;
    const Complex vp_a =
        -jxl / (radius * radius) + k_longitudinal * djxl / radius;
    const Complex vp_b =
        -jxt / (radius * radius) + k_transverse * djxt / radius +
        k_transverse * k_transverse * ddjxt;
    const Complex srt_a = mu * (ur_a / radius + vp_a - va_a / radius);
    const Complex srt_b = mu * (ur_b / radius + vp_b - va_b / radius);
    std::array<std::array<Complex, 3>, 3> matrix{{
        {{fluid_factor * dhx, -ur_a, -ur_b}},
        {{hx, srr_a, srr_b}},
        {{Complex(0.0, 0.0), srt_a, srt_b}},
    }};
    std::array<Complex, 3> rhs{{-fluid_factor * djx, -jx, Complex(0.0, 0.0)}};
    return SolveSmall<3>(matrix, rhs)[0];
}

}  // namespace

ScatteringResult SoftSphereScattering(
    double wavenumber, double radius, const std::vector<double>& points,
    int terms) {
    CheckPositive(wavenumber, "wavenumber");
    CheckPositive(radius, "radius");
    CheckTerms(terms, true);
    const PreparedPoints prepared = PreparePoints(points, radius, false);
    const int count = terms < 0
        ? static_cast<int>(std::ceil(wavenumber * radius)) + 12
        : terms;
    CheckTerms(count, false);
    const auto boundary_j = SphericalJSequence(wavenumber * radius, count);
    const auto boundary_h = SphericalHankelSequence(wavenumber * radius, count);
    const auto legendre = PointLegendre(prepared, count);
    std::vector<std::vector<Complex>> point_h(prepared.count);
    for (std::size_t index = 0; index < prepared.count; ++index)
        point_h[index] = SphericalHankelSequence(
            wavenumber * prepared.radius[index], count);
    std::vector<Complex> scattered(prepared.count);
    std::vector<Complex> last(prepared.count);
    for (int order = 0; order <= count; ++order) {
        const Complex coefficient =
            -PowI(order) * static_cast<double>(2 * order + 1) *
            boundary_j[static_cast<std::size_t>(order)] /
            boundary_h[static_cast<std::size_t>(order)];
        for (std::size_t index = 0; index < prepared.count; ++index) {
            last[index] = coefficient * point_h[index][static_cast<std::size_t>(order)] *
                          legendre[index][static_cast<std::size_t>(order)];
            scattered[index] += last[index];
        }
    }
    ScatteringResult result;
    result.kind = "soft_sphere_plane_wave_scattering_series";
    result.wavenumber = wavenumber;
    result.radius = radius;
    result.terms = count;
    result.truncation_tail = TailMagnitude(last);
    result.scattered = std::move(scattered);
    result.incident = IncidentField(wavenumber, points);
    result.total.resize(prepared.count);
    for (std::size_t index = 0; index < prepared.count; ++index)
        result.total[index] = result.incident[index] + result.scattered[index];
    return result;
}

ScatteringResult RigidSphereScattering(
    double wavenumber, double radius, const std::vector<double>& points,
    int terms) {
    CheckPositive(wavenumber, "wavenumber");
    CheckPositive(radius, "radius");
    CheckTerms(terms, false);
    const PreparedPoints prepared = PreparePoints(points, radius, false);
    const int count = std::max(
        terms, static_cast<int>(std::ceil(wavenumber *
                                         std::max(radius, prepared.maximum_radius))) + 12);
    CheckTerms(count, false);
    const auto boundary_j = SphericalJSequence(wavenumber * radius, count + 1);
    const auto boundary_h = SphericalHankelSequence(wavenumber * radius, count + 1);
    const auto legendre = PointLegendre(prepared, count);
    std::vector<std::vector<Complex>> point_h(prepared.count);
    for (std::size_t index = 0; index < prepared.count; ++index)
        point_h[index] = SphericalHankelSequence(
            wavenumber * prepared.radius[index], count);
    std::vector<Complex> scattered(prepared.count);
    std::vector<Complex> last(prepared.count);
    for (int order = 0; order <= count; ++order) {
        const Complex coefficient =
            -PowI(order) * static_cast<double>(2 * order + 1) *
            FirstDerivative(boundary_j, order, wavenumber * radius) /
            FirstDerivative(boundary_h, order, wavenumber * radius);
        for (std::size_t index = 0; index < prepared.count; ++index) {
            last[index] = coefficient * point_h[index][static_cast<std::size_t>(order)] *
                          legendre[index][static_cast<std::size_t>(order)];
            scattered[index] += last[index];
        }
    }
    ScatteringResult result;
    result.kind = "rigid_sphere_plane_wave_scattering_series";
    result.wavenumber = wavenumber;
    result.radius = radius;
    result.terms = count;
    result.truncation_tail = TailMagnitude(last);
    result.scattered = std::move(scattered);
    result.incident = IncidentField(wavenumber, points);
    result.total.resize(prepared.count);
    for (std::size_t index = 0; index < prepared.count; ++index)
        result.total[index] = result.incident[index] + result.scattered[index];
    return result;
}

ScatteringResult FluidSphereScattering(
    double wavenumber, double radius, const std::vector<double>& points,
    double interior_wavenumber, double density_ratio, int terms) {
    CheckPositive(wavenumber, "wavenumber");
    CheckPositive(radius, "radius");
    CheckPositive(interior_wavenumber, "interior_wavenumber");
    CheckPositive(density_ratio, "density_ratio");
    CheckTerms(terms, true);
    const PreparedPoints prepared = PreparePoints(points, radius, true);
    const int requested = terms < 0 ? 0 : terms;
    const int count = std::max(
        requested,
        static_cast<int>(std::ceil(std::max(
            wavenumber * std::max(radius, prepared.maximum_radius),
            interior_wavenumber * radius))) + 12);
    CheckTerms(count, false);
    const Complex x0 = wavenumber * radius;
    const Complex x1 = interior_wavenumber * radius;
    const auto j0_boundary = SphericalJSequence(x0, count + 1);
    const auto h0_boundary = SphericalHankelSequence(x0, count + 1);
    const auto j1_boundary = SphericalJSequence(x1, count + 1);
    const auto legendre = PointLegendre(prepared, count);
    std::vector<std::vector<Complex>> point_j(prepared.count);
    std::vector<std::vector<Complex>> point_h(prepared.count);
    std::vector<std::uint8_t> inside(prepared.count);
    for (std::size_t index = 0; index < prepared.count; ++index) {
        inside[index] = prepared.radius[index] <= radius * (1.0 + 1.0e-12);
        if (inside[index])
            point_j[index] = SphericalJSequence(
                interior_wavenumber * prepared.safe_radius[index], count);
        else {
            point_j[index] = SphericalJSequence(
                wavenumber * prepared.safe_radius[index], count);
            point_h[index] = SphericalHankelSequence(
                wavenumber * prepared.safe_radius[index], count);
        }
    }
    std::vector<Complex> total(prepared.count);
    std::vector<Complex> last(prepared.count);
    for (int order = 0; order <= count; ++order) {
        const Complex incident_coefficient =
            PowI(order) * static_cast<double>(2 * order + 1);
        const Complex beta =
            (interior_wavenumber / density_ratio) *
            FirstDerivative(j1_boundary, order, x1) /
            j1_boundary[static_cast<std::size_t>(order)];
        const Complex scattered_coefficient =
            -incident_coefficient *
            (wavenumber * FirstDerivative(j0_boundary, order, x0) -
             beta * j0_boundary[static_cast<std::size_t>(order)]) /
            (wavenumber * FirstDerivative(h0_boundary, order, x0) -
             beta * h0_boundary[static_cast<std::size_t>(order)]);
        const Complex interior_coefficient =
            (incident_coefficient * j0_boundary[static_cast<std::size_t>(order)] +
             scattered_coefficient * h0_boundary[static_cast<std::size_t>(order)]) /
            j1_boundary[static_cast<std::size_t>(order)];
        for (std::size_t index = 0; index < prepared.count; ++index) {
            if (inside[index])
                last[index] = interior_coefficient *
                    point_j[index][static_cast<std::size_t>(order)] *
                    legendre[index][static_cast<std::size_t>(order)];
            else
                last[index] =
                    (incident_coefficient * point_j[index][static_cast<std::size_t>(order)] +
                     scattered_coefficient * point_h[index][static_cast<std::size_t>(order)]) *
                    legendre[index][static_cast<std::size_t>(order)];
            total[index] += last[index];
        }
    }
    ScatteringResult result;
    result.kind = "fluid_sphere_transmission_scattering_series";
    result.wavenumber = wavenumber;
    result.interior_wavenumber = interior_wavenumber;
    result.density_ratio = density_ratio;
    result.radius = radius;
    result.terms = count;
    result.truncation_tail = TailMagnitude(last);
    result.incident = IncidentField(wavenumber, points);
    result.total = std::move(total);
    result.inside_mask = std::move(inside);
    return result;
}

ScatteringResult ElasticSphereScattering(
    double wavenumber, double radius, const std::vector<double>& points,
    double longitudinal_speed, double shear_speed, double density_ratio,
    int terms) {
    CheckPositive(wavenumber, "wavenumber");
    CheckPositive(radius, "radius");
    CheckPositive(longitudinal_speed, "longitudinal_speed");
    if (!std::isfinite(shear_speed) || shear_speed < 0.0)
        throw std::invalid_argument("shear_speed must be nonnegative and finite");
    CheckPositive(density_ratio, "density_ratio");
    CheckTerms(terms, false);
    const PreparedPoints prepared = PreparePoints(points, radius, false);
    const int count = terms > 0
        ? terms
        : static_cast<int>(std::ceil(wavenumber * radius)) + 10;
    CheckTerms(count, false);
    const double x = wavenumber * radius;
    const double xl = wavenumber * radius / longitudinal_speed;
    const auto j_fluid = SphericalJSequence(x, count + 1);
    const auto h_fluid = SphericalHankelSequence(x, count + 1);
    const auto j_longitudinal = SphericalJSequence(xl, count + 1);
    std::vector<Complex> j_transverse;
    if (shear_speed > 0.0)
        j_transverse = SphericalJSequence(
            wavenumber * radius / shear_speed, count + 1);
    const auto legendre = PointLegendre(prepared, count);
    std::vector<std::vector<Complex>> point_h(prepared.count);
    for (std::size_t index = 0; index < prepared.count; ++index)
        point_h[index] = SphericalHankelSequence(
            wavenumber * prepared.radius[index], count);
    std::vector<Complex> coefficients(static_cast<std::size_t>(count + 1));
    for (int order = 0; order <= count; ++order)
        coefficients[static_cast<std::size_t>(order)] = ElasticCoefficient(
            order, wavenumber, radius, longitudinal_speed, shear_speed,
            density_ratio, j_fluid, h_fluid, j_longitudinal, j_transverse);
    std::vector<Complex> scattered(prepared.count);
    std::vector<Complex> last(prepared.count);
    for (int order = 0; order <= count; ++order)
        for (std::size_t index = 0; index < prepared.count; ++index) {
            last[index] = PowI(order) * static_cast<double>(2 * order + 1) *
                coefficients[static_cast<std::size_t>(order)] *
                point_h[index][static_cast<std::size_t>(order)] *
                legendre[index][static_cast<std::size_t>(order)];
            scattered[index] += last[index];
        }
    ScatteringResult result;
    result.kind = "elastic_solid_sphere_faran_scattering_series";
    result.wavenumber = wavenumber;
    result.radius = radius;
    result.longitudinal_speed = longitudinal_speed;
    result.shear_speed = shear_speed;
    result.density_ratio = density_ratio;
    result.terms = count;
    result.truncation_tail = TailMagnitude(last);
    result.scattered = std::move(scattered);
    result.incident = IncidentField(wavenumber, points);
    result.total.resize(prepared.count);
    for (std::size_t index = 0; index < prepared.count; ++index)
        result.total[index] = result.incident[index] + result.scattered[index];
    return result;
}

std::vector<Complex> SoftSphereScatteringComplexK(
    Complex wavenumber, double radius, const std::vector<double>& points,
    int terms) {
    if (!Finite(wavenumber) || std::abs(wavenumber) == 0.0)
        throw std::invalid_argument("wavenumber must be finite and nonzero");
    CheckPositive(radius, "radius");
    CheckTerms(terms, false);
    const PreparedPoints prepared = PreparePoints(points, radius, false);
    const auto boundary_j = SphericalJSequence(wavenumber * radius, terms);
    const auto boundary_h = SphericalHankelSequence(wavenumber * radius, terms);
    const auto legendre = PointLegendre(prepared, terms);
    std::vector<std::vector<Complex>> point_h(prepared.count);
    for (std::size_t index = 0; index < prepared.count; ++index)
        point_h[index] = SphericalHankelSequence(
            wavenumber * prepared.radius[index], terms);
    std::vector<Complex> scattered(prepared.count);
    for (int order = 0; order <= terms; ++order) {
        const Complex coefficient =
            -PowI(order) * static_cast<double>(2 * order + 1) *
            boundary_j[static_cast<std::size_t>(order)] /
            boundary_h[static_cast<std::size_t>(order)];
        for (std::size_t index = 0; index < prepared.count; ++index)
            scattered[index] += coefficient *
                point_h[index][static_cast<std::size_t>(order)] *
                legendre[index][static_cast<std::size_t>(order)];
    }
    return scattered;
}

Complex BDFDelta(Complex zeta, const std::string& method) {
    if (!Finite(zeta))
        throw std::invalid_argument("zeta must be finite");
    if (method == "BDF1")
        return 1.0 - zeta;
    if (method == "BDF2")
        return 1.5 - 2.0 * zeta + 0.5 * zeta * zeta;
    throw std::invalid_argument("method must be 'BDF1' or 'BDF2'");
}

CQGridResult BuildCQGrid(
    int sample_count, double time_step, double sound_speed,
    const std::string& method) {
    if (sample_count <= 0 || sample_count > 1048576)
        throw std::invalid_argument("sample_count must be a positive practical integer");
    CheckPositive(time_step, "time_step");
    CheckPositive(sound_speed, "sound_speed");
    BDFDelta(Complex(0.0, 0.0), method);
    CQGridResult result;
    result.radius = std::pow(
        std::numeric_limits<double>::epsilon(), 1.0 / (2.0 * sample_count));
    result.zeta.resize(static_cast<std::size_t>(sample_count));
    result.nodes.resize(static_cast<std::size_t>(sample_count));
    result.wavenumbers.resize(static_cast<std::size_t>(sample_count));
    for (int index = 0; index < sample_count; ++index) {
        const double angle = -2.0 * kPi * index / sample_count;
        const Complex zeta = result.radius * std::exp(Complex(0.0, angle));
        const Complex node = BDFDelta(zeta, method) / time_step;
        result.zeta[static_cast<std::size_t>(index)] = zeta;
        result.nodes[static_cast<std::size_t>(index)] = node;
        result.wavenumbers[static_cast<std::size_t>(index)] =
            Complex(0.0, 1.0) * node / sound_speed;
    }
    return result;
}

}  // namespace radia::acoustics
