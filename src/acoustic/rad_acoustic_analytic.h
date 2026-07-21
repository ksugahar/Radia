#pragma once

#include <complex>
#include <cstdint>
#include <string>
#include <vector>

namespace radia::acoustics {

using Complex = std::complex<double>;

struct ScatteringResult {
    std::string kind;
    double wavenumber = 0.0;
    double radius = 0.0;
    double interior_wavenumber = 0.0;
    double density_ratio = 0.0;
    double longitudinal_speed = 0.0;
    double shear_speed = 0.0;
    int terms = 0;
    double truncation_tail = 0.0;
    std::vector<Complex> scattered;
    std::vector<Complex> incident;
    std::vector<Complex> total;
    std::vector<std::uint8_t> inside_mask;
};

struct CQGridResult {
    double radius = 0.0;
    std::vector<Complex> zeta;
    std::vector<Complex> nodes;
    std::vector<Complex> wavenumbers;
};

ScatteringResult SoftSphereScattering(
    double wavenumber, double radius, const std::vector<double>& points,
    int terms = -1);

ScatteringResult RigidSphereScattering(
    double wavenumber, double radius, const std::vector<double>& points,
    int terms = 1);

ScatteringResult FluidSphereScattering(
    double wavenumber, double radius, const std::vector<double>& points,
    double interior_wavenumber, double density_ratio, int terms = -1);

ScatteringResult ElasticSphereScattering(
    double wavenumber, double radius, const std::vector<double>& points,
    double longitudinal_speed, double shear_speed, double density_ratio,
    int terms = 0);

std::vector<Complex> SoftSphereScatteringComplexK(
    Complex wavenumber, double radius, const std::vector<double>& points,
    int terms = 28);

Complex BDFDelta(Complex zeta, const std::string& method);

CQGridResult BuildCQGrid(
    int sample_count, double time_step, double sound_speed,
    const std::string& method = "BDF2");

}  // namespace radia::acoustics
