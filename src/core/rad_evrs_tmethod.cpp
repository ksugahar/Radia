#include "rad_evrs_tmethod.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace radia {
namespace evrs {
namespace {

void CheckSize(const std::vector<double>& a, int rows, int cols, const char* name) {
    if (rows < 0 || cols < 0)
        throw std::runtime_error(std::string(name) + ": negative dimension");
    const auto expected = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
    if (a.size() != expected)
        throw std::runtime_error(
            std::string(name) + ": expected " + std::to_string(rows) + " x " +
            std::to_string(cols) + " entries, got " + std::to_string(a.size()));
}

std::vector<double> MatMul(
        const std::vector<double>& a, int ar, int ac,
        const std::vector<double>& b, int br, int bc,
        const char* name) {
    if (ac != br)
        throw std::runtime_error(std::string(name) + ": inner dimension mismatch");
    std::vector<double> out(static_cast<std::size_t>(ar) * static_cast<std::size_t>(bc), 0.0);
    for (int i = 0; i < ar; ++i) {
        for (int k = 0; k < ac; ++k) {
            const double aik = a[static_cast<std::size_t>(i) * ac + k];
            if (aik == 0.0) continue;
            for (int j = 0; j < bc; ++j)
                out[static_cast<std::size_t>(i) * bc + j] +=
                    aik * b[static_cast<std::size_t>(k) * bc + j];
        }
    }
    return out;
}

std::vector<double> Transpose(const std::vector<double>& a, int rows, int cols) {
    std::vector<double> out(static_cast<std::size_t>(cols) * static_cast<std::size_t>(rows), 0.0);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            out[static_cast<std::size_t>(j) * rows + i] =
                a[static_cast<std::size_t>(i) * cols + j];
    return out;
}

std::vector<double> Subtract(const std::vector<double>& a, const std::vector<double>& b, const char* name) {
    if (a.size() != b.size())
        throw std::runtime_error(std::string(name) + ": size mismatch");
    std::vector<double> out(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i)
        out[i] = a[i] - b[i];
    return out;
}

double FrobeniusNorm(const std::vector<double>& a) {
    double sum = 0.0;
    for (double x : a)
        sum += x * x;
    return std::sqrt(sum);
}

double SymmetryNorm(const std::vector<double>& a, int n) {
    return FrobeniusNorm(Subtract(a, Transpose(a, n, n), "symmetry"));
}

std::vector<double> TripleProduct(
        const std::vector<double>& left, int lr, int lc,
        const std::vector<double>& middle, int mr, int mc,
        const std::vector<double>& right, int rr, int rc,
        const char* name) {
    auto tmp = MatMul(middle, mr, mc, right, rr, rc, name);
    return MatMul(left, lr, lc, tmp, mr, rc, name);
}

} // namespace

TMethodAlgebraResult BuildTMethodAlgebra(
    const std::vector<double>& curl_map, int n_current, int n_t,
    const std::vector<double>& div_map, int n_rho, int div_cols,
    const std::vector<double>& grad_map, int grad_rows, int n_phi,
    const std::vector<double>& evrs_map, int evrs_rows, int n_evrs,
    const std::vector<double>& resistance_current, int resistance_rows, int resistance_cols,
    const std::vector<double>& inductance_current, int inductance_rows, int inductance_cols,
    const std::vector<double>& port_current, int port_rows, int n_ports) {

    CheckSize(curl_map, n_current, n_t, "curl_map");
    CheckSize(div_map, n_rho, div_cols, "div_map");
    CheckSize(grad_map, grad_rows, n_phi, "grad_map");
    CheckSize(evrs_map, evrs_rows, n_evrs, "evrs_map");
    CheckSize(resistance_current, resistance_rows, resistance_cols, "resistance_current");
    CheckSize(inductance_current, inductance_rows, inductance_cols, "inductance_current");
    CheckSize(port_current, port_rows, n_ports, "port_current");

    if (div_cols != n_current)
        throw std::runtime_error("div_map columns must match curl_map rows");
    if (grad_rows != n_t)
        throw std::runtime_error("grad_map rows must match curl_map columns");
    if (evrs_rows != n_t)
        throw std::runtime_error("evrs_map rows must match curl_map columns");
    if (resistance_rows != n_current || resistance_cols != n_current)
        throw std::runtime_error("resistance_current must be square in current space");
    if (inductance_rows != n_current || inductance_cols != n_current)
        throw std::runtime_error("inductance_current must be square in current space");
    if (port_rows != n_current)
        throw std::runtime_error("port_current rows must match current-space dimension");

    TMethodAlgebraResult result;
    result.n_current = n_current;
    result.n_t = n_t;
    result.n_phi = n_phi;
    result.n_evrs = n_evrs;
    result.n_ports = n_ports;
    result.n_rho = n_rho;

    const auto cT = Transpose(curl_map, n_current, n_t);
    const auto qT = Transpose(evrs_map, n_t, n_evrs);
    const auto gT = Transpose(grad_map, n_t, n_phi);

    result.current_evrs = MatMul(curl_map, n_current, n_t, evrs_map, n_t, n_evrs, "C Q");

    const auto mrC = MatMul(resistance_current, n_current, n_current,
                            curl_map, n_current, n_t, "M_R C");
    const auto mlC = MatMul(inductance_current, n_current, n_current,
                            curl_map, n_current, n_t, "M_L C");
    result.resistance_t = MatMul(cT, n_t, n_current, mrC, n_current, n_t, "C^T M_R C");
    result.inductance_t = MatMul(cT, n_t, n_current, mlC, n_current, n_t, "C^T M_L C");

    result.resistance_evrs = TripleProduct(qT, n_evrs, n_t,
                                           result.resistance_t, n_t, n_t,
                                           evrs_map, n_t, n_evrs,
                                           "Q^T R_T Q");
    result.inductance_evrs = TripleProduct(qT, n_evrs, n_t,
                                           result.inductance_t, n_t, n_t,
                                           evrs_map, n_t, n_evrs,
                                           "Q^T L_T Q");

    result.port_t = MatMul(cT, n_t, n_current,
                           port_current, n_current, n_ports, "C^T P");
    result.port_evrs = MatMul(qT, n_evrs, n_t,
                              result.port_t, n_t, n_ports, "Q^T C^T P");

    const auto dc = MatMul(div_map, n_rho, n_current, curl_map, n_current, n_t, "D C");
    const auto dcq = MatMul(div_map, n_rho, n_current,
                            result.current_evrs, n_current, n_evrs, "D C Q");
    const auto rtG = MatMul(result.resistance_t, n_t, n_t, grad_map, n_t, n_phi, "R_T G");
    const auto ltG = MatMul(result.inductance_t, n_t, n_t, grad_map, n_t, n_phi, "L_T G");
    const auto gtPort = MatMul(gT, n_phi, n_t, result.port_t, n_t, n_ports, "G^T C^T P");

    const auto jqT = Transpose(result.current_evrs, n_current, n_evrs);
    const auto rFromCurrent = TripleProduct(jqT, n_evrs, n_current,
                                            resistance_current, n_current, n_current,
                                            result.current_evrs, n_current, n_evrs,
                                            "(CQ)^T M_R (CQ)");
    const auto lFromCurrent = TripleProduct(jqT, n_evrs, n_current,
                                            inductance_current, n_current, n_current,
                                            result.current_evrs, n_current, n_evrs,
                                            "(CQ)^T M_L (CQ)");

    result.div_curl_norm = FrobeniusNorm(dc);
    result.div_evrs_norm = FrobeniusNorm(dcq);
    result.resistance_gauge_norm = FrobeniusNorm(rtG);
    result.inductance_gauge_norm = FrobeniusNorm(ltG);
    result.port_gauge_norm = FrobeniusNorm(gtPort);
    result.resistance_symmetry_norm = SymmetryNorm(result.resistance_t, n_t);
    result.inductance_symmetry_norm = SymmetryNorm(result.inductance_t, n_t);
    result.evrs_resistance_symmetry_norm = SymmetryNorm(result.resistance_evrs, n_evrs);
    result.evrs_inductance_symmetry_norm = SymmetryNorm(result.inductance_evrs, n_evrs);
    result.evrs_resistance_galerkin_residual =
        FrobeniusNorm(Subtract(result.resistance_evrs, rFromCurrent, "EVRS resistance identity"));
    result.evrs_inductance_galerkin_residual =
        FrobeniusNorm(Subtract(result.inductance_evrs, lFromCurrent, "EVRS inductance identity"));

    return result;
}

} // namespace evrs
} // namespace radia
