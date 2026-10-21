// Host-only config boundary tests. No GPU, domain solve or physics acceptance.
#include "../../native/ih_coupled_config.hpp"
#include <cstdlib>
#include <iostream>
#include <new>
#include <cstring>

namespace {
bool watch_large_allocations=false;
std::size_t large_allocations=0;
constexpr std::size_t scratch_threshold=2048;
}
void* operator new(std::size_t n) {
    if(watch_large_allocations&&n>=scratch_threshold)++large_allocations;
    if(void* p=std::malloc(n?n:1))return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p,std::size_t) noexcept { std::free(p); }

using namespace radia_ih;
void require(bool value) { if(!value)throw std::runtime_error("IH config contract failed"); }
OwnedDense dense(std::uint32_t r,std::uint32_t c,double value) {
    return {r,c,std::vector<double>(std::size_t(r)*c,value)};
}
template<class T> OwnedCSR<T> empty_matrix(std::uint32_t r,std::uint32_t c) {
    return {r,c,std::vector<std::int32_t>(std::size_t(r)+1,0),{}, {}};
}
CoupledConfig fixture() {
    CoupledConfig c{};
    c.identity={std::string(40,'a'),std::string(40,'b'),std::string(64,'c'),
                std::string(64,'d'),std::string(64,'e')};
    c.backend=ConfigBackend::independent_cpu;c.magnetic_mode=MagneticMode::linear_temperature;
    c.budgets={16ULL<<20,16ULL<<20};
    auto& t=c.thermal;t.coefficients=t.volume_samples=t.boundary_samples=t.dimension=1;
    t.temperature_axis_K={300.,500.};t.basis=dense(1,1,1.);
    t.gradient=dense(1,1,0.);t.boundary_basis=dense(1,1,1.);
    t.volume_weights_m3={2.};t.boundary_weights_m2={2.};
    t.density_kg_per_m3={1.};t.convection_W_per_m2K={0.};
    t.constant_coefficients={1.};t.initial_temperature_coefficients_K={300.};
    t.cp_J_per_kgK=dense(1,2,1.);t.kappa_W_per_mK=dense(1,2,1.);
    t.emissivity=dense(1,2,0.);t.volume_element_ids={0};t.boundary_element_ids={0};
    t.volume_points_m={0.};t.boundary_points_m={0.};
    t.nonlinear_relative_limit=1e-10;t.energy_relative_limit=1e-6;
    t.power_floor_W=1e-12;t.max_newton_iterations=30;
    auto& m=c.magnetic;m.elements=2;m.element_ids={10,11};
    m.temperature_evaluation={2,1,{0,1,2},{0,0},{1.,1.}};
    m.temperature_axis_K={300.,500.};m.resistivity_ohm_m=dense(2,2,1.);
    m.linear_permeability_H_per_m.assign(4,{1.,0.});
    auto& e=c.em;e.unknowns=e.field_components=1;
    e.unknown_units={"A"};e.equation_units={"V"};
    for(std::uint32_t h:{1u}) {
        HarmonicOperator op{};op.harmonic_index=h;
        op.base_operator={1,1,{0,1},{0},{{1.,0.}}};
        op.boundary_E_to_rhs=empty_matrix<std::complex<double>>(1,4);
        op.state_to_boundary_H=empty_matrix<std::complex<double>>(4,1);
        op.coil_rhs_per_A={{h==1?1.:0.,0.}};
        op.coil_boundary_H_per_A.resize(4);
        op.port_voltage_projection={{1.,0.}};op.coil_resistance_ohm=0.;
        e.harmonics.push_back(std::move(op));
    }
    e.linear_true_residual_limit=1e-6;e.port_energy_relative_limit=1e-6;
    e.power_floor_W=1e-12;e.max_linear_iterations=100;
    auto& heat=c.heat;heat.measure=HeatMeasure::volume_m3;heat.samples=256;heat.components=1;
    heat.measure_weights.assign(256,1./128);
    heat.to_thermal_load=empty_matrix<double>(1,256);
    heat.to_thermal_load.row_offsets={0,256};
    for(std::int32_t i=0;i<256;++i)heat.to_thermal_load.column_indices.push_back(i);
    heat.to_thermal_load.values=heat.measure_weights;
    e.linear_loss_projection=empty_matrix<std::complex<double>>(256,1);
    e.linear_loss_coefficient.assign(256,1.);
    c.frequency_Hz=1.;c.thermal_macro_dt_s=1.;
    return c;
}
void rejects_without_scratch(CoupledConfig& c,const char* expected) {
    bool rejected=false;large_allocations=0;watch_large_allocations=true;
    try { (void)validate_config(c); }
    catch(const std::invalid_argument& error) {
        watch_large_allocations=false;
        rejected=std::strstr(error.what(),expected)!=nullptr;
    }
    catch(...) { watch_large_allocations=false;throw; }
    watch_large_allocations=false;
    require(rejected&&large_allocations==0);
}
template<class T> std::uint64_t payload(const std::vector<T>& v) {
    return std::uint64_t(v.size())*sizeof(T);
}
int main() {
    // Cancellation scales with the map's physical gradient units, not an
    // absolute inverse-length threshold; a relative constant-field leak fails.
    auto scaled=dense(1,2,0.);scaled.values={1.e8,-1.e8+1.e-5};
    config_detail::preserves(scaled,{1.,1.},0.);
    scaled.values[1]=-1.e8+1.;
    bool leak_rejected=false;
    try { config_detail::preserves(scaled,{1.,1.},0.); }
    catch(const std::invalid_argument&) { leak_rejected=true; }
    require(leak_rejected);
    scaled.values={1.e12,-1.e12};
    bool constant_rejected=false;
    try { config_detail::preserves(scaled,{1.,1.},1.); }
    catch(const std::invalid_argument&) { constant_rejected=true; }
    require(constant_rejected); // Relative gradient rule must not weaken E/B.
    auto valid=fixture();large_allocations=0;watch_large_allocations=true;
    const auto info=validate_config(valid);watch_large_allocations=false;
    require(large_allocations==1); // Sanity: probe sees real projection scratch.
    auto tiny_input=valid;tiny_input.budgets.input_bytes=1;
    rejects_without_scratch(tiny_input,"exceeds declared byte budgets");
    auto tiny_host=valid;tiny_host.budgets.domain_host_bytes=1;
    rejects_without_scratch(tiny_host,"exceeds declared byte budgets");
    auto duplicate=valid;duplicate.magnetic.element_ids[1]=duplicate.magnetic.element_ids[0];
    rejects_without_scratch(duplicate,"Duplicate IH magnetic element identifier");
    auto repeated=valid;auto& t=repeated.thermal;
    t.volume_samples=t.boundary_samples=2;t.basis=dense(2,1,1.);
    t.gradient=dense(2,1,0.);t.boundary_basis=dense(2,1,1.);
    t.volume_weights_m3=t.boundary_weights_m2={1.,1.};
    t.density_kg_per_m3={1.,1.};t.convection_W_per_m2K={0.,0.};
    t.cp_J_per_kgK=t.kappa_W_per_mK=dense(2,2,1.);t.emissivity=dense(2,2,0.);
    t.volume_element_ids=t.boundary_element_ids={0,0};
    t.volume_points_m=t.boundary_points_m={0.,1.};
    (void)validate_config(repeated); // Repeated QUADRATURE provenance is valid.

    auto invalid_mode=valid;invalid_mode.magnetic_mode=MagneticMode::invalid;
    rejects_without_scratch(invalid_mode,"Explicit linear temperature");
    auto invalid_backend=valid;invalid_backend.backend=ConfigBackend::invalid;
    rejects_without_scratch(invalid_backend,"Explicit independent CPU");
    auto one_short=valid;one_short.budgets.domain_host_bytes=info.host_reserve_bytes-1;
    rejects_without_scratch(one_short,"exceeds declared byte budgets");
    std::cout<<"IH CPU config budget/ID contracts PASS (host only; no physics)\n";
}
