// Complete owned CONFIG PROPOSAL; not a declaration of executable physics.
// Both adapters convert language buffers once; this header owns canonical checks.
#pragma once
#include "ih_periodic_config.hpp"
#include <algorithm>
#include <complex>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace radia_ih {
enum class ConfigBackend { invalid, independent_cpu, gpu };
enum class MagneticMode { invalid, linear_temperature, resolved_periodic_bh };
enum class HeatMeasure { invalid, volume_m3, surface_m2 };
struct ConfigIdentity {
    std::string producer_source_sha1, physics_source_sha1;
    std::string mesh_sha256, label_contract_sha256, assembled_operators_sha256;
    std::string expected_native_artifact_sha256; // required GPU, empty CPU
};
struct MemoryBudgets {
    std::uint64_t input_bytes{}, domain_host_bytes{}, device_bytes{};
    // Domain host budget covers simultaneous accepted/trial host snapshots;
    // caller-owned published NumPy/MATLAB copies are separate caller memory.
};
template<class T> struct OwnedCSR {
    std::uint32_t rows{}, columns{};
    std::vector<std::int32_t> row_offsets, column_indices;
    std::vector<T> values;
};
struct OwnedDense {
    std::uint32_t rows{}, columns{};
    std::vector<double> values; // C row-major
};
struct ThermalConfig {
    std::uint32_t coefficients{}, volume_samples{}, boundary_samples{}, dimension{};
    std::vector<double> temperature_axis_K;
    OwnedDense basis, gradient, boundary_basis;
    // gradient rows are [volume sample, Cartesian component], then coefficient.
    std::vector<double> volume_weights_m3, boundary_weights_m2;
    std::vector<double> density_kg_per_m3, convection_W_per_m2K;
    std::vector<double> constant_coefficients, initial_temperature_coefficients_K;
    OwnedDense cp_J_per_kgK, kappa_W_per_mK, emissivity;
    std::vector<std::uint32_t> volume_element_ids, boundary_element_ids;
    std::vector<double> volume_points_m, boundary_points_m; // [sample,dimension]
    double nonlinear_relative_limit{}, energy_relative_limit{}, power_floor_W{};
    std::uint32_t max_newton_iterations{};
};
struct MagneticMaterialConfig {
    std::uint32_t elements{};
    std::vector<std::uint32_t> element_ids;
    OwnedCSR<double> temperature_evaluation; // [element,thermal coefficient]
    std::vector<double> temperature_axis_K;
    OwnedDense resistivity_ohm_m; // [element,T]
    // Linear route: [element,T], exp(+iwt), passive imag(mu)<=0.
    std::vector<std::complex<double>> linear_permeability_H_per_m;
    // Periodic route ONLY: common positive H knots from0; B[element,T,H].
    // Reversible/isotropic/odd/single-valued, no hysteresis or complex mu.
    std::vector<double> H_axis_A_per_m, B_T;
};
struct HarmonicOperator {
    std::uint32_t harmonic_index{};
    OwnedCSR<std::complex<double>> base_operator; // [EM unknown,EM unknown]
    // Equation: A_h*u_h = coil_rhs_per_A*I_h + boundary_E_to_rhs*E_h.
    // Boundary convention uses outward-normal orientation supplied by producer.
    OwnedCSR<std::complex<double>> boundary_E_to_rhs;
    // H_h = state_to_boundary_H*u_h + coil_boundary_H_per_A*I_h.
    OwnedCSR<std::complex<double>> state_to_boundary_H;
    std::vector<std::complex<double>> coil_rhs_per_A, coil_boundary_H_per_A;
    // Independently assembled reaction voltage, NOT reconstructed from heat.
    std::vector<std::complex<double>> port_voltage_projection;
    double coil_resistance_ohm{};
};
struct ElectromagneticConfig {
    std::uint32_t unknowns{}, field_components{};
    // DOF order/units belong to supplied assembly. No node/field inference.
    std::vector<std::string> unknown_units, equation_units;
    std::vector<HarmonicOperator> harmonics;
    // Complete two-face planar diffusion cells: [element,face,component].
    // One or two tangential components. No silent scalarization of vector data.
    std::uint32_t cells_per_element{};
    std::vector<double> thickness_m, face_area_m2;
    // Initial phase-zero B/H, [element,cell,component]; paired validated maps.
    std::vector<double> initial_cell_B_T, initial_cell_H_A_per_m;
    // Linear loss only: q=.5*sum(loss_coefficient*abs(projected_field)^2).
    // Volume: projected J[A/m2], coefficient rho[ohm*m]; surface: H[A/m],
    // coefficient Re(Zs)[ohm]. Periodic mode forbids this surrogate loss route.
    OwnedCSR<std::complex<double>> linear_loss_projection;
    std::vector<double> linear_loss_coefficient;
    double linear_true_residual_limit{}, port_energy_relative_limit{}, power_floor_W{};
    std::uint32_t max_linear_iterations{};
};
struct HeatTransferConfig {
    HeatMeasure measure{HeatMeasure::invalid};
    std::uint32_t samples{}, components{};
    std::vector<double> measure_weights; // m3 or m2, exactly one measure per sample
    OwnedCSR<double> to_thermal_load; // [thermal coefficient,heat sample]
    // constant_coefficients^T * projection[:,j] == measure_weights[j].
};
struct CoupledConfig {
    static constexpr std::uint32_t schema_version = 2;
    ConfigIdentity identity;
    ConfigBackend backend{ConfigBackend::invalid};
    MagneticMode magnetic_mode{MagneticMode::invalid};
    MemoryBudgets budgets;
    ThermalConfig thermal;
    MagneticMaterialConfig magnetic;
    ElectromagneticConfig em;
    HeatTransferConfig heat;
    PeriodicDriveConfig drive{};
    // No user-set "supported"/fallback/accepted flag. Build capabilities are
    // selected by the compiled domain, separately from this owned payload.
};
struct ConfigValidation {
    std::uint64_t input_bytes{}, host_reserve_bytes{}, device_reserve_bytes{};
    std::uint32_t periods_per_macro{}, boundary_components{};
};
namespace config_detail {
[[noreturn]] inline void bad(const char* message) { throw std::invalid_argument(message); }
inline std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
    if (b && a > std::numeric_limits<std::uint64_t>::max()/b) bad("IH config size overflow");
    return a*b;
}
inline void add(std::uint64_t& a, std::uint64_t b) {
    if (a > std::numeric_limits<std::uint64_t>::max()-b) bad("IH config byte overflow");
    a += b;
}
inline bool finite(double x) { return std::isfinite(x); }
inline bool finite(std::complex<double> x) { return finite(x.real()) && finite(x.imag()); }
template<class T> inline void vector(const std::vector<T>& v, std::uint64_t n,
                                    std::uint64_t& bytes) {
    if (v.size()!=n) bad("IH config vector shape mismatch");
    add(bytes,mul(n,sizeof(T)));
    for (const auto& x:v) if (!finite(x)) bad("Nonfinite IH config data");
}
inline void positive(const std::vector<double>& v) {
    for(double x:v) if (!finite(x)||x<=0.) bad("Positive IH physical values required");
}
inline void axis(const std::vector<double>& v,bool origin=false) {
    if(v.size()<2 || v.size()>1024 || (origin?v.front()!=0.:v.front()<=0.)) bad("Invalid IH material axis");
    for(std::size_t i=0;i<v.size();++i)
        if(!finite(v[i])||(i&&v[i]<=v[i-1])) bad("IH material axis must strictly increase");
}
inline void hex(const std::string& s,std::size_t n,std::uint64_t& bytes) {
    if(s.size()!=n || s.find_first_not_of("0123456789abcdef")!=std::string::npos) bad("Invalid IH identity hash");
    add(bytes,n);
}
inline void labels(const std::vector<std::string>& s,std::size_t n,std::uint64_t& bytes) {
    if(s.size()!=n) bad("IH unknown/unit labels mismatch");
    for(const auto& x:s) {
        if(x.empty()||x.size()>128) bad("IH unit labels must be explicit bounded ASCII");
        for(unsigned char c:x) if(c<32||c>126) bad("IH unit labels require printable ASCII");
        add(bytes,x.size());
    }
}
inline void ids(const std::vector<std::uint32_t>& v,std::size_t n,std::uint64_t& bytes) {
    if(v.size()!=n) bad("IH provenance numbering mismatch");
    for(auto x:v) if(x>std::uint32_t(std::numeric_limits<std::int32_t>::max())) bad("IH element identifier overflow");
    add(bytes,mul(n,sizeof(std::uint32_t)));
}
// Magnetic IDs identify whole physical cell blocks, unlike repeated quadrature
// provenance IDs. Bounded allocation-free search keeps validation scratch zero.
inline void unique_magnetic_ids(const std::vector<std::uint32_t>& v) {
    for(std::size_t i=0;i<v.size();++i)
        if(std::find(v.begin(),v.begin()+i,v[i])!=v.begin()+i)
            bad("Duplicate IH magnetic element identifier");
}
inline void dense(const OwnedDense& a,std::uint32_t r,std::uint32_t c,std::uint64_t& bytes) {
    if(a.rows!=r||a.columns!=c) bad("IH dense map shape mismatch");
    vector(a.values,mul(r,c),bytes);
}
template<class T> inline void csr(const OwnedCSR<T>& a,std::uint32_t r,std::uint32_t c,std::uint64_t& bytes) {
    if(!r||!c||a.rows!=r||a.columns!=c||a.row_offsets.size()!=std::size_t(r)+1||
       a.values.size()>std::size_t(std::numeric_limits<std::int32_t>::max())||a.column_indices.size()!=a.values.size()||
       a.row_offsets.front()!=0||a.row_offsets.back()!=std::int32_t(a.values.size())) bad("Invalid IH CSR shape/offsets");
    for(std::size_t i=0;i<r;++i) {
        if(a.row_offsets[i]<0||a.row_offsets[i+1]<a.row_offsets[i]) bad("Invalid IH CSR offsets");
        std::int32_t previous=-1;
        for(auto k=a.row_offsets[i];k<a.row_offsets[i+1];++k) {
            if(k<0||std::size_t(k)>=a.values.size()) bad("IH CSR offset out of range");
            auto col=a.column_indices[k];
            if(col<0||std::uint32_t(col)>=c||col<=previous) bad("IH CSR columns must be sorted unique zero-based");
            previous=col;
        }
    }
    vector(a.values,a.values.size(),bytes);
    add(bytes,mul(a.row_offsets.size()+a.column_indices.size(),sizeof(std::int32_t)));
}
inline void preserves(const OwnedDense& m,const std::vector<double>& one,double expected) {
    for(std::uint32_t i=0;i<m.rows;++i) {
        double value=0.,scale=std::abs(expected);
        for(std::uint32_t j=0;j<m.columns;++j) {
            const double term=m.values[std::size_t(i)*m.columns+j]*one[j];
            value+=term;scale+=std::abs(term);
        }
        // Only gradient rows carry inverse-length units. Basis/trace rows
        // must still preserve the dimensionless constant with absolute error.
        const double tolerance=expected==0. ? 1e-12*scale : 1e-12;
        if(!finite(value)||!finite(scale)||std::abs(value-expected)>tolerance)bad("IH constant-field map mismatch");
    }
}
inline void limit(double value,double maximum) {
    if(!finite(value)||value<=0.||value>maximum) bad("Invalid IH numerical gate");
}
} // namespace config_detail

// STRUCTURAL/PHYSICAL-INPUT validation only. No CUDA allocation/context, solve,
// parity, mesh accuracy, port-independence proof or release acceptance implied.
inline ConfigValidation validate_config(const CoupledConfig& cfg) {
    using namespace config_detail;
    ConfigValidation result;auto& bytes=result.input_bytes;
    const auto& id=cfg.identity;
    hex(id.producer_source_sha1,40,bytes);hex(id.physics_source_sha1,40,bytes);
    hex(id.mesh_sha256,64,bytes);hex(id.label_contract_sha256,64,bytes);hex(id.assembled_operators_sha256,64,bytes);
    if(cfg.backend==ConfigBackend::gpu)hex(id.expected_native_artifact_sha256,64,bytes);
    else if(cfg.backend!=ConfigBackend::independent_cpu||!id.expected_native_artifact_sha256.empty())bad("Invalid explicit IH backend identity");
    if(!cfg.budgets.input_bytes||cfg.budgets.input_bytes>(2ULL<<30)||
       !cfg.budgets.domain_host_bytes||!cfg.budgets.device_bytes)bad("Invalid explicit IH byte budgets");
    result.periods_per_macro=cfg.drive.checked_periods_per_macro();
    const auto& t=cfg.thermal;const auto n=t.coefficients,q=t.volume_samples,b=t.boundary_samples,d=t.dimension;
    if(!n||n>256||!q||q>4096||!b||b>4096||!d||d>3)bad("Bounded IH thermal dimensions exceeded");
    axis(t.temperature_axis_K);vector(t.temperature_axis_K,t.temperature_axis_K.size(),bytes);
    const auto nt=std::uint32_t(t.temperature_axis_K.size());
    dense(t.basis,q,n,bytes);dense(t.gradient,q*d,n,bytes);dense(t.boundary_basis,b,n,bytes);
    vector(t.volume_weights_m3,q,bytes);positive(t.volume_weights_m3);
    vector(t.boundary_weights_m2,b,bytes);positive(t.boundary_weights_m2);
    vector(t.density_kg_per_m3,q,bytes);positive(t.density_kg_per_m3);
    vector(t.convection_W_per_m2K,b,bytes);for(auto h:t.convection_W_per_m2K)if(h<0.)bad("Negative IH convection");
    vector(t.constant_coefficients,n,bytes);vector(t.initial_temperature_coefficients_K,n,bytes);
    preserves(t.basis,t.constant_coefficients,1.);preserves(t.boundary_basis,t.constant_coefficients,1.);preserves(t.gradient,t.constant_coefficients,0.);
    dense(t.cp_J_per_kgK,q,nt,bytes);positive(t.cp_J_per_kgK.values);
    dense(t.kappa_W_per_mK,q,nt,bytes);positive(t.kappa_W_per_mK.values);
    dense(t.emissivity,b,nt,bytes);for(auto x:t.emissivity.values)if(x<0.||x>1.)bad("Invalid IH emissivity");
    ids(t.volume_element_ids,q,bytes);ids(t.boundary_element_ids,b,bytes);
    vector(t.volume_points_m,mul(q,d),bytes);vector(t.boundary_points_m,mul(b,d),bytes);
    // Coefficients need not be nodal Kelvin values; check their evaluated field.
    for(const auto* map:{&t.basis,&t.boundary_basis})for(std::uint32_t i=0;i<map->rows;++i) {
        double value=0.;for(std::uint32_t j=0;j<n;++j)value+=map->values[std::size_t(i)*n+j]*t.initial_temperature_coefficients_K[j];
        if(!finite(value)||value<t.temperature_axis_K.front()||value>t.temperature_axis_K.back())bad("Initial IH temperature outside table");
    }
    limit(t.nonlinear_relative_limit,1e-6);limit(t.energy_relative_limit,1e-6);limit(t.power_floor_W,1e-12);
    if(!t.max_newton_iterations||t.max_newton_iterations>100)bad("Invalid IH Newton iteration budget");
    const auto& m=cfg.magnetic;const auto e=m.elements;
    if(!e||e>4096)bad("Bounded IH magnetic element count exceeded");
    ids(m.element_ids,e,bytes);unique_magnetic_ids(m.element_ids);
    csr(m.temperature_evaluation,e,n,bytes);
    axis(m.temperature_axis_K);vector(m.temperature_axis_K,m.temperature_axis_K.size(),bytes);
    const auto mt=std::uint32_t(m.temperature_axis_K.size());
    dense(m.resistivity_ohm_m,e,mt,bytes);positive(m.resistivity_ohm_m.values);
    for(std::uint32_t i=0;i<e;++i) {
        double c=0.,temp=0.;for(auto k=m.temperature_evaluation.row_offsets[i];k<m.temperature_evaluation.row_offsets[i+1];++k) {
            const auto j=m.temperature_evaluation.column_indices[k];const auto v=m.temperature_evaluation.values[k];
            c+=v*t.constant_coefficients[j];temp+=v*t.initial_temperature_coefficients_K[j];
        }
        if(!finite(c)||std::abs(c-1.)>1e-12||!finite(temp)||temp<m.temperature_axis_K.front()||temp>m.temperature_axis_K.back())bad("Invalid magnetic temperature map/range");
    }
    const bool periodic=cfg.magnetic_mode==MagneticMode::resolved_periodic_bh;
    if(periodic) {
        if(!m.linear_permeability_H_per_m.empty())bad("No complex/frozen/secant mu in resolved nonlinear BH");
        axis(m.H_axis_A_per_m,true);vector(m.H_axis_A_per_m,m.H_axis_A_per_m.size(),bytes);
        const auto nh=m.H_axis_A_per_m.size();vector(m.B_T,mul(mul(e,mt),nh),bytes);
        for(std::size_t row=0;row<std::size_t(e)*mt;++row)for(std::size_t k=0;k<nh;++k) {
            auto x=m.B_T[row*nh+k];
            if((!k&&x!=0.)||(k&&x<=m.B_T[row*nh+k-1]))bad("BH must increase through origin without hysteresis");
        }
    } else if(cfg.magnetic_mode==MagneticMode::linear_temperature) {
        if(!m.H_axis_A_per_m.empty()||!m.B_T.empty())bad("Nonlinear BH cannot enter linear route");
        vector(m.linear_permeability_H_per_m,mul(e,mt),bytes);
        for(auto mu:m.linear_permeability_H_per_m)if(mu.real()<=0.||mu.imag()>0.)bad("Invalid passive linear permeability");
    } else bad("Unknown IH magnetic mode");
    const auto& em=cfg.em;const auto ne=em.unknowns,components=em.field_components;
    if(!ne||ne>4096||!components||components>2||em.harmonics.empty()||em.harmonics.size()>64)bad("Invalid bounded IH electromagnetic dimensions");
    labels(em.unknown_units,ne,bytes);labels(em.equation_units,ne,bytes);
    const auto nb=std::uint32_t(mul(mul(2,e),components));result.boundary_components=nb;
    for(std::size_t k=0;k<em.harmonics.size();++k) {
        const auto& op=em.harmonics[k];
        if(op.harmonic_index!=(periodic?2*k+1:1)||(!periodic&&em.harmonics.size()!=1)||
           op.harmonic_index>=cfg.drive.micro_steps_per_cycle/2)bad("IH harmonic numbering/Nyquist mismatch");
        csr(op.base_operator,ne,ne,bytes);csr(op.boundary_E_to_rhs,ne,nb,bytes);csr(op.state_to_boundary_H,nb,ne,bytes);
        vector(op.coil_rhs_per_A,ne,bytes);vector(op.coil_boundary_H_per_A,nb,bytes);vector(op.port_voltage_projection,ne,bytes);
        if(!finite(op.coil_resistance_ohm)||op.coil_resistance_ohm<0.)bad("Invalid coil resistance");
        if(k)for(auto value:op.coil_rhs_per_A)if(value!=std::complex<double>{})bad("Sinusoidal source has no imposed higher-harmonic current");
        if(k)for(auto value:op.coil_boundary_H_per_A)if(value!=std::complex<double>{})bad("Sinusoidal source has no imposed higher-harmonic boundary field");
    }
    if(periodic&&em.harmonics.size()<3)bad("Resolved nonlinear proposal requires at least harmonics1/3/5");
    limit(em.linear_true_residual_limit,1e-6);limit(em.port_energy_relative_limit,1e-6);limit(em.power_floor_W,1e-12);
    if(!em.max_linear_iterations||em.max_linear_iterations>10000)bad("Invalid IH EM iteration budget");
    const auto& heat=cfg.heat;
    if(!heat.samples||heat.samples>524288||!heat.components||heat.components>2)bad("Invalid IH heat sample dimensions");
    if(heat.measure!=HeatMeasure::volume_m3&&heat.measure!=HeatMeasure::surface_m2)bad("Explicit IH heat measure required");
    vector(heat.measure_weights,heat.samples,bytes);positive(heat.measure_weights);
    csr(heat.to_thermal_load,n,heat.samples,bytes);
    std::uint64_t cell_values=0;
    if(periodic) {
        if(em.cells_per_element<2||em.cells_per_element>128||heat.measure!=HeatMeasure::volume_m3||
           heat.samples!=mul(e,em.cells_per_element)||heat.components!=components)bad("Invalid IH planar diffusion cell/heat layout");
        vector(em.thickness_m,e,bytes);positive(em.thickness_m);vector(em.face_area_m2,e,bytes);positive(em.face_area_m2);
        cell_values=mul(heat.samples,components);
        vector(em.initial_cell_B_T,cell_values,bytes);vector(em.initial_cell_H_A_per_m,cell_values,bytes);
        // Initial schema supports cold-start only. Warm starts are retained
        // accepted snapshots, not independently supplied inconsistent B/H pairs.
        for(auto x:em.initial_cell_B_T)if(x!=0.)bad("Nonzero initial cycle state is not implemented in this schema");
        for(auto x:em.initial_cell_H_A_per_m)if(x!=0.)bad("Nonzero initial cycle state is not implemented in this schema");
        for(std::uint32_t i=0;i<e;++i)for(std::uint32_t c=0;c<em.cells_per_element;++c) {
            const double vol=em.thickness_m[i]*em.face_area_m2[i]/em.cells_per_element;
            if(!finite(vol)||vol<=0.||std::abs(heat.measure_weights[i*em.cells_per_element+c]-vol)>1e-12*vol)bad("IH cell volume inconsistent with declared slab geometry");
        }
        if(em.linear_loss_projection.rows||em.linear_loss_projection.columns||!em.linear_loss_projection.values.empty()||!em.linear_loss_projection.row_offsets.empty()||!em.linear_loss_projection.column_indices.empty()||!em.linear_loss_coefficient.empty())bad("No quadratic amplitude-loss substitute for nonlinear periodic Joule power");
    } else {
        if(em.cells_per_element||!em.thickness_m.empty()||!em.face_area_m2.empty()||!em.initial_cell_B_T.empty()||!em.initial_cell_H_A_per_m.empty())bad("Unused periodic cell data cannot enter linear config");
        csr(em.linear_loss_projection,heat.samples*heat.components,ne,bytes);
        vector(em.linear_loss_coefficient,mul(heat.samples,heat.components),bytes);
        for(auto x:em.linear_loss_coefficient)if(x<0.)bad("Negative IH loss coefficient");
    }
    // Conservative preallocation reserves: owned input plus accepted/trial state.
    // Actual capacities, validation scratch, solver and preconditioner allocations
    // remain charged by the compiled domain before any allocation; these reserves
    // do not certify the complete device peak or the runtime memory gate.
    auto one_snapshot=std::uint64_t(sizeof(CommittedSnapshot));
    add(one_snapshot,mul(n+heat.samples,sizeof(double)));
    add(one_snapshot,mul(mul(ne,em.harmonics.size()),sizeof(std::complex<double>)));
    if(periodic) {
        // Distinct endpoint B/H, full-cycle B/H and harmonic metadata payloads.
        add(one_snapshot,mul(mul(cell_values,2),sizeof(double)));
        add(one_snapshot,mul(mul(mul(cell_values,cfg.drive.micro_steps_per_cycle),2),sizeof(double)));
        add(one_snapshot,mul(em.harmonics.size(),sizeof(std::uint32_t)));
        add(one_snapshot,mul(mul(cell_values,em.harmonics.size()),sizeof(std::complex<double>)));
    }
    result.host_reserve_bytes=mul(one_snapshot,2);
    add(result.host_reserve_bytes,bytes);
    add(result.host_reserve_bytes,mul(heat.samples,sizeof(double)));
    result.device_reserve_bytes=mul(mul(n+heat.samples+cell_values*2,2),sizeof(double));
    add(result.device_reserve_bytes,mul(mul(mul(ne,em.harmonics.size()),2),sizeof(std::complex<double>)));
    add(result.device_reserve_bytes,bytes);
    if(bytes>cfg.budgets.input_bytes||result.host_reserve_bytes>cfg.budgets.domain_host_bytes||
       result.device_reserve_bytes>cfg.budgets.device_bytes)bad("IH config exceeds declared byte budgets");
    // First and only owned validation workspace allocation; all declared input,
    // host (including this scratch) and device reserves have passed above.
    std::vector<double> sums(heat.samples,0.);
    for(std::uint32_t i=0;i<n;++i)for(auto k=heat.to_thermal_load.row_offsets[i];k<heat.to_thermal_load.row_offsets[i+1];++k)
        sums[heat.to_thermal_load.column_indices[k]]+=t.constant_coefficients[i]*heat.to_thermal_load.values[k];
    for(std::size_t i=0;i<sums.size();++i)if(!finite(sums[i])||std::abs(sums[i]-heat.measure_weights[i])>1e-12*heat.measure_weights[i])bad("IH heat projection must apply physical measure exactly once");
    return result;
}

// Capabilities come ONLY from compiled domain implementation, never from user
// config or adapter-provided flags. Structural validity is not execution support.
struct CompiledDomainCapabilities {
    bool independent_cpu{}, gpu{}, linear_temperature{}, resolved_periodic_bh{}, vector_tangential_field{};
};
inline ConfigValidation validate_config_for_execution(const CoupledConfig& cfg,
        const CompiledDomainCapabilities& compiled) {
    auto result=validate_config(cfg);
    if((cfg.backend==ConfigBackend::gpu?!compiled.gpu:!compiled.independent_cpu)||
       (cfg.magnetic_mode==MagneticMode::resolved_periodic_bh?!compiled.resolved_periodic_bh:!compiled.linear_temperature)||
       (cfg.em.field_components==2&&!compiled.vector_tangential_field))
        throw std::runtime_error("IH physical execution route not implemented in this compiled domain; no fallback");
    return result;
}
} // namespace radia_ih
