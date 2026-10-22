// Host transaction orchestration; no MEX, pybind, CUDA or physics substitution.
#include "ih_coupled_backend.hpp"
#include "ih_coupled_transaction.hpp"
#include <algorithm>
#include <new>
#include <type_traits>

namespace radia_ih {
DomainError::DomainError(ErrorCode code,const char* message)
    :std::runtime_error(message),code_(code) {}
ErrorCode DomainError::code() const noexcept { return code_; }
namespace {
template<class T> bool finite_values(const std::vector<T>& values) {
    return std::all_of(values.begin(),values.end(),[](const auto& v) {
        if constexpr(std::is_same_v<T,std::complex<double>>)
            return std::isfinite(v.real())&&std::isfinite(v.imag());
        else return std::isfinite(v);
    });
}
void controls_checked(const Controls& c) {
    if(!std::isfinite(c.current_A)||!std::isfinite(c.ambient_K)||c.ambient_K<=0.||
       (c.drive_mode!=DriveMode::sinusoidal_periodic&&c.drive_mode!=DriveMode::zero_drive_cooling)||
       (c.drive_mode==DriveMode::sinusoidal_periodic?c.current_A<=0.:c.current_A!=0.))
        throw DomainError(ErrorCode::invalid_input,"Invalid IH drive controls");
}
bool same_controls(const Controls& a,const Controls& b) noexcept {
    return a.current_A==b.current_A&&a.ambient_K==b.ambient_K&&a.drive_mode==b.drive_mode;
}
std::uint64_t snapshot_bytes(const CommittedSnapshot& h) {
    std::uint64_t bytes=sizeof(h);
    const auto add_vector=[&](const auto& values) {
        using T=typename std::decay_t<decltype(values)>::value_type;
        config_detail::add(bytes,config_detail::mul(values.capacity(),sizeof(T)));
    };
    add_vector(h.temperature_coefficients_K);add_vector(h.heat_samples);add_vector(h.eddy_solution);
    return bytes;
}
}
struct CoupledDomain::Impl {
    CoupledConfig config;
    std::unique_ptr<detail::PhysicsBackend> backend;
    std::unique_ptr<detail::PhysicalState> physical;
    std::shared_ptr<const CommittedSnapshot> host;
    DomainStatus status{DomainStatus::initialized};
    std::uint64_t config_payload_bytes{};
    explicit Impl(const CoupledConfig& cfg):config(cfg) {
        config_payload_bytes=validate_config(config).input_bytes;
        backend=detail::make_coupled_physics_backend(config);
        if(!backend)throw DomainError(ErrorCode::unsupported_physics,"Missing compiled IH backend");
        try { (void)validate_config_for_execution(config,backend->capabilities()); }
        catch(const std::runtime_error&) {
            throw DomainError(ErrorCode::unsupported_physics,"IH compiled physical route is unavailable; no fallback");
        }
        auto initial=backend->initialize(config);
        backend->synchronize_trial(initial);
        validate(initial,{0,0,0.},nullptr,true);
        commit(initial);
    }
    void usable() const {
        if(status==DomainStatus::closed)throw DomainError(ErrorCode::closed,"IH domain closed");
        if(status==DomainStatus::poisoned)throw DomainError(ErrorCode::poisoned,"IH domain poisoned");
    }
    void validate(const detail::DomainTrial& trial,detail::NextCounters counters,
                  const Controls* controls,bool resetting) const {
        if(!trial.physical||!trial.host)
            throw DomainError(ErrorCode::invalid_input,"Incomplete IH staged state");
        const auto& h=*trial.host;
        // Charge capacities, not just sizes, before validating array shapes.
        std::uint64_t host_bytes=sizeof(Impl);
        config_detail::add(host_bytes,config_payload_bytes);
        config_detail::add(host_bytes,snapshot_bytes(h));
        if(host)config_detail::add(host_bytes,snapshot_bytes(*host));
        config_detail::add(host_bytes,backend->owned_host_workspace_bytes());
        if(host_bytes>config.budgets.domain_host_bytes)
            throw DomainError(ErrorCode::budget,"IH staged host capacities exceed budget");
        const auto expected_status=resetting?DomainStatus::initialized:DomainStatus::committed;
        if(h.revision!=counters.revision||h.steps!=counters.steps||h.time_s!=counters.time_s||
           h.status!=expected_status||h.diagnostics.solved==resetting||
           h.temperature_coefficients_K.size()!=config.thermal.coefficients||
           h.heat_samples.size()!=config.heat.samples||
           h.eddy_solution.size()!=std::size_t(config.em.unknowns)*config.em.harmonics.size()||
           !finite_values(h.temperature_coefficients_K)||!finite_values(h.heat_samples)||
           !finite_values(h.eddy_solution)||(controls&&!same_controls(h.controls,*controls)))
            throw DomainError(ErrorCode::invalid_input,"Invalid IH staged host snapshot");
        if(resetting) {
            if(h.temperature_coefficients_K!=config.thermal.initial_temperature_coefficients_K||
               std::any_of(h.heat_samples.begin(),h.heat_samples.end(),[](double x){return x!=0.;})||
               std::any_of(h.eddy_solution.begin(),h.eddy_solution.end(),[](auto x){return x!=0.;}))
                throw DomainError(ErrorCode::invalid_input,"IH reset is not the configured cold state");
        }
        const auto& d=h.diagnostics;
        for(double v:{d.eddy.true_relative_residual,d.thermal.true_relative_residual,
                d.constitutive_true_residual,d.constitutive_field_residual,
                d.thermal_nonlinear_residual,d.thermal_energy_relative_defect,
                d.thermal_input_power_W,d.thermal_storage_power_W,d.thermal_boundary_loss_W,
                d.thermal_energy_defect_W,d.port_input_power_W,d.coil_loss_power_W,
                d.integrated_heat_power_W,d.port_heat_defect_W,d.port_heat_relative_defect})
            if(!std::isfinite(v))throw DomainError(ErrorCode::invalid_input,"Nonfinite IH staged report");
        if(!resetting) {
            for(auto pair:{std::pair<double,double>{d.eddy.true_relative_residual,config.em.linear_true_residual_limit},
                    {d.thermal_nonlinear_residual,config.thermal.nonlinear_relative_limit},
                    {d.thermal_energy_relative_defect,config.thermal.energy_relative_limit},
                    {d.port_heat_relative_defect,config.em.port_energy_relative_limit}})
                if(pair.first<0.||pair.first>pair.second)
                    throw DomainError(ErrorCode::energy_gate,"IH staged acceptance report failed its declared limit");
        }
        if(backend->owned_device_bytes()!=0)
            throw DomainError(ErrorCode::unsupported_physics,"Host IH domain cannot own device storage");
        // Numerical residual/power gates belong to the compiled
        // physical backend. A finite snapshot is not physical acceptance.
    }
    void commit(detail::DomainTrial& trial) noexcept {
        physical.swap(trial.physical);
        host.swap(trial.host);
        status=host->status;
        // Retired buffers remain in trial; no destructor/callback in commit.
    }
    template<class Stage> std::shared_ptr<const CommittedSnapshot> execute(
            Stage stage,detail::NextCounters counters,const Controls* controls,bool resetting) {
        try {
            auto trial=stage();
            backend->synchronize_trial(trial);
            validate(trial,counters,controls,resetting);
            commit(trial);
        } catch(const DomainError& e) {
            if(e.code()==ErrorCode::device_failure)status=DomainStatus::poisoned;
            throw;
        } catch(const std::bad_alloc&) {
            throw DomainError(ErrorCode::allocation,"IH trial allocation failed");
        }
        return host; // noexcept shared ownership copy, no publication callback
    }
};
CoupledDomain::CoupledDomain(const CoupledConfig& cfg) {
    (void)validate_config(cfg); // gate budgets before copying owned config
    impl_=std::make_unique<Impl>(cfg);
}
CoupledDomain::~CoupledDomain() noexcept { close(); }
std::shared_ptr<const CommittedSnapshot> CoupledDomain::advance(const Controls& controls) {
    impl_->usable();controls_checked(controls);
    detail::NextCounters next;
    try { next=detail::next_counters(*impl_->host,impl_->config.thermal_macro_dt_s,false); }
    catch(const std::overflow_error&) {throw DomainError(ErrorCode::counter_overflow,"IH counter overflow");}
    return impl_->execute([&]{return impl_->backend->advance(impl_->config,*impl_->physical,*impl_->host,controls);},next,&controls,false);
}
std::shared_ptr<const CommittedSnapshot> CoupledDomain::reset() {
    impl_->usable();detail::NextCounters next;
    try { next=detail::next_counters(*impl_->host,impl_->config.thermal_macro_dt_s,true); }
    catch(const std::overflow_error&) {throw DomainError(ErrorCode::counter_overflow,"IH counter overflow");}
    return impl_->execute([&]{return impl_->backend->reset(impl_->config,*impl_->physical,*impl_->host);},next,nullptr,true);
}
std::shared_ptr<const CommittedSnapshot> CoupledDomain::snapshot() const noexcept {return impl_->host;}
RuntimeInfo CoupledDomain::runtime_info() const noexcept {
    return {1,impl_->status==DomainStatus::closed?0:impl_->backend->owned_device_bytes(),impl_->status};
}
void CoupledDomain::close() noexcept {
    if(!impl_||impl_->status==DomainStatus::closed)return;
    impl_->backend->close();impl_->physical.reset();impl_->status=DomainStatus::closed;
}
} // namespace radia_ih
