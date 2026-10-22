// Same production host orchestrator, test-only compiled backend.
// Exercises transactions, NOT numerical physics, CUDA or language adapters.
#include "ih_coupled_test_fixture.hpp"
#include "../../native/ih_coupled_backend.hpp"
#include "../../native/ih_coupled_runtime.hpp"
#include <iostream>
#include <string>
#include <limits>

namespace {
using namespace radia_ih;
std::string injection;
int synchronized_trials=0;
void require(bool value) { if(!value)throw std::runtime_error("IH scenario failed"); }
void fail(const char* where) {
    if(injection!=where)return;
    if(injection=="device_failure")throw DomainError(ErrorCode::device_failure,"Injected device loss");
    if(injection=="host_snapshot_allocation")throw std::bad_alloc();
    throw DomainError(ErrorCode::convergence,"Injected trial failure");
}
struct TestState final:detail::PhysicalState { double marker{}; };
class TestBackend final:public detail::PhysicsBackend {
    detail::DomainTrial capture(const CoupledConfig& c,std::uint64_t revision,
            std::uint64_t steps,double time,const Controls& controls,bool solved) {
        fail("host_snapshot_allocation");
        auto h=std::make_shared<CommittedSnapshot>();
        h->revision=revision;h->steps=steps;h->time_s=time;h->controls=controls;
        h->status=solved?DomainStatus::committed:DomainStatus::initialized;
        h->temperature_coefficients_K=c.thermal.initial_temperature_coefficients_K;
        for(auto& t:h->temperature_coefficients_K)if(solved)t+=double(steps);
        h->heat_samples.assign(c.heat.samples,solved?2.*double(steps):0.);
        h->eddy_solution.assign(c.em.unknowns*c.em.harmonics.size(),solved?std::complex<double>(3.,4.):0.);
        h->diagnostics.solved=solved;
        h->diagnostics.thermal_input_power_W=solved?7.:0.;
        if(injection=="malformed_snapshot")h->revision+=10;
        if(injection=="oversized_snapshot")
            h->eddy_solution.reserve(c.budgets.domain_host_bytes/sizeof(std::complex<double>)+1);
        return {std::make_unique<TestState>(),h};
    }
public:
    CompiledDomainCapabilities capabilities() const noexcept override {return {true,true,false};}
    std::uint64_t owned_device_bytes() const noexcept override {return 0;}
    std::uint64_t owned_host_workspace_bytes() const noexcept override {
        // Test-only storage. Production must use an actual allocation ledger.
        return sizeof(TestBackend)+2*sizeof(TestState);
    }
    detail::DomainTrial initialize(const CoupledConfig& c) override {return capture(c,0,0,0.,{},false);}
    detail::DomainTrial advance(const CoupledConfig& c,const detail::PhysicalState& state,
            const CommittedSnapshot& accepted,const Controls& controls) override {
        require(static_cast<const TestState&>(state).marker==double(accepted.steps));
        // Allocate and modify a private physical candidate before injection;
        // failures are not merely argument checks before any staging work.
        auto candidate=std::make_unique<TestState>();candidate->marker=double(accepted.steps+1);
        fail("eddy_solve");fail("heat_projection");fail("thermal_solve");
        fail("port_gate");fail("thermal_gate");fail("device_failure");
        auto result=capture(c,accepted.revision+1,accepted.steps+1,
            double(accepted.steps+1)*c.thermal_macro_dt_s,controls,true);
        result.physical=std::move(candidate);return result;
    }
    detail::DomainTrial reset(const CoupledConfig& c,const detail::PhysicalState& state,
            const CommittedSnapshot& accepted) override {
        require(static_cast<const TestState&>(state).marker==double(accepted.steps));
        auto candidate=std::make_unique<TestState>();candidate->marker=0.;
        fail("reset_copy_or_sync");
        auto result=capture(c,accepted.revision+1,0,0.,{},false);
        result.physical=std::move(candidate);return result;
    }
    void synchronize_trial(const detail::DomainTrial&) override {fail("trial_sync");++synchronized_trials;}
    void close() noexcept override {}
};
void scenario(const std::string& name) {
    auto c=ih_contract_fixture::fixture();c.backend=ConfigBackend::independent_cpu;
    CoupledDomain domain(c);
    const Controls drive{3.,300.,DriveMode::sinusoidal_periodic};
    const auto accepted=domain.advance(drive);
    require(accepted->revision==1&&accepted->heat_samples[0]!=0.);
    const auto preserved=*accepted; // independent copy to detect in-place mutation
    if(name=="outputs_repeated") {
        for(int i=0;i<3;++i)require(domain.snapshot()==accepted);
    } else if(name=="reset_success") {
        auto h=domain.reset();require(h->revision==2&&h->steps==0&&h->time_s==0.);
        auto next=domain.advance(drive);require(next->revision==3&&next->steps==1);
    } else if(name=="closed_snapshot") {
        domain.close();domain.close();require(domain.snapshot()==accepted);
        bool rejected=false;try{domain.reset();}catch(const DomainError& e){rejected=e.code()==ErrorCode::closed;}
        require(rejected&&domain.runtime_info().status==DomainStatus::closed);
    } else {
        injection=name;
        const int before_sync=synchronized_trials;
        bool rejected=false;
        try {if(name=="reset_copy_or_sync")domain.reset();else domain.advance(drive);}
        catch(const DomainError& e){
            rejected=true;
            if(name=="oversized_snapshot")require(e.code()==ErrorCode::budget);
        }
        injection.clear();require(rejected&&domain.snapshot()==accepted);
        if(name=="trial_sync")require(synchronized_trials==before_sync);
        if(name=="malformed_snapshot"||name=="oversized_snapshot")
            require(synchronized_trials==before_sync+1);
        require(accepted->temperature_coefficients_K==preserved.temperature_coefficients_K&&
            accepted->heat_samples==preserved.heat_samples&&accepted->eddy_solution==preserved.eddy_solution&&
            accepted->diagnostics.thermal_input_power_W==preserved.diagnostics.thermal_input_power_W&&
            accepted->revision==preserved.revision&&accepted->steps==preserved.steps&&accepted->time_s==preserved.time_s);
        if(name=="device_failure") {
            require(domain.runtime_info().status==DomainStatus::poisoned);
            for(bool resetting:{false,true}) {
                bool blocked=false;try{if(resetting)domain.reset();else domain.advance(drive);}
                catch(const DomainError& e){blocked=e.code()==ErrorCode::poisoned;}
                require(blocked);
            }
            domain.close();require(domain.snapshot()==accepted);
        } else {
            auto h=domain.advance(drive);require(h->revision==2&&h->steps==2);
        }
    }
}
void owner_scenario(const std::string& name) {
    auto c=ih_contract_fixture::fixture();c.backend=ConfigBackend::independent_cpu;
    using RuntimeOwner=radia_ih_runtime::CoupledOwner;
    auto owner=std::make_shared<RuntimeOwner>(std::make_unique<CoupledDomain>(c));
    radia_ih_runtime::Registry<RuntimeOwner,0x49485432u> registry;
    const auto token=registry.insert(owner);
    const Controls drive{3.,300.,DriveMode::sinusoidal_periodic};
    const auto first=registry.lookup(token)->advance(0,1,drive);
    require(first.sample_id==1&&first.epoch==1&&first.snapshot->revision==1);
    const auto retained=*first.snapshot;
    if(name=="owner_idempotency") {
        auto again=owner->advance(0,1,drive);
        require(again.snapshot==first.snapshot&&again.epoch==first.epoch);
        bool rejected=false;
        try {owner->advance(0,1,Controls{4.,300.,DriveMode::sinusoidal_periodic});}
        catch(const std::invalid_argument&){rejected=true;}
        require(rejected&&owner->outputs(0).snapshot==first.snapshot);
    } else if(name=="owner_failed_advance") {
        injection="thermal_solve";bool failed=false;
        try {owner->advance(0,2,drive);}catch(const DomainError&){failed=true;}
        const auto after=owner->outputs(0);
        require(failed&&after.snapshot==first.snapshot&&after.sample_id==1&&after.epoch==1);
        injection.clear();auto next=owner->advance(0,2,drive);
        require(next.snapshot->revision==2&&next.sample_id==2&&next.epoch==2);
    } else if(name=="owner_reset_generation") {
        injection="reset_copy_or_sync";bool failed=false;
        try {owner->reset(0);}catch(const DomainError&){failed=true;}
        require(failed&&owner->outputs(0).generation==0&&owner->outputs(0).snapshot==first.snapshot);
        injection.clear();auto reset=owner->reset(0);
        require(reset.generation==1&&reset.sample_id==0&&reset.epoch==2&&reset.snapshot->steps==0);
        bool stale=false;try{owner->outputs(0);}catch(const std::invalid_argument&){stale=true;}
        require(stale&&owner->advance(1,1,drive).snapshot->revision==3);
    } else if(name=="owner_poison") {
        injection="device_failure";bool failed=false;
        try{owner->advance(0,2,drive);}catch(const DomainError& e){failed=e.code()==ErrorCode::device_failure;}
        const auto after=owner->outputs(0);
        require(failed&&after.poisoned&&after.snapshot==first.snapshot&&after.epoch==1);
        bool rejected=false;try{owner->advance(0,2,drive);}catch(const std::runtime_error&){rejected=true;}
        require(rejected);
    } else throw std::invalid_argument("Unknown joint owner/domain scenario");
    require(first.snapshot->temperature_coefficients_K==retained.temperature_coefficients_K&&
            first.snapshot->heat_samples==retained.heat_samples&&
            first.snapshot->eddy_solution==retained.eddy_solution);
    registry.destroy(token);owner->close();
    bool stale=false;try{registry.lookup(token);}catch(const std::invalid_argument&){stale=true;}
    require(stale);
}

}
namespace radia_ih::detail {
// This symbol exists ONLY in the host contract test executable. It is never
// linked into the production MEX/adapter or advertised as a physical backend.
std::unique_ptr<PhysicsBackend> make_coupled_physics_backend(const CoupledConfig&) {
    return std::make_unique<TestBackend>();
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=2)throw std::invalid_argument("Expected one declared scenario injection");
        const std::string name=argv[1];
        const std::vector<std::string> supported={"eddy_solve","heat_projection","thermal_solve",
            "port_gate","thermal_gate",
            "host_snapshot_allocation","reset_copy_or_sync","device_failure",
            "trial_sync","malformed_snapshot","oversized_snapshot",
            "outputs_repeated","reset_success","closed_snapshot",
            "owner_idempotency","owner_failed_advance","owner_reset_generation","owner_poison"};
        if(std::find(supported.begin(),supported.end(),name)==supported.end())
            throw std::invalid_argument("Unsupported scenario; never report PASS");
        if(name.rfind("owner_",0)==0)owner_scenario(name);else scenario(name);std::cout<<"IH host transaction scenario PASS: "<<name<<"\n";return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
}
