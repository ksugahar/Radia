// Internal compiled physics seam. Never exposed as an adapter/config option.
#pragma once
#include "ih_coupled_config.hpp"
#include "ih_coupled_domain.hpp"

namespace radia_ih::detail {
struct PhysicalState {
    virtual ~PhysicalState() noexcept = default;
};
struct DomainTrial {
    std::unique_ptr<PhysicalState> physical;
    std::shared_ptr<const CommittedSnapshot> host;
};
class PhysicsBackend {
public:
    virtual ~PhysicsBackend() noexcept = default;
    virtual CompiledDomainCapabilities capabilities() const noexcept = 0;
    virtual std::uint64_t owned_device_bytes() const noexcept = 0;
    // Host allocations outside canonical config/snapshot objects+vector
    // capacities (including control blocks/physics workspace) must be charged.
    // Real backends must preflight the simultaneous accepted/trial peak BEFORE
    // allocation; the orchestrator's post-stage guard is an additional check.
    virtual std::uint64_t owned_host_workspace_bytes() const noexcept = 0;
    // Every returned trial owns complete physical storage. Stage operations
    // must not mutate the accepted state or publish language objects.
    virtual DomainTrial initialize(const CoupledConfig&) = 0;
    virtual DomainTrial advance(const CoupledConfig&, const PhysicalState&,
        const CommittedSnapshot&, const Controls&) = 0;
    virtual DomainTrial reset(const CoupledConfig&, const PhysicalState&,
        const CommittedSnapshot&) = 0;
    // All device work/errors must be resolved BEFORE host validation/commit.
    virtual void synchronize_trial(const DomainTrial&) = 0;
    virtual void close() noexcept = 0;
};
// Defined by the compiled physical implementation, never by an adapter.
// The production definition rejects unsupported execution until a canonical
// physical provider is implemented and qualified. It never links test physics.
// Transaction scenario tests link a separate test-only definition instead.
std::unique_ptr<PhysicsBackend> make_coupled_physics_backend(const CoupledConfig&);
} // namespace radia_ih::detail
