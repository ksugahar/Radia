// Internal implementation draft; not wired into a gateway or CUDA execution.
// No adapter/registry ownership and no numerical algorithm live here.
#pragma once
#include "ih_coupled_domain.hpp"
#include <cmath>
#include <limits>
#include <utility>

namespace radia_ih::detail {
struct NextCounters {
    std::uint64_t revision, steps;
    double time_s;
};
inline NextCounters next_counters(const CommittedSnapshot& accepted,
                                 double dt, bool resetting) {
    if (!std::isfinite(dt) || dt <= 0.)
        throw std::invalid_argument("IH sample time must be finite and positive");
    const auto max = std::numeric_limits<std::uint64_t>::max();
    if (accepted.revision == max || (!resetting && accepted.steps == max))
        throw std::overflow_error("IH domain counter overflow before staging");
    const std::uint64_t steps = resetting ? 0 : accepted.steps + 1;
    const double time = static_cast<double>(steps) * dt;
    if (!std::isfinite(time))
        throw std::overflow_error("IH domain time overflow before staging");
    return {accepted.revision + 1, steps, time};
}

// DeviceState must own a COMPLETE accepted physical snapshot, not solver-owned
// scratch. Its destructor releases only owned resources and must not throw.
// Material matrices/workspaces may remain private scratch if every new trial
// reconstructs them from immutable config and accepted physical state.
template<class DeviceState>
class AcceptedStorage final {
public:
    struct Staged {
        std::unique_ptr<DeviceState> device;
        std::shared_ptr<const CommittedSnapshot> host;
        // After commit these hold the retired accepted state. Destroy them only
        // OUTSIDE commit, with domain-specific device-error handling in place.
    };
    explicit AcceptedStorage(Staged initial)
        : device_(std::move(initial.device)), host_(std::move(initial.host)) {
        if (!device_ || !host_)
            throw std::invalid_argument("IH initial snapshot is incomplete");
    }
    AcceptedStorage(const AcceptedStorage&) = delete;
    AcceptedStorage& operator=(const AcceptedStorage&) = delete;
    AcceptedStorage(AcceptedStorage&&) = delete;
    AcceptedStorage& operator=(AcceptedStorage&&) = delete;

    const DeviceState& device() const noexcept { return *device_; }
    std::shared_ptr<const CommittedSnapshot> snapshot() const noexcept {
        return host_;
    }
    // Validate BEFORE final commit. The domain fills host vectors/reports and
    // synchronizes ALL CUDA streams before calling this host-only checker.
    void check_staged(const Staged& trial, NextCounters counters) const {
        if (!trial.device || !trial.host)
            throw std::invalid_argument("IH trial snapshot is incomplete");
        const auto& h = *trial.host;
        if (h.revision != counters.revision || h.steps != counters.steps ||
            h.time_s != counters.time_s ||
            h.temperature_coefficients_K.size() != host_->temperature_coefficients_K.size() ||
            h.heat_samples.size() != host_->heat_samples.size() ||
            h.eddy_solution.size() != host_->eddy_solution.size())
            throw std::invalid_argument("IH staged snapshot identity/shape mismatch");
        // Physics gates, finite arrays, units/maps and status/solved semantics
        // belong to the domain staging path, not this pointer-storage helper.
    }
    // Preconditions: lock held; owner usable; check_staged/physics gates and host
    // capture succeeded; no callback/CUDA/allocator may run between here and the
    // runtime's preallocated metadata switch. This routine destroys NOTHING.
    void commit(Staged& trial) noexcept {
        device_.swap(trial.device);
        host_.swap(trial.host);
    }
private:
    std::unique_ptr<DeviceState> device_;
    std::shared_ptr<const CommittedSnapshot> host_;
};
} // namespace radia_ih::detail
