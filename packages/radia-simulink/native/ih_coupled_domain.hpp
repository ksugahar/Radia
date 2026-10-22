// Host transaction interface; canonical physical execution remains unavailable.
#pragma once
#include <complex>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace radia_ih {
// Typed immutable config definition is frozen jointly before implementation.
// It owns canonical SI/row-major/zero-based inputs; never mxArray/PyObject data.
struct CoupledConfig;
enum class DriveMode { sinusoidal_periodic, zero_drive_cooling };
struct Controls {
    double current_A; // positive peak amplitude in sine mode; exactly zero in cooling
    double ambient_K;
    DriveMode drive_mode{DriveMode::sinusoidal_periodic};
};
enum class DomainStatus { initialized, committed, poisoned, closed };
enum class ErrorCode { invalid_input, unsupported_physics, budget,
    convergence, energy_gate, allocation, device_failure, poisoned, closed,
    counter_overflow };
class DomainError : public std::runtime_error {
public:
    DomainError(ErrorCode code, const char* message);
    ErrorCode code() const noexcept;
private:
    ErrorCode code_;
};
struct SolveReport {
    std::uint32_t iterations{};
    double true_relative_residual{};
};
struct Diagnostics {
    SolveReport eddy, thermal;
    std::uint32_t constitutive_iterations{};
    double constitutive_true_residual{}, constitutive_field_residual{};
    double thermal_nonlinear_residual{}, thermal_energy_relative_defect{};
    double thermal_input_power_W{}, thermal_storage_power_W{};
    double thermal_boundary_loss_W{}, thermal_energy_defect_W{};
    double port_input_power_W{}, coil_loss_power_W{}, integrated_heat_power_W{};
    double port_heat_defect_W{}, port_heat_relative_defect{};
    // Reports are valid only after a solved commit, never initial/reset zeros.
    bool solved{};
};
struct CommittedSnapshot {
    // Revision starts at zero; increases on advance AND reset; rejects overflow.
    // Registry generation/sample_id/epoch are owned by the common runtime layer.
    std::uint64_t revision{}, steps{};
    double time_s{};
    Controls controls{};
    DomainStatus status{DomainStatus::initialized};
    std::vector<double> temperature_coefficients_K;
    // Units/numbering of heat and EM unknown blocks are immutable config metadata.
    std::vector<double> heat_samples;
    std::vector<std::complex<double>> eddy_solution;
    Diagnostics diagnostics;
};
struct RuntimeInfo {
    std::uint32_t interface_version{1};
    std::uint64_t owned_device_bytes{};
    DomainStatus status{DomainStatus::initialized};
    // Build/source/backend identity supplied by the checked artifact/config;
    // identity never inferred from successful loading or a backend name alone.
};
class CoupledDomain final {
public:
    explicit CoupledDomain(const CoupledConfig& config);
    ~CoupledDomain() noexcept;
    CoupledDomain(const CoupledDomain&) = delete;
    CoupledDomain& operator=(const CoupledDomain&) = delete;
    CoupledDomain(CoupledDomain&&) = delete;
    CoupledDomain& operator=(CoupledDomain&&) = delete;

    // Caller holds the per-owner mutex and validates generation/sample controls.
    // All allocation, GPU work, synchronization, gates AND host snapshot capture
    // happen on trial storage first. Final accepted pointer/report switch is
    // noexcept. Returned shared_ptr copy allocates nothing after commit.
    std::shared_ptr<const CommittedSnapshot> advance(const Controls& controls);
    std::shared_ptr<const CommittedSnapshot> reset();
    // Retained immutable HOST snapshot: no CUDA operations and no allocation.
    // Nonnull after successful construction, including after close/poison.
    // It remains the last coherent snapshot even when runtime_info is poisoned.
    std::shared_ptr<const CommittedSnapshot> snapshot() const noexcept;
    RuntimeInfo runtime_info() const noexcept;
    // Caller serializes; frees only this owner. Idempotent, no device recovery.
    // Retain the host snapshot; runtime_info becomes closed. Already returned
    // immutable snapshots remain valid and contain no borrowed device memory.
    void close() noexcept;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace radia_ih
