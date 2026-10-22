// Language-neutral binding of the jointly frozen domain interface to ownership.
#pragma once
#include "ih_coupled_domain.hpp"
#include "ih_coupled_registry.hpp"
#include <cmath>

namespace radia_ih_runtime {
struct CoupledPolicy final {
    static void validate(const radia_ih::Controls& controls) {
        using radia_ih::DriveMode;
        if (!std::isfinite(controls.current_A) ||
            !std::isfinite(controls.ambient_K) || controls.ambient_K <= 0.)
            throw std::invalid_argument("IH controls must be finite SI scalars");
        switch (controls.drive_mode) {
        case DriveMode::sinusoidal_periodic:
            if (controls.current_A <= 0.)
                throw std::invalid_argument("IH sinusoidal peak current must be positive");
            break;
        case DriveMode::zero_drive_cooling:
            if (controls.current_A != 0.)
                throw std::invalid_argument("IH cooling current must be zero");
            break;
        default:
            throw std::invalid_argument("Unsupported IH drive mode");
        }
    }
    static bool equal(const radia_ih::Controls& left,
                      const radia_ih::Controls& right) noexcept {
        return left.current_A == right.current_A &&
               left.ambient_K == right.ambient_K &&
               left.drive_mode == right.drive_mode;
    }
    static bool poisoned(const radia_ih::CoupledDomain& domain) noexcept {
        return domain.runtime_info().status == radia_ih::DomainStatus::poisoned;
    }
};
using CoupledOwner = Owner<radia_ih::CoupledDomain, radia_ih::Controls, CoupledPolicy>;
// Each runtime supplies its own registry tag and owns config conversion,
// GIL release, checked token publication, and copies into language-owned arrays.
// This header neither creates CUDA state nor supplies an alternate domain.
} // namespace radia_ih_runtime
