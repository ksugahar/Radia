// rad_particle_trajectory.h - STUB (particle trajectory removed 2026-03-22)
// Use CERN Xsuite/Xtrack for GPU-accelerated beam tracking.
// This stub exists only to satisfy legacy references in rad_transform_impl.cpp.
// TODO: Remove all radTPrtclTrj references from rad_transform_impl.cpp and delete this file.

#ifndef RAD_PARTICLE_TRAJECTORY_H
#define RAD_PARTICLE_TRAJECTORY_H

#include <stdexcept>

class radTPrtclTrj {
public:
    template<typename... Args>
    radTPrtclTrj(Args&&...) {
        throw std::runtime_error("Particle trajectory removed. Use CERN Xsuite/Xtrack.");
    }
    template<typename... Args>
    void Tabulation(Args&&...) { throw std::runtime_error("Removed"); }
    template<typename... Args>
    double FocusingPotential(Args&&...) { throw std::runtime_error("Removed"); return 0; }
    template<typename... Args>
    void ComputeSecondOrderKickPer(Args&&...) { throw std::runtime_error("Removed"); }
    template<typename... Args>
    void ComputeSecondOrderKick(Args&&...) { throw std::runtime_error("Removed"); }
    template<typename... Args>
    void ComposeStrReportSecondOrderKick(Args&&...) { throw std::runtime_error("Removed"); }
    template<typename... Args>
    static void ComposeStrReportSecondOrderKickPer(Args&&...) { throw std::runtime_error("Removed"); }
};

#endif
