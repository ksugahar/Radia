#pragma once

#include <cstddef>
#include <vector>

namespace rad_hdiv {

// Isotropic vector B-input Stop material for HDiv-VIM.  The committed state
// is [s_0, ..., s_{K-1}, B_last], with three components per vector.  Trial
// evaluation is pure; CommitBatch is the only operation that advances state.
class EnergyStopMaterial {
public:
    EnergyStopMaterial(std::vector<double> eta,
                       std::vector<double> table_r,
                       std::vector<double> table_g,
                       std::vector<int> table_offsets,
                       std::vector<double> gamma,
                       double alpha,
                       double b_max);

    std::size_t BranchCount() const noexcept { return tables_.size(); }
    std::size_t StateSize() const noexcept { return 3 * (tables_.size() + 1); }
    double Alpha() const noexcept { return alpha_; }
    double BMax() const noexcept { return b_max_; }
    double NuBound() const noexcept { return nu_bound_; }
    const std::vector<double>& Gamma() const noexcept { return gamma_; }
    const std::vector<double>& Eta() const noexcept { return eta_; }

    void State0(double* output) const;
    void ForwardBatch(const double* B, const double* states, int count,
                      double* H) const;
    void CommitBatch(const double* B, const double* states, int count,
                     double* new_states) const;
    void StoredEnergyBatch(const double* B, const double* states, int count,
                           double* energy) const;

private:
    struct Table {
        std::vector<double> r;
        std::vector<double> g;
        std::vector<double> u;
        double max_slope = 0.0;
        double max_tangent = 0.0;
    };

    struct BranchTrial {
        double s[3];
        double radius;
        double g;
        double u;
    };

    BranchTrial TrialBranch(std::size_t branch, const double* state,
                            const double* delta_B) const;
    void Interpolate(const Table& table, double radius, double& g,
                     double& u) const;

    std::vector<double> eta_;
    std::vector<double> gamma_;
    std::vector<Table> tables_;
    double alpha_ = 0.0;
    double b_max_ = 0.0;
    double nu_bound_ = 0.0;
};

} // namespace rad_hdiv
