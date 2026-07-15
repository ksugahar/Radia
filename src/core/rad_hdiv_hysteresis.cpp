#include "rad_hdiv_hysteresis.h"

#include "rad_parallel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace rad_hdiv {
namespace {

bool finite(double value) { return std::isfinite(value); }

double norm3(const double* value)
{
    return std::sqrt(value[0] * value[0] + value[1] * value[1] +
                     value[2] * value[2]);
}

} // namespace

EnergyStopMaterial::EnergyStopMaterial(
    std::vector<double> eta, std::vector<double> table_r,
    std::vector<double> table_g, std::vector<int> table_offsets,
    std::vector<double> gamma, double alpha, double b_max)
    : eta_(std::move(eta)), gamma_(std::move(gamma)), alpha_(alpha),
      b_max_(b_max)
{
    if (eta_.empty())
        throw std::invalid_argument("EnergyStopMaterial: eta must not be empty");
    if (gamma_.size() != eta_.size())
        throw std::invalid_argument("EnergyStopMaterial: gamma must have one value per branch");
    if (table_offsets.size() != eta_.size() + 1 || table_offsets.front() != 0 ||
        table_offsets.back() != static_cast<int>(table_r.size()) ||
        table_r.size() != table_g.size())
        throw std::invalid_argument("EnergyStopMaterial: invalid flattened table offsets");
    for (std::size_t i = 0; i + 1 < table_offsets.size(); ++i) {
        if (table_offsets[i] < 0 || table_offsets[i + 1] < table_offsets[i] ||
            table_offsets[i + 1] > static_cast<int>(table_r.size()))
            throw std::invalid_argument("EnergyStopMaterial: table offsets must be ordered and in range");
    }
    if (!finite(alpha_) || alpha_ <= 0.0)
        throw std::invalid_argument("EnergyStopMaterial: alpha must be finite and positive");
    if (!(finite(b_max_) && b_max_ > 0.0))
        b_max_ = std::numeric_limits<double>::infinity();

    tables_.reserve(eta_.size());
    nu_bound_ = alpha_;
    for (std::size_t branch = 0; branch < eta_.size(); ++branch) {
        const double limit = eta_[branch];
        if (!finite(limit) || limit <= 0.0)
            throw std::invalid_argument("EnergyStopMaterial: every eta must be finite and positive");
        if (!finite(gamma_[branch]) || gamma_[branch] < 0.0)
            throw std::invalid_argument("EnergyStopMaterial: gamma must be finite and non-negative");

        const int begin = table_offsets[branch];
        const int end = table_offsets[branch + 1];
        if (end - begin < 2)
            throw std::invalid_argument("EnergyStopMaterial: every g table needs at least two samples");

        Table table;
        table.r.assign(table_r.begin() + begin, table_r.begin() + end);
        table.g.assign(table_g.begin() + begin, table_g.begin() + end);
        table.u.assign(table.r.size(), 0.0);
        const double scale = std::max(1.0, std::abs(table.g.back()));
        if (!finite(table.r.front()) || std::abs(table.r.front()) > 1.0e-14 ||
            !finite(table.g.front()) || std::abs(table.g.front()) > 1.0e-12 * scale)
            throw std::invalid_argument("EnergyStopMaterial: every table must start at r=0, g=0");

        for (std::size_t i = 1; i < table.r.size(); ++i) {
            if (!finite(table.r[i]) || !finite(table.g[i]) ||
                table.r[i] <= table.r[i - 1])
                throw std::invalid_argument("EnergyStopMaterial: table radii must be finite and strictly increasing");
            if (table.g[i] < -1.0e-14 ||
                table.g[i] + 1.0e-12 * scale < table.g[i - 1])
                throw std::invalid_argument("EnergyStopMaterial: g tables must be non-negative and monotone");
            const double dr = table.r[i] - table.r[i - 1];
            const double slope = (table.g[i] - table.g[i - 1]) / dr;
            table.max_slope = std::max(table.max_slope, slope);
            table.max_tangent = std::max(table.max_tangent, table.g[i] / table.r[i]);
            table.u[i] = table.u[i - 1] +
                         0.5 * (table.g[i - 1] + table.g[i]) * dr;
        }
        if (table.r.back() + 1.0e-14 * std::max(1.0, limit) < limit)
            throw std::invalid_argument("EnergyStopMaterial: each g table must cover its eta");
        nu_bound_ += std::max(table.max_slope, table.max_tangent);
        tables_.push_back(std::move(table));
    }
}

void EnergyStopMaterial::Interpolate(const Table& table, double radius,
                                     double& g, double& u) const
{
    if (radius <= 0.0) {
        g = 0.0;
        u = 0.0;
        return;
    }
    if (radius >= table.r.back()) {
        g = table.g.back();
        u = table.u.back();
        return;
    }
    const auto upper = std::upper_bound(table.r.begin(), table.r.end(), radius);
    const std::size_t i = static_cast<std::size_t>(upper - table.r.begin());
    const double r0 = table.r[i - 1];
    const double g0 = table.g[i - 1];
    const double slope = (table.g[i] - g0) / (table.r[i] - r0);
    const double dr = radius - r0;
    g = g0 + slope * dr;
    u = table.u[i - 1] + g0 * dr + 0.5 * slope * dr * dr;
}

EnergyStopMaterial::BranchTrial EnergyStopMaterial::TrialBranch(
    std::size_t branch, const double* state, const double* delta_B) const
{
    BranchTrial trial{};
    const double c[3] = {
        state[3 * branch] + delta_B[0],
        state[3 * branch + 1] + delta_B[1],
        state[3 * branch + 2] + delta_B[2],
    };
    const double c_norm = norm3(c);
    if (c_norm == 0.0)
        return trial;

    const double limit = eta_[branch];
    double radius = std::min(c_norm, limit);
    if (gamma_[branch] > 0.0) {
        const double inv_gamma = 1.0 / gamma_[branch];
        double g_limit = 0.0, u_limit = 0.0;
        Interpolate(tables_[branch], limit, g_limit, u_limit);
        if (g_limit + (limit - c_norm) * inv_gamma > 0.0) {
            double lower = 0.0;
            double upper = limit;
            for (int iteration = 0; iteration < 64; ++iteration) {
                const double mid = 0.5 * (lower + upper);
                double g_mid = 0.0, u_mid = 0.0;
                Interpolate(tables_[branch], mid, g_mid, u_mid);
                if (g_mid + (mid - c_norm) * inv_gamma < 0.0)
                    lower = mid;
                else
                    upper = mid;
            }
            radius = 0.5 * (lower + upper);
        }
    }

    const double scale = radius / c_norm;
    trial.s[0] = scale * c[0];
    trial.s[1] = scale * c[1];
    trial.s[2] = scale * c[2];
    trial.radius = radius;
    Interpolate(tables_[branch], radius, trial.g, trial.u);
    return trial;
}

void EnergyStopMaterial::State0(double* output) const
{
    std::fill(output, output + StateSize(), 0.0);
}

void EnergyStopMaterial::ForwardBatch(const double* B, const double* states,
                                      int count, double* H) const
{
    if (count < 0)
        throw std::invalid_argument("EnergyStopMaterial::ForwardBatch: negative count");
    const std::size_t state_size = StateSize();
    const std::size_t b_offset = 3 * BranchCount();
    ngcore::RegionTaskManager task_manager(radia::GetMaxThreads());
    ngcore::ParallelFor(ngcore::IntRange(count), [&](int row) {
        const double* b = B + 3 * row;
        const double* state = states + state_size * row;
        const double delta[3] = {b[0] - state[b_offset],
                                 b[1] - state[b_offset + 1],
                                 b[2] - state[b_offset + 2]};
        double* h = H + 3 * row;
        h[0] = alpha_ * b[0];
        h[1] = alpha_ * b[1];
        h[2] = alpha_ * b[2];
        for (std::size_t branch = 0; branch < BranchCount(); ++branch) {
            const BranchTrial trial = TrialBranch(branch, state, delta);
            if (trial.radius > 0.0) {
                const double factor = trial.g / trial.radius;
                h[0] += factor * trial.s[0];
                h[1] += factor * trial.s[1];
                h[2] += factor * trial.s[2];
            }
        }
    });
}

void EnergyStopMaterial::CommitBatch(const double* B, const double* states,
                                     int count, double* new_states) const
{
    if (count < 0)
        throw std::invalid_argument("EnergyStopMaterial::CommitBatch: negative count");
    const std::size_t state_size = StateSize();
    const std::size_t b_offset = 3 * BranchCount();
    ngcore::RegionTaskManager task_manager(radia::GetMaxThreads());
    ngcore::ParallelFor(ngcore::IntRange(count), [&](int row) {
        const double* b = B + 3 * row;
        const double* state = states + state_size * row;
        double* output = new_states + state_size * row;
        const double delta[3] = {b[0] - state[b_offset],
                                 b[1] - state[b_offset + 1],
                                 b[2] - state[b_offset + 2]};
        for (std::size_t branch = 0; branch < BranchCount(); ++branch) {
            const BranchTrial trial = TrialBranch(branch, state, delta);
            output[3 * branch] = trial.s[0];
            output[3 * branch + 1] = trial.s[1];
            output[3 * branch + 2] = trial.s[2];
        }
        output[b_offset] = b[0];
        output[b_offset + 1] = b[1];
        output[b_offset + 2] = b[2];
    });
}

void EnergyStopMaterial::StoredEnergyBatch(const double* B, const double* states,
                                           int count, double* energy) const
{
    if (count < 0)
        throw std::invalid_argument("EnergyStopMaterial::StoredEnergyBatch: negative count");
    const std::size_t state_size = StateSize();
    const std::size_t b_offset = 3 * BranchCount();
    ngcore::RegionTaskManager task_manager(radia::GetMaxThreads());
    ngcore::ParallelFor(ngcore::IntRange(count), [&](int row) {
        const double* b = B + 3 * row;
        const double* state = states + state_size * row;
        const double delta[3] = {b[0] - state[b_offset],
                                 b[1] - state[b_offset + 1],
                                 b[2] - state[b_offset + 2]};
        double value = 0.5 * alpha_ *
                       (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
        for (std::size_t branch = 0; branch < BranchCount(); ++branch)
            value += TrialBranch(branch, state, delta).u;
        energy[row] = value;
    });
}

} // namespace rad_hdiv
