/* rad_hacapk_hdiv_entry.cpp -- element-specific ChargeGram entry strategies. */

#include "rad_hacapk_hdiv.h"

#include <cmath>
#include <memory>
#include <stdexcept>

namespace {

constexpr double kInvFourPi = 0.07957747154594766788;

}  // namespace

struct RadHACApKChargeGram::SampledLaplaceEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override;
};

struct RadHACApKChargeGram::SampledPlanarEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override;
};

struct RadHACApKChargeGram::PlanarEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override
    { return EvaluateHostBlock(owner, row, col); }
};

struct RadHACApKChargeGram::HexEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override
    { return EvaluateHostBlock(owner, row, col); }
};

struct RadHACApKChargeGram::WedgeEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override
    { return EvaluateHostBlock(owner, row, col); }
};

struct RadHACApKChargeGram::HighOrderTetEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override;
};

struct RadHACApKChargeGram::AnalyticEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override;
};

struct RadHACApKChargeGram::MonopoleEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override;
};

const RadHACApKChargeGram::EntryStrategy& RadHACApKChargeGram::GetEntryStrategy() const
{
    std::call_once(m_entryStrategyOnce, [this] {
        if (m_sampledLaplace)
            m_entryStrategy = std::make_unique<SampledLaplaceEntryStrategy>();
        else if (m_sampledPlanarLog)
            m_entryStrategy = std::make_unique<SampledPlanarEntryStrategy>();
        else if (m_d2)
            m_entryStrategy = std::make_unique<PlanarEntryStrategy>();
        else if (m_hexmode)
            m_entryStrategy = std::make_unique<HexEntryStrategy>();
        else if (m_wedgemode)
            m_entryStrategy = std::make_unique<WedgeEntryStrategy>();
        else if (m_highorder)
            m_entryStrategy = std::make_unique<HighOrderTetEntryStrategy>();
        else if (m_analytic)
            m_entryStrategy = std::make_unique<AnalyticEntryStrategy>();
        else
            m_entryStrategy = std::make_unique<MonopoleEntryStrategy>();
    });
    return *m_entryStrategy;
}

double RadHACApKChargeGram::SampledLaplaceEntryStrategy::Evaluate(
    const RadHACApKChargeGram& owner, int a, int b) const
{
    const double dx = owner.m_cent[3 * a] - owner.m_cent[3 * b];
    const double dy = owner.m_cent[3 * a + 1] - owner.m_cent[3 * b + 1];
    const double dz = owner.m_cent[3 * a + 2] - owner.m_cent[3 * b + 2];
    const double eps2 = owner.m_sampledKernelEpsilon * owner.m_sampledKernelEpsilon;
    return owner.m_meas[a] * owner.m_meas[b] * kInvFourPi /
           std::sqrt(dx * dx + dy * dy + dz * dz + eps2);
}

double RadHACApKChargeGram::SampledPlanarEntryStrategy::Evaluate(
    const RadHACApKChargeGram& owner, int a, int b) const
{
    const double dx = owner.m_cent[3 * a] - owner.m_cent[3 * b];
    const double dy = owner.m_cent[3 * a + 1] - owner.m_cent[3 * b + 1];
    const double eps2 = owner.m_sampledKernelEpsilon * owner.m_sampledKernelEpsilon;
    const double distance = std::sqrt(dx * dx + dy * dy + eps2);
    return -2.0 * kInvFourPi * owner.m_meas[a] * owner.m_meas[b] *
           std::log(distance / owner.m_sampledReferenceLength);
}

double RadHACApKChargeGram::EntryStrategy::EvaluateHostBlock(
    const RadHACApKChargeGram& owner, int a, int b)
{
    const int kind_a = owner.m_kind[a], host_a = owner.m_host[a];
    const int kind_b = owner.m_kind[b], host_b = owner.m_host[b];
    const int local_a = owner.m_hexLocalOf[a], local_b = owner.m_hexLocalOf[b];
    const int count_b = kind_b == 0 ? (int)owner.m_cellCharges[host_b].size()
                                    : (int)owner.m_faceCharges[host_b].size();
    double value = owner.GetHexSymBlock(kind_a, host_a, kind_b, host_b)
        [(size_t)local_a * count_b + local_b];
    for (size_t image = 0; image < owner.m_image_masks.size(); ++image)
        value += owner.m_image_signs[image] *
            owner.GetHexSymBlock(kind_a, host_a, kind_b, host_b, (int)image + 1)
                [(size_t)local_a * count_b + local_b];
    return value;
}

double RadHACApKChargeGram::HighOrderTetEntryStrategy::Evaluate(
    const RadHACApKChargeGram& owner, int a, int b) const
{
    const bool far_pair = a != b && owner.m_ho_far_factor < 1e29 && [&] {
        const double dx = owner.m_cent[3*a] - owner.m_cent[3*b];
        const double dy = owner.m_cent[3*a+1] - owner.m_cent[3*b+1];
        const double dz = owner.m_cent[3*a+2] - owner.m_cent[3*b+2];
        return std::sqrt(dx*dx + dy*dy + dz*dz) >
               owner.m_ho_far_factor * (owner.m_size[a] + owner.m_size[b]);
    }();
    double value;
    const bool use_host_block = owner.m_polyCombo || owner.m_curved ||
        (owner.m_hoAnalyticBlock && HOAnalyticBlockEnabled());
    if (use_host_block && !far_pair) {
        const int kind_a = owner.m_kind[a], host_a = owner.m_host[a];
        const int kind_b = owner.m_kind[b], host_b = owner.m_host[b];
        const int local_a = owner.m_hoLocalOf[a], local_b = owner.m_hoLocalOf[b];
        const int count_b = kind_b == 0 ? (int)owner.m_hoCellCharges[host_b].size()
                                        : (int)owner.m_hoFaceCharges[host_b].size();
        if (!owner.m_curved ||
            !owner.CurvedTouchBlockValue(kind_a, host_a, local_a, kind_b, host_b, local_b, value))
            value = owner.GetHOTetSymBlock(kind_a, host_a, kind_b, host_b)
                [(size_t)local_a * count_b + local_b];
    }
    else if (a == b)
        value = owner.QuadDot(a, a);
    else if (far_pair)
        value = HOFarOneSidedEnabled()
            ? owner.QuadDotFar(a, b)
            : 0.5 * (owner.QuadDotFar(a, b) + owner.QuadDotFar(b, a));
    else
        value = 0.5 * (owner.QuadDot(a, b) + owner.QuadDot(b, a));

    for (size_t image = 0; image < owner.m_image_masks.size(); ++image) {
        const int image_id = (int)image + 1;
        const bool rotation_image = image < owner.m_image_rot_angle.size() &&
                                    owner.m_image_rot_angle[image] != 0.0;
        if (rotation_image && !owner.m_curved && !owner.m_polyCombo && !owner.m_hexmode &&
            !owner.m_wedgemode && owner.m_hoAnalyticBlock && HOAnalyticBlockEnabled() &&
            HOTetImageBlockEnabled()) {
            const int kind_a = owner.m_kind[a], host_a = owner.m_host[a];
            const int kind_b = owner.m_kind[b], host_b = owner.m_host[b];
            const int local_a = owner.m_hoLocalOf[a], local_b = owner.m_hoLocalOf[b];
            const int count_b = kind_b == 0 ? (int)owner.m_hoCellCharges[host_b].size()
                                            : (int)owner.m_hoFaceCharges[host_b].size();
            value += owner.m_image_signs[image] *
                owner.GetHOTetSymBlock(kind_a, host_a, kind_b, host_b, image_id)
                    [(size_t)local_a * count_b + local_b];
        }
        else
            value += owner.m_image_signs[image] * 0.5 *
                (owner.QuadDotRefl(a, b, image_id) + owner.QuadDotRefl(b, a, image_id));
    }
    return value;
}

double RadHACApKChargeGram::AnalyticEntryStrategy::Evaluate(
    const RadHACApKChargeGram& owner, int a, int b) const
{
    double value;
    if (a == b)
        value = owner.QuadDot(a, a);
    else {
        const double dx = owner.m_cent[3*a] - owner.m_cent[3*b];
        const double dy = owner.m_cent[3*a+1] - owner.m_cent[3*b+1];
        const double dz = owner.m_cent[3*a+2] - owner.m_cent[3*b+2];
        const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (distance <= owner.m_near_factor * (owner.m_size[a] + owner.m_size[b]))
            value = 0.5 * (owner.QuadDot(a, b) + owner.QuadDot(b, a));
        else if (owner.m_far_quad > 0)
            value = owner.QuadDotFarLow(a, b);
        else
            value = owner.m_meas[a] * owner.m_meas[b] * kInvFourPi / distance;
    }
    for (size_t image = 0; image < owner.m_image_masks.size(); ++image)
        value += owner.m_image_signs[image] * 0.5 *
            (owner.QuadDotRefl(a, b, (int)image + 1) +
             owner.QuadDotRefl(b, a, (int)image + 1));
    return value;
}

double RadHACApKChargeGram::MonopoleEntryStrategy::Evaluate(
    const RadHACApKChargeGram& owner, int a, int b) const
{
    if (a == b) return owner.m_self[a];
    const double dx = owner.m_cent[3*a] - owner.m_cent[3*b];
    const double dy = owner.m_cent[3*a+1] - owner.m_cent[3*b+1];
    const double dz = owner.m_cent[3*a+2] - owner.m_cent[3*b+2];
    return owner.m_meas[a] * owner.m_meas[b] * kInvFourPi /
           std::sqrt(dx*dx + dy*dy + dz*dz);
}

double RadHACApKChargeGram::EvaluateDirectSelfEntry(int a) const
{
    if (m_sampledLaplace || m_sampledPlanarLog)
        return GetEntryStrategy().Evaluate(*this, a, a);
    if (m_d2 || m_hexmode || m_wedgemode) {
        const int kind = m_kind[a], host = m_host[a];
        const int local = m_hexLocalOf[a];
        const int count = kind == 0 ? (int)m_cellCharges[host].size()
                                    : (int)m_faceCharges[host].size();
        return GetHexSymBlock(kind, host, kind, host)
            [(size_t)local * count + local];
    }
    if (m_highorder) {
        const bool use_host_block = m_polyCombo || m_curved ||
            (m_hoAnalyticBlock && HOAnalyticBlockEnabled());
        if (!use_host_block)
            return QuadDot(a, a);

        const int kind = m_kind[a], host = m_host[a];
        const int local = m_hoLocalOf[a];
        const int count = kind == 0 ? (int)m_hoCellCharges[host].size()
                                    : (int)m_hoFaceCharges[host].size();
        double value;
        if (!m_curved ||
            !CurvedTouchBlockValue(kind, host, local, kind, host, local, value))
            value = GetHOTetSymBlock(kind, host, kind, host)
                [(size_t)local * count + local];
        return value;
    }
    if (m_analytic)
        return QuadDot(a, a);
    return m_self[a];
}
