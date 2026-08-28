/* rad_hacapk_hdiv_entry.cpp -- element-specific ChargeGram entry strategies. */

#include "rad_hacapk_hdiv.h"

#include <cmath>
#include <memory>
#include <stdexcept>

namespace {

constexpr double kInvFourPi = 0.07957747154594766788;

}  // namespace

struct RadHACApKChargeGram::SampledLaplaceEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override
    { return owner.EvaluateSampledLaplaceEntry(row, col); }
};

struct RadHACApKChargeGram::SampledPlanarEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override
    { return owner.EvaluateSampledPlanarEntry(row, col); }
};

struct RadHACApKChargeGram::PlanarEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override
    { return owner.EvaluatePlanarEntry(row, col); }
};

struct RadHACApKChargeGram::HexEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override
    { return owner.EvaluateHexEntry(row, col); }
};

struct RadHACApKChargeGram::WedgeEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override
    { return owner.EvaluateWedgeEntry(row, col); }
};

struct RadHACApKChargeGram::HighOrderTetEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override
    { return owner.EvaluateHighOrderTetEntry(row, col); }
};

struct RadHACApKChargeGram::AnalyticEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override
    { return owner.EvaluateAnalyticEntry(row, col); }
};

struct RadHACApKChargeGram::MonopoleEntryStrategy final : EntryStrategy {
    double Evaluate(const RadHACApKChargeGram& owner, int row, int col) const override
    { return owner.EvaluateMonopoleEntry(row, col); }
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

double RadHACApKChargeGram::EvaluateSampledLaplaceEntry(int a, int b) const
{
    const double dx = m_cent[3 * a] - m_cent[3 * b];
    const double dy = m_cent[3 * a + 1] - m_cent[3 * b + 1];
    const double dz = m_cent[3 * a + 2] - m_cent[3 * b + 2];
    const double eps2 = m_sampledKernelEpsilon * m_sampledKernelEpsilon;
    return m_meas[a] * m_meas[b] * kInvFourPi /
           std::sqrt(dx * dx + dy * dy + dz * dz + eps2);
}

double RadHACApKChargeGram::EvaluateSampledPlanarEntry(int a, int b) const
{
    const double dx = m_cent[3 * a] - m_cent[3 * b];
    const double dy = m_cent[3 * a + 1] - m_cent[3 * b + 1];
    const double eps2 = m_sampledKernelEpsilon * m_sampledKernelEpsilon;
    const double distance = std::sqrt(dx * dx + dy * dy + eps2);
    return -2.0 * kInvFourPi * m_meas[a] * m_meas[b] *
           std::log(distance / m_sampledReferenceLength);
}

double RadHACApKChargeGram::EvaluatePlanarEntry(int a, int b) const
{
    return EvaluateHostBlockEntry(a, b);
}

double RadHACApKChargeGram::EvaluateHexEntry(int a, int b) const
{
    return EvaluateHostBlockEntry(a, b);
}

double RadHACApKChargeGram::EvaluateWedgeEntry(int a, int b) const
{
    return EvaluateHostBlockEntry(a, b);
}

double RadHACApKChargeGram::EvaluateHostBlockEntry(int a, int b) const
{
    const int kind_a = m_kind[a], host_a = m_host[a];
    const int kind_b = m_kind[b], host_b = m_host[b];
    const int local_a = m_hexLocalOf[a], local_b = m_hexLocalOf[b];
    const int count_b = kind_b == 0 ? (int)m_cellCharges[host_b].size()
                                    : (int)m_faceCharges[host_b].size();
    double value = GetHexSymBlock(kind_a, host_a, kind_b, host_b)
        [(size_t)local_a * count_b + local_b];
    for (size_t image = 0; image < m_image_masks.size(); ++image)
        value += m_image_signs[image] *
            GetHexSymBlock(kind_a, host_a, kind_b, host_b, (int)image + 1)
                [(size_t)local_a * count_b + local_b];
    return value;
}

double RadHACApKChargeGram::EvaluateHighOrderTetEntry(int a, int b) const
{
    const bool far_pair = a != b && m_ho_far_factor < 1e29 && [&] {
        const double dx = m_cent[3*a] - m_cent[3*b];
        const double dy = m_cent[3*a+1] - m_cent[3*b+1];
        const double dz = m_cent[3*a+2] - m_cent[3*b+2];
        return std::sqrt(dx*dx + dy*dy + dz*dz) >
               m_ho_far_factor * (m_size[a] + m_size[b]);
    }();
    double value;
    const bool use_host_block = m_polyCombo || m_curved ||
        (m_hoAnalyticBlock && HOAnalyticBlockEnabled());
    if (use_host_block && !far_pair) {
        const int kind_a = m_kind[a], host_a = m_host[a];
        const int kind_b = m_kind[b], host_b = m_host[b];
        const int local_a = m_hoLocalOf[a], local_b = m_hoLocalOf[b];
        const int count_b = kind_b == 0 ? (int)m_hoCellCharges[host_b].size()
                                        : (int)m_hoFaceCharges[host_b].size();
        if (!m_curved ||
            !CurvedTouchBlockValue(kind_a, host_a, local_a, kind_b, host_b, local_b, value))
            value = GetHOTetSymBlock(kind_a, host_a, kind_b, host_b)
                [(size_t)local_a * count_b + local_b];
    }
    else if (a == b)
        value = QuadDot(a, a);
    else if (far_pair)
        value = HOFarOneSidedEnabled()
            ? QuadDotFar(a, b)
            : 0.5 * (QuadDotFar(a, b) + QuadDotFar(b, a));
    else
        value = 0.5 * (QuadDot(a, b) + QuadDot(b, a));

    for (size_t image = 0; image < m_image_masks.size(); ++image) {
        const int image_id = (int)image + 1;
        const bool rotation_image = image < m_image_rot_angle.size() &&
                                    m_image_rot_angle[image] != 0.0;
        if (rotation_image && !m_curved && !m_polyCombo && !m_hexmode &&
            !m_wedgemode && m_hoAnalyticBlock && HOAnalyticBlockEnabled() &&
            HOTetImageBlockEnabled()) {
            const int kind_a = m_kind[a], host_a = m_host[a];
            const int kind_b = m_kind[b], host_b = m_host[b];
            const int local_a = m_hoLocalOf[a], local_b = m_hoLocalOf[b];
            const int count_b = kind_b == 0 ? (int)m_hoCellCharges[host_b].size()
                                            : (int)m_hoFaceCharges[host_b].size();
            value += m_image_signs[image] *
                GetHOTetSymBlock(kind_a, host_a, kind_b, host_b, image_id)
                    [(size_t)local_a * count_b + local_b];
        }
        else
            value += m_image_signs[image] * 0.5 *
                (QuadDotRefl(a, b, image_id) + QuadDotRefl(b, a, image_id));
    }
    return value;
}

double RadHACApKChargeGram::EvaluateAnalyticEntry(int a, int b) const
{
    double value;
    if (a == b)
        value = QuadDot(a, a);
    else {
        const double dx = m_cent[3*a] - m_cent[3*b];
        const double dy = m_cent[3*a+1] - m_cent[3*b+1];
        const double dz = m_cent[3*a+2] - m_cent[3*b+2];
        const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (distance <= m_near_factor * (m_size[a] + m_size[b]))
            value = 0.5 * (QuadDot(a, b) + QuadDot(b, a));
        else if (m_far_quad > 0)
            value = QuadDotFarLow(a, b);
        else
            value = m_meas[a] * m_meas[b] * kInvFourPi / distance;
    }
    for (size_t image = 0; image < m_image_masks.size(); ++image)
        value += m_image_signs[image] * 0.5 *
            (QuadDotRefl(a, b, (int)image + 1) + QuadDotRefl(b, a, (int)image + 1));
    return value;
}

double RadHACApKChargeGram::EvaluateMonopoleEntry(int a, int b) const
{
    if (a == b) return m_self[a];
    const double dx = m_cent[3*a] - m_cent[3*b];
    const double dy = m_cent[3*a+1] - m_cent[3*b+1];
    const double dz = m_cent[3*a+2] - m_cent[3*b+2];
    return m_meas[a] * m_meas[b] * kInvFourPi /
           std::sqrt(dx*dx + dy*dy + dz*dz);
}

double RadHACApKChargeGram::EvaluateDirectSelfEntry(int a) const
{
    if (m_sampledLaplace)
        return EvaluateSampledLaplaceEntry(a, a);
    if (m_sampledPlanarLog)
        return EvaluateSampledPlanarEntry(a, a);
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
