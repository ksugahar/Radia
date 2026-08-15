#ifndef RAD_NGSOLVE_FIELD_COEFFICIENTS_H
#define RAD_NGSOLVE_FIELD_COEFFICIENTS_H

#include "rad_hdiv_field_evaluator.h"
#include "rad_planar_charges.h"

#include <coefficient.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace radia::ngsolve_bridge {

class HDivFieldCoefficient final : public ngfem::CoefficientFunction {
public:
    HDivFieldCoefficient(
        std::shared_ptr<rad_hdiv::HDivFieldEvaluator> evaluator,
        const std::string& algorithm = "direct")
        : CoefficientFunction(3), evaluator_(std::move(evaluator)),
          algorithm_(rad_hdiv::HDivFieldEvaluator::ParseAlgorithm(algorithm)) {
        if (!evaluator_)
            throw std::invalid_argument(
                "HDivFieldCoefficient: evaluator must not be null");
    }

    double Evaluate(const ngfem::BaseMappedIntegrationPoint&) const override {
        return 0.0;
    }

    void Evaluate(const ngfem::BaseMappedIntegrationPoint& mip,
                  ngbla::FlatVector<> result) const override {
        const auto point = mip.GetPoint();
        const int dimension = point.Size();
        double observation[3] = {
            point[0], dimension >= 2 ? point[1] : 0.0,
            dimension >= 3 ? point[2] : 0.0};
        double value[3];
        evaluator_->EvaluateSerial(observation, 1, value, algorithm_);
        result(0) = scale_ * value[0];
        result(1) = scale_ * value[1];
        result(2) = scale_ * value[2];
    }

    void Evaluate(const ngfem::BaseMappedIntegrationRule& rule,
                  ngbla::BareSliceMatrix<> result) const override {
        const std::size_t count = rule.Size();
        thread_local std::vector<double> observations;
        thread_local std::vector<double> values;
        observations.resize(3 * count);
        values.resize(3 * count);
        for (std::size_t i = 0; i < count; ++i) {
            const auto point = rule[i].GetPoint();
            const int dimension = point.Size();
            observations[3 * i] = point[0];
            observations[3 * i + 1] = dimension >= 2 ? point[1] : 0.0;
            observations[3 * i + 2] = dimension >= 3 ? point[2] : 0.0;
        }
        evaluator_->EvaluateSerial(
            observations.data(), count, values.data(), algorithm_);
        for (std::size_t i = 0; i < count; ++i) {
            result(i, 0) = scale_ * values[3 * i];
            result(i, 1) = scale_ * values[3 * i + 1];
            result(i, 2) = scale_ * values[3 * i + 2];
        }
    }

    const std::shared_ptr<rad_hdiv::HDivFieldEvaluator>& Evaluator() const {
        return evaluator_;
    }
    const char* AlgorithmName() const {
        return rad_hdiv::HDivFieldEvaluator::AlgorithmName(algorithm_);
    }

private:
    static constexpr double scale_ =
        0.079577471545947667884441881686257181;
    std::shared_ptr<rad_hdiv::HDivFieldEvaluator> evaluator_;
    rad_hdiv::HDivFieldEvaluator::Algorithm algorithm_;
};

class HDivVectorPotentialCoefficient final
    : public ngfem::CoefficientFunction {
public:
    HDivVectorPotentialCoefficient(
        std::vector<double> source_points,
        std::vector<double> integrated_magnetization)
        : CoefficientFunction(3), source_points_(std::move(source_points)),
          integrated_magnetization_(std::move(integrated_magnetization)) {
        if (source_points_.empty() || source_points_.size() % 3 != 0)
            throw std::invalid_argument(
                "HDivVectorPotentialCoefficient: source_points must have "
                "shape (n,3) with n > 0");
        if (integrated_magnetization_.size() != source_points_.size())
            throw std::invalid_argument(
                "HDivVectorPotentialCoefficient: integrated_magnetization "
                "must have the same shape as source_points");
        for (double value : source_points_)
            if (!std::isfinite(value))
                throw std::invalid_argument(
                    "HDivVectorPotentialCoefficient: source_points must be finite");
        for (double value : integrated_magnetization_)
            if (!std::isfinite(value))
                throw std::invalid_argument(
                    "HDivVectorPotentialCoefficient: integrated_magnetization "
                    "must be finite");
    }

    double Evaluate(const ngfem::BaseMappedIntegrationPoint&) const override {
        return 0.0;
    }

    void Evaluate(const ngfem::BaseMappedIntegrationPoint& mip,
                  ngbla::FlatVector<> result) const override {
        const auto point = mip.GetPoint();
        const int dimension = point.Size();
        const double observation[3] = {
            point[0], dimension >= 2 ? point[1] : 0.0,
            dimension >= 3 ? point[2] : 0.0};
        double value[3];
        EvaluatePoint(observation, value);
        for (int component = 0; component < 3; ++component)
            result(component) = value[component];
    }

    void Evaluate(const ngfem::BaseMappedIntegrationRule& rule,
                  ngbla::BareSliceMatrix<> result) const override {
        for (std::size_t index = 0; index < rule.Size(); ++index) {
            const auto point = rule[index].GetPoint();
            const int dimension = point.Size();
            const double observation[3] = {
                point[0], dimension >= 2 ? point[1] : 0.0,
                dimension >= 3 ? point[2] : 0.0};
            double value[3];
            EvaluatePoint(observation, value);
            for (int component = 0; component < 3; ++component)
                result(index, component) = value[component];
        }
    }

    std::size_t SourceCount() const { return source_points_.size() / 3; }

private:
    void EvaluatePoint(const double observation[3], double output[3]) const {
        // A(r) = mu0/(4*pi) int M(r') x (r-r') / |r-r'|^3 dV'.
        // The source vectors already include the mapped quadrature weights,
        // so this immutable coefficient can be evaluated safely by NGSolve's
        // worker threads without reconstructing FE basis/orientation data.
        double sum[3] = {0.0, 0.0, 0.0};
        double correction[3] = {0.0, 0.0, 0.0};
        const std::size_t count = SourceCount();
        for (std::size_t source = 0; source < count; ++source) {
            const double displacement[3] = {
                observation[0] - source_points_[3 * source],
                observation[1] - source_points_[3 * source + 1],
                observation[2] - source_points_[3 * source + 2]};
            const double radius_squared =
                displacement[0] * displacement[0]
                + displacement[1] * displacement[1]
                + displacement[2] * displacement[2];
            if (!(radius_squared > singular_radius_squared_))
                throw std::runtime_error(
                    "HDivVectorPotentialCoefficient: an observation point "
                    "coincides with a magnetization quadrature source; this "
                    "exterior-field coefficient must not be evaluated inside iron");
            const double inverse_radius_cubed =
                1.0 / (radius_squared * std::sqrt(radius_squared));
            const double mx = integrated_magnetization_[3 * source];
            const double my = integrated_magnetization_[3 * source + 1];
            const double mz = integrated_magnetization_[3 * source + 2];
            const double contribution[3] = {
                scale_ * (my * displacement[2] - mz * displacement[1])
                    * inverse_radius_cubed,
                scale_ * (mz * displacement[0] - mx * displacement[2])
                    * inverse_radius_cubed,
                scale_ * (mx * displacement[1] - my * displacement[0])
                    * inverse_radius_cubed};
            for (int component = 0; component < 3; ++component) {
                const double next = sum[component] + contribution[component];
                correction[component] +=
                    std::fabs(sum[component]) >= std::fabs(contribution[component])
                    ? (sum[component] - next) + contribution[component]
                    : (contribution[component] - next) + sum[component];
                sum[component] = next;
            }
        }
        for (int component = 0; component < 3; ++component)
            output[component] = sum[component] + correction[component];
    }

    static constexpr double scale_ = 1.0e-7; // mu0/(4*pi), SI
    static constexpr double singular_radius_squared_ =
        64.0 * std::numeric_limits<double>::epsilon()
        * 64.0 * std::numeric_limits<double>::epsilon();
    std::vector<double> source_points_;
    std::vector<double> integrated_magnetization_;
};

class PlanarHDivFieldCoefficient final : public ngfem::CoefficientFunction {
public:
    PlanarHDivFieldCoefficient(
        std::shared_ptr<rad_planar_charges::PlanarFieldEvaluator> evaluator,
        double source_angle = 0.0, double target_angle = 0.0,
        double center_x = 0.0, double center_y = 0.0)
        : CoefficientFunction(2), evaluator_(std::move(evaluator)),
          center_x_(center_x), center_y_(center_y),
          cos_delta_(std::cos(source_angle - target_angle)),
          sin_delta_(std::sin(source_angle - target_angle)),
          source_angle_(source_angle), target_angle_(target_angle) {
        if (!evaluator_)
            throw std::invalid_argument(
                "PlanarHDivFieldCoefficient: evaluator must not be null");
        if (!std::isfinite(source_angle_) || !std::isfinite(target_angle_) ||
            !std::isfinite(center_x_) || !std::isfinite(center_y_))
            throw std::invalid_argument(
                "PlanarHDivFieldCoefficient: angles and center must be finite");
    }

    double Evaluate(const ngfem::BaseMappedIntegrationPoint&) const override {
        return 0.0;
    }

    void Evaluate(const ngfem::BaseMappedIntegrationPoint& mip,
                  ngbla::FlatVector<> result) const override {
        const auto point = mip.GetPoint();
        double source_point[2];
        TargetToSource(point[0], point[1], source_point[0], source_point[1]);
        double source_field[2];
        evaluator_->EvaluateFieldSerial(source_point, 1, source_field);
        SourceFieldToTarget(
            source_field[0], source_field[1], result(0), result(1));
    }

    void Evaluate(const ngfem::BaseMappedIntegrationRule& rule,
                  ngbla::BareSliceMatrix<> result) const override {
        const std::size_t count = rule.Size();
        thread_local std::vector<double> source_points;
        thread_local std::vector<double> source_fields;
        source_points.resize(2 * count);
        source_fields.resize(2 * count);
        for (std::size_t i = 0; i < count; ++i) {
            const auto point = rule[i].GetPoint();
            TargetToSource(point[0], point[1], source_points[2 * i],
                           source_points[2 * i + 1]);
        }
        evaluator_->EvaluateFieldSerial(
            source_points.data(), count, source_fields.data());
        for (std::size_t i = 0; i < count; ++i)
            SourceFieldToTarget(
                source_fields[2 * i], source_fields[2 * i + 1],
                result(i, 0), result(i, 1));
    }

    double SourceAngle() const { return source_angle_; }
    double TargetAngle() const { return target_angle_; }

private:
    void TargetToSource(double target_x, double target_y,
                        double& source_x, double& source_y) const {
        const double dx = target_x - center_x_;
        const double dy = target_y - center_y_;
        source_x = center_x_ + cos_delta_ * dx + sin_delta_ * dy;
        source_y = center_y_ - sin_delta_ * dx + cos_delta_ * dy;
    }

    void SourceFieldToTarget(double source_x, double source_y,
                             double& target_x, double& target_y) const {
        target_x = cos_delta_ * source_x - sin_delta_ * source_y;
        target_y = sin_delta_ * source_x + cos_delta_ * source_y;
    }

    std::shared_ptr<rad_planar_charges::PlanarFieldEvaluator> evaluator_;
    double center_x_ = 0.0;
    double center_y_ = 0.0;
    double cos_delta_ = 1.0;
    double sin_delta_ = 0.0;
    double source_angle_ = 0.0;
    double target_angle_ = 0.0;
};

} // namespace radia::ngsolve_bridge

#endif
