#ifndef RAD_NGSOLVE_FIELD_COEFFICIENTS_H
#define RAD_NGSOLVE_FIELD_COEFFICIENTS_H

#include "rad_hdiv_field_evaluator.h"
#include "rad_planar_charges.h"

#include <coefficient.hpp>

#include <cmath>
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
