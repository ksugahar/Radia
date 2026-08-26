#ifndef RAD_NGSOLVE_RADIA_FIELD_H
#define RAD_NGSOLVE_RADIA_FIELD_H

#include <coefficient.hpp>

#include <array>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#if defined(_WIN32) && (defined(ALPHA__DLL__) || defined(MATLAB_MEX_FILE))
#define RADIA_FIELD_C_API __declspec(dllexport)
#define RADIA_FIELD_C_CALL __cdecl
#else
#define RADIA_FIELD_C_API
#define RADIA_FIELD_C_CALL
#endif
extern "C" {
RADIA_FIELD_C_API int RADIA_FIELD_C_CALL RadFldPhiSerial(
    double* phi, int n_points, double* points, int object);
RADIA_FIELD_C_API int RADIA_FIELD_C_CALL RadFldASerial(
    double* vector_potential, int n_points, double* points, int object);
RADIA_FIELD_C_API int RADIA_FIELD_C_CALL RadFldBatchSerial(
    double* flux_density, double* field_strength, int n_points,
    double* points, int object);
RADIA_FIELD_C_API int RADIA_FIELD_C_CALL RadFldCmpPrc(
    int* result, char* options);
}
#undef RADIA_FIELD_C_API
#undef RADIA_FIELD_C_CALL

namespace radia::ngsolve_bridge {

struct RadiaFieldCacheStats {
    bool enabled = false;
    std::size_t size = 0;
    std::size_t hits = 0;
    std::size_t misses = 0;
    double hit_rate = 0.0;
};

class RadiaVoxelCoefficient final : public ngfem::CoefficientFunction {
public:
    RadiaVoxelCoefficient(std::array<double, 3> lower,
                          std::array<double, 3> upper, int resolution,
                          int components, std::vector<double> values)
        : CoefficientFunction(components), lower_(lower), upper_(upper),
          resolution_(resolution), components_(components),
          values_(std::move(values)) {
        const std::size_t side = static_cast<std::size_t>(resolution_);
        if (resolution_ < 2 || (components_ != 1 && components_ != 3) ||
            values_.size() != side * side * side * components_)
            throw std::invalid_argument("invalid Radia voxel field shape");
    }

    double Evaluate(
        const ngfem::BaseMappedIntegrationPoint& point) const override {
        if (components_ != 1) return 0.0;
        return Interpolate(point, 0);
    }

    void Evaluate(const ngfem::BaseMappedIntegrationPoint& point,
                  ngbla::FlatVector<> result) const override {
        for (int component = 0; component < components_; ++component)
            result(component) = Interpolate(point, component);
    }

    void Evaluate(const ngfem::BaseMappedIntegrationRule& rule,
                  ngbla::BareSliceMatrix<> result) const override {
        for (std::size_t i = 0; i < rule.Size(); ++i)
            for (int component = 0; component < components_; ++component)
                result(i, component) = Interpolate(rule[i], component);
    }

private:
    double Interpolate(const ngfem::BaseMappedIntegrationPoint& point,
                       int component) const {
        const auto mapped = point.GetPoint();
        std::array<int, 3> lower_index{};
        std::array<double, 3> fraction{};
        for (int axis = 0; axis < 3; ++axis) {
            const double coordinate = axis < mapped.Size() ? mapped[axis] : 0.0;
            const double normalized = std::clamp(
                (coordinate - lower_[axis]) /
                    (upper_[axis] - lower_[axis]) * (resolution_ - 1),
                0.0, static_cast<double>(resolution_ - 1));
            lower_index[axis] = std::min(
                static_cast<int>(std::floor(normalized)), resolution_ - 2);
            fraction[axis] = normalized - lower_index[axis];
        }
        double result = 0.0;
        for (int dz = 0; dz < 2; ++dz)
            for (int dy = 0; dy < 2; ++dy)
                for (int dx = 0; dx < 2; ++dx) {
                    const double weight =
                        (dx ? fraction[0] : 1.0 - fraction[0]) *
                        (dy ? fraction[1] : 1.0 - fraction[1]) *
                        (dz ? fraction[2] : 1.0 - fraction[2]);
                    const std::size_t x = lower_index[0] + dx;
                    const std::size_t y = lower_index[1] + dy;
                    const std::size_t z = lower_index[2] + dz;
                    const std::size_t index =
                        (((z * resolution_ + y) * resolution_ + x) *
                         components_) + component;
                    result += weight * values_[index];
                }
        return result;
    }

    std::array<double, 3> lower_{};
    std::array<double, 3> upper_{};
    int resolution_ = 0;
    int components_ = 0;
    std::vector<double> values_;
};

class RadiaFieldCoefficient final
    : public ngfem::CoefficientFunctionNoDerivative {
public:
    RadiaFieldCoefficient(
        int object, const std::string& field_type = "b",
        std::optional<std::vector<double>> origin = std::nullopt,
        std::optional<std::vector<double>> u_axis = std::nullopt,
        std::optional<std::vector<double>> v_axis = std::nullopt,
        std::optional<std::vector<double>> w_axis = std::nullopt,
        std::optional<double> precision = std::nullopt,
        const std::string& units = "m")
        : CoefficientFunctionNoDerivative(field_type == "phi" ? 1 : 3),
          object_(object),
          field_type_(field_type), precision_(precision) {
        if (object_ <= 0)
            throw std::invalid_argument(
                "RadiaField requires a positive Radia object handle");
        if (field_type_ != "b" && field_type_ != "h" &&
            field_type_ != "a" && field_type_ != "m" &&
            field_type_ != "phi")
            throw std::invalid_argument(
                "field_type must be 'b', 'h', 'a', 'm', or 'phi'");
        if (units != "m")
            throw std::invalid_argument(
                "RadiaField requires units='m'. Radia always uses meters.");

        ApplyVector(origin, origin_, false);
        ApplyVector(u_axis, u_axis_, true);
        ApplyVector(v_axis, v_axis_, true);
        ApplyVector(w_axis, w_axis_, true);

        if (precision_) {
            if (!std::isfinite(*precision_) || *precision_ <= 0.0)
                throw std::invalid_argument("precision must be positive and finite");
            std::string options =
                "PrcB->" + std::to_string(*precision_) +
                ",PrcA->" + std::to_string(*precision_) +
                ",PrcH->" + std::to_string(*precision_) +
                ",PrcM->" + std::to_string(*precision_);
            std::vector<char> mutable_options(options.begin(), options.end());
            mutable_options.push_back('\0');
            int result = 0;
            CheckRadiaError(RadFldCmpPrc(&result, mutable_options.data()),
                            "precision configuration");
        }
    }

    int Object() const { return object_; }
    const std::string& FieldType() const { return field_type_; }
    bool UsesTransform() const { return use_transform_; }
    const std::optional<double>& Precision() const { return precision_; }

    void PrepareCache(const std::vector<double>& global_points) {
        if (global_points.size() % 3 != 0)
            throw std::invalid_argument(
                "RadiaField cache points must contain coordinate triples");
        const std::size_t count = global_points.size() / 3;
        point_cache_.clear();
        cache_hits_ = 0;
        cache_misses_ = 0;
        if (count == 0) {
            use_cache_ = false;
            return;
        }

        std::vector<double> local_points(global_points.size());
        for (std::size_t i = 0; i < count; ++i)
            TransformToLocal(&global_points[3 * i], &local_points[3 * i]);
        std::vector<double> local_values;
        ComputeLocalField(local_points, count, local_values);
        for (std::size_t i = 0; i < count; ++i) {
            std::array<double, 3> cached{};
            if (field_type_ == "phi") {
                cached[0] = local_values[i];
            } else {
                TransformToGlobal(&local_values[3 * i], cached.data());
            }
            point_cache_[HashPoint(&global_points[3 * i])] = cached;
        }
        use_cache_ = true;
    }

    void ClearCache() {
        point_cache_.clear();
        use_cache_ = false;
        cache_hits_ = 0;
        cache_misses_ = 0;
    }

    RadiaFieldCacheStats CacheStats() const {
        RadiaFieldCacheStats result;
        result.enabled = use_cache_;
        result.size = point_cache_.size();
        result.hits = cache_hits_.load();
        result.misses = cache_misses_.load();
        const double total = static_cast<double>(result.hits + result.misses);
        result.hit_rate = total > 0.0 ? result.hits / total : 0.0;
        return result;
    }

    std::shared_ptr<ngfem::CoefficientFunction> AsVoxelCoefficient(
        const std::array<double, 3>& lower,
        const std::array<double, 3>& upper, int resolution) const {
        if (resolution < 2)
            throw std::invalid_argument("voxel resolution must be at least 2");
        for (int component = 0; component < 3; ++component)
            if (!std::isfinite(lower[component]) ||
                !std::isfinite(upper[component]) ||
                upper[component] <= lower[component])
                throw std::invalid_argument(
                    "voxel bounds must be finite and strictly increasing");

        const std::size_t side = static_cast<std::size_t>(resolution);
        const std::size_t total = side * side * side;
        std::vector<double> global_points(3 * total);
        std::size_t index = 0;
        for (int z = 0; z < resolution; ++z)
            for (int y = 0; y < resolution; ++y)
                for (int x = 0; x < resolution; ++x) {
                    global_points[3 * index] = lower[0] +
                        (upper[0] - lower[0]) * x / (resolution - 1);
                    global_points[3 * index + 1] = lower[1] +
                        (upper[1] - lower[1]) * y / (resolution - 1);
                    global_points[3 * index + 2] = lower[2] +
                        (upper[2] - lower[2]) * z / (resolution - 1);
                    ++index;
                }

        std::vector<double> global_values;
        EvaluateGlobalPoints(global_points, global_values);
        return std::make_shared<RadiaVoxelCoefficient>(
            lower, upper, resolution, field_type_ == "phi" ? 1 : 3,
            std::move(global_values));
    }

    double Evaluate(
        const ngfem::BaseMappedIntegrationPoint& point) const override {
        if (field_type_ != "phi") return 0.0;
        const auto mapped = point.GetPoint();
        const int dimension = mapped.Size();
        double global[3] = {mapped[0], dimension >= 2 ? mapped[1] : 0.0,
                            dimension >= 3 ? mapped[2] : 0.0};
        std::vector<double> local(3);
        TransformToLocal(global, local.data());
        std::vector<double> values;
        ComputeLocalField(local, 1, values);
        return values[0];
    }

    void Evaluate(const ngfem::BaseMappedIntegrationPoint& point,
                  ngbla::FlatVector<> result) const override {
        const auto mapped = point.GetPoint();
        const int dimension = mapped.Size();
        double global[3] = {mapped[0], dimension >= 2 ? mapped[1] : 0.0,
                            dimension >= 3 ? mapped[2] : 0.0};
        if (field_type_ == "phi") {
            result(0) = Evaluate(point);
            return;
        }
        if (ReadCache(global, result)) return;
        std::vector<double> local(3);
        TransformToLocal(global, local.data());
        std::vector<double> values;
        ComputeLocalField(local, 1, values);
        double transformed[3];
        TransformToGlobal(values.data(), transformed);
        for (int component = 0; component < 3; ++component)
            result(component) = transformed[component];
    }

    void Evaluate(const ngfem::BaseMappedIntegrationRule& rule,
                  ngbla::BareSliceMatrix<> result) const override {
        const std::size_t count = rule.Size();
        if (use_cache_ && ReadCache(rule, result)) return;

        std::vector<double> local_points(3 * count);
        for (std::size_t i = 0; i < count; ++i) {
            const auto mapped = rule[i].GetPoint();
            const int dimension = mapped.Size();
            double global[3] = {
                mapped[0], dimension >= 2 ? mapped[1] : 0.0,
                dimension >= 3 ? mapped[2] : 0.0};
            TransformToLocal(global, &local_points[3 * i]);
        }
        std::vector<double> values;
        ComputeLocalField(local_points, count, values);
        if (field_type_ == "phi") {
            for (std::size_t i = 0; i < count; ++i) result(i, 0) = values[i];
            return;
        }
        for (std::size_t i = 0; i < count; ++i) {
            double transformed[3];
            TransformToGlobal(&values[3 * i], transformed);
            for (int component = 0; component < 3; ++component)
                result(i, component) = transformed[component];
        }
    }

private:
    static void CheckRadiaError(int error, const char* operation) {
        if (error != 0)
            throw std::runtime_error(
                std::string("RadiaField: Radia error ") +
                std::to_string(error) + " during " + operation);
    }

    static void Normalize(std::array<double, 3>& vector) {
        const double norm = std::sqrt(
            vector[0] * vector[0] + vector[1] * vector[1] +
            vector[2] * vector[2]);
        if (norm < 1.0e-12)
            throw std::invalid_argument("cannot normalize a zero axis");
        for (double& value : vector) value /= norm;
    }

    void ApplyVector(const std::optional<std::vector<double>>& value,
                     std::array<double, 3>& destination, bool normalize) {
        if (!value) return;
        if (value->size() != 3)
            throw std::invalid_argument(
                "RadiaField transformation vectors must have three components");
        std::copy(value->begin(), value->end(), destination.begin());
        if (normalize) Normalize(destination);
        use_transform_ = true;
    }

    static double Dot(const std::array<double, 3>& left,
                      const double right[3]) {
        return left[0] * right[0] + left[1] * right[1] +
               left[2] * right[2];
    }

    void TransformToLocal(const double global[3], double local[3]) const {
        if (!use_transform_) {
            std::copy(global, global + 3, local);
            return;
        }
        double shifted[3] = {global[0] - origin_[0], global[1] - origin_[1],
                             global[2] - origin_[2]};
        local[0] = Dot(u_axis_, shifted);
        local[1] = Dot(v_axis_, shifted);
        local[2] = Dot(w_axis_, shifted);
    }

    void TransformToGlobal(const double local[3], double global[3]) const {
        if (!use_transform_) {
            std::copy(local, local + 3, global);
            return;
        }
        for (int row = 0; row < 3; ++row)
            global[row] = u_axis_[row] * local[0] +
                          v_axis_[row] * local[1] +
                          w_axis_[row] * local[2];
    }

    std::uint64_t HashPoint(const double point[3]) const {
        const std::int64_t ix =
            static_cast<std::int64_t>(point[0] / cache_tolerance_);
        const std::int64_t iy =
            static_cast<std::int64_t>(point[1] / cache_tolerance_);
        const std::int64_t iz =
            static_cast<std::int64_t>(point[2] / cache_tolerance_);
        std::uint64_t hash = 14695981039346656037ULL;
        hash ^= static_cast<std::uint64_t>(ix); hash *= 1099511628211ULL;
        hash ^= static_cast<std::uint64_t>(iy); hash *= 1099511628211ULL;
        hash ^= static_cast<std::uint64_t>(iz); hash *= 1099511628211ULL;
        return hash;
    }

    void ComputeLocalField(std::vector<double>& points,
                           std::size_t count,
                           std::vector<double>& values) const {
        const int native_count = static_cast<int>(count);
        if (field_type_ == "phi") {
            values.assign(count, 0.0);
            CheckRadiaError(RadFldPhiSerial(
                values.data(), native_count, points.data(), object_),
                "scalar-potential evaluation");
            return;
        }
        if (field_type_ == "a") {
            values.assign(3 * count, 0.0);
            CheckRadiaError(RadFldASerial(
                values.data(), native_count, points.data(), object_),
                "vector-potential evaluation");
            return;
        }
        std::vector<double> flux_density(3 * count, 0.0);
        std::vector<double> field_strength(3 * count, 0.0);
        CheckRadiaError(RadFldBatchSerial(
            flux_density.data(), field_strength.data(), native_count,
            points.data(), object_), "field evaluation");
        if (field_type_ == "b") values = std::move(flux_density);
        else if (field_type_ == "h") values = std::move(field_strength);
        else {
            constexpr double inverse_mu0 =
                1.0 / (4.0 * 3.14159265358979323846e-7);
            values.resize(3 * count);
            for (std::size_t i = 0; i < values.size(); ++i)
                values[i] = flux_density[i] * inverse_mu0 - field_strength[i];
        }
    }

    void EvaluateGlobalPoints(const std::vector<double>& global_points,
                              std::vector<double>& global_values) const {
        const std::size_t count = global_points.size() / 3;
        std::vector<double> local_points(global_points.size());
        for (std::size_t i = 0; i < count; ++i)
            TransformToLocal(&global_points[3 * i], &local_points[3 * i]);
        std::vector<double> local_values;
        ComputeLocalField(local_points, count, local_values);
        if (field_type_ == "phi") {
            global_values = std::move(local_values);
            return;
        }
        global_values.resize(3 * count);
        for (std::size_t i = 0; i < count; ++i)
            TransformToGlobal(&local_values[3 * i], &global_values[3 * i]);
    }

    bool ReadCache(const double global[3], ngbla::FlatVector<> result) const {
        if (!use_cache_) return false;
        const auto found = point_cache_.find(HashPoint(global));
        if (found == point_cache_.end()) {
            ++cache_misses_;
            return false;
        }
        ++cache_hits_;
        for (int component = 0; component < Dimension(); ++component)
            result(component) = found->second[component];
        return true;
    }

    bool ReadCache(const ngfem::BaseMappedIntegrationRule& rule,
                   ngbla::BareSliceMatrix<> result) const {
        for (std::size_t i = 0; i < rule.Size(); ++i) {
            const auto mapped = rule[i].GetPoint();
            const int dimension = mapped.Size();
            double global[3] = {
                mapped[0], dimension >= 2 ? mapped[1] : 0.0,
                dimension >= 3 ? mapped[2] : 0.0};
            const auto found = point_cache_.find(HashPoint(global));
            if (found == point_cache_.end()) {
                ++cache_misses_;
                return false;
            }
            ++cache_hits_;
            for (int component = 0; component < Dimension(); ++component)
                result(i, component) = found->second[component];
        }
        return true;
    }

    int object_ = 0;
    std::string field_type_;
    std::array<double, 3> origin_{0.0, 0.0, 0.0};
    std::array<double, 3> u_axis_{1.0, 0.0, 0.0};
    std::array<double, 3> v_axis_{0.0, 1.0, 0.0};
    std::array<double, 3> w_axis_{0.0, 0.0, 1.0};
    bool use_transform_ = false;
    std::optional<double> precision_;
    std::unordered_map<std::uint64_t, std::array<double, 3>> point_cache_;
    bool use_cache_ = false;
    double cache_tolerance_ = 1.0e-10;
    mutable std::atomic<std::size_t> cache_hits_{0};
    mutable std::atomic<std::size_t> cache_misses_{0};
};

} // namespace radia::ngsolve_bridge

#endif
