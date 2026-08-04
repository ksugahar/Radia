/* rad_hdiv_field_evaluator.h -- persistent BDM1 HDiv-VIM field evaluator.
 *
 * The solve owns one immutable evaluator.  Source packing and the source tree
 * are therefore paid once, while arbitrary observation batches can be applied
 * without rebuilding Python lists or charge geometry.
 */
#ifndef RAD_HDIV_FIELD_EVALUATOR_H
#define RAD_HDIV_FIELD_EVALUATOR_H

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace rad_hdiv {

struct FieldEvaluatorOptions {
    int leaf_size = 32;
    double theta = 0.05;
    std::size_t tree_min_sources = 256;
    std::size_t auto_min_work = 500000000;
    double tree_relative_tolerance = 1.0e-5;
    int probe_count = 16;
};

class HDivFieldEvaluator {
public:
    enum class Algorithm { Auto, Direct, Tree };

    static std::shared_ptr<HDivFieldEvaluator> FromTet(
        std::vector<double> volume,
        std::vector<double> surface,
        std::vector<int> image_masks,
        std::vector<double> image_signs,
        const FieldEvaluatorOptions& options = {});
    static std::shared_ptr<HDivFieldEvaluator> FromPolynomialTet(
        std::vector<double> volume,
        std::vector<double> surface,
        std::vector<int> image_masks,
        std::vector<double> image_signs,
        const FieldEvaluatorOptions& options = {});

    static std::shared_ptr<HDivFieldEvaluator> FromCloud(
        std::vector<double> xyz,
        std::vector<double> strength,
        std::vector<int> image_masks,
        std::vector<double> image_signs,
        const FieldEvaluatorOptions& options = {});

    // Curved P2 tetrahedra retain their geometry and BDM1/BDM2 reference-charge
    // polynomials.  Direct leaves integrate the actual element instead of a
    // permanently sampled point cloud; tree nodes use moments built once.
    // volume: [30 P2-node coordinates, 4 coefficients (1,xi,eta,zeta)]
    // surface: [18 P2-node coordinates, 6 coefficients
    //           (1,eta,eta^2,xi,xi*eta,xi^2)]
    static std::shared_ptr<HDivFieldEvaluator> FromCurvedTet(
        std::vector<double> volume,
        std::vector<double> surface,
        std::vector<double> gauss_points,
        std::vector<double> gauss_weights,
        std::vector<int> image_masks,
        std::vector<double> image_signs,
        const FieldEvaluatorOptions& options = {});

    ~HDivFieldEvaluator();
    HDivFieldEvaluator(const HDivFieldEvaluator&) = delete;
    HDivFieldEvaluator& operator=(const HDivFieldEvaluator&) = delete;

    void Evaluate(const double* observations, std::size_t n_observations,
                  double* output, Algorithm algorithm = Algorithm::Auto) const;

    // Evaluate without opening a TaskManager region.  NGSolve calls
    // CoefficientFunction::Evaluate from its existing worker threads, so the
    // coefficient-function adapter must not create a nested parallel region.
    void EvaluateSerial(const double* observations, std::size_t n_observations,
                        double* output, Algorithm algorithm = Algorithm::Auto) const;

    Algorithm AlgorithmFor(std::size_t n_observations) const;
    static Algorithm ParseAlgorithm(const std::string& name);
    static const char* AlgorithmName(Algorithm algorithm);

    std::size_t SourceCount() const;
    std::size_t ImageCount() const;
    std::size_t TreeNodeCount() const;
    int LeafSize() const;
    double Theta() const;
    std::size_t TreeMinSources() const;
    std::size_t AutoMinWork() const;
    double TreeRelativeTolerance() const;
    int ProbeCount() const;
    const char* SourceRepresentation() const;
    Algorithm LastAlgorithm() const;
    void Bounds(double lower[3], double upper[3]) const;

private:
    struct Impl;
    explicit HDivFieldEvaluator(std::unique_ptr<Impl> impl);
    std::unique_ptr<Impl> m_impl;
};

} // namespace rad_hdiv

#endif
