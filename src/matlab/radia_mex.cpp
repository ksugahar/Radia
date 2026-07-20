#include "mex.h"

#include "rad_evrs_tmethod.h"
#include "rad_biot_savart_filaments.h"
#include "rad_biot_savart_surface.h"
#include "rad_bem_galerkin.h"
#include "rad_cln.h"
#include "rad_hacapk_bem.h"
#include "rad_hacapk_hdiv.h"
#include "rad_hacapk_peec.h"
#include "rad_hdiv_hysteresis.h"
#include "rad_hdiv_field_evaluator.h"
#include "rad_hdiv_vim.h"
#include "rad_hybrid_vim_schur.h"
#include "rad_ngsolve_field_coefficients.h"
#include "rad_ngsolve_operators.h"
#include "rad_ngsolve_radia_field.h"
#include "rad_planar_charges.h"
#include "rad_average_field.h"
#include "rad_equivalence_source.h"
#include "rad_stream_function.h"

#include <core/taskmanager.hpp>
#include <bdbequations.hpp>
#include <bilinearform.hpp>
#include <cg.hpp>
#include <coefficient.hpp>
#include <flags.hpp>
#include <h1hofespace.hpp>
#include <hcurl_equations.hpp>
#include <hcurlhofespace.hpp>
#include <hdiv_equations.hpp>
#include <hdivhofespace.hpp>
#include <gridfunction.hpp>
#include <linearform.hpp>
#include <meshaccess.hpp>
#include <postproc.hpp>
#include <sparsematrix.hpp>
#include <symbolicintegrator.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <climits>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// The legacy C API header defines short macros (EXP, CALL, OK), so include it
// only after the NGSolve and standard-library headers have been parsed.
#include "radentry.h"

namespace {

using Complex = std::complex<double>;
using EnergyStopMaterial = rad_hdiv::EnergyStopMaterial;
using HACApKBEMManager = RadHACApKBEMManager;
using HACApKPEECManager = RadHACApKPEECManager;
using HACApKChargeGram = RadHACApKChargeGram;
using HACApKChargeGramDerivative = RadHACApKChargeGramDerivative;
using HDivFieldEvaluator = rad_hdiv::HDivFieldEvaluator;
using PlanarFieldEvaluator = rad_planar_charges::PlanarFieldEvaluator;

struct BEMHandle {
    std::vector<double> coordinates;
    std::vector<double> entries;
    std::unique_ptr<HACApKBEMManager> manager;
};

struct PEECHandle {
    std::unique_ptr<radia::PEECMatrixBuilder> builder;
    std::unique_ptr<HACApKPEECManager> manager;
};

struct ChargeGramHandle {
    std::shared_ptr<HACApKChargeGram> manager;
    int n_dof = 0;
};

struct ChargeGramDerivativeHandle {
    std::shared_ptr<HACApKChargeGram> parent;
    std::unique_ptr<HACApKChargeGramDerivative> manager;
    int n_dof = 0;
};

struct NGSolveCoefficientHandle {
    std::shared_ptr<ngfem::CoefficientFunction> coefficient;
};

struct NGSolveGridFunctionHandle {
    std::shared_ptr<ngcomp::MeshAccess> mesh;
    std::shared_ptr<ngcomp::FESpace> fespace;
    std::shared_ptr<ngcomp::GridFunction> gridfunction;
    std::string space;
    int order = 0;
    bool nograds = true;
};

struct NGSolveLinearFormHandle {
    std::shared_ptr<ngcomp::MeshAccess> mesh;
    std::shared_ptr<ngcomp::FESpace> fespace;
    std::shared_ptr<ngcomp::LinearForm> form;
    std::shared_ptr<ngla::BaseVector> vector;
    std::shared_ptr<ngfem::CoefficientFunction> coefficient;
    std::string space;
    std::string source_name;
    std::string label;
    Complex source_value = Complex(1.0, 0.0);
};

struct NGSolveVectorHandle {
    std::shared_ptr<ngla::BaseVector> vector;
    // A GridFunction view must retain its owner while MATLAB holds the vector.
    std::shared_ptr<ngcomp::GridFunction> parent_gridfunction;
    // A LinearForm vector must retain its assembled form while MATLAB holds it.
    std::shared_ptr<ngcomp::LinearForm> parent_linear_form;
    // Vectors created from a native matrix retain that matrix for the same
    // lifetime guarantee when MATLAB keeps only the vector handle.
    std::shared_ptr<ngla::BaseMatrix> parent_matrix;
    // A solver result retains the solver and its matrix/preconditioner graph.
    std::shared_ptr<ngla::KrylovSpaceSolver> parent_solver;
    int component = 0;
    bool is_view = false;
};

struct NGSolveMeshHandle {
    std::shared_ptr<ngcomp::MeshAccess> mesh;
    std::string path;
};

struct NGSolveFESpaceHandle {
    std::shared_ptr<ngcomp::MeshAccess> mesh;
    std::shared_ptr<ngcomp::FESpace> fespace;
    std::string space;
    int order = 0;
    bool nograds = true;
    bool is_complex = false;
};

struct NGSolveBilinearFormHandle {
    std::shared_ptr<ngcomp::MeshAccess> mesh;
    std::shared_ptr<ngcomp::FESpace> fespace;
    std::shared_ptr<ngcomp::BilinearForm> form;
    std::shared_ptr<ngla::BaseMatrix> matrix;
    std::shared_ptr<ngfem::CoefficientFunction> coefficient;
    std::string space;
    std::string form_name;
    std::string label;
};

struct NGSolveMatrixHandle {
    std::shared_ptr<ngla::BaseMatrix> matrix;
    std::shared_ptr<ngcomp::FESpace> fespace;
    std::string kind;
};

struct NGSolveSolverHandle {
    std::shared_ptr<ngla::KrylovSpaceSolver> solver;
    std::shared_ptr<ngla::BaseMatrix> matrix;
    std::shared_ptr<ngla::BaseMatrix> preconditioner;
    std::string method;
    double tolerance = 1e-8;
    int max_steps = 1000;
};

struct NativeStateSpaceHandle {
    std::size_t state_size = 0;
    std::size_t input_size = 0;
    std::size_t output_size = 0;
    std::vector<double> A;
    std::vector<double> B;
    std::vector<double> C;
    std::vector<double> D;
    std::vector<double> initial_state;
    std::vector<double> state;
    std::size_t step_count = 0;
};

std::mutex registry_mutex;
std::unordered_map<std::uint64_t, std::unique_ptr<EnergyStopMaterial>> energy_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<BEMHandle>> bem_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<PEECHandle>> peec_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<ChargeGramHandle>> charge_gram_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<ChargeGramDerivativeHandle>>
    charge_gram_derivative_registry;
std::unordered_map<std::uint64_t, std::shared_ptr<HDivFieldEvaluator>> field_registry;
std::unordered_map<std::uint64_t, std::shared_ptr<PlanarFieldEvaluator>> planar_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<NGSolveCoefficientHandle>>
    coefficient_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<NGSolveGridFunctionHandle>>
    gridfunction_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<NGSolveLinearFormHandle>>
    linear_form_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<NGSolveVectorHandle>>
    vector_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<NGSolveMeshHandle>>
    mesh_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<NGSolveFESpaceHandle>>
    fespace_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<NGSolveBilinearFormHandle>>
    bilinear_form_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<NGSolveMatrixHandle>>
    matrix_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<NGSolveSolverHandle>>
    solver_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<NativeStateSpaceHandle>>
    state_space_registry;
std::uint64_t next_handle = 1;
std::size_t lock_count = 0;
bool exit_handler_registered = false;

extern "C" {
    void cHACApK_hlu_get_timings(double* out_t_decomp, double* out_t_solve,
                                 long* out_n_dense_lu, long* out_n_dense_gemm);
    void cHACApK_hlu_set_trunc_tol(double tol);
    double cHACApK_hlu_get_trunc_tol(void);
    void cHACApK_hlu_get_materialize_stats(long* out_n_calls, long* out_n_elems);
    void cHACApK_hlu_get_materialize_split(long* out_internal, long* out_leaf);
    void cHACApK_hlu_get_mixed_breakdown(long* out_addmul9, long* out_lln9,
                                         long* out_run9);
    void cHACApK_hlu_set_parallel(int on);
    int cHACApK_hlu_get_parallel(void);
    void cHACApK_hlu_set_par_cutoff(long cutoff);
    int chacapk_max_threads(void);
    void cHACApK_hlu_set_accum_cap(int cap);
    int cHACApK_hlu_get_accum_cap(void);
    int cHACApK_get_cluster_strategy(void);
    double cHACApK_harith_self_test(int depth, int n_per_block);
    double cHACApK_harith_self_test_rk(int n_per_block, int rk_rank);
    double cHACApK_harith_self_test_addmul_rkrk(int m, int n, int inner,
                                                int kA, int kB, int kC);
    double cHACApK_harith_self_test_radia_exact_with_matrix(
        const double* A_full, const double* b);
    double cHACApK_harith_self_test_radia_exact_diag(double diag_boost);
    double cHACApK_harith_self_test_radia_exact(void);
    double cHACApK_harith_self_test_depth3_asymmetric(int nb_tiny);
    double cHACApK_harith_self_test_mixed_sibling_via_conversion(int nb_small);
    double cHACApK_harith_self_test_mixed_sibling_nonuniform(
        int n1, int n2, int m1, int m3);
    double cHACApK_harith_self_test_mixed_sibling(int nb_small);
    double cHACApK_harith_self_test_rk_deep(int n_per_block, int rk_rank);
}

[[noreturn]] void BadArgument(const std::string& message) {
    throw std::invalid_argument(message);
}

void CheckArity(int nrhs, int expected_rhs, int nlhs, int expected_lhs,
                const char* usage) {
    if (nrhs != expected_rhs || nlhs != expected_lhs)
        BadArgument(std::string("usage: ") + usage);
}

std::string Text(const mxArray* value, const char* name) {
    if (!mxIsChar(value))
        BadArgument(std::string(name) + " must be a character vector");
    char* raw = mxArrayToString(value);
    if (raw == nullptr)
        throw std::runtime_error(std::string("could not decode ") + name);
    std::string result(raw);
    mxFree(raw);
    return result;
}

double Scalar(const mxArray* value, const char* name) {
    if (!mxIsDouble(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) != 1)
        BadArgument(std::string(name) + " must be a real double scalar");
    const double result = mxGetScalar(value);
    if (!std::isfinite(result))
        BadArgument(std::string(name) + " must be finite");
    return result;
}

int PositiveInteger(const mxArray* value, const char* name) {
    const double raw = Scalar(value, name);
    if (raw < 1.0 || raw != std::floor(raw) ||
        raw > static_cast<double>(std::numeric_limits<int>::max()))
        BadArgument(std::string(name) + " must be a positive integer");
    return static_cast<int>(raw);
}

int IntegerScalar(const mxArray* value, const char* name) {
    const double raw = Scalar(value, name);
    if (raw != std::floor(raw) ||
        raw < static_cast<double>(std::numeric_limits<int>::min()) ||
        raw > static_cast<double>(std::numeric_limits<int>::max()))
        BadArgument(std::string(name) + " must be an integer scalar");
    return static_cast<int>(raw);
}

long NonnegativeLong(const mxArray* value, const char* name) {
    const double raw = Scalar(value, name);
    if (raw < 0.0 || raw != std::floor(raw) ||
        raw > static_cast<double>(std::numeric_limits<long>::max()))
        BadArgument(std::string(name) + " must be a nonnegative integer scalar");
    return static_cast<long>(raw);
}

int NonnegativeInteger(const mxArray* value, const char* name) {
    const long raw = NonnegativeLong(value, name);
    if (raw > static_cast<long>(std::numeric_limits<int>::max()))
        BadArgument(std::string(name) + " is too large for the native integer contract");
    return static_cast<int>(raw);
}

bool Boolean(const mxArray* value, const char* name) {
    if (mxIsLogicalScalar(value))
        return mxIsLogicalScalarTrue(value);
    const double raw = Scalar(value, name);
    if (raw != 0.0 && raw != 1.0)
        BadArgument(std::string(name) + " must be logical or 0/1");
    return raw != 0.0;
}

std::vector<double> RealVector(const mxArray* value, const char* name) {
    if (!mxIsDouble(value) || mxIsComplex(value))
        BadArgument(std::string(name) + " must be a real double array");
    const std::size_t count = mxGetNumberOfElements(value);
    const double* data = mxGetDoubles(value);
    return std::vector<double>(data, data + count);
}

std::vector<int> IntegerVector(const mxArray* value, const char* name) {
    const std::size_t count = mxGetNumberOfElements(value);
    std::vector<int> result(count);
    if (mxIsInt32(value) && !mxIsComplex(value)) {
        const auto* data = static_cast<const std::int32_t*>(mxGetData(value));
        std::copy(data, data + count, result.begin());
        return result;
    }
    if (!mxIsDouble(value) || mxIsComplex(value))
        BadArgument(std::string(name) + " must be int32 or real double");
    const double* data = mxGetDoubles(value);
    for (std::size_t i = 0; i < count; ++i) {
        if (!std::isfinite(data[i]) || data[i] != std::floor(data[i]) ||
            data[i] < static_cast<double>(std::numeric_limits<int>::min()) ||
            data[i] > static_cast<double>(std::numeric_limits<int>::max()))
            BadArgument(std::string(name) + " must contain integer values");
        result[i] = static_cast<int>(data[i]);
    }
    return result;
}

std::vector<int> IntegerMatrix(const mxArray* value, std::size_t& rows,
                               std::size_t& cols, const char* name) {
    if (mxGetNumberOfDimensions(value) != 2)
        BadArgument(std::string(name) + " must be a two-dimensional integer matrix");
    rows = mxGetM(value);
    cols = mxGetN(value);
    std::vector<int> result(rows * cols);
    if (mxIsInt32(value) && !mxIsComplex(value)) {
        const auto* data = static_cast<const std::int32_t*>(mxGetData(value));
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result[i * cols + j] = data[i + j * rows];
        return result;
    }
    if (!mxIsDouble(value) || mxIsComplex(value))
        BadArgument(std::string(name) + " must be int32 or real double");
    const double* data = mxGetDoubles(value);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < cols; ++j) {
            const double item = data[i + j * rows];
            if (!std::isfinite(item) || item != std::floor(item) ||
                item < static_cast<double>(std::numeric_limits<int>::min()) ||
                item > static_cast<double>(std::numeric_limits<int>::max()))
                BadArgument(std::string(name) + " must contain integer values");
            result[i * cols + j] = static_cast<int>(item);
        }
    return result;
}

std::vector<std::int64_t> Integer64Matrix(const mxArray* value,
                                          std::size_t& rows,
                                          std::size_t& cols,
                                          const char* name) {
    if (mxGetNumberOfDimensions(value) != 2)
        BadArgument(std::string(name) + " must be a two-dimensional integer matrix");
    rows = mxGetM(value);
    cols = mxGetN(value);
    std::vector<std::int64_t> result(rows * cols);
    if (mxIsInt64(value) && !mxIsComplex(value)) {
        const auto* data = static_cast<const std::int64_t*>(mxGetData(value));
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result[i * cols + j] = data[i + j * rows];
        return result;
    }
    if (mxIsInt32(value) && !mxIsComplex(value)) {
        const auto* data = static_cast<const std::int32_t*>(mxGetData(value));
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result[i * cols + j] = data[i + j * rows];
        return result;
    }
    if (!mxIsDouble(value) || mxIsComplex(value))
        BadArgument(std::string(name) + " must be int64, int32, or real double");
    const double* data = mxGetDoubles(value);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < cols; ++j) {
            const double item = data[i + j * rows];
            if (!std::isfinite(item) || item != std::floor(item) ||
                item < static_cast<double>(std::numeric_limits<std::int64_t>::min()) ||
                item > static_cast<double>(std::numeric_limits<std::int64_t>::max()))
                BadArgument(std::string(name) + " must contain integer values");
            result[i * cols + j] = static_cast<std::int64_t>(item);
        }
    return result;
}

std::vector<double> FixedRealVector(const mxArray* value, std::size_t size,
                                    const char* name) {
    auto result = RealVector(value, name);
    if (result.size() != size)
        BadArgument(std::string(name) + " must have " + std::to_string(size) + " entries");
    return result;
}

std::vector<double> RealMatrix(const mxArray* value, std::size_t& rows,
                               std::size_t& cols, const char* name) {
    if (!mxIsDouble(value) || mxIsComplex(value) || mxGetNumberOfDimensions(value) != 2)
        BadArgument(std::string(name) + " must be a real double matrix");
    rows = mxGetM(value);
    cols = mxGetN(value);
    const double* data = mxGetDoubles(value);
    std::vector<double> result(rows * cols);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < cols; ++j)
            result[i * cols + j] = data[i + j * rows];
    return result;
}

std::vector<double> RealTensor3(const mxArray* value, std::size_t& dim0,
                                std::size_t& dim1, std::size_t& dim2,
                                const char* name) {
    if (!mxIsDouble(value) || mxIsComplex(value) ||
        mxGetNumberOfDimensions(value) != 3)
        BadArgument(std::string(name) + " must be a real three-dimensional array");
    const mwSize* dims = mxGetDimensions(value);
    dim0 = static_cast<std::size_t>(dims[0]);
    dim1 = static_cast<std::size_t>(dims[1]);
    dim2 = static_cast<std::size_t>(dims[2]);
    const double* data = mxGetDoubles(value);
    std::vector<double> result(dim0 * dim1 * dim2);
    for (std::size_t i = 0; i < dim0; ++i)
        for (std::size_t j = 0; j < dim1; ++j)
            for (std::size_t k = 0; k < dim2; ++k)
                result[(i * dim1 + j) * dim2 + k] =
                    data[i + dim0 * (j + dim1 * k)];
    return result;
}

std::vector<double> RealTensor4(const mxArray* value, std::size_t& dim0,
                                std::size_t& dim1, std::size_t& dim2,
                                std::size_t& dim3, const char* name) {
    if (!mxIsDouble(value) || mxIsComplex(value) ||
        mxGetNumberOfDimensions(value) != 4)
        BadArgument(std::string(name) + " must be a real four-dimensional array");
    const mwSize* dims = mxGetDimensions(value);
    dim0 = static_cast<std::size_t>(dims[0]);
    dim1 = static_cast<std::size_t>(dims[1]);
    dim2 = static_cast<std::size_t>(dims[2]);
    dim3 = static_cast<std::size_t>(dims[3]);
    const double* data = mxGetDoubles(value);
    std::vector<double> result(dim0 * dim1 * dim2 * dim3);
    for (std::size_t i = 0; i < dim0; ++i)
        for (std::size_t j = 0; j < dim1; ++j)
            for (std::size_t k = 0; k < dim2; ++k)
                for (std::size_t l = 0; l < dim3; ++l)
                    result[((i * dim1 + j) * dim2 + k) * dim3 + l] =
                        data[i + dim0 * (j + dim1 * (k + dim2 * l))];
    return result;
}

void ReadQuadrature(const mxArray* point_value, const mxArray* weight_value,
                    std::size_t dimension, bool required,
                    const char* point_name, const char* weight_name,
                    std::vector<double>& points,
                    std::vector<double>& weights) {
    std::size_t rows = 0, cols = 0;
    points = RealMatrix(point_value, rows, cols, point_name);
    weights = RealVector(weight_value, weight_name);
    if (points.empty() && weights.empty()) {
        if (required)
            BadArgument(std::string(point_name) + " and " + weight_name + " may not be empty");
        return;
    }
    if (rows == 0 || cols != dimension || weights.size() != rows)
        BadArgument(std::string(point_name) + " must have one row per " + weight_name + " entry");
    for (double point : points)
        if (!std::isfinite(point))
            BadArgument(std::string(point_name) + " must contain finite values");
    for (double weight : weights)
        if (!std::isfinite(weight) || weight <= 0.0)
            BadArgument(std::string(weight_name) + " must contain positive finite values");
}

void ReadRule1D(const mxArray* point_value, const mxArray* weight_value,
                const char* point_name, const char* weight_name,
                std::vector<double>& points, std::vector<double>& weights) {
    points = RealVector(point_value, point_name);
    weights = RealVector(weight_value, weight_name);
    if (points.empty() || points.size() != weights.size())
        BadArgument(std::string(point_name) + " and " + weight_name +
                    " must be nonempty vectors of equal length");
    for (double point : points)
        if (!std::isfinite(point) || point < 0.0 || point > 1.0)
            BadArgument(std::string(point_name) + " points must be finite and lie in [0,1]");
    for (double weight : weights)
        if (!std::isfinite(weight) || weight <= 0.0)
            BadArgument(std::string(weight_name) + " weights must be positive and finite");
}

void ValidateChargeDescriptors(const std::vector<int>& charge_host,
                               const std::vector<int>& charge_kind,
                               const std::vector<int>& charge_expo,
                               std::size_t exponent_rows,
                               std::size_t exponent_cols,
                               int n_cell, int n_face) {
    const std::size_t n_charge = charge_host.size();
    if (n_charge == 0 || charge_kind.size() != n_charge ||
        exponent_rows != n_charge || exponent_cols != 3)
        BadArgument("charge_host/kind must have n_charge entries and charge_expo must be n_charge-by-3");
    if (n_charge > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("n_charge exceeds the native integer range");
    for (std::size_t charge = 0; charge < n_charge; ++charge) {
        const int kind = charge_kind[charge];
        const int host = charge_host[charge];
        if (kind != 0 && kind != 1)
            BadArgument("charge_kind entries must be 0 (cell) or 1 (face)");
        const int host_count = kind == 0 ? n_cell : n_face;
        if (host < 0 || host >= host_count)
            BadArgument("charge_host contains an index outside its cell/face host range");
        const int ex = charge_expo[3 * charge];
        const int ey = charge_expo[3 * charge + 1];
        const int ez = charge_expo[3 * charge + 2];
        const long long total_degree = static_cast<long long>(ex) + ey + ez;
        if (ex < 0 || ey < 0 || ez < 0 || total_degree > 18)
            BadArgument("charge_expo degrees must be nonnegative with total degree at most 18");
        if (kind == 1 && ez != 0)
            BadArgument("face charge_expo must have a zero third exponent");
    }
}

Complex ComplexScalar(const mxArray* value, const char* name) {
    if (!mxIsDouble(value) || mxGetNumberOfElements(value) != 1)
        BadArgument(std::string(name) + " must be a double scalar");
    if (!mxIsComplex(value))
        return Complex(mxGetScalar(value), 0.0);
    const mxComplexDouble* data = mxGetComplexDoubles(value);
    return Complex(data[0].real, data[0].imag);
}

std::vector<Complex> ComplexMatrix(const mxArray* value, std::size_t& rows,
                                   std::size_t& cols, const char* name) {
    if (!mxIsDouble(value) || mxGetNumberOfDimensions(value) != 2)
        BadArgument(std::string(name) + " must be a double matrix");
    rows = mxGetM(value);
    cols = mxGetN(value);
    std::vector<Complex> result(rows * cols);
    if (mxIsComplex(value)) {
        const mxComplexDouble* data = mxGetComplexDoubles(value);
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j) {
                const auto& item = data[i + j * rows];
                result[i * cols + j] = Complex(item.real, item.imag);
            }
    } else {
        const double* data = mxGetDoubles(value);
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result[i * cols + j] = Complex(data[i + j * rows], 0.0);
    }
    return result;
}

mxArray* RealRow(const std::vector<double>& values) {
    mxArray* result = mxCreateDoubleMatrix(1, values.size(), mxREAL);
    std::copy(values.begin(), values.end(), mxGetDoubles(result));
    return result;
}

mxArray* ComplexMatrixOutput(const std::vector<Complex>& values,
                             std::size_t rows, std::size_t cols);

mxArray* ComplexRow(const std::vector<Complex>& values) {
    return ComplexMatrixOutput(values, 1, values.size());
}

mxArray* RealColumn(const std::vector<double>& values) {
    mxArray* result = mxCreateDoubleMatrix(values.size(), 1, mxREAL);
    std::copy(values.begin(), values.end(), mxGetDoubles(result));
    return result;
}

mxArray* RealMatrixOutput(const std::vector<double>& values,
                          std::size_t rows, std::size_t cols) {
    mxArray* result = mxCreateDoubleMatrix(rows, cols, mxREAL);
    double* data = mxGetDoubles(result);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < cols; ++j)
            data[i + j * rows] = values[i * cols + j];
    return result;
}

mxArray* TextOutput(const char* value) {
    return mxCreateString(value ? value : "");
}

mxArray* PairStructOutput(const std::vector<std::pair<std::string, double>>& values) {
    std::vector<const char*> fields;
    fields.reserve(values.size());
    for (const auto& item : values) fields.push_back(item.first.c_str());
    mxArray* result = mxCreateStructMatrix(1, 1, fields.size(), fields.data());
    for (const auto& item : values)
        mxSetField(result, 0, item.first.c_str(), mxCreateDoubleScalar(item.second));
    return result;
}

rad_hdiv::FieldEvaluatorOptions FieldOptions(int nrhs, const mxArray* prhs[],
                                             int first_option, const char* usage) {
    const int count = nrhs - first_option;
    if (count != 0 && count != 6) BadArgument(std::string("usage: ") + usage);
    rad_hdiv::FieldEvaluatorOptions options;
    if (count == 6) {
        options.leaf_size = PositiveInteger(prhs[first_option], "leaf_size");
        options.theta = Scalar(prhs[first_option + 1], "theta");
        options.tree_min_sources = static_cast<std::size_t>(
            NonnegativeLong(prhs[first_option + 2], "tree_min_sources"));
        options.auto_min_work = static_cast<std::size_t>(
            NonnegativeLong(prhs[first_option + 3], "auto_min_work"));
        options.tree_relative_tolerance = Scalar(
            prhs[first_option + 4], "tree_relative_tolerance");
        options.probe_count = NonnegativeInteger(prhs[first_option + 5], "probe_count");
    }
    return options;
}

mxArray* ComplexScalarOutput(Complex value) {
    mxArray* result = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    mxComplexDouble* data = mxGetComplexDoubles(result);
    data[0].real = value.real();
    data[0].imag = value.imag();
    return result;
}

mxArray* ComplexMatrixOutput(const std::vector<Complex>& values,
                             std::size_t rows, std::size_t cols) {
    mxArray* result = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
    mxComplexDouble* data = mxGetComplexDoubles(result);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < cols; ++j) {
            const Complex value = values[i * cols + j];
            data[i + j * rows].real = value.real();
            data[i + j * rows].imag = value.imag();
        }
    return result;
}

mxArray* Uint64Output(std::uint64_t value) {
    mxArray* result = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *static_cast<std::uint64_t*>(mxGetData(result)) = value;
    return result;
}

int MatrixDimension(std::size_t value, const char* name) {
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument(std::string(name) + " is too large");
    return static_cast<int>(value);
}

std::uint64_t Handle(const mxArray* value) {
    if (!mxIsUint64(value) || mxIsComplex(value) || mxGetNumberOfElements(value) != 1)
        BadArgument("handle must be a uint64 scalar");
    return *static_cast<const std::uint64_t*>(mxGetData(value));
}

std::vector<std::uint64_t> HandleVector(const mxArray* value,
                                        const char* name) {
    if (!mxIsUint64(value) || mxIsComplex(value))
        BadArgument(std::string(name) + " must be a uint64 array");
    const std::size_t count = mxGetNumberOfElements(value);
    const auto* data = static_cast<const std::uint64_t*>(mxGetData(value));
    return std::vector<std::uint64_t>(data, data + count);
}

void CheckRadia(int error_code) {
    if (error_code == 0)
        return;
    char text[2048] = {0};
    RadErrGetText(text, error_code);
    throw std::runtime_error(std::string("Radia error: ") + text);
}

void Cleanup() {
    std::lock_guard<std::mutex> guard(registry_mutex);
    energy_registry.clear();
    bem_registry.clear();
    peec_registry.clear();
    charge_gram_registry.clear();
    charge_gram_derivative_registry.clear();
    field_registry.clear();
    planar_registry.clear();
    coefficient_registry.clear();
    gridfunction_registry.clear();
    linear_form_registry.clear();
    vector_registry.clear();
    mesh_registry.clear();
    fespace_registry.clear();
    bilinear_form_registry.clear();
    matrix_registry.clear();
    solver_registry.clear();
    state_space_registry.clear();
    while (lock_count > 0) {
        mexUnlock();
        --lock_count;
    }
}

bool HandleInUse(std::uint64_t handle) {
    return handle == 0 || energy_registry.count(handle) != 0 ||
           bem_registry.count(handle) != 0 || peec_registry.count(handle) != 0 ||
           charge_gram_registry.count(handle) != 0 ||
           charge_gram_derivative_registry.count(handle) != 0 ||
           field_registry.count(handle) != 0 ||
           planar_registry.count(handle) != 0 || coefficient_registry.count(handle) != 0 ||
           gridfunction_registry.count(handle) != 0 || vector_registry.count(handle) != 0 ||
           linear_form_registry.count(handle) != 0 ||
           mesh_registry.count(handle) != 0 || fespace_registry.count(handle) != 0 ||
           bilinear_form_registry.count(handle) != 0 || matrix_registry.count(handle) != 0 ||
           solver_registry.count(handle) != 0 || state_space_registry.count(handle) != 0;
}

std::uint64_t RegisterBEM(std::unique_ptr<BEMHandle> bem) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    bem_registry.emplace(handle, std::move(bem));
    mexLock();
    ++lock_count;
    return handle;
}

void EnsureExitHandler() {
    if (!exit_handler_registered) {
        mexAtExit(Cleanup);
        exit_handler_registered = true;
    }
}

std::uint64_t Register(std::unique_ptr<EnergyStopMaterial> material) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    energy_registry.emplace(handle, std::move(material));
    mexLock();
    ++lock_count;
    return handle;
}

EnergyStopMaterial& Energy(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = energy_registry.find(handle);
    if (found == energy_registry.end())
        BadArgument("invalid or stale EnergyStopMaterial handle");
    return *found->second;
}

void Destroy(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (energy_registry.erase(handle) == 0)
        BadArgument("invalid or stale EnergyStopMaterial handle");
    mexUnlock();
    --lock_count;
}

BEMHandle& BEM(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = bem_registry.find(handle);
    if (found == bem_registry.end())
        BadArgument("invalid or stale HACApK BEM handle");
    return *found->second;
}

void DestroyBEM(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (bem_registry.erase(handle) == 0)
        BadArgument("invalid or stale HACApK BEM handle");
    mexUnlock();
    --lock_count;
}

PEECHandle& PEEC(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = peec_registry.find(handle);
    if (found == peec_registry.end())
        BadArgument("invalid or stale HACApK PEEC handle");
    return *found->second;
}

std::uint64_t RegisterPEEC(std::unique_ptr<PEECHandle> peec) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    peec_registry.emplace(handle, std::move(peec));
    mexLock();
    ++lock_count;
    return handle;
}

void DestroyPEEC(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (peec_registry.erase(handle) == 0)
        BadArgument("invalid or stale HACApK PEEC handle");
    mexUnlock();
    --lock_count;
}

ChargeGramHandle& ChargeGram(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = charge_gram_registry.find(handle);
    if (found == charge_gram_registry.end())
        BadArgument("invalid or stale HACApK charge-Gram handle");
    return *found->second;
}

std::uint64_t RegisterChargeGram(std::unique_ptr<ChargeGramHandle> gram) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    charge_gram_registry.emplace(handle, std::move(gram));
    mexLock();
    ++lock_count;
    return handle;
}

void DestroyChargeGram(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (charge_gram_registry.erase(handle) == 0)
        BadArgument("invalid or stale HACApK charge-Gram handle");
    mexUnlock();
    --lock_count;
}

ChargeGramDerivativeHandle& ChargeGramDerivative(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = charge_gram_derivative_registry.find(handle);
    if (found == charge_gram_derivative_registry.end())
        BadArgument("invalid or stale HACApK charge-Gram derivative handle");
    return *found->second;
}

std::uint64_t RegisterChargeGramDerivative(
    std::unique_ptr<ChargeGramDerivativeHandle> derivative) {
    if (!derivative || !derivative->parent || !derivative->manager)
        BadArgument("charge-Gram derivative construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    charge_gram_derivative_registry.emplace(handle, std::move(derivative));
    mexLock();
    ++lock_count;
    return handle;
}

void DestroyChargeGramDerivative(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (charge_gram_derivative_registry.erase(handle) == 0)
        BadArgument("invalid or stale HACApK charge-Gram derivative handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterField(std::shared_ptr<HDivFieldEvaluator> evaluator) {
    if (!evaluator) BadArgument("field evaluator construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle)) ++next_handle;
    const std::uint64_t handle = next_handle++;
    field_registry.emplace(handle, std::move(evaluator));
    mexLock();
    ++lock_count;
    return handle;
}

std::shared_ptr<HDivFieldEvaluator> Field(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = field_registry.find(handle);
    if (found == field_registry.end())
        BadArgument("invalid or stale HDiv field evaluator handle");
    return found->second;
}

void DestroyField(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (field_registry.erase(handle) == 0)
        BadArgument("invalid or stale HDiv field evaluator handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterPlanar(std::shared_ptr<PlanarFieldEvaluator> evaluator) {
    if (!evaluator) BadArgument("planar field evaluator construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle)) ++next_handle;
    const std::uint64_t handle = next_handle++;
    planar_registry.emplace(handle, std::move(evaluator));
    mexLock();
    ++lock_count;
    return handle;
}

std::shared_ptr<PlanarFieldEvaluator> Planar(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = planar_registry.find(handle);
    if (found == planar_registry.end())
        BadArgument("invalid or stale planar field evaluator handle");
    return found->second;
}

void DestroyPlanar(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (planar_registry.erase(handle) == 0)
        BadArgument("invalid or stale planar field evaluator handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterCoefficient(
    std::shared_ptr<ngfem::CoefficientFunction> coefficient) {
    if (!coefficient)
        BadArgument("NGSolve CoefficientFunction construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    auto holder = std::make_unique<NGSolveCoefficientHandle>();
    holder->coefficient = std::move(coefficient);
    coefficient_registry.emplace(handle, std::move(holder));
    mexLock();
    ++lock_count;
    return handle;
}

std::shared_ptr<ngfem::CoefficientFunction> Coefficient(
    std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = coefficient_registry.find(handle);
    if (found == coefficient_registry.end())
        BadArgument("invalid or stale NGSolve CoefficientFunction handle");
    return found->second->coefficient;
}

void DestroyCoefficient(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (coefficient_registry.erase(handle) == 0)
        BadArgument("invalid or stale NGSolve CoefficientFunction handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterGridFunction(
    std::unique_ptr<NGSolveGridFunctionHandle> gridfunction) {
    if (!gridfunction || !gridfunction->gridfunction)
        BadArgument("NGSolve GridFunction construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    gridfunction_registry.emplace(handle, std::move(gridfunction));
    mexLock();
    ++lock_count;
    return handle;
}

NGSolveGridFunctionHandle& GridFunction(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = gridfunction_registry.find(handle);
    if (found == gridfunction_registry.end())
        BadArgument("invalid or stale NGSolve GridFunction handle");
    return *found->second;
}

void DestroyGridFunction(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (gridfunction_registry.erase(handle) == 0)
        BadArgument("invalid or stale NGSolve GridFunction handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterLinearForm(
    std::unique_ptr<NGSolveLinearFormHandle> linear_form) {
    if (!linear_form || !linear_form->form || !linear_form->vector)
        BadArgument("NGSolve LinearForm construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    linear_form_registry.emplace(handle, std::move(linear_form));
    mexLock();
    ++lock_count;
    return handle;
}

NGSolveLinearFormHandle& LinearForm(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = linear_form_registry.find(handle);
    if (found == linear_form_registry.end())
        BadArgument("invalid or stale NGSolve LinearForm handle");
    return *found->second;
}

void DestroyLinearForm(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (linear_form_registry.erase(handle) == 0)
        BadArgument("invalid or stale NGSolve LinearForm handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterVector(std::unique_ptr<NGSolveVectorHandle> vector) {
    if (!vector || !vector->vector)
        BadArgument("NGSolve vector construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    vector_registry.emplace(handle, std::move(vector));
    mexLock();
    ++lock_count;
    return handle;
}

NGSolveVectorHandle& Vector(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = vector_registry.find(handle);
    if (found == vector_registry.end())
        BadArgument("invalid or stale NGSolve vector handle");
    return *found->second;
}

void DestroyVector(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (vector_registry.erase(handle) == 0)
        BadArgument("invalid or stale NGSolve vector handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterMesh(std::unique_ptr<NGSolveMeshHandle> mesh) {
    if (!mesh || !mesh->mesh)
        BadArgument("NGSolve MeshAccess construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    mesh_registry.emplace(handle, std::move(mesh));
    mexLock();
    ++lock_count;
    return handle;
}

NGSolveMeshHandle& Mesh(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = mesh_registry.find(handle);
    if (found == mesh_registry.end())
        BadArgument("invalid or stale NGSolve MeshAccess handle");
    return *found->second;
}

void DestroyMesh(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (mesh_registry.erase(handle) == 0)
        BadArgument("invalid or stale NGSolve MeshAccess handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterFESpace(std::unique_ptr<NGSolveFESpaceHandle> fespace) {
    if (!fespace || !fespace->fespace)
        BadArgument("NGSolve FESpace construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    fespace_registry.emplace(handle, std::move(fespace));
    mexLock();
    ++lock_count;
    return handle;
}

NGSolveFESpaceHandle& FESpace(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = fespace_registry.find(handle);
    if (found == fespace_registry.end())
        BadArgument("invalid or stale NGSolve FESpace handle");
    return *found->second;
}

void DestroyFESpace(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (fespace_registry.erase(handle) == 0)
        BadArgument("invalid or stale NGSolve FESpace handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterBilinearForm(
    std::unique_ptr<NGSolveBilinearFormHandle> form) {
    if (!form || !form->form || !form->matrix)
        BadArgument("NGSolve BilinearForm construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    bilinear_form_registry.emplace(handle, std::move(form));
    mexLock();
    ++lock_count;
    return handle;
}

NGSolveBilinearFormHandle& BilinearForm(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = bilinear_form_registry.find(handle);
    if (found == bilinear_form_registry.end())
        BadArgument("invalid or stale NGSolve BilinearForm handle");
    return *found->second;
}

void DestroyBilinearForm(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (bilinear_form_registry.erase(handle) == 0)
        BadArgument("invalid or stale NGSolve BilinearForm handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterMatrix(std::unique_ptr<NGSolveMatrixHandle> matrix) {
    if (!matrix || !matrix->matrix)
        BadArgument("NGSolve BaseMatrix construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    matrix_registry.emplace(handle, std::move(matrix));
    mexLock();
    ++lock_count;
    return handle;
}

NGSolveMatrixHandle& Matrix(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = matrix_registry.find(handle);
    if (found == matrix_registry.end())
        BadArgument("invalid or stale NGSolve BaseMatrix handle");
    return *found->second;
}

void DestroyMatrix(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (matrix_registry.erase(handle) == 0)
        BadArgument("invalid or stale NGSolve BaseMatrix handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterSolver(std::unique_ptr<NGSolveSolverHandle> solver) {
    if (!solver || !solver->solver || !solver->matrix)
        BadArgument("NGSolve Krylov solver construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    solver_registry.emplace(handle, std::move(solver));
    mexLock();
    ++lock_count;
    return handle;
}

NGSolveSolverHandle& Solver(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = solver_registry.find(handle);
    if (found == solver_registry.end())
        BadArgument("invalid or stale NGSolve Solver handle");
    return *found->second;
}

void DestroySolver(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (solver_registry.erase(handle) == 0)
        BadArgument("invalid or stale NGSolve Solver handle");
    mexUnlock();
    --lock_count;
}

std::uint64_t RegisterStateSpace(
    std::unique_ptr<NativeStateSpaceHandle> state_space) {
    if (!state_space)
        BadArgument("native Simulink state-space construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    state_space_registry.emplace(handle, std::move(state_space));
    mexLock();
    ++lock_count;
    return handle;
}

NativeStateSpaceHandle& StateSpace(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = state_space_registry.find(handle);
    if (found == state_space_registry.end())
        BadArgument("invalid or stale native Simulink state-space handle");
    return *found->second;
}

void DestroyStateSpace(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (state_space_registry.erase(handle) == 0)
        BadArgument("invalid or stale native Simulink state-space handle");
    mexUnlock();
    --lock_count;
}

mxArray* Commands() {
    static const char* names[] = {
        "api.info", "api.commands", "taskmanager.probe", "ngsolve.space_info",
        "ngsolve.matrix_dump",
        "ngsolve.mesh.create", "ngsolve.mesh.info",
        "ngsolve.mesh.set_deformation", "ngsolve.mesh.unset_deformation",
        "ngsolve.mesh.trafo_quality", "ngsolve.mesh.destroy",
        "ngsolve.fespace.create", "ngsolve.fespace.info",
        "ngsolve.fespace.destroy",
        "ngsolve.bilinear_form.create", "ngsolve.bilinear_form.info",
        "ngsolve.bilinear_form.create_from_coefficient",
        "ngsolve.bilinear_form.create_boundary_from_coefficient",
        "ngsolve.bilinear_form.matrix", "ngsolve.bilinear_form.destroy",
        "ngsolve.matrix.info", "ngsolve.matrix.values", "ngsolve.matrix.vector",
        "ngsolve.matrix.matvec", "ngsolve.matrix.inverse",
        "ngsolve.matrix.projected_create",
        "ngsolve.matrix.reduced_block_create",
        "ngsolve.matrix.diagonal_preconditioner",
        "ngsolve.matrix.term_count",
        "ngsolve.matrix.destroy",
        "ngsolve.solver.create", "ngsolve.solver.info",
        "ngsolve.solver.solve", "ngsolve.solver.destroy",
        "ngsolve.coefficient_function.constant_create",
        "ngsolve.coefficient_function.add",
        "ngsolve.coefficient_function.subtract",
        "ngsolve.coefficient_function.multiply",
        "ngsolve.coefficient_function.scale",
        "ngsolve.coefficient_function.info",
        "ngsolve.coefficient_function.evaluate",
        "ngsolve.coefficient_function.destroy",
        "ngsolve.radia_field.create",
        "ngsolve.radia_field.info",
        "ngsolve.radia_field.prepare_cache",
        "ngsolve.radia_field.clear_cache",
        "ngsolve.radia_field.cache_stats",
        "ngsolve.radia_field.as_voxel_coefficient",
        "ngsolve.grid_function.create",
        "ngsolve.grid_function.info",
        "ngsolve.grid_function.vector",
        "ngsolve.grid_function.set_vector",
        "ngsolve.grid_function.interpolate",
        "ngsolve.grid_function.as_coefficient",
        "ngsolve.grid_function.destroy",
        "ngsolve.grid_function.vector_handle",
        "ngsolve.grid_function.from_fespace",
        "ngsolve.linear_form.create", "ngsolve.linear_form.info",
        "ngsolve.linear_form.create_from_coefficient",
        "ngsolve.linear_form.create_boundary_from_coefficient",
        "ngsolve.linear_form.vector", "ngsolve.linear_form.destroy",
        "ngsolve.vector.info",
        "ngsolve.vector.copy",
        "ngsolve.vector.set_zero",
        "ngsolve.vector.scale",
        "ngsolve.vector.axpy",
        "ngsolve.vector.dot",
        "ngsolve.vector.norm",
        "ngsolve.vector.values",
        "ngsolve.vector.set_values",
        "ngsolve.vector.destroy",
        "simulink.state_space.create",
        "simulink.state_space.info",
        "simulink.state_space.step",
        "simulink.state_space.reset",
        "simulink.state_space.destroy",
        "hcurl.eddy_cln.native_basis",
        "energy_stop.create", "energy_stop.destroy", "energy_stop.info",
        "energy_stop.state0", "energy_stop.forward", "energy_stop.commit",
        "energy_stop.stored_energy", "hybrid_vim.solve", "hybrid_vim.schur",
        "hybrid_vim.skin_impedance", "hybrid_vim.sibc_admittance_tail",
        "hybrid_vim.sibc_termination_impedance",
        "hybrid_vim.sibc_termination_admittance", "cln.lanczos",
        "cln.build_tridiagonal", "cln.impedance", "cln.impedance_sweep",
        "cln.transform_coupling", "cln.transform_port", "cln.aca_compress",
        "evrs.tmethod",
        "hcurl.tet_reduced_gram",
        "hdiv.affine_cell_self_energy_shape_derivative",
        "stream.aca_tsvd",
        "biot_savart.h_segments_complex", "biot_savart.a_segments_complex",
        "biot_savart.a_triangles_complex", "biot_savart.b_triangles_complex",
        "bem.assemble_sldl", "bem.assemble_sldl_p2",
        "hacapk.bem.create", "hacapk.bem.destroy", "hacapk.bem.build",
        "hacapk.bem.matvec", "hacapk.bem.info", "hacapk.peec.create",
        "hacapk.peec.destroy", "hacapk.peec.build", "hacapk.peec.matvec",
        "hacapk.peec.info", "hacapk.charge_gram.create_monopole",
        "hacapk.charge_gram.create_sampled_laplace",
        "hacapk.charge_gram.create_sampled_planar_log",
        "hacapk.charge_gram.create_local_polynomials",
        "hacapk.charge_gram.create_analytic_tet",
        "hacapk.charge_gram.create_analytic_polytope",
        "hacapk.charge_gram.create_high_order_tet",
        "hacapk.charge_gram.create_curved_high_order_tet",
        "hacapk.charge_gram.create_hex",
        "hacapk.charge_gram.create_wedge",
        "hacapk.charge_gram.create_planar_2d",
        "hacapk.charge_gram.create_curved_polytope",
        "hacapk.charge_gram.destroy", "hacapk.charge_gram.build",
        "hacapk.charge_gram.matvec", "hacapk.charge_gram.matvec_transpose",
        "hacapk.charge_gram.matvec_sym", "hacapk.charge_gram.entry",
        "hacapk.charge_gram.hex_volume_self_block_directional_derivative",
        "hacapk.charge_gram.hex_face_self_block_directional_derivative",
        "hacapk.charge_gram.hex_directional_derivative",
        "hacapk.charge_gram.tet_volume_self_block_directional_derivative",
        "hacapk.charge_gram.tet_face_self_block_directional_derivative",
        "hacapk.charge_gram.tet_directional_derivative",
        "hacapk.charge_gram.tet_charge_map_row_directional_rates",
        "hacapk.charge_gram.wedge_volume_self_block_directional_derivative",
        "hacapk.charge_gram.wedge_face_self_block_directional_derivative",
        "hacapk.charge_gram.wedge_directional_derivative",
        "hacapk.charge_gram.directional_derivative_operator",
        "hacapk.charge_gram.directional_derivative_contractions",
        "hacapk.charge_gram_derivative.destroy",
        "hacapk.charge_gram_derivative.info",
        "hacapk.charge_gram_derivative.entry",
        "hacapk.charge_gram_derivative.matvec_sym",
        "hacapk.charge_gram.info",
        "hacapk.charge_gram.hex_state_check", "hacapk.charge_gram.hex_stored_nodes",
        "hacapk.charge_gram.hex_state_breakdown",
        "hacapk.charge_gram.configure_charge_map",
        "hacapk.charge_gram.configure_vector_charge_map",
        "hacapk.charge_gram.configure_mass_matrix",
        "hacapk.charge_gram.configure_geometry_mass_matrix",
        "hacapk.charge_gram.configure_mass_matrix_ngsolve",
        "hacapk.charge_gram.configure_geometry_mass_matrix_ngsolve",
        "hacapk.charge_gram.restore_geometry_mass_matrix",
        "hacapk.charge_gram.operator_info", "hacapk.charge_gram.demag_matrix",
        "hacapk.charge_gram.demag_apply",
        "hacapk.charge_gram.geometry_mass_apply", "hacapk.charge_gram.mass_riesz",
        "hacapk.charge_gram.solve_configured_linear_material",
        "hacapk.charge_gram.solve_configured_linear_material_auto_prec",
        "hacapk.charge_gram.create_field_evaluator",
        "hacapk.charge_gram.create_planar_field_evaluator",
        "hacapk.charge_gram.stats",
        "hdiv.field_evaluator.from_tet", "hdiv.field_evaluator.from_cloud",
        "hdiv.field_evaluator.from_curved_tet", "hdiv.field_evaluator.destroy",
        "hdiv.field_evaluator.field", "hdiv.field_evaluator.candidate_algorithm",
        "hdiv.field_evaluator.last_algorithm", "hdiv.field_evaluator.stats",
        "hdiv.field_evaluator.as_coefficient",
        "hdiv.planar_evaluator.create", "hdiv.planar_evaluator.destroy",
        "hdiv.planar_evaluator.field", "hdiv.planar_evaluator.az",
        "hdiv.planar_evaluator.stats", "hdiv.planar_evaluator.as_coefficient",
        "radia.ObjHexahedron", "radia.ObjTetrahedron", "radia.ObjWedge",
        "radia.ObjPyramid", "radia.ObjThckPgn", "radia.ObjCylMag",
        "radia.ObjRecCur", "radia.ObjArcCur", "radia.ObjRaceTrk",
        "radia.ObjFlmCur", "radia.ObjArcPgnMag", "radia.ObjBckg",
        "radia.ObjCnt", "radia.MatLin",
        "radia.MatSatIsoTab", "radia.MatSatIsoFrm", "radia.MatSatAniso",
        "radia.MatSatLamTab", "radia.MatSatLamFrm", "radia.MatMvsH",
        "radia.MatEnergyHysteresis", "radia.MatPlayHysteresis",
        "radia.MatHysSaveState", "radia.MatHysRestoreState",
        "radia.MatHysCommitState", "radia.MatHysGetNuRev",
        "radia.MatHysIrreversible", "radia.MatHysForwardBatch",
        "radia.MatHysCommitBatch",
        "radia.MatApl", "radia.Solve", "radia.SolveNonl", "radia.GetSolveStats",
        "radia.SolverConfig", "radia.GetSolverConfig",
        "radia.BuildMatrix", "radia.GetInteractMatrix", "radia.GetFaceGeom",
        "radia.PlanarChargeField", "radia.PlanarChargeAz",
        "radia.PlanarMaxwellTorqueCircle", "radia.PlanarMaxwellForceCircle",
        "radia.AverageBInBox", "radia.AverageDemagTensor",
        "equivalence.static_h", "equivalence.harmonic",
        "hlu.set_trunc_tol", "hlu.get_trunc_tol", "hlu.last_timings",
        "hlu.materialize_stats", "hlu.set_parallel", "hlu.get_parallel",
        "hlu.set_par_cutoff", "hlu.max_threads", "hlu.set_accum_cap",
        "hlu.get_accum_cap", "hlu.mixed_breakdown", "hlu.cluster_strategy",
        "radia.GetClusterStrategy",
        "hlu.self_test", "hlu.self_test_rk", "hlu.self_test_addmul_rkrk",
        "hlu.self_test_radia_exact_with_matrix", "hlu.self_test_radia_exact_diag",
        "hlu.self_test_radia_exact", "hlu.self_test_depth3_asymmetric",
        "hlu.self_test_mixed_sibling_via_conversion",
        "hlu.self_test_mixed_sibling_nonuniform", "hlu.self_test_mixed_sibling",
        "hlu.self_test_rk_deep",
        "radia.Fld",
        "radia.FldFrcShpRtg", "radia.FldFrc", "radia.FldLst",
        "radia.FldInt", "radia.ObjCenFld", "radia.FldCmpCrt",
        "radia.FldCmpPrc", "radia.FldLenRndSw", "radia.FldLenTol",
        "radia.ObjGeoVol", "radia.ObjDegFre", "radia.UtiDel", "radia.UtiDelAll",
        "radia.ObjAddToCnt", "radia.ObjCntSize", "radia.ObjCntStuf",
        "radia.ObjDpl", "radia.ObjM", "radia.ObjSetM", "radia.ObjScaleCur",
        "radia.TrfTrsl", "radia.TrfRot", "radia.TrfInv", "radia.TrfCmbL",
        "radia.TrfCmbR", "radia.TrfOrnt", "radia.MatPM", "radia.UtiVer"
    };
    constexpr std::size_t count = sizeof(names) / sizeof(names[0]);
    mxArray* result = mxCreateCellMatrix(1, count);
    for (std::size_t i = 0; i < count; ++i)
        mxSetCell(result, i, mxCreateString(names[i]));
    return result;
}

void ApiInfo(int nlhs, mxArray* plhs[], int nrhs) {
    CheckArity(nrhs, 1, nlhs, 1, "info = radia_mex('api.info')");
    const char* fields[] = {"api_version", "handle_count", "taskmanager_max_threads"};
    plhs[0] = mxCreateStructMatrix(1, 1, 3, fields);
    mxSetField(plhs[0], 0, "api_version", mxCreateDoubleScalar(1.0));
    {
        std::lock_guard<std::mutex> guard(registry_mutex);
        mxSetField(plhs[0], 0, "handle_count",
                   mxCreateDoubleScalar(static_cast<double>(energy_registry.size() +
                                                           bem_registry.size() +
                                                           peec_registry.size() +
                                                           charge_gram_registry.size() +
                                                           field_registry.size() +
                                                           planar_registry.size() +
                                                           coefficient_registry.size() +
                                                           gridfunction_registry.size() +
                                                           vector_registry.size() +
                                                           state_space_registry.size())));
    }
    mxSetField(plhs[0], 0, "taskmanager_max_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetMaxThreads()));
}

void TaskProbe(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if ((nrhs != 1 && nrhs != 2) || nlhs != 1)
        BadArgument("usage: info = radia_mex('taskmanager.probe' [, n])");
    const std::size_t n = nrhs == 2
        ? static_cast<std::size_t>(PositiveInteger(prhs[1], "n"))
        : std::size_t(100000);
    std::vector<double> values(n);
    std::vector<std::atomic<int>> touched(256);
    for (auto& item : touched)
        item.store(0);
    int active_threads = 1;
    {
        ngcore::RegionTaskManager task_manager;
        active_threads = std::max(1, ngcore::TaskManager::GetNumThreads());
        ngcore::ParallelFor(ngcore::IntRange(n), [&](std::size_t i) {
            const int id = ngcore::TaskManager::GetThreadId();
            if (id >= 0 && id < static_cast<int>(touched.size()))
                touched[static_cast<std::size_t>(id)].store(1);
            values[i] = std::sin(0.0001 * static_cast<double>(i)) +
                        std::cos(0.00007 * static_cast<double>(i));
        });
    }
    int used_threads = 0;
    int max_thread_id = 0;
    for (std::size_t i = 0; i < touched.size(); ++i) {
        if (touched[i].load() != 0) {
            ++used_threads;
            max_thread_id = static_cast<int>(i);
        }
    }
    double checksum = 0.0;
    for (double value : values)
        checksum += value;
    const char* fields[] = {"active_threads", "used_threads", "max_thread_id", "checksum"};
    plhs[0] = mxCreateStructMatrix(1, 1, 4, fields);
    mxSetField(plhs[0], 0, "active_threads", mxCreateDoubleScalar(active_threads));
    mxSetField(plhs[0], 0, "used_threads", mxCreateDoubleScalar(used_threads));
    mxSetField(plhs[0], 0, "max_thread_id", mxCreateDoubleScalar(max_thread_id));
    mxSetField(plhs[0], 0, "checksum", mxCreateDoubleScalar(checksum));
}

void SpaceInfo(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if ((nrhs != 3 && nrhs != 4) || nlhs != 1)
        BadArgument("usage: info = radia_mex('ngsolve.space_info', vol_path, order [, nograds])");
    const std::string path = Text(prhs[1], "vol_path");
    const int order = PositiveInteger(prhs[2], "order");
    const bool nograds = nrhs == 4 ? Boolean(prhs[3], "nograds") : true;

    ngcore::RegionTaskManager task_manager;
    auto mesh = std::make_shared<ngcomp::MeshAccess>(path);
    ngcore::Flags hcurl_flags;
    hcurl_flags.SetFlag("order", order);
    hcurl_flags.SetFlag("nograds", nograds);
    auto hcurl = std::make_shared<ngcomp::HCurlHighOrderFESpace>(mesh, hcurl_flags);
    hcurl->Update();
    hcurl->FinalizeUpdate();
    ngcore::Flags hdiv_flags;
    hdiv_flags.SetFlag("order", order);
    auto hdiv = std::make_shared<ngcomp::HDivHighOrderFESpace>(mesh, hdiv_flags);
    hdiv->Update();
    hdiv->FinalizeUpdate();

    const char* fields[] = {"dimension", "vertices", "elements", "order",
        "hcurl_nograds", "hcurl_ndof", "hdiv_ndof", "taskmanager_threads"};
    plhs[0] = mxCreateStructMatrix(1, 1, 8, fields);
    mxSetField(plhs[0], 0, "dimension", mxCreateDoubleScalar(mesh->GetDimension()));
    mxSetField(plhs[0], 0, "vertices", mxCreateDoubleScalar(mesh->GetNV()));
    mxSetField(plhs[0], 0, "elements", mxCreateDoubleScalar(mesh->GetNE()));
    mxSetField(plhs[0], 0, "order", mxCreateDoubleScalar(order));
    mxSetField(plhs[0], 0, "hcurl_nograds", mxCreateLogicalScalar(nograds));
    mxSetField(plhs[0], 0, "hcurl_ndof", mxCreateDoubleScalar(hcurl->GetNDof()));
    mxSetField(plhs[0], 0, "hdiv_ndof", mxCreateDoubleScalar(hdiv->GetNDof()));
    mxSetField(plhs[0], 0, "taskmanager_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetNumThreads()));
}

std::string Lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

std::shared_ptr<ngcomp::FESpace> MakeNGSolveSpace(
    const std::shared_ptr<ngcomp::MeshAccess>& mesh,
    const std::string& space_name, int order, bool nograds,
    bool complex_space = false) {
    ngcore::Flags flags;
    flags.SetFlag("order", order);
    if (complex_space)
        flags.SetFlag("complex", true);
    if (space_name == "h1")
        return std::make_shared<ngcomp::H1HighOrderFESpace>(mesh, flags);
    if (space_name == "vectorh1") {
        auto scalar = std::make_shared<ngcomp::H1HighOrderFESpace>(mesh, flags);
        return std::make_shared<ngcomp::CompoundFESpaceAllSame>(
            scalar, mesh->GetDimension(), flags);
    }
    if (space_name == "hcurl") {
        flags.SetFlag("nograds", nograds);
        return std::make_shared<ngcomp::HCurlHighOrderFESpace>(mesh, flags);
    }
    if (space_name == "hdiv")
        return std::make_shared<ngcomp::HDivHighOrderFESpace>(mesh, flags);
    BadArgument("space must be 'h1', 'vectorh1', 'hcurl', or 'hdiv'");
}

std::shared_ptr<ngfem::BilinearFormIntegrator> MakeNGSolveIntegrator(
    const std::string& space_name, const std::string& form_name, int dimension,
    std::shared_ptr<ngfem::CoefficientFunction> coefficient = nullptr,
    bool allow_complex = false) {
    if (dimension != 2 && dimension != 3)
        BadArgument("ngsolve.matrix_dump supports only 2D and 3D volume meshes");

    if (!coefficient)
        coefficient = ngfem::ConstantCF(1.0);
    if (coefficient->IsComplex() && !allow_complex)
        BadArgument("real NGSolve BilinearForms do not accept a complex CoefficientFunction");
    if (coefficient->Dimension() != 1)
        BadArgument("bilinear-form coefficients must be scalar CoefficientFunctions");
    if (space_name == "h1") {
        if (form_name == "mass") {
            if (dimension == 2)
                return std::make_shared<ngfem::MassIntegrator<2>>(coefficient);
            return std::make_shared<ngfem::MassIntegrator<3>>(coefficient);
        }
        if (form_name == "stiffness") {
            if (dimension == 2)
                return std::make_shared<ngfem::LaplaceIntegrator<2>>(coefficient);
            return std::make_shared<ngfem::LaplaceIntegrator<3>>(coefficient);
        }
    } else if (space_name == "hcurl") {
        if (form_name == "mass") {
            if (dimension == 2)
                return std::make_shared<ngfem::MassEdgeIntegrator<2>>(coefficient);
            return std::make_shared<ngfem::MassEdgeIntegrator<3>>(coefficient);
        }
        if (form_name == "stiffness" || form_name == "curlcurl") {
            if (dimension == 2)
                return std::make_shared<ngfem::CurlCurlEdgeIntegrator<2>>(coefficient);
            return std::make_shared<ngfem::CurlCurlEdgeIntegrator<3>>(coefficient);
        }
    } else if (space_name == "hdiv") {
        if (form_name == "mass") {
            if (dimension == 2)
                return std::make_shared<ngfem::MassHDivIntegrator<2>>(coefficient);
            return std::make_shared<ngfem::MassHDivIntegrator<3>>(coefficient);
        }
        if (form_name == "stiffness" || form_name == "divdiv") {
            if (dimension == 2)
                return std::make_shared<ngfem::DivDivHDivIntegrator<2>>(coefficient);
            return std::make_shared<ngfem::DivDivHDivIntegrator<3>>(coefficient);
        }
    }
    BadArgument("form is incompatible with the selected NGSolve space");
}

std::shared_ptr<ngfem::LinearFormIntegrator> MakeNGSolveLinearIntegrator(
    const std::string& space_name, const std::string& source_name,
    int dimension, Complex source_value, bool allow_complex = false) {
    if (space_name != "h1")
        BadArgument("native persistent LinearForm currently supports the h1 space only");
    if (source_name != "constant")
        BadArgument("source must be 'constant' for the native persistent LinearForm");
    if (source_value.imag() != 0.0 && !allow_complex)
        BadArgument("real NGSolve LinearForms do not accept a complex source value");
    const auto coefficient = allow_complex
        ? ngfem::MakeConstantCoefficientFunction(source_value)
        : ngfem::ConstantCF(source_value.real());
    if (dimension == 2)
        return std::make_shared<ngfem::SourceIntegrator<2>>(
            coefficient);
    if (dimension == 3)
        return std::make_shared<ngfem::SourceIntegrator<3>>(
            coefficient);
    BadArgument("native persistent LinearForm supports only 2D and 3D meshes");
}

std::shared_ptr<ngfem::LinearFormIntegrator>
MakeNGSolveCoefficientLinearIntegrator(
    const std::shared_ptr<ngcomp::FESpace>& fespace,
    const std::shared_ptr<ngfem::CoefficientFunction>& coefficient,
    int dimension, bool allow_complex = false) {
    if (!fespace || !coefficient)
        BadArgument("NGSolve coefficient LinearForm requires a valid space and coefficient");
    if (dimension != 2 && dimension != 3)
        BadArgument("native persistent LinearForm supports only 2D and 3D meshes");
    if (coefficient->IsComplex() && !allow_complex)
        BadArgument("native persistent real LinearForm does not accept a complex CoefficientFunction");

    // Keep the test-function expression inside NGSolve.  This preserves the
    // element-specific Piola maps and local orientation transforms for HCurl
    // and HDiv spaces instead of reconstructing them at the MATLAB boundary.
    const auto test_function = fespace->GetTestFunction();
    const auto integrand = ngfem::InnerProduct(coefficient, *test_function);
    return std::make_shared<ngfem::SymbolicLinearFormIntegrator>(
        integrand, ngfem::VOL, ngfem::VOL);
}

std::shared_ptr<ngfem::LinearFormIntegrator>
MakeNGSolveBoundaryCoefficientLinearIntegrator(
    const std::shared_ptr<ngcomp::FESpace>& fespace,
    const std::shared_ptr<ngfem::CoefficientFunction>& coefficient,
    bool allow_complex = false) {
    if (!fespace || !coefficient)
        BadArgument("NGSolve boundary LinearForm requires a valid space and coefficient");
    if (coefficient->IsComplex() && !allow_complex)
        BadArgument("native persistent real boundary LinearForm does not accept a complex CoefficientFunction");

    // Trace() is deliberately evaluated by NGSolve.  Its meaning is
    // element-family specific (scalar, tangential, or normal trace), so the
    // MEX boundary must not manufacture a surface basis in MATLAB.
    const auto test_function = fespace->GetTestFunction();
    const auto test_trace = (*test_function)->Trace();
    const auto integrand = ngfem::InnerProduct(coefficient, test_trace);
    return std::make_shared<ngfem::SymbolicFacetLinearFormIntegrator>(
        integrand, ngfem::BND);
}

std::shared_ptr<ngfem::BilinearFormIntegrator>
MakeNGSolveBoundaryCoefficientBilinearIntegrator(
    const std::shared_ptr<ngcomp::FESpace>& fespace,
    const std::shared_ptr<ngfem::CoefficientFunction>& coefficient,
    bool allow_complex = false) {
    if (!fespace || !coefficient)
        BadArgument("NGSolve boundary BilinearForm requires a valid space and coefficient");
    if (coefficient->IsComplex() && !allow_complex)
        BadArgument("real NGSolve boundary BilinearForms do not accept a complex CoefficientFunction");
    if (coefficient->Dimension() != 1)
        BadArgument("boundary bilinear-form coefficients must be scalar CoefficientFunctions");

    const auto trial_function = fespace->GetTrialFunction();
    const auto test_function = fespace->GetTestFunction();
    const auto trial_trace = (*trial_function)->Trace();
    const auto test_trace = (*test_function)->Trace();
    const auto integrand = ngfem::InnerProduct(
        ngfem::operator*(coefficient, trial_trace), test_trace);
    return std::make_shared<ngfem::SymbolicFacetBilinearFormIntegrator>(
        integrand, ngfem::BND, false);
}

void NGSolveMatrixDump(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    if ((nrhs != 5 && nrhs != 6) || nlhs != 1)
        BadArgument("usage: out = radia_mex('ngsolve.matrix_dump', "
                    "vol_path, space, order, form [, nograds])");

    const std::string path = Text(prhs[1], "vol_path");
    const std::string space_name = Lowercase(Text(prhs[2], "space"));
    const int order = PositiveInteger(prhs[3], "order");
    const std::string requested_form = Lowercase(Text(prhs[4], "form"));
    const bool nograds = nrhs == 6 ? Boolean(prhs[5], "nograds") : true;
    if (space_name != "h1" && space_name != "hcurl" && space_name != "hdiv")
        BadArgument("space must be 'h1', 'hcurl', or 'hdiv'");

    std::string form_name = requested_form;
    if (space_name == "hcurl" && form_name == "curl_curl")
        form_name = "curlcurl";
    if (space_name == "hdiv" && form_name == "div_div")
        form_name = "divdiv";

    ngcore::RegionTaskManager task_manager;
    auto mesh = std::make_shared<ngcomp::MeshAccess>(path);
    auto fespace = MakeNGSolveSpace(mesh, space_name, order, nograds);
    fespace->Update();
    fespace->FinalizeUpdate();

    ngcore::Flags biform_flags;
    auto biform = std::make_shared<ngcomp::T_BilinearForm<double>>(
        fespace, "radia_matlab_dump", biform_flags);
    biform->AddIntegrator(
        MakeNGSolveIntegrator(space_name, form_name, mesh->GetDimension()));
    // High-order HCurl/HDiv spaces (for example p=6 on a tetrahedron) need
    // more local scratch space during mapped element integration than the
    // old p<=2 default.  Keep this aligned with the native HCurl reduction
    // path so matrix_dump remains a valid independent verification route.
    ngstd::LocalHeap local_heap(1 << 26, "radia_matlab_matrix_dump");
    biform->Assemble(local_heap);

    auto matrix = biform->GetMatrixPtr();
    auto sparse = std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(matrix);
    if (!sparse)
        BadArgument("NGSolve returned a non-real sparse matrix");

    std::vector<double> rows;
    std::vector<double> cols;
    std::vector<double> values;
    rows.reserve(sparse->NZE());
    cols.reserve(sparse->NZE());
    values.reserve(sparse->NZE());
    std::size_t structural_nze = 0;
    for (int row = 0; row < sparse->VHeight(); ++row) {
        auto column_indices = sparse->GetRowIndices(row);
        auto row_values = sparse->GetRowValues(row);
        structural_nze += static_cast<std::size_t>(column_indices.Size());
        for (int entry = 0; entry < column_indices.Size(); ++entry) {
            const double value = row_values[entry];
            if (value == 0.0)
                continue;
            rows.push_back(static_cast<double>(row + 1));
            cols.push_back(static_cast<double>(column_indices[entry] + 1));
            values.push_back(value);
        }
    }

    const char* fields[] = {
        "row", "col", "values", "shape", "space", "form", "dimension",
        "order", "dof_count", "structural_nze", "nonzero_nze", "nograds",
        "symmetric", "taskmanager_threads"};
    plhs[0] = mxCreateStructMatrix(1, 1, 14, fields);
    mxArray* row_array = mxCreateDoubleMatrix(1, rows.size(), mxREAL);
    mxArray* col_array = mxCreateDoubleMatrix(1, cols.size(), mxREAL);
    mxArray* value_array = mxCreateDoubleMatrix(1, values.size(), mxREAL);
    std::copy(rows.begin(), rows.end(), mxGetDoubles(row_array));
    std::copy(cols.begin(), cols.end(), mxGetDoubles(col_array));
    std::copy(values.begin(), values.end(), mxGetDoubles(value_array));
    mxSetField(plhs[0], 0, "row", row_array);
    mxSetField(plhs[0], 0, "col", col_array);
    mxSetField(plhs[0], 0, "values", value_array);
    mxArray* shape = mxCreateDoubleMatrix(1, 2, mxREAL);
    mxGetDoubles(shape)[0] = static_cast<double>(sparse->VHeight());
    mxGetDoubles(shape)[1] = static_cast<double>(sparse->VWidth());
    mxSetField(plhs[0], 0, "shape", shape);
    mxSetField(plhs[0], 0, "space", mxCreateString(space_name.c_str()));
    mxSetField(plhs[0], 0, "form", mxCreateString(form_name.c_str()));
    mxSetField(plhs[0], 0, "dimension",
               mxCreateDoubleScalar(mesh->GetDimension()));
    mxSetField(plhs[0], 0, "order", mxCreateDoubleScalar(order));
    mxSetField(plhs[0], 0, "dof_count",
               mxCreateDoubleScalar(fespace->GetNDof()));
    mxSetField(plhs[0], 0, "structural_nze",
               mxCreateDoubleScalar(structural_nze));
    mxSetField(plhs[0], 0, "nonzero_nze",
               mxCreateDoubleScalar(values.size()));
    mxSetField(plhs[0], 0, "nograds", mxCreateLogicalScalar(nograds));
    mxSetField(plhs[0], 0, "symmetric",
               mxCreateLogicalScalar(matrix->IsSymmetric().IsTrue()));
    mxSetField(plhs[0], 0, "taskmanager_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetNumThreads()));
}

struct NGSolveAssembledSparse {
    std::shared_ptr<ngcomp::T_BilinearForm<double>> form;
    std::shared_ptr<ngla::SparseMatrix<double>> matrix;
};

NGSolveAssembledSparse AssembleNGSolveSparse(
    const std::shared_ptr<ngcomp::FESpace>& fespace,
    const std::string& space_name, const std::string& form_name,
    const std::string& label, ngstd::LocalHeap& local_heap) {
    ngcore::Flags biform_flags;
    auto biform = std::make_shared<ngcomp::T_BilinearForm<double>>(
        fespace, label.c_str(), biform_flags);
    biform->AddIntegrator(MakeNGSolveIntegrator(
        space_name, form_name, fespace->GetMeshAccess()->GetDimension()));
    biform->Assemble(local_heap);
    auto sparse = std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(
        biform->GetMatrixPtr());
    if (!sparse)
        BadArgument("NGSolve returned a non-real sparse matrix");
    return {std::move(biform), std::move(sparse)};
}

NGSolveAssembledSparse AssembleNGSolveEddySparse(
    const std::shared_ptr<ngcomp::FESpace>& fespace,
    const std::string& label, ngstd::LocalHeap& local_heap) {
    ngcore::Flags biform_flags;
    auto biform = std::make_shared<ngcomp::T_BilinearForm<double>>(
        fespace, label.c_str(), biform_flags);
    const int dimension = fespace->GetMeshAccess()->GetDimension();
    biform->AddIntegrator(MakeNGSolveIntegrator("hcurl", "mass", dimension));
    biform->AddIntegrator(
        MakeNGSolveIntegrator("hcurl", "curlcurl", dimension));
    biform->Assemble(local_heap);
    auto sparse = std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(
        biform->GetMatrixPtr());
    if (!sparse)
        BadArgument("NGSolve returned a non-real sparse eddy matrix");
    return {std::move(biform), std::move(sparse)};
}

std::shared_ptr<ngla::BaseVector> NewNGSolveVector(
    const std::shared_ptr<ngla::BaseMatrix>& matrix) {
    auto vector = matrix->CreateColVector();
    return std::move(vector);
}

std::vector<double> NGSolveVectorValues(
    const std::shared_ptr<ngla::BaseVector>& vector, int size) {
    auto values = vector->FVDouble();
    if (values.Size() != size)
        BadArgument("NGSolve vector export size " + std::to_string(values.Size()) +
                    " does not match expected size " + std::to_string(size));
    std::vector<double> result(static_cast<std::size_t>(size));
    for (int i = 0; i < size; ++i)
        result[static_cast<std::size_t>(i)] = values[i];
    return result;
}

void SetNGSolveVectorValues(const std::shared_ptr<ngla::BaseVector>& vector,
                            const std::vector<double>& values) {
    auto target = vector->FVDouble();
    if (target.Size() != static_cast<int>(values.size()))
        BadArgument("NGSolve vector input size " + std::to_string(target.Size()) +
                    " does not match supplied size " +
                    std::to_string(values.size()));
    for (std::size_t i = 0; i < values.size(); ++i)
        target[static_cast<int>(i)] = values[i];
}

std::vector<double> ApplyNGSolveMatrixValues(
    const std::shared_ptr<ngla::BaseMatrix>& matrix,
    const std::vector<double>& values, int size) {
    auto input = NewNGSolveVector(matrix);
    SetNGSolveVectorValues(input, values);
    auto output = NewNGSolveVector(matrix);
    matrix->Mult(*input, *output);
    return NGSolveVectorValues(output, size);
}

std::vector<double> ApplyNGSolveInverseValues(
    const std::shared_ptr<ngla::BaseMatrix>& inverse,
    const std::shared_ptr<ngla::BaseMatrix>& vector_prototype,
    const std::vector<double>& values, int size) {
    auto input = NewNGSolveVector(vector_prototype);
    SetNGSolveVectorValues(input, values);
    auto output = NewNGSolveVector(inverse);
    inverse->Mult(*input, *output);
    return NGSolveVectorValues(output, size);
}

void NGSolveMeshCreate(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "h = radia_mex('ngsolve.mesh.create', vol_path)");
    auto holder = std::make_unique<NGSolveMeshHandle>();
    holder->path = Text(prhs[1], "vol_path");
    {
        ngcore::RegionTaskManager task_manager;
        holder->mesh = std::make_shared<ngcomp::MeshAccess>(holder->path);
    }
    plhs[0] = Uint64Output(RegisterMesh(std::move(holder)));
}

void NGSolveMeshInfo(int nlhs, mxArray* plhs[], int nrhs,
                     const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "info = radia_mex('ngsolve.mesh.info', handle)");
    const auto& holder = Mesh(Handle(prhs[1]));
    const char* fields[] = {"path", "dimension", "vertices", "elements",
                            "has_deformation", "taskmanager_threads"};
    plhs[0] = mxCreateStructMatrix(1, 1, 6, fields);
    mxSetField(plhs[0], 0, "path", TextOutput(holder.path.c_str()));
    mxSetField(plhs[0], 0, "dimension",
               mxCreateDoubleScalar(holder.mesh->GetDimension()));
    mxSetField(plhs[0], 0, "vertices",
               mxCreateDoubleScalar(holder.mesh->GetNV()));
    mxSetField(plhs[0], 0, "elements",
               mxCreateDoubleScalar(holder.mesh->GetNE()));
    mxSetField(plhs[0], 0, "has_deformation",
               mxCreateLogicalScalar(holder.mesh->GetDeformation() != nullptr));
    mxSetField(plhs[0], 0, "taskmanager_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetNumThreads()));
}

template <int D>
double JacobianCondition2(const ngbla::Mat<D, D>& jacobian) {
    double gram[D][D]{};
    for (int row = 0; row < D; ++row)
        for (int col = 0; col < D; ++col)
            for (int component = 0; component < D; ++component)
                gram[row][col] += jacobian(component, row) *
                                  jacobian(component, col);

    // Symmetric Jacobi rotations give the singular-value extrema without
    // introducing a second dense-linear-algebra dependency into the gateway.
    for (int sweep = 0; sweep < 16; ++sweep) {
        int p = 0;
        int q = 1;
        double maximum = std::abs(gram[p][q]);
        for (int row = 0; row < D; ++row)
            for (int col = row + 1; col < D; ++col)
                if (std::abs(gram[row][col]) > maximum) {
                    p = row;
                    q = col;
                    maximum = std::abs(gram[row][col]);
                }
        if (maximum <= 1e-15)
            break;
        const double tau = (gram[q][q] - gram[p][p]) /
                           (2.0 * gram[p][q]);
        const double tangent = (tau >= 0.0 ? 1.0 : -1.0) /
            (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double cosine = 1.0 / std::sqrt(1.0 + tangent * tangent);
        const double sine = tangent * cosine;
        const double app = gram[p][p];
        const double aqq = gram[q][q];
        const double apq = gram[p][q];
        gram[p][p] = app - tangent * apq;
        gram[q][q] = aqq + tangent * apq;
        gram[p][q] = 0.0;
        gram[q][p] = 0.0;
        for (int index = 0; index < D; ++index) {
            if (index == p || index == q)
                continue;
            const double aip = gram[index][p];
            const double aiq = gram[index][q];
            gram[index][p] = cosine * aip - sine * aiq;
            gram[p][index] = gram[index][p];
            gram[index][q] = sine * aip + cosine * aiq;
            gram[q][index] = gram[index][q];
        }
    }
    double minimum = gram[0][0];
    double maximum = gram[0][0];
    for (int index = 1; index < D; ++index) {
        minimum = std::min(minimum, gram[index][index]);
        maximum = std::max(maximum, gram[index][index]);
    }
    if (!(minimum > std::numeric_limits<double>::epsilon()) || maximum < 0.0)
        return std::numeric_limits<double>::infinity();
    return std::sqrt(maximum / minimum);
}

template <int D>
void SampleElementTrafoQuality(ngfem::ElementTransformation& transformation,
                               ngfem::ELEMENT_TYPE element_type,
                               int integration_order, double& minimum_determinant,
                               double& maximum_condition) {
    const auto& rule = ngfem::SelectIntegrationRule(element_type, integration_order);
    minimum_determinant = std::numeric_limits<double>::infinity();
    maximum_condition = 0.0;
    for (int index = 0; index < rule.Size(); ++index) {
        ngfem::MappedIntegrationPoint<D, D> mapped(rule[index], transformation);
        const auto& jacobian = mapped.GetJacobian();
        minimum_determinant = std::min(
            minimum_determinant, static_cast<double>(ngbla::Det(jacobian)));
        maximum_condition = std::max(
            maximum_condition, JacobianCondition2<D>(jacobian));
    }
}

void NGSolveMeshSetDeformation(int nlhs, mxArray* plhs[], int nrhs,
                               const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 0,
               "radia_mex('ngsolve.mesh.set_deformation', mesh, gridfunction)");
    auto& mesh = Mesh(Handle(prhs[1]));
    const auto& deformation = GridFunction(Handle(prhs[2]));
    if (mesh.mesh != deformation.mesh)
        BadArgument("deformation GridFunction must share the target Mesh handle");
    if (deformation.space != "vectorh1")
        BadArgument("mesh deformation requires a real VectorH1 GridFunction");
    if (deformation.gridfunction->GetVectorPtr(0)->IsComplex())
        BadArgument("mesh deformation requires a real VectorH1 GridFunction");
    ngcore::RegionTaskManager task_manager;
    mesh.mesh->SetDeformation(deformation.gridfunction);
}

void NGSolveMeshUnsetDeformation(int nlhs, mxArray* plhs[], int nrhs,
                                 const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 0,
               "radia_mex('ngsolve.mesh.unset_deformation', mesh)");
    auto& mesh = Mesh(Handle(prhs[1]));
    ngcore::RegionTaskManager task_manager;
    mesh.mesh->SetDeformation(nullptr);
}

void NGSolveMeshTrafoQuality(int nlhs, mxArray* plhs[], int nrhs,
                             const mxArray* prhs[]) {
    if ((nrhs != 3 && nrhs != 4) || nlhs != 1)
        BadArgument(
            "usage: q = radia_mex('ngsolve.mesh.trafo_quality', mesh, order "
            "[, reference_determinants])");
    const auto& holder = Mesh(Handle(prhs[1]));
    const int integration_order = PositiveInteger(prhs[2], "integration_order");
    const int element_count = holder.mesh->GetNE();
    if (element_count <= 0)
        BadArgument("Trafo quality requires at least one volume element");
    std::vector<double> reference;
    if (nrhs == 4) {
        reference = RealVector(prhs[3], "reference_determinants");
        if (reference.size() != static_cast<std::size_t>(element_count))
            BadArgument("reference_determinants must match the volume element count");
        for (double value : reference)
            if (value == 0.0 || !std::isfinite(value))
                BadArgument("reference_determinants must be finite and nonzero");
    }

    std::vector<double> raw(static_cast<std::size_t>(element_count));
    std::vector<double> determinants(static_cast<std::size_t>(element_count));
    std::vector<double> conditions(static_cast<std::size_t>(element_count));
    ngcore::RegionTaskManager task_manager;
    ngstd::LocalHeap local_heap(1 << 20, "radia_matlab_mesh_trafo_quality");
    const int dimension = holder.mesh->GetDimension();
    for (int element = 0; element < element_count; ++element) {
        local_heap.CleanUp();
        const ngfem::ElementId element_id(ngfem::VOL, element);
        const auto element_type = holder.mesh->GetElement(element_id).GetType();
        auto& transformation = holder.mesh->GetTrafo(element_id, local_heap);
        if (dimension == 2)
            SampleElementTrafoQuality<2>(transformation, element_type,
                                         integration_order, raw[element],
                                         conditions[element]);
        else if (dimension == 3)
            SampleElementTrafoQuality<3>(transformation, element_type,
                                         integration_order, raw[element],
                                         conditions[element]);
        else
            BadArgument("Trafo quality supports only 2D and 3D volume meshes");
        determinants[element] = reference.empty()
            ? raw[element]
            : raw[element] / reference[static_cast<std::size_t>(element)];
    }
    const double minimum = *std::min_element(determinants.begin(), determinants.end());
    const double maximum = *std::max_element(conditions.begin(), conditions.end());
    const char* fields[] = {"schema", "jacobian_determinants",
                            "raw_jacobian_determinants", "jacobian_conditions",
                            "minimum_jacobian", "maximum_condition",
                            "integration_order", "element_count",
                            "reference_applied", "has_deformation"};
    plhs[0] = mxCreateStructMatrix(1, 1, 10, fields);
    mxSetField(plhs[0], 0, "schema",
               TextOutput("radia.ngsolve.trafo-quality/v1"));
    mxSetField(plhs[0], 0, "jacobian_determinants", RealColumn(determinants));
    mxSetField(plhs[0], 0, "raw_jacobian_determinants", RealColumn(raw));
    mxSetField(plhs[0], 0, "jacobian_conditions", RealColumn(conditions));
    mxSetField(plhs[0], 0, "minimum_jacobian", mxCreateDoubleScalar(minimum));
    mxSetField(plhs[0], 0, "maximum_condition", mxCreateDoubleScalar(maximum));
    mxSetField(plhs[0], 0, "integration_order",
               mxCreateDoubleScalar(integration_order));
    mxSetField(plhs[0], 0, "element_count", mxCreateDoubleScalar(element_count));
    mxSetField(plhs[0], 0, "reference_applied",
               mxCreateLogicalScalar(!reference.empty()));
    mxSetField(plhs[0], 0, "has_deformation",
               mxCreateLogicalScalar(holder.mesh->GetDeformation() != nullptr));
}

void NGSolveFESpaceCreate(int nlhs, mxArray* plhs[], int nrhs,
                          const mxArray* prhs[]) {
    if ((nrhs < 5 || nrhs > 7) || nlhs != 1)
        BadArgument(
            "usage: h = radia_mex('ngsolve.fespace.create', "
            "mesh_handle, space, order [, nograds [, complex]])");
    const auto& mesh_holder = Mesh(Handle(prhs[1]));
    const std::string space_name = Lowercase(Text(prhs[2], "space"));
    const int order = PositiveInteger(prhs[3], "order");
    const bool nograds = nrhs >= 5 ? Boolean(prhs[4], "nograds") : true;
    const bool complex_space = nrhs >= 6 ? Boolean(prhs[5], "complex") : false;
    auto holder = std::make_unique<NGSolveFESpaceHandle>();
    holder->mesh = mesh_holder.mesh;
    holder->space = space_name;
    holder->order = order;
    holder->nograds = nograds;
    holder->is_complex = complex_space;
    {
        ngcore::RegionTaskManager task_manager;
        holder->fespace = MakeNGSolveSpace(
            holder->mesh, space_name, order, nograds, complex_space);
        holder->fespace->Update();
        holder->fespace->FinalizeUpdate();
    }
    plhs[0] = Uint64Output(RegisterFESpace(std::move(holder)));
}

void NGSolveFESpaceInfo(int nlhs, mxArray* plhs[], int nrhs,
                        const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "info = radia_mex('ngsolve.fespace.info', handle)");
    const auto& holder = FESpace(Handle(prhs[1]));
    auto free_dofs = holder.fespace->GetFreeDofs(false);
    std::size_t free_count = 0;
    for (int i = 0; i < holder.fespace->GetNDof(); ++i)
        if (free_dofs->Test(i))
            ++free_count;
    const char* fields[] = {"space", "order", "nograds", "is_complex",
                            "dimension", "vertices", "elements", "dof_count",
                            "free_dof_count", "taskmanager_threads"};
    plhs[0] = mxCreateStructMatrix(1, 1, 10, fields);
    mxSetField(plhs[0], 0, "space", TextOutput(holder.space.c_str()));
    mxSetField(plhs[0], 0, "order", mxCreateDoubleScalar(holder.order));
    mxSetField(plhs[0], 0, "nograds", mxCreateLogicalScalar(holder.nograds));
    mxSetField(plhs[0], 0, "is_complex",
               mxCreateLogicalScalar(holder.is_complex));
    mxSetField(plhs[0], 0, "dimension",
               mxCreateDoubleScalar(holder.mesh->GetDimension()));
    mxSetField(plhs[0], 0, "vertices",
               mxCreateDoubleScalar(holder.mesh->GetNV()));
    mxSetField(plhs[0], 0, "elements",
               mxCreateDoubleScalar(holder.mesh->GetNE()));
    mxSetField(plhs[0], 0, "dof_count",
               mxCreateDoubleScalar(holder.fespace->GetNDof()));
    mxSetField(plhs[0], 0, "free_dof_count",
               mxCreateDoubleScalar(free_count));
    mxSetField(plhs[0], 0, "taskmanager_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetNumThreads()));
}

void NGSolveBilinearFormCreate(int nlhs, mxArray* plhs[], int nrhs,
                               const mxArray* prhs[]) {
    if ((nrhs < 3 || nrhs > 4) || nlhs != 1)
        BadArgument(
            "usage: h = radia_mex('ngsolve.bilinear_form.create', "
            "fespace_handle, form [, label])");
    const auto& space_holder = FESpace(Handle(prhs[1]));
    const std::string form_name = Lowercase(Text(prhs[2], "form"));
    const std::string label = nrhs == 4
        ? Text(prhs[3], "label")
        : "radia_matlab_bilinear_form";
    if (label.empty())
        BadArgument("label must not be empty");

    auto holder = std::make_unique<NGSolveBilinearFormHandle>();
    holder->mesh = space_holder.mesh;
    holder->fespace = space_holder.fespace;
    holder->space = space_holder.space;
    holder->form_name = form_name;
    holder->label = label;
    {
        ngcore::RegionTaskManager task_manager;
        ngcore::Flags flags;
        if (space_holder.is_complex) {
            auto form = std::make_shared<ngcomp::T_BilinearForm<Complex>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveIntegrator(
                holder->space, form_name, holder->mesh->GetDimension(), nullptr, true));
            ngstd::LocalHeap local_heap(1 << 20, "radia_matlab_persistent_complex_form");
            form->Assemble(local_heap);
            holder->matrix = form->GetMatrixPtr();
            holder->form = std::move(form);
        } else {
            auto form = std::make_shared<ngcomp::T_BilinearForm<double>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveIntegrator(
                holder->space, form_name, holder->mesh->GetDimension()));
            ngstd::LocalHeap local_heap(1 << 20, "radia_matlab_persistent_form");
            form->Assemble(local_heap);
            holder->matrix = form->GetMatrixPtr();
            holder->form = std::move(form);
        }
    }
    if (!holder->matrix)
        BadArgument("NGSolve BilinearForm did not produce a matrix");
    plhs[0] = Uint64Output(RegisterBilinearForm(std::move(holder)));
}

void NGSolveBilinearFormCreateFromCoefficient(int nlhs, mxArray* plhs[], int nrhs,
                                              const mxArray* prhs[]) {
    if ((nrhs < 4 || nrhs > 5) || nlhs != 1)
        BadArgument(
            "usage: h = radia_mex('ngsolve.bilinear_form.create_from_coefficient', "
            "fespace_handle, form, coefficient_handle [, label])");
    const auto& space_holder = FESpace(Handle(prhs[1]));
    const std::string form_name = Lowercase(Text(prhs[2], "form"));
    const auto coefficient = Coefficient(Handle(prhs[3]));
    const std::string label = nrhs == 5
        ? Text(prhs[4], "label")
        : "radia_matlab_coefficient_bilinear_form";
    if (label.empty())
        BadArgument("label must not be empty");

    auto holder = std::make_unique<NGSolveBilinearFormHandle>();
    holder->mesh = space_holder.mesh;
    holder->fespace = space_holder.fespace;
    holder->space = space_holder.space;
    holder->form_name = form_name;
    holder->label = label;
    holder->coefficient = coefficient;
    {
        ngcore::RegionTaskManager task_manager;
        ngcore::Flags flags;
        if (!space_holder.is_complex && coefficient->IsComplex())
            BadArgument("real NGSolve BilinearForms do not accept a complex CoefficientFunction");
        if (space_holder.is_complex) {
            auto form = std::make_shared<ngcomp::T_BilinearForm<Complex>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveIntegrator(
                holder->space, form_name, holder->mesh->GetDimension(), coefficient, true));
            ngstd::LocalHeap local_heap(1 << 20,
                                        "radia_matlab_coefficient_complex_bilinear_form");
            form->Assemble(local_heap);
            holder->matrix = form->GetMatrixPtr();
            holder->form = std::move(form);
        } else {
            auto form = std::make_shared<ngcomp::T_BilinearForm<double>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveIntegrator(
                holder->space, form_name, holder->mesh->GetDimension(), coefficient));
            ngstd::LocalHeap local_heap(1 << 20,
                                        "radia_matlab_coefficient_bilinear_form");
            form->Assemble(local_heap);
            holder->matrix = form->GetMatrixPtr();
            holder->form = std::move(form);
        }
    }
    if (!holder->matrix)
        BadArgument("NGSolve BilinearForm did not produce a matrix");
    plhs[0] = Uint64Output(RegisterBilinearForm(std::move(holder)));
}

void NGSolveBilinearFormCreateBoundaryFromCoefficient(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if ((nrhs < 3 || nrhs > 4) || nlhs != 1)
        BadArgument(
            "usage: h = radia_mex('ngsolve.bilinear_form.create_boundary_from_coefficient', "
            "fespace_handle, coefficient_handle [, label])");
    const auto& space_holder = FESpace(Handle(prhs[1]));
    const auto coefficient = Coefficient(Handle(prhs[2]));
    const std::string label = nrhs == 4
        ? Text(prhs[3], "label")
        : "radia_matlab_boundary_coefficient_bilinear_form";
    if (label.empty())
        BadArgument("label must not be empty");

    auto holder = std::make_unique<NGSolveBilinearFormHandle>();
    holder->mesh = space_holder.mesh;
    holder->fespace = space_holder.fespace;
    holder->space = space_holder.space;
    holder->form_name = "boundary_coefficient";
    holder->label = label;
    holder->coefficient = coefficient;
    {
        ngcore::RegionTaskManager task_manager;
        ngcore::Flags flags;
        if (!space_holder.is_complex && coefficient->IsComplex())
            BadArgument("real NGSolve boundary BilinearForms do not accept a complex CoefficientFunction");
        if (space_holder.is_complex) {
            auto form = std::make_shared<ngcomp::T_BilinearForm<Complex>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveBoundaryCoefficientBilinearIntegrator(
                holder->fespace, coefficient, true));
            ngstd::LocalHeap local_heap(1 << 20,
                                        "radia_matlab_boundary_complex_bilinear_form");
            form->Assemble(local_heap);
            holder->matrix = form->GetMatrixPtr();
            holder->form = std::move(form);
        } else {
            auto form = std::make_shared<ngcomp::T_BilinearForm<double>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveBoundaryCoefficientBilinearIntegrator(
                holder->fespace, coefficient));
            ngstd::LocalHeap local_heap(1 << 20,
                                        "radia_matlab_boundary_bilinear_form");
            form->Assemble(local_heap);
            holder->matrix = form->GetMatrixPtr();
            holder->form = std::move(form);
        }
    }
    if (!holder->matrix)
        BadArgument("NGSolve boundary BilinearForm did not produce a matrix");
    plhs[0] = Uint64Output(RegisterBilinearForm(std::move(holder)));
}

void NGSolveBilinearFormInfo(int nlhs, mxArray* plhs[], int nrhs,
                             const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "info = radia_mex('ngsolve.bilinear_form.info', handle)");
    const auto& holder = BilinearForm(Handle(prhs[1]));
    const char* fields[] = {"space", "form", "label", "rows", "cols",
                            "is_complex", "symmetric", "taskmanager_threads"};
    plhs[0] = mxCreateStructMatrix(1, 1, 8, fields);
    mxSetField(plhs[0], 0, "space", TextOutput(holder.space.c_str()));
    mxSetField(plhs[0], 0, "form", TextOutput(holder.form_name.c_str()));
    mxSetField(plhs[0], 0, "label", TextOutput(holder.label.c_str()));
    mxSetField(plhs[0], 0, "rows",
               mxCreateDoubleScalar(holder.matrix->VHeight()));
    mxSetField(plhs[0], 0, "cols",
               mxCreateDoubleScalar(holder.matrix->VWidth()));
    mxSetField(plhs[0], 0, "is_complex",
               mxCreateLogicalScalar(holder.matrix->IsComplex()));
    mxSetField(plhs[0], 0, "symmetric",
               mxCreateLogicalScalar(holder.matrix->IsSymmetric().IsTrue()));
    mxSetField(plhs[0], 0, "taskmanager_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetNumThreads()));
}

std::unique_ptr<NGSolveMatrixHandle> MakeNGSolveMatrixHandle(
    std::shared_ptr<ngla::BaseMatrix> matrix,
    std::shared_ptr<ngcomp::FESpace> fespace,
    const std::string& kind) {
    auto holder = std::make_unique<NGSolveMatrixHandle>();
    holder->matrix = std::move(matrix);
    holder->fespace = std::move(fespace);
    holder->kind = kind;
    return holder;
}

void NGSolveBilinearFormMatrix(int nlhs, mxArray* plhs[], int nrhs,
                               const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "h = radia_mex('ngsolve.bilinear_form.matrix', handle)");
    const auto& form = BilinearForm(Handle(prhs[1]));
    plhs[0] = Uint64Output(RegisterMatrix(MakeNGSolveMatrixHandle(
        form.matrix, form.fespace, "bilinear_form:" + form.form_name)));
}

void NGSolveMatrixInfo(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "info = radia_mex('ngsolve.matrix.info', handle)");
    const auto& holder = Matrix(Handle(prhs[1]));
    const char* fields[] = {"kind", "rows", "cols", "is_complex",
                            "symmetric", "is_sparse", "nonzero_count",
                            "taskmanager_threads"};
    plhs[0] = mxCreateStructMatrix(1, 1, 8, fields);
    const auto real_sparse = std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(
        holder.matrix);
    const auto complex_sparse = std::dynamic_pointer_cast<ngla::SparseMatrix<Complex>>(
        holder.matrix);
    mxSetField(plhs[0], 0, "kind", TextOutput(holder.kind.c_str()));
    mxSetField(plhs[0], 0, "rows",
               mxCreateDoubleScalar(holder.matrix->VHeight()));
    mxSetField(plhs[0], 0, "cols",
               mxCreateDoubleScalar(holder.matrix->VWidth()));
    mxSetField(plhs[0], 0, "is_complex",
               mxCreateLogicalScalar(holder.matrix->IsComplex()));
    mxSetField(plhs[0], 0, "symmetric",
               mxCreateLogicalScalar(holder.matrix->IsSymmetric().IsTrue()));
    mxSetField(plhs[0], 0, "is_sparse",
               mxCreateLogicalScalar(real_sparse != nullptr || complex_sparse != nullptr));
    mxSetField(plhs[0], 0, "nonzero_count", mxCreateDoubleScalar(
        real_sparse ? real_sparse->NZE() :
        (complex_sparse ? complex_sparse->NZE() : 0)));
    mxSetField(plhs[0], 0, "taskmanager_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetNumThreads()));
}

void NGSolveProjectedMatrixCreate(int nlhs, mxArray* plhs[], int nrhs,
                                  const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
        "handle = radia_mex('ngsolve.matrix.projected_create', "
        "parent_handle, projection)");
    const auto& parent = Matrix(Handle(prhs[1]));
    std::size_t rows = 0, cols = 0;
    auto projection = ComplexMatrix(prhs[2], rows, cols, "projection");
    if (rows != static_cast<std::size_t>(parent.matrix->VHeight()) ||
        parent.matrix->VHeight() != parent.matrix->VWidth() || cols == 0)
        BadArgument(
            "projection must have one row per square parent-matrix row");
    auto matrix = std::make_shared<radia::ngsolve_bridge::ProjectedBaseMatrix>(
        parent.matrix, std::move(projection), MatrixDimension(rows, "parent rows"),
        MatrixDimension(cols, "reduced columns"));
    plhs[0] = Uint64Output(RegisterMatrix(MakeNGSolveMatrixHandle(
        std::move(matrix), nullptr, "projected_base_matrix")));
}

void NGSolveReducedBlockMatrixCreate(int nlhs, mxArray* plhs[], int nrhs,
                                     const mxArray* prhs[]) {
    CheckArity(nrhs, 6, nlhs, 1,
        "handle = radia_mex('ngsolve.matrix.reduced_block_create', "
        "dense, matrix_handles, starts, stops, scales)");
    std::size_t dense_rows = 0, dense_cols = 0;
    auto dense = ComplexMatrix(prhs[1], dense_rows, dense_cols, "dense");
    if (dense_rows == 0 || dense_rows != dense_cols)
        BadArgument("dense must be a nonempty square matrix");
    const auto handles = HandleVector(prhs[2], "matrix_handles");
    const auto starts = IntegerVector(prhs[3], "starts");
    const auto stops = IntegerVector(prhs[4], "stops");
    std::size_t scale_rows = 0, scale_cols = 0;
    const auto scales = ComplexMatrix(
        prhs[5], scale_rows, scale_cols, "scales");
    if (handles.size() != starts.size() || starts.size() != stops.size() ||
        stops.size() != scales.size())
        BadArgument(
            "matrix_handles, starts, stops, and scales must have equal lengths");
    std::vector<radia::ngsolve_bridge::ReducedBlockMatrix::Term> terms;
    terms.reserve(handles.size());
    for (std::size_t i = 0; i < handles.size(); ++i) {
        const auto& matrix = Matrix(handles[i]);
        terms.push_back({matrix.matrix, starts[i], stops[i], scales[i]});
    }
    auto matrix = std::make_shared<radia::ngsolve_bridge::ReducedBlockMatrix>(
        std::move(dense), MatrixDimension(dense_rows, "dense rows"),
        std::move(terms));
    plhs[0] = Uint64Output(RegisterMatrix(MakeNGSolveMatrixHandle(
        std::move(matrix), nullptr, "reduced_block_matrix")));
}

void NGSolveMatrixDiagonalPreconditioner(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if ((nrhs != 2 && nrhs != 3) || nlhs != 1)
        BadArgument(
            "handle = radia_mex('ngsolve.matrix.diagonal_preconditioner', "
            "matrix_handle [, relative_floor])");
    const auto& holder = Matrix(Handle(prhs[1]));
    const auto reduced = std::dynamic_pointer_cast<
        radia::ngsolve_bridge::ReducedBlockMatrix>(holder.matrix);
    if (!reduced)
        BadArgument(
            "diagonal_preconditioner requires a reduced block matrix");
    const double relative_floor =
        nrhs == 3 ? Scalar(prhs[2], "relative_floor") : 1.0e-14;
    auto matrix = reduced->DiagonalPreconditioner(relative_floor);
    plhs[0] = Uint64Output(RegisterMatrix(MakeNGSolveMatrixHandle(
        std::move(matrix), nullptr, "complex_diagonal_inverse")));
}

void NGSolveMatrixTermCount(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
        "count = radia_mex('ngsolve.matrix.term_count', matrix_handle)");
    const auto& holder = Matrix(Handle(prhs[1]));
    const auto reduced = std::dynamic_pointer_cast<
        radia::ngsolve_bridge::ReducedBlockMatrix>(holder.matrix);
    if (!reduced)
        BadArgument("term_count requires a reduced block matrix");
    plhs[0] = mxCreateDoubleScalar(reduced->TermCount());
}

mxArray* NGSolveSparseOutput(
    const std::shared_ptr<ngla::SparseMatrix<double>>& sparse) {
    std::vector<double> rows;
    std::vector<double> cols;
    std::vector<double> values;
    std::size_t structural_nze = 0;
    for (int row = 0; row < sparse->VHeight(); ++row) {
        auto column_indices = sparse->GetRowIndices(row);
        auto row_values = sparse->GetRowValues(row);
        structural_nze += static_cast<std::size_t>(column_indices.Size());
        for (int entry = 0; entry < column_indices.Size(); ++entry) {
            const double value = row_values[entry];
            if (value == 0.0)
                continue;
            rows.push_back(static_cast<double>(row + 1));
            cols.push_back(static_cast<double>(column_indices[entry] + 1));
            values.push_back(value);
        }
    }
    const char* fields[] = {"row", "col", "values", "shape",
                            "structural_nze", "nonzero_nze", "symmetric"};
    mxArray* result = mxCreateStructMatrix(1, 1, 7, fields);
    mxSetField(result, 0, "row", RealRow(rows));
    mxSetField(result, 0, "col", RealRow(cols));
    mxSetField(result, 0, "values", RealRow(values));
    mxArray* shape = mxCreateDoubleMatrix(1, 2, mxREAL);
    mxGetDoubles(shape)[0] = static_cast<double>(sparse->VHeight());
    mxGetDoubles(shape)[1] = static_cast<double>(sparse->VWidth());
    mxSetField(result, 0, "shape", shape);
    mxSetField(result, 0, "structural_nze",
               mxCreateDoubleScalar(structural_nze));
    mxSetField(result, 0, "nonzero_nze",
               mxCreateDoubleScalar(values.size()));
    mxSetField(result, 0, "symmetric",
               mxCreateLogicalScalar(sparse->IsSymmetric().IsTrue()));
    return result;
}

struct ScalarSparseCOO {
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> values;
    int size = 0;
};

ScalarSparseCOO ExtractNGSolveScalarSparse(
    const std::shared_ptr<ngla::BaseMatrix>& matrix, const char* caller) {
    if (!matrix)
        BadArgument(std::string(caller) + ": null NGSolve matrix");
    if (matrix->VHeight() != matrix->VWidth())
        BadArgument(std::string(caller) + ": matrix must be square");
    const auto sparse = std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(matrix);
    if (!sparse)
        BadArgument(
            std::string(caller) +
            ": expected an assembled scalar ngla::SparseMatrix<double>; "
            "assemble the NGSolve BilinearForm before configuring HDiv-VIM");

    ScalarSparseCOO result;
    result.size = matrix->VHeight();
    result.rows.reserve(sparse->NZE());
    result.cols.reserve(sparse->NZE());
    result.values.reserve(sparse->NZE());
    for (int row = 0; row < result.size; ++row) {
        const auto indices = sparse->GetRowIndices(row);
        const auto values = sparse->GetRowValues(row);
        if (indices.Size() != values.Size())
            BadArgument(std::string(caller) + ": inconsistent sparse row");
        for (int local = 0; local < indices.Size(); ++local) {
            result.rows.push_back(row);
            result.cols.push_back(indices[local]);
            result.values.push_back(values[local]);
        }
    }
    return result;
}

mxArray* NGSolveSparseOutput(
    const std::shared_ptr<ngla::SparseMatrix<Complex>>& sparse) {
    std::vector<double> rows;
    std::vector<double> cols;
    std::vector<Complex> values;
    std::size_t structural_nze = 0;
    for (int row = 0; row < sparse->VHeight(); ++row) {
        auto column_indices = sparse->GetRowIndices(row);
        auto row_values = sparse->GetRowValues(row);
        structural_nze += static_cast<std::size_t>(column_indices.Size());
        for (int entry = 0; entry < column_indices.Size(); ++entry) {
            const Complex value = row_values[entry];
            if (value == Complex(0.0, 0.0))
                continue;
            rows.push_back(static_cast<double>(row + 1));
            cols.push_back(static_cast<double>(column_indices[entry] + 1));
            values.push_back(value);
        }
    }
    const char* fields[] = {"row", "col", "values", "shape",
                            "structural_nze", "nonzero_nze", "symmetric"};
    mxArray* result = mxCreateStructMatrix(1, 1, 7, fields);
    mxSetField(result, 0, "row", RealRow(rows));
    mxSetField(result, 0, "col", RealRow(cols));
    mxSetField(result, 0, "values", ComplexRow(values));
    mxArray* shape = mxCreateDoubleMatrix(1, 2, mxREAL);
    mxGetDoubles(shape)[0] = static_cast<double>(sparse->VHeight());
    mxGetDoubles(shape)[1] = static_cast<double>(sparse->VWidth());
    mxSetField(result, 0, "shape", shape);
    mxSetField(result, 0, "structural_nze",
               mxCreateDoubleScalar(structural_nze));
    mxSetField(result, 0, "nonzero_nze",
               mxCreateDoubleScalar(values.size()));
    mxSetField(result, 0, "symmetric",
               mxCreateLogicalScalar(sparse->IsSymmetric().IsTrue()));
    return result;
}

void NGSolveMatrixValues(int nlhs, mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "out = radia_mex('ngsolve.matrix.values', handle)");
    const auto& holder = Matrix(Handle(prhs[1]));
    const auto real_sparse = std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(
        holder.matrix);
    if (real_sparse) {
        plhs[0] = NGSolveSparseOutput(real_sparse);
        return;
    }
    const auto complex_sparse = std::dynamic_pointer_cast<ngla::SparseMatrix<Complex>>(
        holder.matrix);
    if (complex_sparse) {
        plhs[0] = NGSolveSparseOutput(complex_sparse);
        return;
    }
    BadArgument("ngsolve.matrix.values requires an assembled sparse matrix");
}

void NGSolveMatrixVector(int nlhs, mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "h = radia_mex('ngsolve.matrix.vector', handle)");
    const auto& matrix = Matrix(Handle(prhs[1]));
    if (matrix.matrix->VHeight() != matrix.matrix->VWidth())
        BadArgument("ngsolve.matrix.vector currently requires a square matrix");
    auto holder = std::make_unique<NGSolveVectorHandle>();
    holder->vector = matrix.matrix->CreateColVector();
    holder->parent_matrix = matrix.matrix;
    holder->is_view = false;
    plhs[0] = Uint64Output(RegisterVector(std::move(holder)));
}

void NGSolveMatrixMatVec(int nlhs, mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]) {
    if ((nrhs != 3 && nrhs != 4) || nlhs != 1)
        BadArgument(
            "out = radia_mex('ngsolve.matrix.matvec', matrix, vector "
            "[, transpose])");
    const auto& matrix = Matrix(Handle(prhs[1]));
    const auto& input = Vector(Handle(prhs[2]));
    const bool transpose = nrhs == 4 ? Boolean(prhs[3], "transpose") : false;
    const int expected = transpose ? matrix.matrix->VHeight()
                                   : matrix.matrix->VWidth();
    if (input.vector->Size() != expected)
        BadArgument("vector size does not match the NGSolve matrix dimension");
    if (input.vector->IsComplex() != matrix.matrix->IsComplex())
        BadArgument("NGSolve matrix and vector scalar types must match for matvec");
    auto output = transpose ? matrix.matrix->CreateRowVector()
                            : matrix.matrix->CreateColVector();
    if (transpose)
        matrix.matrix->MultTrans(*input.vector, *output);
    else
        matrix.matrix->Mult(*input.vector, *output);
    auto holder = std::make_unique<NGSolveVectorHandle>();
    holder->vector = std::move(output);
    holder->parent_matrix = matrix.matrix;
    holder->is_view = false;
    plhs[0] = Uint64Output(RegisterVector(std::move(holder)));
}

std::shared_ptr<ngla::KrylovSpaceSolver> MakeNGSolveSolver(
    const std::string& method,
    const std::shared_ptr<ngla::BaseMatrix>& matrix,
    const std::shared_ptr<ngla::BaseMatrix>& preconditioner) {
    if (matrix->IsComplex()) {
        if (method != "gmres")
            BadArgument("complex native NGSolve Solver currently supports 'gmres'");
        if (preconditioner)
            return std::make_shared<ngla::GMRESSolver<ngla::ComplexConjugate>>(
                matrix, preconditioner);
        return std::make_shared<ngla::GMRESSolver<ngla::ComplexConjugate>>(matrix);
    }
    if (method == "cg") {
        if (preconditioner)
            return std::make_shared<ngla::CGSolver<double>>(
                matrix, preconditioner);
        return std::make_shared<ngla::CGSolver<double>>(matrix);
    }
    if (method == "gmres") {
        if (preconditioner)
            return std::make_shared<ngla::GMRESSolver<double>>(
                matrix, preconditioner);
        return std::make_shared<ngla::GMRESSolver<double>>(matrix);
    }
    if (method == "bicgstab") {
        if (preconditioner)
            return std::make_shared<ngla::BiCGStabSolver<double>>(
                matrix, preconditioner);
        return std::make_shared<ngla::BiCGStabSolver<double>>(matrix);
    }
    BadArgument("method must be 'cg', 'gmres', or 'bicgstab'");
}

void NGSolveSolverCreate(int nlhs, mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]) {
    if ((nrhs < 5 || nrhs > 6) || nlhs != 1)
        BadArgument(
            "usage: h = radia_mex('ngsolve.solver.create', matrix_handle, "
            "method, tolerance, max_steps [, preconditioner_handle])");
    const auto& matrix_holder = Matrix(Handle(prhs[1]));
    const std::string method = Lowercase(Text(prhs[2], "method"));
    const double tolerance = Scalar(prhs[3], "tolerance");
    const int max_steps = PositiveInteger(prhs[4], "max_steps");
    if (!(tolerance > 0.0) || !std::isfinite(tolerance))
        BadArgument("tolerance must be finite and positive");
    if (matrix_holder.matrix->VHeight() != matrix_holder.matrix->VWidth())
        BadArgument("ngsolve.solver.create requires a square matrix");

    std::shared_ptr<ngla::BaseMatrix> preconditioner;
    if (nrhs == 6) {
        const auto& preconditioner_holder = Matrix(Handle(prhs[5]));
        if (preconditioner_holder.matrix->VHeight() != matrix_holder.matrix->VHeight() ||
            preconditioner_holder.matrix->VWidth() != matrix_holder.matrix->VWidth())
            BadArgument("preconditioner dimensions must match the solver matrix");
        if (preconditioner_holder.matrix->IsComplex() != matrix_holder.matrix->IsComplex())
            BadArgument("preconditioner scalar type must match the solver matrix");
        preconditioner = preconditioner_holder.matrix;
    }

    auto holder = std::make_unique<NGSolveSolverHandle>();
    holder->matrix = matrix_holder.matrix;
    holder->preconditioner = preconditioner;
    holder->method = method;
    holder->tolerance = tolerance;
    holder->max_steps = max_steps;
    holder->solver = MakeNGSolveSolver(method, holder->matrix, preconditioner);
    holder->solver->SetPrecision(tolerance);
    holder->solver->SetMaxSteps(max_steps);
    plhs[0] = Uint64Output(RegisterSolver(std::move(holder)));
}

void NGSolveSolverInfo(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "info = radia_mex('ngsolve.solver.info', handle)");
    const auto& holder = Solver(Handle(prhs[1]));
    const char* fields[] = {"method", "rows", "cols", "is_complex",
                            "tolerance", "max_steps", "steps",
                            "has_preconditioner", "taskmanager_threads"};
    plhs[0] = mxCreateStructMatrix(1, 1, 9, fields);
    mxSetField(plhs[0], 0, "method", TextOutput(holder.method.c_str()));
    mxSetField(plhs[0], 0, "rows",
               mxCreateDoubleScalar(holder.matrix->VHeight()));
    mxSetField(plhs[0], 0, "cols",
               mxCreateDoubleScalar(holder.matrix->VWidth()));
    mxSetField(plhs[0], 0, "is_complex",
               mxCreateLogicalScalar(holder.matrix->IsComplex()));
    mxSetField(plhs[0], 0, "tolerance",
               mxCreateDoubleScalar(holder.tolerance));
    mxSetField(plhs[0], 0, "max_steps",
               mxCreateDoubleScalar(holder.max_steps));
    mxSetField(plhs[0], 0, "steps",
               mxCreateDoubleScalar(holder.solver->GetSteps()));
    mxSetField(plhs[0], 0, "has_preconditioner",
               mxCreateLogicalScalar(holder.preconditioner != nullptr));
    mxSetField(plhs[0], 0, "taskmanager_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetNumThreads()));
}

void NGSolveSolverSolve(int nlhs, mxArray* plhs[], int nrhs,
                        const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
               "x = radia_mex('ngsolve.solver.solve', solver_handle, rhs_handle)");
    const auto& holder = Solver(Handle(prhs[1]));
    const auto& rhs = Vector(Handle(prhs[2]));
    if (rhs.vector->Size() != holder.matrix->VHeight())
        BadArgument("right-hand-side size does not match the solver matrix");
    if (rhs.vector->IsComplex() != holder.matrix->IsComplex())
        BadArgument("right-hand-side scalar type must match the solver matrix");
    ngcore::RegionTaskManager task_manager;
    auto solution = holder.solver->CreateColVector();
    holder.solver->Mult(*rhs.vector, *solution);
    auto result = std::make_unique<NGSolveVectorHandle>();
    result->vector = std::move(solution);
    result->parent_solver = holder.solver;
    result->is_view = false;
    plhs[0] = Uint64Output(RegisterVector(std::move(result)));
}

void NGSolveMatrixInverse(int nlhs, mxArray* plhs[], int nrhs,
                          const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "h = radia_mex('ngsolve.matrix.inverse', handle)");
    const auto& matrix = Matrix(Handle(prhs[1]));
    if (matrix.matrix->VHeight() != matrix.matrix->VWidth())
        BadArgument("ngsolve.matrix.inverse requires a square matrix");
    if (!matrix.fespace)
        BadArgument("NGSolve matrix has no FESpace for free-DoF inversion");
    ngcore::RegionTaskManager task_manager;
    auto free_dofs = matrix.fespace->GetFreeDofs(false);
    auto inverse = matrix.matrix->InverseMatrix(free_dofs);
    plhs[0] = Uint64Output(RegisterMatrix(MakeNGSolveMatrixHandle(
        std::move(inverse), matrix.fespace, "inverse(" + matrix.kind + ")")));
}

void HCurlEddyCLNNativeBasis(int nlhs, mxArray* plhs[], int nrhs,
                             const mxArray* prhs[]) {
    if ((nrhs < 5 || nrhs > 7) || nlhs != 1)
        BadArgument(
            "usage: out = radia_mex('hcurl.eddy_cln.native_basis', "
            "vol_path, order, ports, steps [, nograds [, rtol]])");

    const std::string path = Text(prhs[1], "vol_path");
    const int order = PositiveInteger(prhs[2], "order");
    std::size_t port_rows = 0;
    std::size_t port_cols = 0;
    const auto ports = RealMatrix(prhs[3], port_rows, port_cols, "ports");
    const int steps = PositiveInteger(prhs[4], "steps");
    const bool nograds = nrhs >= 6 ? Boolean(prhs[5], "nograds") : true;
    const double rtol = nrhs >= 7 ? Scalar(prhs[6], "rtol") : 1.0e-12;
    if (port_rows == 0 || port_cols == 0)
        BadArgument("ports must be a non-empty ndof-by-nports matrix");
    if (!std::isfinite(rtol) || rtol <= 0.0)
        BadArgument("rtol must be positive and finite");
    if (!std::all_of(ports.begin(), ports.end(),
                     [](double value) { return std::isfinite(value); }))
        BadArgument("ports must contain only finite values");

    ngcore::RegionTaskManager task_manager;
    auto mesh = std::make_shared<ngcomp::MeshAccess>(path);
    auto fespace = MakeNGSolveSpace(mesh, "hcurl", order, nograds);
    fespace->Update();
    fespace->FinalizeUpdate();
    const int ndof = MatrixDimension(fespace->GetNDof(), "HCurl ndof");
    if (port_rows != static_cast<std::size_t>(ndof))
        BadArgument("ports row count must equal the HCurl DoF count");

    ngstd::LocalHeap local_heap(1 << 26, "radia_matlab_hcurl_eddy_cln");
    auto mass_assembly = AssembleNGSolveSparse(
        fespace, "hcurl", "mass", "radia_hcurl_eddy_mass", local_heap);
    auto mass = mass_assembly.matrix;
    auto system_assembly = AssembleNGSolveEddySparse(
        fespace, "radia_hcurl_eddy_system", local_heap);
    auto system = system_assembly.matrix;
    auto free_dofs = fespace->GetFreeDofs(false);
    std::size_t active_count = 0;
    for (int i = 0; i < ndof; ++i)
        if (free_dofs->Test(i))
            ++active_count;
    if (active_count == 0)
        BadArgument("HCurl space has no free DoFs");

    // The unit-shifted eddy operator is positive definite: mass + curl-curl.
    // This keeps the HCurl gradient kernel controlled while leaving the
    // physical frequency/material scaling as an explicit future contract.
    std::shared_ptr<ngla::BaseMatrix> inverse =
        system->InverseMatrix(free_dofs);
    std::vector<std::vector<double>> current;
    current.reserve(port_cols);
    for (std::size_t port = 0; port < port_cols; ++port) {
        std::vector<double> rhs(static_cast<std::size_t>(ndof), 0.0);
        for (int row = 0; row < ndof; ++row)
            if (free_dofs->Test(row))
                rhs[static_cast<std::size_t>(row)] =
                    ports[static_cast<std::size_t>(row) * port_cols + port];
        current.push_back(std::move(rhs));
    }

    std::vector<std::vector<double>> modes;
    std::vector<std::vector<double>> metric_modes;
    modes.reserve(static_cast<std::size_t>(steps) * port_cols);
    metric_modes.reserve(static_cast<std::size_t>(steps) * port_cols);
    for (int step = 0; step < steps; ++step) {
        std::vector<std::vector<double>> next_current;
        next_current.reserve(port_cols);
        for (std::size_t port = 0; port < port_cols; ++port) {
            auto solved = ApplyNGSolveInverseValues(
                inverse, system, current[port], ndof);
            auto metric_solved = ApplyNGSolveMatrixValues(mass, solved, ndof);
            const double initial_norm2 = std::inner_product(
                solved.begin(), solved.end(), metric_solved.begin(), 0.0);
            if (!std::isfinite(initial_norm2) || initial_norm2 <= 0.0) {
                next_current.emplace_back(static_cast<std::size_t>(ndof), 0.0);
                continue;
            }
            auto candidate = std::move(solved);
            auto metric_candidate = std::move(metric_solved);
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t mode = 0; mode < modes.size(); ++mode) {
                    const double coefficient = std::inner_product(
                        modes[mode].begin(), modes[mode].end(),
                        metric_candidate.begin(), 0.0);
                    for (int row = 0; row < ndof; ++row) {
                        candidate[static_cast<std::size_t>(row)] -=
                            coefficient * modes[mode][static_cast<std::size_t>(row)];
                        metric_candidate[static_cast<std::size_t>(row)] -=
                            coefficient * metric_modes[mode][static_cast<std::size_t>(row)];
                    }
                }
            }
            const double norm2 = std::inner_product(
                candidate.begin(), candidate.end(), metric_candidate.begin(), 0.0);
            const double threshold =
                rtol * rtol * std::max(initial_norm2, 1.0e-300);
            if (!std::isfinite(norm2) || norm2 <= threshold) {
                next_current.emplace_back(static_cast<std::size_t>(ndof), 0.0);
                continue;
            }
            const double norm = std::sqrt(norm2);
            for (int row = 0; row < ndof; ++row) {
                candidate[static_cast<std::size_t>(row)] /= norm;
                metric_candidate[static_cast<std::size_t>(row)] /= norm;
            }
            modes.push_back(std::move(candidate));
            metric_modes.push_back(std::move(metric_candidate));
            next_current.push_back(
                ApplyNGSolveMatrixValues(mass, modes.back(), ndof));
        }
        current = std::move(next_current);
        if (current.empty())
            break;
    }

    const std::size_t rank = modes.size();
    std::vector<double> vectors(static_cast<std::size_t>(ndof) * rank, 0.0);
    for (std::size_t mode = 0; mode < rank; ++mode) {
        for (int row = 0; row < ndof; ++row)
            vectors[static_cast<std::size_t>(row) * rank + mode] =
                modes[mode][static_cast<std::size_t>(row)];
    }

    // Export the reduced FE operators together with the response basis.  The
    // MATLAB side can therefore form a CLN model without rebuilding the
    // high-order NGSolve space or copying the parent sparse matrices through
    // Python.  The system operator is M + K, so subtracting the separately
    // applied mass matrix gives the curl-curl block K without another global
    // assembly.
    std::vector<double> mass_gram(rank * rank, 0.0);
    std::vector<double> curlcurl_gram(rank * rank, 0.0);
    for (std::size_t column = 0; column < rank; ++column) {
        const auto mass_applied = ApplyNGSolveMatrixValues(
            mass, modes[column], ndof);
        const auto system_applied = ApplyNGSolveMatrixValues(
            system, modes[column], ndof);
        for (std::size_t row = 0; row < rank; ++row) {
            const auto row_index = row * rank + column;
            mass_gram[row_index] = std::inner_product(
                modes[row].begin(), modes[row].end(),
                mass_applied.begin(), 0.0);
            curlcurl_gram[row_index] = std::inner_product(
                modes[row].begin(), modes[row].end(),
                system_applied.begin(), 0.0) - mass_gram[row_index];
        }
    }
    std::vector<double> port_rhs(rank * port_cols, 0.0);
    for (std::size_t mode = 0; mode < rank; ++mode)
        for (std::size_t port = 0; port < port_cols; ++port)
            for (int row = 0; row < ndof; ++row)
                port_rhs[mode * port_cols + port] +=
                    modes[mode][static_cast<std::size_t>(row)] *
                    ports[static_cast<std::size_t>(row) * port_cols + port];

    double orthogonality_error = 0.0;
    for (std::size_t i = 0; i < rank; ++i)
        for (std::size_t j = 0; j < rank; ++j) {
            const double expected = i == j ? 1.0 : 0.0;
            orthogonality_error = std::max(
                orthogonality_error,
                std::fabs(std::inner_product(
                              modes[i].begin(), modes[i].end(),
                              metric_modes[j].begin(), 0.0) - expected));
        }

    const char* fields[] = {
        "vectors", "free_dofs", "dof_count",
        "active_count", "rank", "order", "steps", "nograds", "rtol",
        "orthogonality_error", "taskmanager_threads", "operator",
        "mass_gram", "curlcurl_gram", "port_rhs", "projection"};
    plhs[0] = mxCreateStructMatrix(1, 1, 16, fields);
    mxSetField(plhs[0], 0, "vectors",
               RealMatrixOutput(vectors, static_cast<std::size_t>(ndof), rank));
    mxArray* free_array = mxCreateLogicalMatrix(1, static_cast<std::size_t>(ndof));
    auto* free_values = mxGetLogicals(free_array);
    for (int row = 0; row < ndof; ++row)
        free_values[row] = free_dofs->Test(row);
    mxSetField(plhs[0], 0, "free_dofs", free_array);
    mxSetField(plhs[0], 0, "dof_count", mxCreateDoubleScalar(ndof));
    mxSetField(plhs[0], 0, "active_count",
               mxCreateDoubleScalar(static_cast<double>(active_count)));
    mxSetField(plhs[0], 0, "rank", mxCreateDoubleScalar(static_cast<double>(rank)));
    mxSetField(plhs[0], 0, "order", mxCreateDoubleScalar(order));
    mxSetField(plhs[0], 0, "steps", mxCreateDoubleScalar(steps));
    mxSetField(plhs[0], 0, "nograds", mxCreateLogicalScalar(nograds));
    mxSetField(plhs[0], 0, "rtol", mxCreateDoubleScalar(rtol));
    mxSetField(plhs[0], 0, "orthogonality_error",
               mxCreateDoubleScalar(orthogonality_error));
    mxSetField(plhs[0], 0, "taskmanager_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetNumThreads()));
    mxSetField(plhs[0], 0, "operator", mxCreateString("mass+curlcurl"));
    mxSetField(plhs[0], 0, "mass_gram",
               RealMatrixOutput(mass_gram, rank, rank));
    mxSetField(plhs[0], 0, "curlcurl_gram",
               RealMatrixOutput(curlcurl_gram, rank, rank));
    mxSetField(plhs[0], 0, "port_rhs",
               RealMatrixOutput(port_rhs, rank, port_cols));
    mxSetField(plhs[0], 0, "projection",
               mxCreateString("mass_gram=V'*M*V; curlcurl_gram=V'*K*V; port_rhs=V'*ports"));
}

std::shared_ptr<ngfem::CoefficientFunction> MakeNGSolveConstantCoefficient(
    const mxArray* value, std::size_t& rows, std::size_t& cols) {
    const auto values = ComplexMatrix(value, rows, cols, "values");
    if (values.empty())
        BadArgument("values must be non-empty");

    const bool is_complex = mxIsComplex(value);
    auto make_scalar = [is_complex](Complex scalar) {
        if (is_complex)
            return ngfem::MakeConstantCoefficientFunction(scalar);
        return ngfem::ConstantCF(scalar.real());
    };
    if (values.size() == 1)
        return make_scalar(values.front());

    const int component_count = MatrixDimension(values.size(), "values");
    ngcore::Array<std::shared_ptr<ngfem::CoefficientFunction>> components(
        component_count);
    for (int i = 0; i < component_count; ++i)
        components[i] = make_scalar(values[static_cast<std::size_t>(i)]);
    auto coefficient = ngfem::MakeVectorialCoefficientFunction(std::move(components));
    if (rows == 1 || cols == 1) {
        coefficient->SetDimensions(ngcore::Array<int>{component_count});
    } else {
        coefficient->SetDimensions(ngcore::Array<int>{
            MatrixDimension(rows, "values rows"),
            MatrixDimension(cols, "values columns")});
    }
    return coefficient;
}

std::vector<double> NGSolveCoefficientDimensions(
    const ngfem::CoefficientFunction& coefficient) {
    std::vector<double> dimensions;
    const auto shape = coefficient.Dimensions();
    dimensions.reserve(shape.Size());
    for (int value : shape)
        dimensions.push_back(static_cast<double>(value));
    return dimensions;
}

void NGSolveCoefficientConstantCreate(int nlhs, mxArray* plhs[], int nrhs,
                                      const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "h = radia_mex('ngsolve.coefficient_function.constant_create', values)");
    std::size_t rows = 0;
    std::size_t cols = 0;
    auto coefficient = MakeNGSolveConstantCoefficient(
        prhs[1], rows, cols);
    plhs[0] = Uint64Output(RegisterCoefficient(std::move(coefficient)));
}

void NGSolveCoefficientBinary(const std::string& command, int nlhs,
                              mxArray* plhs[], int nrhs,
                              const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
               "h = radia_mex('ngsolve.coefficient_function.<op>', left, right)");
    const auto left = Coefficient(Handle(prhs[1]));
    const auto right = Coefficient(Handle(prhs[2]));
    std::shared_ptr<ngfem::CoefficientFunction> result;
    if (command == "ngsolve.coefficient_function.add")
        result = ngfem::operator+(left, right);
    else if (command == "ngsolve.coefficient_function.subtract")
        result = ngfem::operator-(left, right);
    else if (command == "ngsolve.coefficient_function.multiply")
        result = ngfem::operator*(left, right);
    else
        BadArgument("unknown NGSolve CoefficientFunction binary operation");
    plhs[0] = Uint64Output(RegisterCoefficient(std::move(result)));
}

void NGSolveCoefficientScale(int nlhs, mxArray* plhs[], int nrhs,
                             const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
               "h = radia_mex('ngsolve.coefficient_function.scale', scalar, coefficient)");
    const Complex scalar = ComplexScalar(prhs[1], "scalar");
    const auto coefficient = Coefficient(Handle(prhs[2]));
    plhs[0] = Uint64Output(RegisterCoefficient(
        ngfem::operator*(scalar, coefficient)));
}

void NGSolveCoefficientInfo(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "info = radia_mex('ngsolve.coefficient_function.info', handle)");
    const auto coefficient = Coefficient(Handle(prhs[1]));
    const auto hdiv = std::dynamic_pointer_cast<
        radia::ngsolve_bridge::HDivFieldCoefficient>(coefficient);
    const auto planar = std::dynamic_pointer_cast<
        radia::ngsolve_bridge::PlanarHDivFieldCoefficient>(coefficient);
    const char* fields[] = {"dimension", "dimensions", "is_complex",
                            "description", "equivalence_key", "kind",
                            "algorithm", "source_angle", "target_angle"};
    plhs[0] = mxCreateStructMatrix(1, 1, 9, fields);
    mxSetField(plhs[0], 0, "dimension",
               mxCreateDoubleScalar(static_cast<double>(coefficient->Dimension())));
    mxSetField(plhs[0], 0, "dimensions",
               RealRow(NGSolveCoefficientDimensions(*coefficient)));
    mxSetField(plhs[0], 0, "is_complex",
               mxCreateLogicalScalar(coefficient->IsComplex()));
    mxSetField(plhs[0], 0, "description",
               TextOutput(coefficient->GetDescription().c_str()));
    mxSetField(plhs[0], 0, "equivalence_key",
               TextOutput(coefficient->EquivalenceKey().c_str()));
    mxSetField(plhs[0], 0, "kind", TextOutput(
        hdiv ? "hdiv_field" : (planar ? "planar_hdiv_field" : "generic")));
    mxSetField(plhs[0], 0, "algorithm",
               TextOutput(hdiv ? hdiv->AlgorithmName() : ""));
    mxSetField(plhs[0], 0, "source_angle", mxCreateDoubleScalar(
        planar ? planar->SourceAngle() : std::numeric_limits<double>::quiet_NaN()));
    mxSetField(plhs[0], 0, "target_angle", mxCreateDoubleScalar(
        planar ? planar->TargetAngle() : std::numeric_limits<double>::quiet_NaN()));
}

void NGSolveCoefficientEvaluate(int nlhs, mxArray* plhs[], int nrhs,
                                const mxArray* prhs[]) {
    if (nrhs != 4 || nlhs != 1)
        BadArgument(
            "values = radia_mex('ngsolve.coefficient_function.evaluate', "
            "vol_path, handle, points)");
    const std::string path = Text(prhs[1], "vol_path");
    const auto coefficient = Coefficient(Handle(prhs[2]));
    std::size_t point_count = 0;
    std::size_t point_dimension = 0;
    const auto points = RealMatrix(prhs[3], point_count, point_dimension, "points");
    if (point_count == 0)
        BadArgument("points must be non-empty");
    if (point_dimension != 2 && point_dimension != 3)
        BadArgument("points must have two or three columns");

    ngcore::RegionTaskManager task_manager;
    auto mesh = std::make_shared<ngcomp::MeshAccess>(path);
    if (mesh->GetDimension() != static_cast<int>(point_dimension))
        BadArgument("points dimension must match the volume mesh dimension");
    const std::size_t value_dimension = coefficient->Dimension();
    if (value_dimension == 0)
        BadArgument("CoefficientFunction has an invalid zero dimension");
    std::vector<Complex> values(point_count * value_dimension, Complex(0.0, 0.0));
    ngstd::LocalHeap local_heap(1 << 20, "radia_matlab_coefficient_evaluate");

    for (std::size_t point = 0; point < point_count; ++point) {
        local_heap.CleanUp();
        double* point_data = const_cast<double*>(
            points.data() + point * point_dimension);
        ngbla::FlatVector<double> physical_point(
            point_dimension, point_data);
        ngfem::IntegrationPoint integration_point;
        const ngfem::ElementId element = mesh->FindElementOfPoint(
            physical_point, integration_point, true);
        if (element.IsInvalid() || !element.IsVolume())
            BadArgument("a requested evaluation point is outside the volume mesh");
        auto& transformation = mesh->GetTrafo(element, local_heap);

        auto store_real = [&](auto& mapped_point) {
            ngbla::FlatVector<double> result(value_dimension, local_heap);
            coefficient->Evaluate(mapped_point, result);
            for (std::size_t component = 0; component < value_dimension; ++component)
                values[point * value_dimension + component] =
                    Complex(result[component], 0.0);
        };
        auto store_complex = [&](auto& mapped_point) {
            ngbla::FlatVector<Complex> result(value_dimension, local_heap);
            coefficient->Evaluate(mapped_point, result);
            for (std::size_t component = 0; component < value_dimension; ++component)
                values[point * value_dimension + component] = result[component];
        };
        if (point_dimension == 2) {
            ngfem::MappedIntegrationPoint<2, 2> mapped_point(
                integration_point, transformation);
            if (coefficient->IsComplex())
                store_complex(mapped_point);
            else
                store_real(mapped_point);
        } else {
            ngfem::MappedIntegrationPoint<3, 3> mapped_point(
                integration_point, transformation);
            if (coefficient->IsComplex())
                store_complex(mapped_point);
            else
                store_real(mapped_point);
        }
    }
    if (coefficient->IsComplex()) {
        plhs[0] = ComplexMatrixOutput(values, point_count, value_dimension);
    } else {
        std::vector<double> real_values(values.size());
        for (std::size_t i = 0; i < values.size(); ++i)
            real_values[i] = values[i].real();
        plhs[0] = RealMatrixOutput(real_values, point_count, value_dimension);
    }
}

std::optional<std::vector<double>> OptionalVector3(
    const mxArray* value, const char* name) {
    if (mxIsEmpty(value)) return std::nullopt;
    return FixedRealVector(value, 3, name);
}

std::shared_ptr<radia::ngsolve_bridge::RadiaFieldCoefficient>
NGSolveRadiaField(std::uint64_t handle) {
    auto field = std::dynamic_pointer_cast<
        radia::ngsolve_bridge::RadiaFieldCoefficient>(Coefficient(handle));
    if (!field)
        BadArgument("handle does not refer to an NGSolve RadiaField");
    return field;
}

void NGSolveRadiaFieldCreate(int nlhs, mxArray* plhs[], int nrhs,
                             const mxArray* prhs[]) {
    CheckArity(nrhs, 9, nlhs, 1,
        "handle = radia_mex('ngsolve.radia_field.create', object, field_type, "
        "origin, u_axis, v_axis, w_axis, precision, units)");
    std::optional<double> precision;
    if (!mxIsEmpty(prhs[7])) precision = Scalar(prhs[7], "precision");
    auto field = std::make_shared<
        radia::ngsolve_bridge::RadiaFieldCoefficient>(
            PositiveInteger(prhs[1], "object"), Text(prhs[2], "field_type"),
            OptionalVector3(prhs[3], "origin"),
            OptionalVector3(prhs[4], "u_axis"),
            OptionalVector3(prhs[5], "v_axis"),
            OptionalVector3(prhs[6], "w_axis"), precision,
            Text(prhs[8], "units"));
    plhs[0] = Uint64Output(RegisterCoefficient(std::move(field)));
}

void NGSolveRadiaFieldInfo(int nlhs, mxArray* plhs[], int nrhs,
                           const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
        "info = radia_mex('ngsolve.radia_field.info', handle)");
    const auto field = NGSolveRadiaField(Handle(prhs[1]));
    const char* fields[] = {"radia_obj", "field_type", "use_transform",
                            "has_precision", "precision"};
    plhs[0] = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(plhs[0], 0, "radia_obj",
               mxCreateDoubleScalar(field->Object()));
    mxSetField(plhs[0], 0, "field_type",
               TextOutput(field->FieldType().c_str()));
    mxSetField(plhs[0], 0, "use_transform",
               mxCreateLogicalScalar(field->UsesTransform()));
    mxSetField(plhs[0], 0, "has_precision",
               mxCreateLogicalScalar(field->Precision().has_value()));
    mxSetField(plhs[0], 0, "precision", mxCreateDoubleScalar(
        field->Precision().value_or(std::numeric_limits<double>::quiet_NaN())));
}

void NGSolveRadiaFieldPrepareCache(int nlhs, mxArray* plhs[], int nrhs,
                                   const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 0,
        "radia_mex('ngsolve.radia_field.prepare_cache', handle, points)");
    const auto field = NGSolveRadiaField(Handle(prhs[1]));
    std::size_t rows = 0, cols = 0;
    auto points = RealMatrix(prhs[2], rows, cols, "points");
    if (cols != 3)
        BadArgument("RadiaField cache points must have shape N-by-3");
    field->PrepareCache(points);
}

void NGSolveRadiaFieldClearCache(int nlhs, mxArray* plhs[], int nrhs,
                                 const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 0,
        "radia_mex('ngsolve.radia_field.clear_cache', handle)");
    NGSolveRadiaField(Handle(prhs[1]))->ClearCache();
}

void NGSolveRadiaFieldCacheStats(int nlhs, mxArray* plhs[], int nrhs,
                                 const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
        "stats = radia_mex('ngsolve.radia_field.cache_stats', handle)");
    const auto stats = NGSolveRadiaField(Handle(prhs[1]))->CacheStats();
    const char* fields[] = {"enabled", "size", "hits", "misses", "hit_rate"};
    plhs[0] = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(plhs[0], 0, "enabled", mxCreateLogicalScalar(stats.enabled));
    mxSetField(plhs[0], 0, "size", mxCreateDoubleScalar(stats.size));
    mxSetField(plhs[0], 0, "hits", mxCreateDoubleScalar(stats.hits));
    mxSetField(plhs[0], 0, "misses", mxCreateDoubleScalar(stats.misses));
    mxSetField(plhs[0], 0, "hit_rate", mxCreateDoubleScalar(stats.hit_rate));
}

void NGSolveRadiaFieldAsVoxelCoefficient(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "coefficient = radia_mex('ngsolve.radia_field.as_voxel_coefficient', "
        "field_handle, mesh_handle, resolution)");
    const auto field = NGSolveRadiaField(Handle(prhs[1]));
    const auto& mesh = Mesh(Handle(prhs[2]));
    if (mesh.mesh->GetDimension() != 3 || mesh.mesh->GetNV() == 0)
        BadArgument("RadiaField voxelization requires a nonempty 3D mesh");
    std::array<double, 3> lower{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<double, 3> upper{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()};
    for (std::size_t vertex = 0; vertex < mesh.mesh->GetNV(); ++vertex) {
        const auto point = mesh.mesh->GetPoint<3>(vertex);
        for (int component = 0; component < 3; ++component) {
            lower[component] = std::min(lower[component], point[component]);
            upper[component] = std::max(upper[component], point[component]);
        }
    }
    double maximum_extent = 0.0;
    for (int component = 0; component < 3; ++component)
        maximum_extent = std::max(
            maximum_extent, upper[component] - lower[component]);
    if (!(maximum_extent > 0.0))
        BadArgument("RadiaField voxel mesh has a degenerate bounding box");
    const double margin = 0.01 * maximum_extent;
    for (int component = 0; component < 3; ++component) {
        lower[component] -= margin;
        upper[component] += margin;
    }
    auto voxel = field->AsVoxelCoefficient(
        lower, upper, PositiveInteger(prhs[3], "resolution"));
    plhs[0] = Uint64Output(RegisterCoefficient(std::move(voxel)));
}

void NGSolveGridFunctionVector(int nlhs, mxArray* plhs[], int nrhs,
                               const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "values = radia_mex('ngsolve.grid_function.vector', handle)");
    const auto& holder = GridFunction(Handle(prhs[1]));
    const int multidim = std::max(1, holder.gridfunction->GetMultiDim());
    const auto prototype = holder.gridfunction->GetVectorPtr(0);
    const std::size_t rows = prototype->Size();
    if (prototype->IsComplex()) {
        std::vector<Complex> values(rows * static_cast<std::size_t>(multidim));
        for (int component = 0; component < multidim; ++component) {
            const auto vector = holder.gridfunction->GetVectorPtr(component);
            const auto data = vector->FVComplex();
            if (data.Size() != rows)
                BadArgument("NGSolve GridFunction component sizes differ");
            for (std::size_t row = 0; row < rows; ++row)
                values[row * static_cast<std::size_t>(multidim) +
                       static_cast<std::size_t>(component)] = data[row];
        }
        plhs[0] = ComplexMatrixOutput(values, rows,
                                      static_cast<std::size_t>(multidim));
    } else {
        std::vector<double> values(rows * static_cast<std::size_t>(multidim));
        for (int component = 0; component < multidim; ++component) {
            const auto vector = holder.gridfunction->GetVectorPtr(component);
            const auto data = vector->FVDouble();
            if (data.Size() != rows)
                BadArgument("NGSolve GridFunction component sizes differ");
            for (std::size_t row = 0; row < rows; ++row)
                values[row * static_cast<std::size_t>(multidim) +
                       static_cast<std::size_t>(component)] = data[row];
        }
        plhs[0] = RealMatrixOutput(values, rows,
                                   static_cast<std::size_t>(multidim));
    }
}

void SetNGSolveGridFunctionVector(NGSolveGridFunctionHandle& holder,
                                  const mxArray* value) {
    std::size_t rows = 0;
    std::size_t cols = 0;
    const auto values = ComplexMatrix(value, rows, cols, "values");
    const int multidim = std::max(1, holder.gridfunction->GetMultiDim());
    const auto prototype = holder.gridfunction->GetVectorPtr(0);
    const std::size_t dof_count = prototype->Size();
    if (multidim == 1) {
        if (values.size() != dof_count)
            BadArgument("values must contain exactly the GridFunction DoF count");
    } else if (rows != dof_count || cols != static_cast<std::size_t>(multidim)) {
        BadArgument("values must have shape [dof_count, multidim]");
    }
    if (!prototype->IsComplex()) {
        for (const auto& value_item : values)
            if (value_item.imag() != 0.0)
                BadArgument("complex values cannot be assigned to a real GridFunction");
        for (int component = 0; component < multidim; ++component) {
            auto target = holder.gridfunction->GetVectorPtr(component)->FVDouble();
            for (std::size_t row = 0; row < dof_count; ++row) {
                const std::size_t index = multidim == 1
                    ? row
                    : row * static_cast<std::size_t>(multidim) +
                      static_cast<std::size_t>(component);
                target[row] = values[index].real();
            }
        }
    } else {
        for (int component = 0; component < multidim; ++component) {
            auto target = holder.gridfunction->GetVectorPtr(component)->FVComplex();
            for (std::size_t row = 0; row < dof_count; ++row) {
                const std::size_t index = multidim == 1
                    ? row
                    : row * static_cast<std::size_t>(multidim) +
                      static_cast<std::size_t>(component);
                target[row] = values[index];
            }
        }
    }
}

void NGSolveGridFunctionCreate(int nlhs, mxArray* plhs[], int nrhs,
                               const mxArray* prhs[]) {
    if ((nrhs < 4 || nrhs > 7) || nlhs != 1)
        BadArgument(
            "usage: h = radia_mex('ngsolve.grid_function.create', "
            "vol_path, space, order [, nograds [, name [, complex]]])");
    const std::string path = Text(prhs[1], "vol_path");
    const std::string space_name = Lowercase(Text(prhs[2], "space"));
    const int order = PositiveInteger(prhs[3], "order");
    const bool nograds = nrhs >= 5 ? Boolean(prhs[4], "nograds") : true;
    const std::string name = nrhs >= 6 ? Text(prhs[5], "name") : "gfu";
    const bool complex_space = nrhs >= 7 ? Boolean(prhs[6], "complex") : false;
    if (name.empty())
        BadArgument("name must not be empty");

    ngcore::RegionTaskManager task_manager;
    auto mesh = std::make_shared<ngcomp::MeshAccess>(path);
    auto fespace = MakeNGSolveSpace(
        mesh, space_name, order, nograds, complex_space);
    fespace->Update();
    fespace->FinalizeUpdate();
    auto gridfunction = ngcomp::CreateGridFunction(
        fespace, name, ngcore::Flags());
    if (!gridfunction)
        BadArgument("NGSolve could not create the GridFunction");
    gridfunction->Update();
    auto holder = std::make_unique<NGSolveGridFunctionHandle>();
    holder->mesh = std::move(mesh);
    holder->fespace = std::move(fespace);
    holder->gridfunction = std::move(gridfunction);
    holder->space = space_name;
    holder->order = order;
    holder->nograds = nograds;
    plhs[0] = Uint64Output(RegisterGridFunction(std::move(holder)));
}

void NGSolveGridFunctionFromFESpace(int nlhs, mxArray* plhs[], int nrhs,
                                    const mxArray* prhs[]) {
    if ((nrhs < 2 || nrhs > 3) || nlhs != 1)
        BadArgument(
            "usage: h = radia_mex('ngsolve.grid_function.from_fespace', "
            "fespace_handle [, name])");
    const auto& space_holder = FESpace(Handle(prhs[1]));
    const std::string name = nrhs == 3 ? Text(prhs[2], "name") : "gfu";
    if (name.empty())
        BadArgument("name must not be empty");

    ngcore::RegionTaskManager task_manager;
    auto gridfunction = ngcomp::CreateGridFunction(
        space_holder.fespace, name, ngcore::Flags());
    if (!gridfunction)
        BadArgument("NGSolve could not create the GridFunction");
    gridfunction->Update();
    auto holder = std::make_unique<NGSolveGridFunctionHandle>();
    holder->mesh = space_holder.mesh;
    holder->fespace = space_holder.fespace;
    holder->gridfunction = std::move(gridfunction);
    holder->space = space_holder.space;
    holder->order = space_holder.order;
    holder->nograds = space_holder.nograds;
    plhs[0] = Uint64Output(RegisterGridFunction(std::move(holder)));
}

void NGSolveLinearFormCreate(int nlhs, mxArray* plhs[], int nrhs,
                             const mxArray* prhs[]) {
    if ((nrhs < 3 || nrhs > 5) || nlhs != 1)
        BadArgument(
            "usage: h = radia_mex('ngsolve.linear_form.create', "
            "fespace_handle, source [, value [, label]])");
    const auto& space_holder = FESpace(Handle(prhs[1]));
    const std::string source_name = Lowercase(Text(prhs[2], "source"));
    const Complex source_value = nrhs >= 4
        ? ComplexScalar(prhs[3], "value")
        : Complex(1.0, 0.0);
    if (!space_holder.is_complex && source_value.imag() != 0.0)
        BadArgument("real NGSolve LinearForms do not accept a complex source value");
    const std::string label = nrhs == 5
        ? Text(prhs[4], "label")
        : "radia_matlab_linear_form";
    if (label.empty())
        BadArgument("label must not be empty");

    auto holder = std::make_unique<NGSolveLinearFormHandle>();
    holder->mesh = space_holder.mesh;
    holder->fespace = space_holder.fespace;
    holder->space = space_holder.space;
    holder->source_name = source_name;
    holder->source_value = source_value;
    holder->label = label;
    {
        ngcore::RegionTaskManager task_manager;
        ngcore::Flags flags;
        if (space_holder.is_complex) {
            auto form = std::make_shared<ngcomp::T_LinearForm<Complex>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveLinearIntegrator(
                holder->space, source_name, holder->mesh->GetDimension(),
                source_value, true));
            ngstd::LocalHeap local_heap(1 << 20,
                                        "radia_matlab_persistent_complex_linear_form");
            form->Assemble(local_heap);
            holder->vector = form->GetVectorPtr();
            holder->form = std::move(form);
        } else {
            auto form = std::make_shared<ngcomp::T_LinearForm<double>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveLinearIntegrator(
                holder->space, source_name, holder->mesh->GetDimension(),
                source_value));
            ngstd::LocalHeap local_heap(1 << 20, "radia_matlab_persistent_linear_form");
            form->Assemble(local_heap);
            holder->vector = form->GetVectorPtr();
            holder->form = std::move(form);
        }
    }
    if (!holder->vector)
        BadArgument("NGSolve LinearForm did not produce a vector");
    plhs[0] = Uint64Output(RegisterLinearForm(std::move(holder)));
}

void NGSolveLinearFormCreateFromCoefficient(int nlhs, mxArray* plhs[], int nrhs,
                                            const mxArray* prhs[]) {
    if ((nrhs < 3 || nrhs > 4) || nlhs != 1)
        BadArgument(
            "usage: h = radia_mex('ngsolve.linear_form.create_from_coefficient', "
            "fespace_handle, coefficient_handle [, label])");
    const auto& space_holder = FESpace(Handle(prhs[1]));
    const auto coefficient = Coefficient(Handle(prhs[2]));
    const std::string label = nrhs == 4
        ? Text(prhs[3], "label")
        : "radia_matlab_coefficient_linear_form";
    if (label.empty())
        BadArgument("label must not be empty");

    auto holder = std::make_unique<NGSolveLinearFormHandle>();
    holder->mesh = space_holder.mesh;
    holder->fespace = space_holder.fespace;
    holder->space = space_holder.space;
    holder->source_name = "coefficient";
    holder->label = label;
    holder->coefficient = coefficient;
    {
        ngcore::RegionTaskManager task_manager;
        ngcore::Flags flags;
        if (space_holder.is_complex) {
            auto form = std::make_shared<ngcomp::T_LinearForm<Complex>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveCoefficientLinearIntegrator(
                holder->fespace, coefficient, holder->mesh->GetDimension(), true));
            ngstd::LocalHeap local_heap(1 << 20,
                                        "radia_matlab_coefficient_complex_linear_form");
            form->Assemble(local_heap);
            holder->vector = form->GetVectorPtr();
            holder->form = std::move(form);
        } else {
            auto form = std::make_shared<ngcomp::T_LinearForm<double>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveCoefficientLinearIntegrator(
                holder->fespace, coefficient, holder->mesh->GetDimension()));
            ngstd::LocalHeap local_heap(1 << 20,
                                        "radia_matlab_coefficient_linear_form");
            form->Assemble(local_heap);
            holder->vector = form->GetVectorPtr();
            holder->form = std::move(form);
        }
    }
    if (!holder->vector)
        BadArgument("NGSolve LinearForm did not produce a vector");
    plhs[0] = Uint64Output(RegisterLinearForm(std::move(holder)));
}

void NGSolveLinearFormCreateBoundaryFromCoefficient(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if ((nrhs < 3 || nrhs > 4) || nlhs != 1)
        BadArgument(
            "usage: h = radia_mex('ngsolve.linear_form.create_boundary_from_coefficient', "
            "fespace_handle, coefficient_handle [, label])");
    const auto& space_holder = FESpace(Handle(prhs[1]));
    const auto coefficient = Coefficient(Handle(prhs[2]));
    const std::string label = nrhs == 4
        ? Text(prhs[3], "label")
        : "radia_matlab_boundary_coefficient_linear_form";
    if (label.empty())
        BadArgument("label must not be empty");

    auto holder = std::make_unique<NGSolveLinearFormHandle>();
    holder->mesh = space_holder.mesh;
    holder->fespace = space_holder.fespace;
    holder->space = space_holder.space;
    holder->source_name = "boundary_coefficient";
    holder->label = label;
    holder->coefficient = coefficient;
    {
        ngcore::RegionTaskManager task_manager;
        ngcore::Flags flags;
        if (space_holder.is_complex) {
            auto form = std::make_shared<ngcomp::T_LinearForm<Complex>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveBoundaryCoefficientLinearIntegrator(
                holder->fespace, coefficient, true));
            ngstd::LocalHeap local_heap(1 << 20,
                                        "radia_matlab_boundary_complex_linear_form");
            form->Assemble(local_heap);
            holder->vector = form->GetVectorPtr();
            holder->form = std::move(form);
        } else {
            auto form = std::make_shared<ngcomp::T_LinearForm<double>>(
                holder->fespace, label.c_str(), flags);
            form->AddIntegrator(MakeNGSolveBoundaryCoefficientLinearIntegrator(
                holder->fespace, coefficient));
            ngstd::LocalHeap local_heap(1 << 20,
                                        "radia_matlab_boundary_linear_form");
            form->Assemble(local_heap);
            holder->vector = form->GetVectorPtr();
            holder->form = std::move(form);
        }
    }
    if (!holder->vector)
        BadArgument("NGSolve boundary LinearForm did not produce a vector");
    plhs[0] = Uint64Output(RegisterLinearForm(std::move(holder)));
}

void NGSolveLinearFormInfo(int nlhs, mxArray* plhs[], int nrhs,
                           const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "info = radia_mex('ngsolve.linear_form.info', handle)");
    const auto& holder = LinearForm(Handle(prhs[1]));
    const char* fields[] = {"space", "source", "label", "size",
                            "is_complex", "taskmanager_threads"};
    plhs[0] = mxCreateStructMatrix(1, 1, 6, fields);
    mxSetField(plhs[0], 0, "space", TextOutput(holder.space.c_str()));
    mxSetField(plhs[0], 0, "source", TextOutput(holder.source_name.c_str()));
    mxSetField(plhs[0], 0, "label", TextOutput(holder.label.c_str()));
    mxSetField(plhs[0], 0, "size",
               mxCreateDoubleScalar(holder.vector->Size()));
    mxSetField(plhs[0], 0, "is_complex",
               mxCreateLogicalScalar(holder.vector->IsComplex()));
    mxSetField(plhs[0], 0, "taskmanager_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetNumThreads()));
}

void NGSolveLinearFormVector(int nlhs, mxArray* plhs[], int nrhs,
                             const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "h = radia_mex('ngsolve.linear_form.vector', handle)");
    const auto& linear_form = LinearForm(Handle(prhs[1]));
    auto holder = std::make_unique<NGSolveVectorHandle>();
    holder->vector = linear_form.vector;
    holder->parent_linear_form = linear_form.form;
    holder->is_view = true;
    plhs[0] = Uint64Output(RegisterVector(std::move(holder)));
}

void NGSolveGridFunctionInfo(int nlhs, mxArray* plhs[], int nrhs,
                             const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "info = radia_mex('ngsolve.grid_function.info', handle)");
    const auto& holder = GridFunction(Handle(prhs[1]));
    const auto vector = holder.gridfunction->GetVectorPtr(0);
    const char* fields[] = {"name", "space", "order", "nograds",
                            "dimension", "vertices", "elements",
                            "dof_count", "multidim", "is_complex",
                            "taskmanager_threads"};
    plhs[0] = mxCreateStructMatrix(1, 1, 11, fields);
    mxSetField(plhs[0], 0, "name",
               TextOutput(holder.gridfunction->GetName().c_str()));
    mxSetField(plhs[0], 0, "space", TextOutput(holder.space.c_str()));
    mxSetField(plhs[0], 0, "order", mxCreateDoubleScalar(holder.order));
    mxSetField(plhs[0], 0, "nograds", mxCreateLogicalScalar(holder.nograds));
    mxSetField(plhs[0], 0, "dimension",
               mxCreateDoubleScalar(holder.mesh->GetDimension()));
    mxSetField(plhs[0], 0, "vertices",
               mxCreateDoubleScalar(holder.mesh->GetNV()));
    mxSetField(plhs[0], 0, "elements",
               mxCreateDoubleScalar(holder.mesh->GetNE()));
    mxSetField(plhs[0], 0, "dof_count",
               mxCreateDoubleScalar(static_cast<double>(vector->Size())));
    mxSetField(plhs[0], 0, "multidim",
               mxCreateDoubleScalar(std::max(1, holder.gridfunction->GetMultiDim())));
    mxSetField(plhs[0], 0, "is_complex",
               mxCreateLogicalScalar(vector->IsComplex()));
    mxSetField(plhs[0], 0, "taskmanager_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetNumThreads()));
}

void NGSolveGridFunctionSetVector(int nlhs, mxArray* plhs[], int nrhs,
                                  const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 0,
               "radia_mex('ngsolve.grid_function.set_vector', handle, values)");
    auto& holder = GridFunction(Handle(prhs[1]));
    SetNGSolveGridFunctionVector(holder, prhs[2]);
}

void NGSolveGridFunctionInterpolate(int nlhs, mxArray* plhs[], int nrhs,
                                    const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 0,
               "radia_mex('ngsolve.grid_function.interpolate', gridfunction, coefficient)");
    auto& holder = GridFunction(Handle(prhs[1]));
    const auto coefficient = Coefficient(Handle(prhs[2]));
    ngcore::RegionTaskManager task_manager;
    ngstd::LocalHeap local_heap(1 << 20, "radia_matlab_gridfunction_interpolate");
    ngcomp::SetValues(coefficient, *holder.gridfunction, ngfem::VOL, nullptr,
                      local_heap);
}

void NGSolveGridFunctionAsCoefficient(int nlhs, mxArray* plhs[], int nrhs,
                                      const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "h = radia_mex('ngsolve.grid_function.as_coefficient', handle)");
    const auto& holder = GridFunction(Handle(prhs[1]));
    auto coefficient = std::make_shared<ngcomp::GridFunctionCoefficientFunction>(
        holder.gridfunction);
    plhs[0] = Uint64Output(RegisterCoefficient(std::move(coefficient)));
}

void NGSolveGridFunctionVectorHandle(int nlhs, mxArray* plhs[], int nrhs,
                                     const mxArray* prhs[]) {
    if ((nrhs < 2 || nrhs > 3) || nlhs != 1)
        BadArgument(
            "usage: h = radia_mex('ngsolve.grid_function.vector_handle', "
            "gridfunction [, component])");
    const auto& gridfunction = GridFunction(Handle(prhs[1]));
    const int multidim = std::max(1, gridfunction.gridfunction->GetMultiDim());
    const int component = nrhs == 3
        ? PositiveInteger(prhs[2], "component") - 1
        : 0;
    if (component >= multidim)
        BadArgument("component must be between 1 and the GridFunction multidim count");
    auto holder = std::make_unique<NGSolveVectorHandle>();
    holder->vector = gridfunction.gridfunction->GetVectorPtr(component);
    holder->parent_gridfunction = gridfunction.gridfunction;
    holder->component = component;
    holder->is_view = true;
    plhs[0] = Uint64Output(RegisterVector(std::move(holder)));
}

void NGSolveVectorInfo(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "info = radia_mex('ngsolve.vector.info', handle)");
    const auto& holder = Vector(Handle(prhs[1]));
    const char* fields[] = {"size", "is_complex", "is_view", "component",
                            "taskmanager_threads"};
    plhs[0] = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(plhs[0], 0, "size",
               mxCreateDoubleScalar(static_cast<double>(holder.vector->Size())));
    mxSetField(plhs[0], 0, "is_complex",
               mxCreateLogicalScalar(holder.vector->IsComplex()));
    mxSetField(plhs[0], 0, "is_view",
               mxCreateLogicalScalar(holder.is_view));
    mxSetField(plhs[0], 0, "component",
               mxCreateDoubleScalar(holder.is_view
                   ? static_cast<double>(holder.component + 1)
                   : 0.0));
    mxSetField(plhs[0], 0, "taskmanager_threads",
               mxCreateDoubleScalar(ngcore::TaskManager::GetNumThreads()));
}

void NGSolveVectorCopy(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "h = radia_mex('ngsolve.vector.copy', handle)");
    const auto& source = Vector(Handle(prhs[1]));
    ngcore::RegionTaskManager task_manager;
    auto copy = source.vector->CreateVector();
    std::shared_ptr<ngla::BaseVector> copied = std::move(copy);
    copied->Set(1.0, *source.vector);
    auto holder = std::make_unique<NGSolveVectorHandle>();
    holder->vector = std::move(copied);
    holder->is_view = false;
    plhs[0] = Uint64Output(RegisterVector(std::move(holder)));
}

void NGSolveVectorSetZero(int nlhs, mxArray* plhs[], int nrhs,
                          const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 0,
               "radia_mex('ngsolve.vector.set_zero', handle)");
    auto& holder = Vector(Handle(prhs[1]));
    ngcore::RegionTaskManager task_manager;
    holder.vector->SetZero();
}

void NGSolveVectorScale(int nlhs, mxArray* plhs[], int nrhs,
                        const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 0,
               "radia_mex('ngsolve.vector.scale', scalar, handle)");
    const Complex scalar = ComplexScalar(prhs[1], "scalar");
    auto& holder = Vector(Handle(prhs[2]));
    if (!holder.vector->IsComplex() && scalar.imag() != 0.0)
        BadArgument("a real NGSolve vector cannot be scaled by a complex scalar");
    ngcore::RegionTaskManager task_manager;
    if (holder.vector->IsComplex())
        holder.vector->Scale(scalar);
    else
        holder.vector->Scale(scalar.real());
}

void NGSolveVectorAxpy(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 0,
               "radia_mex('ngsolve.vector.axpy', alpha, y, x)");
    const Complex alpha = ComplexScalar(prhs[1], "alpha");
    auto& y = Vector(Handle(prhs[2]));
    const auto& x = Vector(Handle(prhs[3]));
    if (y.vector->Size() != x.vector->Size())
        BadArgument("NGSolve vector sizes must match for axpy");
    if (y.vector->IsComplex() != x.vector->IsComplex())
        BadArgument("NGSolve vector scalar types must match for axpy");
    if (!y.vector->IsComplex() && alpha.imag() != 0.0)
        BadArgument("a real NGSolve vector cannot use a complex axpy scalar");
    ngcore::RegionTaskManager task_manager;
    if (y.vector->IsComplex())
        y.vector->Add(alpha, *x.vector);
    else
        y.vector->Add(alpha.real(), *x.vector);
}

void NGSolveVectorDot(int nlhs, mxArray* plhs[], int nrhs,
                      const mxArray* prhs[]) {
    if ((nrhs < 3 || nrhs > 4) || nlhs != 1)
        BadArgument(
            "usage: value = radia_mex('ngsolve.vector.dot', left, right "
            "[, conjugate])");
    const auto& left = Vector(Handle(prhs[1]));
    const auto& right = Vector(Handle(prhs[2]));
    if (left.vector->Size() != right.vector->Size())
        BadArgument("NGSolve vector sizes must match for dot");
    if (left.vector->IsComplex() != right.vector->IsComplex())
        BadArgument("NGSolve vector scalar types must match for dot");
    const bool conjugate = nrhs == 4 ? Boolean(prhs[3], "conjugate") : false;
    ngcore::RegionTaskManager task_manager;
    if (left.vector->IsComplex())
        plhs[0] = ComplexScalarOutput(
            left.vector->InnerProductC(*right.vector, conjugate));
    else
        plhs[0] = mxCreateDoubleScalar(
            left.vector->InnerProductD(*right.vector));
}

void NGSolveVectorNorm(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "value = radia_mex('ngsolve.vector.norm', handle)");
    const auto& holder = Vector(Handle(prhs[1]));
    ngcore::RegionTaskManager task_manager;
    plhs[0] = mxCreateDoubleScalar(holder.vector->L2Norm());
}

void NGSolveVectorValues(int nlhs, mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "values = radia_mex('ngsolve.vector.values', handle)");
    const auto& holder = Vector(Handle(prhs[1]));
    const std::size_t size = holder.vector->Size();
    if (holder.vector->IsComplex()) {
        const auto data = holder.vector->FVComplex();
        std::vector<Complex> values(size);
        for (std::size_t i = 0; i < size; ++i)
            values[i] = data[i];
        plhs[0] = ComplexMatrixOutput(values, size, 1);
    } else {
        const auto data = holder.vector->FVDouble();
        std::vector<double> values(size);
        for (std::size_t i = 0; i < size; ++i)
            values[i] = data[i];
        plhs[0] = RealColumn(values);
    }
}

void NGSolveVectorSetValues(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 0,
               "radia_mex('ngsolve.vector.set_values', handle, values)");
    auto& holder = Vector(Handle(prhs[1]));
    std::size_t rows = 0;
    std::size_t cols = 0;
    const auto values = ComplexMatrix(prhs[2], rows, cols, "values");
    if (values.size() != holder.vector->Size())
        BadArgument("values must contain exactly the NGSolve vector size");
    if (!holder.vector->IsComplex()) {
        for (const auto& value : values)
            if (value.imag() != 0.0)
                BadArgument("complex values cannot be assigned to a real NGSolve vector");
        auto data = holder.vector->FVDouble();
        for (std::size_t i = 0; i < values.size(); ++i)
            data[i] = values[i].real();
    } else {
        auto data = holder.vector->FVComplex();
        for (std::size_t i = 0; i < values.size(); ++i)
            data[i] = values[i];
    }
}

bool AllFinite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) { return std::isfinite(value); });
}

void SimulinkStateSpaceCreate(int nlhs, mxArray* plhs[], int nrhs,
                              const mxArray* prhs[]) {
    CheckArity(nrhs, 6, nlhs, 1,
               "h = radia_mex('simulink.state_space.create', A, B, C, D, x0)");
    std::size_t a_rows = 0;
    std::size_t a_cols = 0;
    std::size_t b_rows = 0;
    std::size_t b_cols = 0;
    std::size_t c_rows = 0;
    std::size_t c_cols = 0;
    std::size_t d_rows = 0;
    std::size_t d_cols = 0;
    auto A = RealMatrix(prhs[1], a_rows, a_cols, "A");
    auto B = RealMatrix(prhs[2], b_rows, b_cols, "B");
    auto C = RealMatrix(prhs[3], c_rows, c_cols, "C");
    auto D = RealMatrix(prhs[4], d_rows, d_cols, "D");
    auto initial = RealVector(prhs[5], "x0");
    if (a_rows == 0 || a_rows != a_cols)
        BadArgument("A must be a non-empty square matrix");
    if (b_rows != a_rows)
        BadArgument("B row count must equal the A state dimension");
    if (c_cols != a_rows)
        BadArgument("C column count must equal the A state dimension");
    if (d_rows != c_rows || d_cols != b_cols)
        BadArgument("D must have shape [size(C,1), size(B,2)]");
    if (initial.size() != a_rows)
        BadArgument("x0 must contain exactly the A state dimension");
    if (!AllFinite(A) || !AllFinite(B) || !AllFinite(C) || !AllFinite(D) ||
        !AllFinite(initial))
        BadArgument("state-space matrices and x0 must contain finite values");

    auto holder = std::make_unique<NativeStateSpaceHandle>();
    holder->state_size = a_rows;
    holder->input_size = b_cols;
    holder->output_size = c_rows;
    holder->A = std::move(A);
    holder->B = std::move(B);
    holder->C = std::move(C);
    holder->D = std::move(D);
    holder->initial_state = std::move(initial);
    holder->state = holder->initial_state;
    plhs[0] = Uint64Output(RegisterStateSpace(std::move(holder)));
}

void SimulinkStateSpaceInfo(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "info = radia_mex('simulink.state_space.info', handle)");
    const auto& holder = StateSpace(Handle(prhs[1]));
    const char* fields[] = {"state_size", "input_size", "output_size",
                            "step_count"};
    plhs[0] = mxCreateStructMatrix(1, 1, 4, fields);
    mxSetField(plhs[0], 0, "state_size",
               mxCreateDoubleScalar(static_cast<double>(holder.state_size)));
    mxSetField(plhs[0], 0, "input_size",
               mxCreateDoubleScalar(static_cast<double>(holder.input_size)));
    mxSetField(plhs[0], 0, "output_size",
               mxCreateDoubleScalar(static_cast<double>(holder.output_size)));
    mxSetField(plhs[0], 0, "step_count",
               mxCreateDoubleScalar(static_cast<double>(holder.step_count)));
}

void SimulinkStateSpaceStep(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
               "y = radia_mex('simulink.state_space.step', handle, u)");
    auto& holder = StateSpace(Handle(prhs[1]));
    const auto input = RealVector(prhs[2], "u");
    if (input.size() != holder.input_size)
        BadArgument("u must contain exactly the state-space input dimension");
    if (!AllFinite(input))
        BadArgument("u must contain finite values");

    std::vector<double> output(holder.output_size, 0.0);
    for (std::size_t row = 0; row < holder.output_size; ++row) {
        for (std::size_t col = 0; col < holder.state_size; ++col)
            output[row] += holder.C[row * holder.state_size + col] *
                           holder.state[col];
        for (std::size_t col = 0; col < holder.input_size; ++col)
            output[row] += holder.D[row * holder.input_size + col] *
                           input[col];
    }

    std::vector<double> next_state(holder.state_size, 0.0);
    for (std::size_t row = 0; row < holder.state_size; ++row) {
        for (std::size_t col = 0; col < holder.state_size; ++col)
            next_state[row] += holder.A[row * holder.state_size + col] *
                               holder.state[col];
        for (std::size_t col = 0; col < holder.input_size; ++col)
            next_state[row] += holder.B[row * holder.input_size + col] *
                               input[col];
    }
    holder.state = std::move(next_state);
    ++holder.step_count;
    plhs[0] = RealColumn(output);
}

void SimulinkStateSpaceReset(int nlhs, mxArray* plhs[], int nrhs,
                             const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 0,
               "radia_mex('simulink.state_space.reset', handle)");
    auto& holder = StateSpace(Handle(prhs[1]));
    holder.state = holder.initial_state;
    holder.step_count = 0;
}

void EnergyCreate(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 8, nlhs, 1,
        "h = radia_mex('energy_stop.create', eta, table_r, table_g, offsets, gamma, alpha, b_max)");
    auto material = std::make_unique<EnergyStopMaterial>(
        RealVector(prhs[1], "eta"), RealVector(prhs[2], "table_r"),
        RealVector(prhs[3], "table_g"), IntegerVector(prhs[4], "offsets"),
        RealVector(prhs[5], "gamma"), Scalar(prhs[6], "alpha"),
        Scalar(prhs[7], "b_max"));
    plhs[0] = Uint64Output(Register(std::move(material)));
}

void EnergyInfo(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1, "info = radia_mex('energy_stop.info', handle)");
    EnergyStopMaterial& material = Energy(Handle(prhs[1]));
    const char* fields[] = {"alpha", "b_max", "nu_bound", "branch_count",
                            "state_size", "eta", "gamma"};
    plhs[0] = mxCreateStructMatrix(1, 1, 7, fields);
    mxSetField(plhs[0], 0, "alpha", mxCreateDoubleScalar(material.Alpha()));
    mxSetField(plhs[0], 0, "b_max", mxCreateDoubleScalar(material.BMax()));
    mxSetField(plhs[0], 0, "nu_bound", mxCreateDoubleScalar(material.NuBound()));
    mxSetField(plhs[0], 0, "branch_count",
               mxCreateDoubleScalar(static_cast<double>(material.BranchCount())));
    mxSetField(plhs[0], 0, "state_size",
               mxCreateDoubleScalar(static_cast<double>(material.StateSize())));
    mxSetField(plhs[0], 0, "eta", RealRow(material.Eta()));
    mxSetField(plhs[0], 0, "gamma", RealRow(material.Gamma()));
}

void EnergyState0(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1, "state = radia_mex('energy_stop.state0', handle)");
    EnergyStopMaterial& material = Energy(Handle(prhs[1]));
    std::vector<double> result(material.StateSize());
    material.State0(result.data());
    plhs[0] = RealRow(result);
}

enum class BatchOperation { Forward, Commit, StoredEnergy };

void EnergyBatch(BatchOperation operation, int nlhs, mxArray* plhs[], int nrhs,
                 const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "out = radia_mex('energy_stop.<forward|commit|stored_energy>', handle, B, states)");
    EnergyStopMaterial& material = Energy(Handle(prhs[1]));
    std::size_t b_rows = 0, b_cols = 0, s_rows = 0, s_cols = 0;
    const auto b = RealMatrix(prhs[2], b_rows, b_cols, "B");
    const auto states = RealMatrix(prhs[3], s_rows, s_cols, "states");
    if (b_cols != 3)
        BadArgument("B must have shape count-by-3");
    if (s_rows != b_rows || s_cols != material.StateSize())
        BadArgument("states must have shape count-by-state_size");
    if (b_rows > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("batch is too large");
    const int count = static_cast<int>(b_rows);
    if (operation == BatchOperation::Forward) {
        std::vector<double> result(b_rows * 3);
        material.ForwardBatch(b.data(), states.data(), count, result.data());
        plhs[0] = RealMatrixOutput(result, b_rows, 3);
    } else if (operation == BatchOperation::Commit) {
        std::vector<double> result(b_rows * material.StateSize());
        material.CommitBatch(b.data(), states.data(), count, result.data());
        plhs[0] = RealMatrixOutput(result, b_rows, material.StateSize());
    } else {
        std::vector<double> result(b_rows);
        material.StoredEnergyBatch(b.data(), states.data(), count, result.data());
        plhs[0] = RealMatrixOutput(result, b_rows, 1);
    }
}

void HybridSolve(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1, "x = radia_mex('hybrid_vim.solve', A, b)");
    std::size_t ar = 0, ac = 0, br = 0, bc = 0;
    const auto a = ComplexMatrix(prhs[1], ar, ac, "A");
    const auto b = ComplexMatrix(prhs[2], br, bc, "b");
    const auto result = radia::hybrid_vim::DenseSolve(
        a, static_cast<int>(ar), static_cast<int>(ac),
        b, static_cast<int>(br), static_cast<int>(bc));
    plhs[0] = ComplexMatrixOutput(result, ar, bc);
}

void HybridSchur(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 5, nlhs, 1,
        "S = radia_mex('hybrid_vim.schur', Kkk, Kke, Kek, Kee)");
    std::size_t kkr = 0, kkc = 0, ker = 0, kec = 0;
    std::size_t ekr = 0, ekc = 0, eer = 0, eec = 0;
    const auto kk = ComplexMatrix(prhs[1], kkr, kkc, "Kkk");
    const auto ke = ComplexMatrix(prhs[2], ker, kec, "Kke");
    const auto ek = ComplexMatrix(prhs[3], ekr, ekc, "Kek");
    const auto ee = ComplexMatrix(prhs[4], eer, eec, "Kee");
    const auto result = radia::hybrid_vim::DenseSchurComplement(
        kk, static_cast<int>(kkr), static_cast<int>(kkc),
        ke, static_cast<int>(ker), static_cast<int>(kec),
        ek, static_cast<int>(ekr), static_cast<int>(ekc),
        ee, static_cast<int>(eer), static_cast<int>(eec));
    plhs[0] = ComplexMatrixOutput(result, kkr, kkc);
}

void HybridScalar(const std::string& command, int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[]) {
    Complex result;
    if (command == "hybrid_vim.skin_impedance") {
        CheckArity(nrhs, 4, nlhs, 1,
            "z = radia_mex('hybrid_vim.skin_impedance', s, sigma, mu)");
        result = radia::hybrid_vim::SkinImpedance(
            ComplexScalar(prhs[1], "s"), Scalar(prhs[2], "sigma"),
            Scalar(prhs[3], "mu"));
    } else if (command == "hybrid_vim.sibc_admittance_tail") {
        CheckArity(nrhs, 5, nlhs, 1,
            "y = radia_mex('hybrid_vim.sibc_admittance_tail', s, measure, sigma, mu)");
        result = radia::hybrid_vim::SIBCAdmittanceTail(
            ComplexScalar(prhs[1], "s"), Scalar(prhs[2], "surface_measure"),
            Scalar(prhs[3], "sigma"), Scalar(prhs[4], "mu"));
    } else if (command == "hybrid_vim.sibc_termination_impedance") {
        CheckArity(nrhs, 4, nlhs, 1,
            "z = radia_mex('hybrid_vim.sibc_termination_impedance', s, k_sibc, d)");
        result = radia::hybrid_vim::SIBCSchurTerminationImpedance(
            ComplexScalar(prhs[1], "s"), Scalar(prhs[2], "k_sibc"),
            Scalar(prhs[3], "d"));
    } else {
        CheckArity(nrhs, 4, nlhs, 1,
            "y = radia_mex('hybrid_vim.sibc_termination_admittance', s, k_sibc, d)");
        result = radia::hybrid_vim::SIBCSchurTerminationAdmittance(
            ComplexScalar(prhs[1], "s"), Scalar(prhs[2], "k_sibc"),
            Scalar(prhs[3], "d"));
    }
    plhs[0] = ComplexScalarOutput(result);
}

void CLNLanczos(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if ((nrhs < 3 || nrhs > 5) || nlhs != 1)
        BadArgument("usage: result = radia_mex('cln.lanczos', K, N [, n_iter [, tol]])");
    std::size_t kr = 0, kc = 0, nr = 0, nc = 0;
    const auto K = RealMatrix(prhs[1], kr, kc, "K");
    const auto N = RealMatrix(prhs[2], nr, nc, "N");
    if (kr == 0 || kr != kc || nr != kr || nc != kc)
        BadArgument("K and N must be non-empty square matrices of the same size");
    const int n = MatrixDimension(kr, "K");
    const int n_iter = nrhs >= 4 ? IntegerScalar(prhs[3], "n_iter") : -1;
    const double tol = nrhs >= 5 ? Scalar(prhs[4], "tol") : 1e-30;
    const auto result = radia::cln::lanczos(K.data(), N.data(), n, n_iter, tol);

    const char* fields[] = {
        "L", "R", "Q", "R_diag", "L_tridiag", "n_input", "n_output", "converged"};
    plhs[0] = mxCreateStructMatrix(1, 1, 8, fields);
    mxSetField(plhs[0], 0, "L", RealColumn(result.L));
    mxSetField(plhs[0], 0, "R", RealColumn(result.R));
    mxSetField(plhs[0], 0, "Q", RealMatrixOutput(
        result.Q, static_cast<std::size_t>(result.n_input),
        static_cast<std::size_t>(result.n_output)));
    mxSetField(plhs[0], 0, "R_diag", RealMatrixOutput(
        result.R_diag, static_cast<std::size_t>(result.n_output),
        static_cast<std::size_t>(result.n_output)));
    mxSetField(plhs[0], 0, "L_tridiag", RealMatrixOutput(
        result.L_tridiag, static_cast<std::size_t>(result.n_output),
        static_cast<std::size_t>(result.n_output)));
    mxSetField(plhs[0], 0, "n_input", mxCreateDoubleScalar(result.n_input));
    mxSetField(plhs[0], 0, "n_output", mxCreateDoubleScalar(result.n_output));
    mxSetField(plhs[0], 0, "converged", mxCreateLogicalScalar(result.converged));
}

void CLNBuildTridiagonal(int nlhs, mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "T = radia_mex('cln.build_tridiagonal', diag)");
    const auto diag = RealVector(prhs[1], "diag");
    if (diag.empty())
        BadArgument("diag must be non-empty");
    const int n = MatrixDimension(diag.size(), "diag");
    plhs[0] = RealMatrixOutput(radia::cln::build_tridiagonal(diag.data(), n), n, n);
}

void CLNImpedance(int nlhs, mxArray* plhs[], int nrhs,
                  const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
               "Z = radia_mex('cln.impedance', R_diag, L_tridiag, freq)");
    std::size_t rr = 0, rc = 0, lr = 0, lc = 0;
    const auto R = RealMatrix(prhs[1], rr, rc, "R_diag");
    const auto L = RealMatrix(prhs[2], lr, lc, "L_tridiag");
    if (rr == 0 || rr != rc || lr != rr || lc != rc)
        BadArgument("R_diag and L_tridiag must be non-empty square matrices of the same size");
    plhs[0] = ComplexScalarOutput(radia::cln::compute_cln_impedance(
        R.data(), L.data(), MatrixDimension(rr, "R_diag"),
        Scalar(prhs[3], "freq")));
}

void CLNImpedanceSweep(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
               "Z = radia_mex('cln.impedance_sweep', R_diag, L_tridiag, freqs)");
    std::size_t rr = 0, rc = 0, lr = 0, lc = 0;
    const auto R = RealMatrix(prhs[1], rr, rc, "R_diag");
    const auto L = RealMatrix(prhs[2], lr, lc, "L_tridiag");
    const auto freqs = RealVector(prhs[3], "freqs");
    if (rr == 0 || rr != rc || lr != rr || lc != rc)
        BadArgument("R_diag and L_tridiag must be non-empty square matrices of the same size");
    if (freqs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("freqs is too large");
    std::vector<Complex> result(freqs.size());
    radia::cln::compute_cln_impedance_sweep(
        R.data(), L.data(), MatrixDimension(rr, "R_diag"), freqs.data(),
        static_cast<int>(freqs.size()), result.data());
    plhs[0] = ComplexMatrixOutput(result, freqs.size(), 1);
}

void CLNTransformCoupling(int nlhs, mxArray* plhs[], int nrhs,
                           const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
               "M = radia_mex('cln.transform_coupling', Q, M_LS)");
    std::size_t qr = 0, qc = 0, mr = 0, mc = 0;
    const auto Q = RealMatrix(prhs[1], qr, qc, "Q");
    const auto M = RealMatrix(prhs[2], mr, mc, "M_LS");
    if (qr == 0 || qc == 0 || mr != qr)
        BadArgument("Q and M_LS have incompatible dimensions");
    const int n_loop = MatrixDimension(qr, "Q");
    const int n_reduced = MatrixDimension(qc, "Q columns");
    const int n_star = MatrixDimension(mc, "M_LS columns");
    plhs[0] = RealMatrixOutput(radia::cln::transform_coupling(
        Q.data(), M.data(), n_loop, n_reduced, n_star), qc, mc);
}

void CLNTransformPort(int nlhs, mxArray* plhs[], int nrhs,
                      const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
               "v = radia_mex('cln.transform_port', Q, port)");
    std::size_t qr = 0, qc = 0;
    const auto Q = RealMatrix(prhs[1], qr, qc, "Q");
    const auto port = RealVector(prhs[2], "port");
    if (qr == 0 || qc == 0 || port.size() != qr)
        BadArgument("Q and port have incompatible dimensions");
    plhs[0] = RealColumn(radia::cln::transform_port_vector(
        Q.data(), port.data(), MatrixDimension(qr, "Q"),
        MatrixDimension(qc, "Q columns")));
}

void CLNACACompress(int nlhs, mxArray* plhs[], int nrhs,
                    const mxArray* prhs[]) {
    if ((nrhs < 2 || nrhs > 4) || nlhs != 1)
        BadArgument("usage: result = radia_mex('cln.aca_compress', P [, eps [, kmax]])");
    std::size_t pr = 0, pc = 0;
    const auto P = RealMatrix(prhs[1], pr, pc, "P");
    if (pr == 0 || pr != pc)
        BadArgument("P must be a non-empty square matrix");
    const double eps = nrhs >= 3 ? Scalar(prhs[2], "eps") : 1e-4;
    const int kmax = nrhs >= 4 ? IntegerScalar(prhs[3], "kmax") : -1;
    const auto result = radia::cln::aca_compress(
        P.data(), MatrixDimension(pr, "P"), eps, kmax);
    const char* fields[] = {"U", "V", "n", "k", "compression_ratio", "converged"};
    plhs[0] = mxCreateStructMatrix(1, 1, 6, fields);
    mxSetField(plhs[0], 0, "U", RealMatrixOutput(
        result.U, static_cast<std::size_t>(result.n), static_cast<std::size_t>(result.k)));
    mxSetField(plhs[0], 0, "V", RealMatrixOutput(
        result.V, static_cast<std::size_t>(result.n), static_cast<std::size_t>(result.k)));
    mxSetField(plhs[0], 0, "n", mxCreateDoubleScalar(result.n));
    mxSetField(plhs[0], 0, "k", mxCreateDoubleScalar(result.k));
    mxSetField(plhs[0], 0, "compression_ratio",
               mxCreateDoubleScalar(result.compression_ratio));
    mxSetField(plhs[0], 0, "converged", mxCreateLogicalScalar(result.converged));
}

void EVRSTMethodAlgebra(int nlhs, mxArray* plhs[], int nrhs,
                        const mxArray* prhs[]) {
    CheckArity(nrhs, 8, nlhs, 1,
        "out = radia_mex('evrs.tmethod', curl_map, div_map, grad_map, "
        "evrs_map, resistance_current, inductance_current, port_current)");

    std::size_t curl_rows = 0, curl_cols = 0;
    std::size_t div_rows = 0, div_cols = 0;
    std::size_t grad_rows = 0, grad_cols = 0;
    std::size_t evrs_rows = 0, evrs_cols = 0;
    std::size_t resistance_rows = 0, resistance_cols = 0;
    std::size_t inductance_rows = 0, inductance_cols = 0;
    std::size_t port_rows = 0, port_cols = 0;

    const auto curl = RealMatrix(prhs[1], curl_rows, curl_cols, "curl_map");
    const auto div = RealMatrix(prhs[2], div_rows, div_cols, "div_map");
    const auto grad = RealMatrix(prhs[3], grad_rows, grad_cols, "grad_map");
    const auto evrs = RealMatrix(prhs[4], evrs_rows, evrs_cols, "evrs_map");
    const auto resistance = RealMatrix(prhs[5], resistance_rows, resistance_cols,
                                        "resistance_current");
    const auto inductance = RealMatrix(prhs[6], inductance_rows, inductance_cols,
                                       "inductance_current");
    const auto port = RealMatrix(prhs[7], port_rows, port_cols, "port_current");

    const auto r = radia::evrs::BuildTMethodAlgebra(
        curl, MatrixDimension(curl_rows, "curl_map rows"),
        MatrixDimension(curl_cols, "curl_map columns"),
        div, MatrixDimension(div_rows, "div_map rows"),
        MatrixDimension(div_cols, "div_map columns"),
        grad, MatrixDimension(grad_rows, "grad_map rows"),
        MatrixDimension(grad_cols, "grad_map columns"),
        evrs, MatrixDimension(evrs_rows, "evrs_map rows"),
        MatrixDimension(evrs_cols, "evrs_map columns"),
        resistance, MatrixDimension(resistance_rows, "resistance rows"),
        MatrixDimension(resistance_cols, "resistance columns"),
        inductance, MatrixDimension(inductance_rows, "inductance rows"),
        MatrixDimension(inductance_cols, "inductance columns"),
        port, MatrixDimension(port_rows, "port rows"),
        MatrixDimension(port_cols, "port columns"));

    const char* fields[] = {
        "current_evrs", "resistance_t", "inductance_t", "resistance_evrs",
        "inductance_evrs", "port_t", "port_evrs", "diagnostics"};
    plhs[0] = mxCreateStructMatrix(1, 1, 8, fields);
    mxSetField(plhs[0], 0, "current_evrs",
               RealMatrixOutput(r.current_evrs, r.n_current, r.n_evrs));
    mxSetField(plhs[0], 0, "resistance_t",
               RealMatrixOutput(r.resistance_t, r.n_t, r.n_t));
    mxSetField(plhs[0], 0, "inductance_t",
               RealMatrixOutput(r.inductance_t, r.n_t, r.n_t));
    mxSetField(plhs[0], 0, "resistance_evrs",
               RealMatrixOutput(r.resistance_evrs, r.n_evrs, r.n_evrs));
    mxSetField(plhs[0], 0, "inductance_evrs",
               RealMatrixOutput(r.inductance_evrs, r.n_evrs, r.n_evrs));
    mxSetField(plhs[0], 0, "port_t",
               RealMatrixOutput(r.port_t, r.n_t, r.n_ports));
    mxSetField(plhs[0], 0, "port_evrs",
               RealMatrixOutput(r.port_evrs, r.n_evrs, r.n_ports));

    const char* diagnostic_fields[] = {
        "n_current", "n_t", "n_phi", "n_evrs", "n_ports", "n_rho",
        "div_curl_norm", "div_evrs_norm", "resistance_gauge_norm",
        "inductance_gauge_norm", "port_gauge_norm", "resistance_symmetry_norm",
        "inductance_symmetry_norm", "evrs_resistance_symmetry_norm",
        "evrs_inductance_symmetry_norm", "evrs_resistance_galerkin_residual",
        "evrs_inductance_galerkin_residual"};
    mxArray* diagnostics = mxCreateStructMatrix(1, 1, 17, diagnostic_fields);
    mxSetField(diagnostics, 0, "n_current", mxCreateDoubleScalar(r.n_current));
    mxSetField(diagnostics, 0, "n_t", mxCreateDoubleScalar(r.n_t));
    mxSetField(diagnostics, 0, "n_phi", mxCreateDoubleScalar(r.n_phi));
    mxSetField(diagnostics, 0, "n_evrs", mxCreateDoubleScalar(r.n_evrs));
    mxSetField(diagnostics, 0, "n_ports", mxCreateDoubleScalar(r.n_ports));
    mxSetField(diagnostics, 0, "n_rho", mxCreateDoubleScalar(r.n_rho));
    mxSetField(diagnostics, 0, "div_curl_norm",
               mxCreateDoubleScalar(r.div_curl_norm));
    mxSetField(diagnostics, 0, "div_evrs_norm",
               mxCreateDoubleScalar(r.div_evrs_norm));
    mxSetField(diagnostics, 0, "resistance_gauge_norm",
               mxCreateDoubleScalar(r.resistance_gauge_norm));
    mxSetField(diagnostics, 0, "inductance_gauge_norm",
               mxCreateDoubleScalar(r.inductance_gauge_norm));
    mxSetField(diagnostics, 0, "port_gauge_norm",
               mxCreateDoubleScalar(r.port_gauge_norm));
    mxSetField(diagnostics, 0, "resistance_symmetry_norm",
               mxCreateDoubleScalar(r.resistance_symmetry_norm));
    mxSetField(diagnostics, 0, "inductance_symmetry_norm",
               mxCreateDoubleScalar(r.inductance_symmetry_norm));
    mxSetField(diagnostics, 0, "evrs_resistance_symmetry_norm",
               mxCreateDoubleScalar(r.evrs_resistance_symmetry_norm));
    mxSetField(diagnostics, 0, "evrs_inductance_symmetry_norm",
               mxCreateDoubleScalar(r.evrs_inductance_symmetry_norm));
    mxSetField(diagnostics, 0, "evrs_resistance_galerkin_residual",
               mxCreateDoubleScalar(r.evrs_resistance_galerkin_residual));
    mxSetField(diagnostics, 0, "evrs_inductance_galerkin_residual",
               mxCreateDoubleScalar(r.evrs_inductance_galerkin_residual));
    mxSetField(plhs[0], 0, "diagnostics", diagnostics);
}

void TetHCurlReducedGram(int nlhs, mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]) {
    CheckArity(nrhs, 7, nlhs, 1,
        "gram = radia_mex('hcurl.tet_reduced_gram', cell_verts, exponents, "
        "coefficients, n_modes, ref_points, ref_weights)");

    const auto cell_verts = RealVector(prhs[1], "cell_verts");
    const auto exponents_flat = IntegerVector(prhs[2], "exponents");
    const auto coefficients = RealVector(prhs[3], "coefficients");
    const int n_modes = PositiveInteger(prhs[4], "n_modes");
    const auto ref_points = RealVector(prhs[5], "ref_points");
    const auto ref_weights = RealVector(prhs[6], "ref_weights");

    if (cell_verts.empty() || cell_verts.size() % 12 != 0)
        BadArgument("cell_verts must contain 12 values per tetrahedron");
    if (exponents_flat.empty() || exponents_flat.size() % 3 != 0)
        BadArgument("exponents must contain triples of non-negative integers");
    if (ref_points.empty() || ref_points.size() % 3 != 0 ||
        ref_weights.size() != ref_points.size() / 3)
        BadArgument("ref_points and ref_weights have incompatible sizes");

    const std::size_t n_cells = cell_verts.size() / 12;
    const std::size_t n_monomials = exponents_flat.size() / 3;
    const std::size_t expected_coefficients =
        static_cast<std::size_t>(n_modes) * n_cells * n_monomials * 3;
    if (coefficients.size() != expected_coefficients)
        BadArgument("coefficients must have shape (n_modes,n_cells,n_monomials,3)");
    if (n_cells > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        n_monomials > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("too many cells or monomials");

    std::vector<std::array<int, 3>> exponents(n_monomials);
    for (std::size_t i = 0; i < n_monomials; ++i)
        for (int component = 0; component < 3; ++component)
            exponents[i][static_cast<std::size_t>(component)] =
                exponents_flat[3 * i + static_cast<std::size_t>(component)];

    std::vector<double> gram;
    {
        ngcore::RegionTaskManager task_manager;
        gram = rad_hdiv::TetHCurlReducedGram(
            cell_verts, exponents, coefficients, n_modes,
            ref_points, ref_weights);
    }
    plhs[0] = RealMatrixOutput(gram,
                               static_cast<std::size_t>(n_modes),
                               static_cast<std::size_t>(n_modes));
}

void AffineCellSelfEnergyShapeDerivative(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "result = radia_mex('hdiv.affine_cell_self_energy_shape_derivative', "
        "cell_type, nodes, velocities)");
    const std::string cell_type = Text(prhs[1], "cell_type");
    const int kind = cell_type == "tet" ? 0 :
                     cell_type == "hex" ? 1 :
                     cell_type == "wedge" ? 2 : -1;
    if (kind < 0)
        BadArgument("cell_type must be 'tet', 'hex', or 'wedge'");
    const std::size_t node_count = kind == 0 ? 4u : (kind == 1 ? 8u : 6u);

    std::size_t node_rows = 0, node_cols = 0;
    auto nodes = RealMatrix(prhs[2], node_rows, node_cols, "nodes");
    if (node_rows != node_count || node_cols != 3)
        BadArgument("nodes must have shape element_nodes-by-3");

    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
        mxGetNumberOfDimensions(prhs[3]) != 3)
        BadArgument("velocities must be a real modes-by-element_nodes-by-3 array");
    const mwSize* dims = mxGetDimensions(prhs[3]);
    const std::size_t mode_count = static_cast<std::size_t>(dims[0]);
    if (mode_count == 0 || static_cast<std::size_t>(dims[1]) != node_count ||
        dims[2] != 3)
        BadArgument("velocities must have shape modes-by-element_nodes-by-3");
    if (mode_count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("velocities has too many modes");

    const double* matlab_velocity = mxGetDoubles(prhs[3]);
    std::vector<double> velocities(mode_count * node_count * 3);
    for (std::size_t mode = 0; mode < mode_count; ++mode)
        for (std::size_t node = 0; node < node_count; ++node)
            for (std::size_t component = 0; component < 3; ++component)
                velocities[(mode * node_count + node) * 3 + component] =
                    matlab_velocity[mode + mode_count *
                        (node + node_count * component)];
    const auto finite = [](double value) { return std::isfinite(value); };
    if (!std::all_of(nodes.begin(), nodes.end(), finite) ||
        !std::all_of(velocities.begin(), velocities.end(), finite))
        BadArgument("nodes and velocities must contain finite values");

    const auto result = rad_hdiv::AffineCellSelfEnergyShapeDerivative(
        kind, nodes, velocities, static_cast<int>(mode_count));
    if (result.size() != mode_count + 1)
        throw std::runtime_error(
            "affine self-energy derivative returned an invalid result size");
    const char* fields[] = {"value", "derivative"};
    plhs[0] = mxCreateStructMatrix(1, 1, 2, fields);
    mxSetField(plhs[0], 0, "value", mxCreateDoubleScalar(result[0]));
    mxSetField(plhs[0], 0, "derivative", RealColumn(
        std::vector<double>(result.begin() + 1, result.end())));
}

void StreamACATSVD(int nlhs, mxArray* plhs[], int nrhs,
                   const mxArray* prhs[]) {
    CheckArity(nrhs, 5, nlhs, 4,
        "[U,S,V,k_aca] = radia_mex('stream.aca_tsvd', A, modes, kmax, aca_eps)");
    std::size_t rows = 0, cols = 0;
    const auto matrix = RealMatrix(prhs[1], rows, cols, "A");
    if (rows == 0 || cols == 0 ||
        rows > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        cols > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("A must be a nonempty matrix with native integer dimensions");
    const int modes = PositiveInteger(prhs[2], "modes");
    const int kmax = NonnegativeInteger(prhs[3], "kmax");
    const double aca_eps = Scalar(prhs[4], "aca_eps");
    if (aca_eps <= 0.0)
        BadArgument("aca_eps must be positive");
    if (!std::all_of(matrix.begin(), matrix.end(),
                     [](double value) { return std::isfinite(value); }))
        BadArgument("A must contain finite values");

    const int row_count = static_cast<int>(rows);
    const int col_count = static_cast<int>(cols);
    const auto result = radia::stream_function::ACATSVD(
        row_count, col_count,
        [&matrix, col_count](int row, int col) {
            return matrix[static_cast<std::size_t>(row) * col_count + col];
        }, modes, kmax, aca_eps);
    plhs[0] = RealMatrixOutput(result.U, result.M, result.modes);
    plhs[1] = RealColumn(result.S);
    plhs[2] = RealMatrixOutput(result.V, result.N, result.modes);
    plhs[3] = mxCreateDoubleScalar(result.k_aca);
}

void BiotSavartSegments(const std::string& command, int nlhs,
                        mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if ((nrhs != 5 && nrhs != 6) || nlhs != 2)
        BadArgument("usage: [out_re,out_im] = radia_mex('biot_savart.<h|a>_segments_complex', "
                    "segments, obs, I_re, I_im [, n_threads])");
    const auto segments = RealVector(prhs[1], "segments");
    const auto obs = RealVector(prhs[2], "obs");
    const auto current_re = RealVector(prhs[3], "I_re");
    const auto current_im = RealVector(prhs[4], "I_im");
    if (segments.empty() || segments.size() % 6 != 0)
        BadArgument("segments must contain 6 values per segment");
    if (obs.empty() || obs.size() % 3 != 0)
        BadArgument("obs must contain triples of coordinates");
    const std::size_t n_segments = segments.size() / 6;
    const std::size_t n_obs = obs.size() / 3;
    if (current_re.size() != n_segments || current_im.size() != n_segments)
        BadArgument("I_re and I_im must have one entry per segment");
    if (n_segments > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        n_obs > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("too many segments or observation points");
    const int n_threads = nrhs == 6 ? IntegerScalar(prhs[5], "n_threads") : 0;
    if (n_threads < 0)
        BadArgument("n_threads must be non-negative");

    std::vector<double> out_re(n_obs * 3, 0.0);
    std::vector<double> out_im(n_obs * 3, 0.0);
    if (command == "biot_savart.h_segments_complex") {
        radia::bs::HFromSegmentsComplex(
            segments.data(), static_cast<int>(n_segments), obs.data(),
            static_cast<int>(n_obs), current_re.data(), current_im.data(),
            out_re.data(), out_im.data(), n_threads);
    } else {
        radia::bs::AFromSegmentsComplex(
            segments.data(), static_cast<int>(n_segments), obs.data(),
            static_cast<int>(n_obs), current_re.data(), current_im.data(),
            out_re.data(), out_im.data(), n_threads);
    }
    plhs[0] = RealMatrixOutput(out_re, n_obs, 3);
    plhs[1] = RealMatrixOutput(out_im, n_obs, 3);
}

void BiotSavartTriangles(const std::string& command, int nlhs,
                         mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if ((nrhs != 5 && nrhs != 6) || nlhs != 2)
        BadArgument("usage: [out_re,out_im] = radia_mex('biot_savart.<a|b>_triangles_complex', "
                    "verts, J_re, J_im, obs [, n_threads])");
    const auto vertices = RealVector(prhs[1], "verts");
    const auto current_re = RealVector(prhs[2], "J_re");
    const auto current_im = RealVector(prhs[3], "J_im");
    const auto obs = RealVector(prhs[4], "obs");
    if (vertices.empty() || vertices.size() % 9 != 0)
        BadArgument("verts must contain 9 values per triangle");
    if (obs.empty() || obs.size() % 3 != 0)
        BadArgument("obs must contain triples of coordinates");
    const std::size_t n_triangles = vertices.size() / 9;
    const std::size_t n_obs = obs.size() / 3;
    if (current_re.size() != n_triangles * 3 ||
        current_im.size() != n_triangles * 3)
        BadArgument("J_re and J_im must have shape (n_triangles,3)");
    if (n_triangles > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        n_obs > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("too many triangles or observation points");
    const int n_threads = nrhs == 6 ? IntegerScalar(prhs[5], "n_threads") : 0;
    if (n_threads < 0)
        BadArgument("n_threads must be non-negative");

    std::vector<double> out_re(n_obs * 3, 0.0);
    std::vector<double> out_im(n_obs * 3, 0.0);
    if (command == "biot_savart.a_triangles_complex") {
        radia::bs::AFromTrianglesComplex(
            vertices.data(), static_cast<int>(n_triangles), current_re.data(),
            current_im.data(), obs.data(), static_cast<int>(n_obs),
            out_re.data(), out_im.data(), n_threads);
    } else {
        radia::bs::BFromTrianglesComplex(
            vertices.data(), static_cast<int>(n_triangles), current_re.data(),
            current_im.data(), obs.data(), static_cast<int>(n_obs),
            out_re.data(), out_im.data(), n_threads);
    }
    plhs[0] = RealMatrixOutput(out_re, n_obs, 3);
    plhs[1] = RealMatrixOutput(out_im, n_obs, 3);
}

void BemGalerkin(const std::string& command, int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[]) {
    const bool p2 = command == "bem.assemble_sldl_p2";
    const int required = p2 ? 6 : 4;
    const int maximum = required + 3;
    if (nlhs != 2 || nrhs < required || nrhs > maximum)
        BadArgument(p2
            ? "usage: [SL,DL] = radia_mex('bem.assemble_sldl_p2', verts, tris, "
              "p2_nodes, dofs_per_tri, n_dof [, regular_degree, singular_n_q, n_threads])"
            : "usage: [SL,DL] = radia_mex('bem.assemble_sldl', verts, tris, "
              "p2_nodes [, regular_degree, singular_n_q, n_threads])");

    std::size_t vertex_rows = 0, vertex_cols = 0;
    const auto vertices = RealMatrix(prhs[1], vertex_rows, vertex_cols, "verts");
    if (vertex_rows == 0 || vertex_cols != 3)
        BadArgument("verts must have shape (n_vertices,3)");
    std::size_t triangle_rows = 0, triangle_cols = 0;
    const auto triangles = Integer64Matrix(prhs[2], triangle_rows, triangle_cols,
                                            "tris");
    if (triangle_rows == 0 || triangle_cols != 3)
        BadArgument("tris must have shape (n_triangles,3)");
    std::size_t p2_rows = 0, p2_cols = 0;
    const auto p2_nodes = RealMatrix(prhs[3], p2_rows, p2_cols, "p2_nodes");
    if (p2_rows != triangle_rows || p2_cols != 18)
        BadArgument("p2_nodes must have shape (n_triangles,18)");
    if (vertex_rows > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        triangle_rows > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("too many vertices or triangles");

    int n_dof = static_cast<int>(vertex_rows);
    std::vector<std::int64_t> dofs;
    if (p2) {
        std::size_t dof_rows = 0, dof_cols = 0;
        dofs = Integer64Matrix(prhs[4], dof_rows, dof_cols, "dofs_per_tri");
        if (dof_rows != triangle_rows || dof_cols != 6)
            BadArgument("dofs_per_tri must have shape (n_triangles,6)");
        n_dof = PositiveInteger(prhs[5], "n_dof");
    }
    const int regular_degree = nrhs > (p2 ? 6 : 3)
        ? IntegerScalar(prhs[p2 ? 6 : 4], "regular_quad_degree") : 11;
    const int singular_order = nrhs > (p2 ? 7 : 4)
        ? IntegerScalar(prhs[p2 ? 7 : 5], "singular_n_q") : 8;
    const int n_threads = nrhs > (p2 ? 8 : 5)
        ? IntegerScalar(prhs[p2 ? 8 : 6], "n_threads") : 0;
    if (regular_degree <= 0 || singular_order <= 0 || n_threads < 0)
        BadArgument("quadrature orders must be positive and n_threads non-negative");

    std::vector<double> sl(static_cast<std::size_t>(n_dof) * n_dof, 0.0);
    std::vector<double> dl(static_cast<std::size_t>(n_dof) * n_dof, 0.0);
    if (p2) {
        radia::bem::AssembleSLDL_P2(
            vertices.data(), static_cast<int>(vertex_rows), triangles.data(),
            static_cast<int>(triangle_rows), p2_nodes.data(), dofs.data(), n_dof,
            regular_degree, singular_order, n_threads, sl.data(), dl.data());
    } else {
        radia::bem::AssembleSLDL(
            vertices.data(), static_cast<int>(vertex_rows), triangles.data(),
            static_cast<int>(triangle_rows), p2_nodes.data(), regular_degree,
            singular_order, n_threads, sl.data(), dl.data());
    }
    plhs[0] = RealMatrixOutput(sl, static_cast<std::size_t>(n_dof),
                               static_cast<std::size_t>(n_dof));
    plhs[1] = RealMatrixOutput(dl, static_cast<std::size_t>(n_dof),
                               static_cast<std::size_t>(n_dof));
}

void BEMCreate(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
        "handle = radia_mex('hacapk.bem.create', coordinates, entries)");
    std::size_t coordinate_rows = 0, coordinate_cols = 0;
    std::size_t entry_rows = 0, entry_cols = 0;
    auto coordinates = RealMatrix(prhs[1], coordinate_rows, coordinate_cols,
                                  "coordinates");
    auto entries = RealMatrix(prhs[2], entry_rows, entry_cols, "entries");
    if (coordinate_rows == 0 || coordinate_cols != 3)
        BadArgument("coordinates must have shape (n,3)");
    if (entry_rows != coordinate_rows || entry_cols != coordinate_rows)
        BadArgument("entries must have shape (n,n) matching coordinates");
    if (coordinate_rows > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("too many BEM degrees of freedom");

    auto holder = std::make_unique<BEMHandle>();
    holder->coordinates = std::move(coordinates);
    holder->entries = std::move(entries);
    holder->manager = std::make_unique<HACApKBEMManager>(
        holder->coordinates.data(), holder->entries.data(),
        static_cast<int>(coordinate_rows));
    plhs[0] = Uint64Output(RegisterBEM(std::move(holder)));
}

void BEMBuild(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 7, nlhs, 1,
        "ok = radia_mex('hacapk.bem.build', handle, aca_eps, leaf_size, "
        "eta, max_rank, print_level)");
    BEMHandle& holder = BEM(Handle(prhs[1]));
    RadHACApKParams params = RadHACApKBEMDefaultParams();
    const double aca_eps = Scalar(prhs[2], "aca_eps");
    const int leaf_size = IntegerScalar(prhs[3], "leaf_size");
    const double eta = Scalar(prhs[4], "eta");
    const int max_rank = IntegerScalar(prhs[5], "max_rank");
    const int print_level = IntegerScalar(prhs[6], "print_level");
    if (aca_eps > 0.0) params.aca_eps = aca_eps;
    if (leaf_size > 0) params.leaf_size = leaf_size;
    if (eta > 0.0) params.eta = eta;
    if (max_rank > 0) params.max_rank = max_rank;
    params.print_level = print_level;
    plhs[0] = mxCreateLogicalScalar(holder.manager->BuildHMatrix(params));
}

void BEMMatVec(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
        "y = radia_mex('hacapk.bem.matvec', handle, x)");
    BEMHandle& holder = BEM(Handle(prhs[1]));
    const auto x = RealVector(prhs[2], "x");
    const int n = holder.manager->GetNDOF();
    if (x.size() != static_cast<std::size_t>(n))
        BadArgument("x must have one entry per BEM degree of freedom");
    std::vector<double> y(static_cast<std::size_t>(n));
    holder.manager->MatVec(x, y);
    plhs[0] = RealColumn(y);
}

void BEMInfo(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
        "info = radia_mex('hacapk.bem.info', handle)");
    BEMHandle& holder = BEM(Handle(prhs[1]));
    const auto& stats = holder.manager->GetStats();
    const char* fields[] = {
        "n_dof", "valid", "n_lowrank", "n_dense", "max_rank", "n_leaves",
        "compression", "build_time", "memory_mb", "dense_memory_mb"};
    plhs[0] = mxCreateStructMatrix(1, 1, 10, fields);
    mxSetField(plhs[0], 0, "n_dof",
               mxCreateDoubleScalar(holder.manager->GetNDOF()));
    mxSetField(plhs[0], 0, "valid",
               mxCreateLogicalScalar(holder.manager->IsValid()));
    mxSetField(plhs[0], 0, "n_lowrank", mxCreateDoubleScalar(stats.n_lowrank));
    mxSetField(plhs[0], 0, "n_dense", mxCreateDoubleScalar(stats.n_dense));
    mxSetField(plhs[0], 0, "max_rank", mxCreateDoubleScalar(stats.max_rank));
    mxSetField(plhs[0], 0, "n_leaves", mxCreateDoubleScalar(stats.n_leaves));
    mxSetField(plhs[0], 0, "compression", mxCreateDoubleScalar(stats.compression));
    mxSetField(plhs[0], 0, "build_time", mxCreateDoubleScalar(stats.build_time));
    mxSetField(plhs[0], 0, "memory_mb", mxCreateDoubleScalar(stats.memory_mb));
    mxSetField(plhs[0], 0, "dense_memory_mb",
               mxCreateDoubleScalar(stats.dense_memory_mb));
}

void PEECCreate(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 7, nlhs, 1,
        "handle = radia_mex('hacapk.peec.create', centers, directions, "
        "lengths, widths, heights, sigmas)");
    std::size_t center_rows = 0, center_cols = 0;
    std::size_t direction_rows = 0, direction_cols = 0;
    const auto centers = RealMatrix(prhs[1], center_rows, center_cols, "centers");
    const auto directions = RealMatrix(prhs[2], direction_rows, direction_cols,
                                       "directions");
    const auto lengths = RealVector(prhs[3], "lengths");
    const auto widths = RealVector(prhs[4], "widths");
    const auto heights = RealVector(prhs[5], "heights");
    const auto sigmas = RealVector(prhs[6], "sigmas");
    if (center_rows == 0 || center_cols != 3 || direction_rows != center_rows ||
        direction_cols != 3 || lengths.size() != center_rows ||
        widths.size() != center_rows || heights.size() != center_rows ||
        sigmas.size() != center_rows)
        BadArgument("PEEC geometry arrays must have matching N-by-3/N shapes");
    if (center_rows > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("too many PEEC filaments");

    auto holder = std::make_unique<PEECHandle>();
    holder->builder = std::make_unique<radia::PEECMatrixBuilder>();
    for (std::size_t i = 0; i < center_rows; ++i) {
        radia::PEECSegment segment(
            TVector3d(centers[i * 3 + 0], centers[i * 3 + 1], centers[i * 3 + 2]),
            TVector3d(directions[i * 3 + 0], directions[i * 3 + 1],
                      directions[i * 3 + 2]),
            lengths[i], widths[i], heights[i], sigmas[i]);
        holder->builder->AddSegment(segment);
    }
    holder->manager = std::make_unique<HACApKPEECManager>(*holder->builder);
    plhs[0] = Uint64Output(RegisterPEEC(std::move(holder)));
}

void PEECBuild(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 7, nlhs, 1,
        "ok = radia_mex('hacapk.peec.build', handle, aca_eps, leaf_size, "
        "eta, max_rank, print_level)");
    PEECHandle& holder = PEEC(Handle(prhs[1]));
    RadHACApKParams params = RadHACApKPEECDefaultParams();
    const double aca_eps = Scalar(prhs[2], "aca_eps");
    const int leaf_size = IntegerScalar(prhs[3], "leaf_size");
    const double eta = Scalar(prhs[4], "eta");
    const int max_rank = IntegerScalar(prhs[5], "max_rank");
    params.print_level = IntegerScalar(prhs[6], "print_level");
    if (aca_eps > 0.0) params.aca_eps = aca_eps;
    if (leaf_size > 0) params.leaf_size = leaf_size;
    if (eta > 0.0) params.eta = eta;
    if (max_rank > 0) params.max_rank = max_rank;
    plhs[0] = mxCreateLogicalScalar(holder.manager->BuildHMatrix(params));
}

void PEECMatVec(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
        "y = radia_mex('hacapk.peec.matvec', handle, x)");
    PEECHandle& holder = PEEC(Handle(prhs[1]));
    const auto x = RealVector(prhs[2], "x");
    const int n = holder.manager->GetNDOF();
    if (n <= 0 || x.size() != static_cast<std::size_t>(n))
        BadArgument("x must have one entry per PEEC filament after build");
    std::vector<double> y(static_cast<std::size_t>(n));
    holder.manager->MatVec(x, y);
    plhs[0] = RealColumn(y);
}

void PEECInfo(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
        "info = radia_mex('hacapk.peec.info', handle)");
    PEECHandle& holder = PEEC(Handle(prhs[1]));
    const auto& stats = holder.manager->GetStats();
    const char* fields[] = {
        "n_dof", "valid", "n_lowrank", "n_dense", "max_rank", "n_leaves",
        "compression", "build_time", "memory_mb", "dense_memory_mb"};
    plhs[0] = mxCreateStructMatrix(1, 1, 10, fields);
    mxSetField(plhs[0], 0, "n_dof",
               mxCreateDoubleScalar(holder.manager->GetNDOF()));
    mxSetField(plhs[0], 0, "valid",
               mxCreateLogicalScalar(holder.manager->IsValid()));
    mxSetField(plhs[0], 0, "n_lowrank", mxCreateDoubleScalar(stats.n_lowrank));
    mxSetField(plhs[0], 0, "n_dense", mxCreateDoubleScalar(stats.n_dense));
    mxSetField(plhs[0], 0, "max_rank", mxCreateDoubleScalar(stats.max_rank));
    mxSetField(plhs[0], 0, "n_leaves", mxCreateDoubleScalar(stats.n_leaves));
    mxSetField(plhs[0], 0, "compression", mxCreateDoubleScalar(stats.compression));
    mxSetField(plhs[0], 0, "build_time", mxCreateDoubleScalar(stats.build_time));
    mxSetField(plhs[0], 0, "memory_mb", mxCreateDoubleScalar(stats.memory_mb));
    mxSetField(plhs[0], 0, "dense_memory_mb",
               mxCreateDoubleScalar(stats.dense_memory_mb));
}

void ChargeGramCreateMonopole(int nlhs, mxArray* plhs[], int nrhs,
                              const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "handle = radia_mex('hacapk.charge_gram.create_monopole', "
        "centroids, measures, self_energy)");
    std::size_t rows = 0, cols = 0;
    auto centroids = RealMatrix(prhs[1], rows, cols, "centroids");
    auto measures = RealVector(prhs[2], "measures");
    auto self_energy = RealVector(prhs[3], "self_energy");
    if (rows == 0 || cols != 3 || measures.size() != rows ||
        self_energy.size() != rows)
        BadArgument("centroids must be N-by-3 and measures/self_energy must have N entries");
    auto holder = std::make_unique<ChargeGramHandle>();
    holder->n_dof = static_cast<int>(rows);
    holder->manager = std::make_unique<HACApKChargeGram>(
        std::move(centroids), std::move(measures), std::move(self_energy));
    plhs[0] = Uint64Output(RegisterChargeGram(std::move(holder)));
}

void ChargeGramCreateSampled(const std::string& command, int nlhs,
                             mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    const bool planar = command == "hacapk.charge_gram.create_sampled_planar_log";
    CheckArity(nrhs, planar ? 5 : 4, nlhs, 1,
        planar
            ? "handle = radia_mex('hacapk.charge_gram.create_sampled_planar_log', points, weights, kernel_epsilon, reference_length)"
            : "handle = radia_mex('hacapk.charge_gram.create_sampled_laplace', points, weights, kernel_epsilon)");
    std::size_t rows = 0, cols = 0;
    auto points = RealMatrix(prhs[1], rows, cols, "points");
    auto weights = RealVector(prhs[2], "weights");
    if (rows == 0 || cols != 3 || weights.size() != rows)
        BadArgument("points must be N-by-3 and weights must have N entries");
    const double kernel_epsilon = Scalar(prhs[3], "kernel_epsilon");
    auto holder = std::make_unique<ChargeGramHandle>();
    holder->n_dof = static_cast<int>(rows);
    if (planar) {
        holder->manager = std::make_unique<HACApKChargeGram>(
            std::move(points), std::move(weights), kernel_epsilon,
            Scalar(prhs[4], "reference_length"));
    } else {
        holder->manager = std::make_unique<HACApKChargeGram>(
            std::move(points), std::move(weights), kernel_epsilon);
    }
    plhs[0] = Uint64Output(RegisterChargeGram(std::move(holder)));
}

void ChargeGramCreateAnalyticTet(int nlhs, mxArray* plhs[], int nrhs,
                                 const mxArray* prhs[]) {
    CheckArity(nrhs, 8, nlhs, 1,
        "handle = radia_mex('hacapk.charge_gram.create_analytic_tet', "
        "cell_verts, face_verts, n_el, near_factor, image_masks, "
        "image_signs, far_quad)");
    const int n_el = PositiveInteger(prhs[3], "n_el");
    std::size_t cell_rows = 0, cell_cols = 0;
    auto cell_verts = RealMatrix(prhs[1], cell_rows, cell_cols, "cell_verts");
    if (cell_rows != 4u * static_cast<std::size_t>(n_el) || cell_cols != 3)
        BadArgument("cell_verts must have shape (4*n_el)-by-3");
    std::size_t face_rows = 0, face_cols = 0;
    auto face_verts = RealMatrix(prhs[2], face_rows, face_cols, "face_verts");
    if (face_rows != 0 && (face_cols != 3 || face_rows % 3 != 0))
        BadArgument("face_verts must have shape (3*n_face)-by-3");
    auto image_masks = IntegerVector(prhs[5], "image_masks");
    auto image_signs = RealVector(prhs[6], "image_signs");
    if (image_masks.size() != image_signs.size())
        BadArgument("image_masks and image_signs must have equal lengths");
    const int far_quad = IntegerScalar(prhs[7], "far_quad");
    if (far_quad < 0)
        BadArgument("far_quad must be nonnegative");
    auto holder = std::make_unique<ChargeGramHandle>();
    holder->n_dof = n_el + static_cast<int>(face_rows / 3);
    holder->manager = std::make_unique<HACApKChargeGram>(
        std::move(cell_verts), std::move(face_verts), n_el,
        Scalar(prhs[4], "near_factor"), std::move(image_masks),
        std::move(image_signs), far_quad);
    plhs[0] = Uint64Output(RegisterChargeGram(std::move(holder)));
}

void ChargeGramCreateAnalyticPolytope(int nlhs, mxArray* plhs[], int nrhs,
                                      const mxArray* prhs[]) {
    CheckArity(nrhs, 14, nlhs, 1,
        "handle = radia_mex('hacapk.charge_gram.create_analytic_polytope', "
        "cell_tris, cell_troff, cell_cent, cell_meas, face_tris, "
        "face_troff, face_cent, face_meas, n_el, near_factor, "
        "image_masks, image_signs, far_quad)");
    const int n_el = PositiveInteger(prhs[9], "n_el");
    std::size_t cell_tri_rows = 0, cell_tri_cols = 0;
    auto cell_tris = RealMatrix(prhs[1], cell_tri_rows, cell_tri_cols, "cell_tris");
    if (cell_tri_rows == 0 || cell_tri_cols != 3 || cell_tri_rows % 3 != 0)
        BadArgument("cell_tris must have shape (3*n_triangle)-by-3");
    auto cell_troff = IntegerVector(prhs[2], "cell_troff");
    std::size_t cell_cent_rows = 0, cell_cent_cols = 0;
    auto cell_cent = RealMatrix(prhs[3], cell_cent_rows, cell_cent_cols, "cell_cent");
    auto cell_meas = RealVector(prhs[4], "cell_meas");
    if (cell_cent_rows != static_cast<std::size_t>(n_el) || cell_cent_cols != 3 ||
        cell_meas.size() != static_cast<std::size_t>(n_el))
        BadArgument("cell_cent must be n_el-by-3 with n_el cell_meas entries");

    std::size_t face_tri_rows = 0, face_tri_cols = 0;
    auto face_tris = RealMatrix(prhs[5], face_tri_rows, face_tri_cols, "face_tris");
    if (face_tri_rows != 0 && (face_tri_cols != 3 || face_tri_rows % 3 != 0))
        BadArgument("face_tris must have shape (3*n_triangle)-by-3");
    auto face_troff = IntegerVector(prhs[6], "face_troff");
    if (face_troff.empty())
        BadArgument("face_troff must contain at least the initial zero");
    const std::size_t n_face = face_troff.size() - 1;
    std::size_t face_cent_rows = 0, face_cent_cols = 0;
    auto face_cent = RealMatrix(prhs[7], face_cent_rows, face_cent_cols, "face_cent");
    auto face_meas = RealVector(prhs[8], "face_meas");
    if (face_cent_rows != n_face || (n_face != 0 && face_cent_cols != 3) ||
        face_meas.size() != n_face)
        BadArgument("face_cent must be n_face-by-3 with n_face face_meas entries");

    auto validate_offsets = [](const std::vector<int>& offsets,
                               std::size_t entities, std::size_t triangles,
                               const char* name) {
        if (offsets.size() != entities + 1 || offsets.front() != 0 ||
            offsets.back() != static_cast<int>(triangles))
            BadArgument(std::string(name) +
                        " must be a zero-based CSR offset vector ending at n_triangle");
        for (std::size_t i = 1; i < offsets.size(); ++i)
            if (offsets[i] < offsets[i - 1])
                BadArgument(std::string(name) + " must be nondecreasing");
    };
    validate_offsets(cell_troff, static_cast<std::size_t>(n_el),
                     cell_tri_rows / 3, "cell_troff");
    validate_offsets(face_troff, n_face, face_tri_rows / 3, "face_troff");

    auto image_masks = IntegerVector(prhs[11], "image_masks");
    auto image_signs = RealVector(prhs[12], "image_signs");
    if (image_masks.size() != image_signs.size())
        BadArgument("image_masks and image_signs must have equal lengths");
    const int far_quad = IntegerScalar(prhs[13], "far_quad");
    if (far_quad < 0)
        BadArgument("far_quad must be nonnegative");
    auto holder = std::make_unique<ChargeGramHandle>();
    holder->n_dof = n_el + static_cast<int>(n_face);
    holder->manager = std::make_unique<HACApKChargeGram>(
        std::move(cell_tris), std::move(cell_troff), std::move(cell_cent),
        std::move(cell_meas), std::move(face_tris), std::move(face_troff),
        std::move(face_cent), std::move(face_meas), n_el,
        Scalar(prhs[10], "near_factor"), std::move(image_masks),
        std::move(image_signs), far_quad);
    plhs[0] = Uint64Output(RegisterChargeGram(std::move(holder)));
}

void ChargeGramCreateHighOrderTet(int nlhs, mxArray* plhs[], int nrhs,
                                  const mxArray* prhs[]) {
    CheckArity(nrhs, 22, nlhs, 1,
        "handle = radia_mex('hacapk.charge_gram.create_high_order_tet', "
        "cell_verts, face_verts, n_el, charge_host, charge_kind, charge_expo, "
        "ref_tet_pts, ref_tet_w, ref_tri_pts, ref_tri_w, ref_tet_pts_lo, "
        "ref_tet_w_lo, ref_tri_pts_lo, ref_tri_w_lo, ho_far_factor, "
        "ref_tet_pts_in, ref_tet_w_in, ref_tri_pts_in, ref_tri_w_in, "
        "image_masks, image_signs)");
    const int n_el = PositiveInteger(prhs[3], "n_el");
    std::size_t cell_rows = 0, cell_cols = 0;
    auto cell_verts = RealMatrix(prhs[1], cell_rows, cell_cols, "cell_verts");
    if (cell_rows != 4u * static_cast<std::size_t>(n_el) || cell_cols != 3)
        BadArgument("cell_verts must have shape (4*n_el)-by-3");
    std::size_t face_rows = 0, face_cols = 0;
    auto face_verts = RealMatrix(prhs[2], face_rows, face_cols, "face_verts");
    if (face_rows != 0 && (face_cols != 3 || face_rows % 3 != 0))
        BadArgument("face_verts must have shape (3*n_face)-by-3");
    const int n_face = static_cast<int>(face_rows / 3);

    auto charge_host = IntegerVector(prhs[4], "charge_host");
    auto charge_kind = IntegerVector(prhs[5], "charge_kind");
    std::size_t exponent_rows = 0, exponent_cols = 0;
    auto charge_expo = IntegerMatrix(
        prhs[6], exponent_rows, exponent_cols, "charge_expo");
    const std::size_t n_charge = charge_host.size();
    ValidateChargeDescriptors(charge_host, charge_kind, charge_expo,
                              exponent_rows, exponent_cols, n_el, n_face);
    std::vector<double> ref_tet_pts, ref_tet_w, ref_tri_pts, ref_tri_w;
    std::vector<double> ref_tet_pts_lo, ref_tet_w_lo, ref_tri_pts_lo, ref_tri_w_lo;
    std::vector<double> ref_tet_pts_in, ref_tet_w_in, ref_tri_pts_in, ref_tri_w_in;
    ReadQuadrature(prhs[7], prhs[8], 3, true, "ref_tet_pts", "ref_tet_w",
                   ref_tet_pts, ref_tet_w);
    ReadQuadrature(prhs[9], prhs[10], 2, true, "ref_tri_pts", "ref_tri_w",
                   ref_tri_pts, ref_tri_w);
    ReadQuadrature(prhs[11], prhs[12], 3, false, "ref_tet_pts_lo", "ref_tet_w_lo",
                   ref_tet_pts_lo, ref_tet_w_lo);
    ReadQuadrature(prhs[13], prhs[14], 2, false, "ref_tri_pts_lo", "ref_tri_w_lo",
                   ref_tri_pts_lo, ref_tri_w_lo);
    if (ref_tet_w_lo.empty() != ref_tri_w_lo.empty())
        BadArgument("low tet and triangle quadrature rules must be supplied together");
    const double ho_far_factor = Scalar(prhs[15], "ho_far_factor");
    if (ho_far_factor <= 0.0)
        BadArgument("ho_far_factor must be positive");
    ReadQuadrature(prhs[16], prhs[17], 3, false, "ref_tet_pts_in", "ref_tet_w_in",
                   ref_tet_pts_in, ref_tet_w_in);
    ReadQuadrature(prhs[18], prhs[19], 2, false, "ref_tri_pts_in", "ref_tri_w_in",
                   ref_tri_pts_in, ref_tri_w_in);
    if (ref_tet_w_in.empty() != ref_tri_w_in.empty())
        BadArgument("inner tet and triangle quadrature rules must be supplied together");
    auto image_masks = IntegerVector(prhs[20], "image_masks");
    auto image_signs = RealVector(prhs[21], "image_signs");
    if (image_masks.size() != image_signs.size())
        BadArgument("image_masks and image_signs must have equal lengths");
    for (double sign : image_signs)
        if (!std::isfinite(sign))
            BadArgument("image_signs must contain finite values");

    auto holder = std::make_unique<ChargeGramHandle>();
    holder->n_dof = static_cast<int>(n_charge);
    holder->manager = std::make_unique<HACApKChargeGram>(
        std::move(cell_verts), std::move(face_verts), n_el,
        std::move(charge_host), std::move(charge_kind), std::move(charge_expo),
        std::move(ref_tet_pts), std::move(ref_tet_w),
        std::move(ref_tri_pts), std::move(ref_tri_w),
        std::move(ref_tet_pts_lo), std::move(ref_tet_w_lo),
        std::move(ref_tri_pts_lo), std::move(ref_tri_w_lo), ho_far_factor,
        std::move(ref_tet_pts_in), std::move(ref_tet_w_in),
        std::move(ref_tri_pts_in), std::move(ref_tri_w_in),
        std::move(image_masks), std::move(image_signs));
    plhs[0] = Uint64Output(RegisterChargeGram(std::move(holder)));
}

void ChargeGramCreateCurvedHighOrderTet(int nlhs, mxArray* plhs[], int nrhs,
                                        const mxArray* prhs[]) {
    CheckArity(nrhs, 24, nlhs, 1,
        "handle = radia_mex('hacapk.charge_gram.create_curved_high_order_tet', "
        "cell_nodes, face_nodes, cell_vertices, face_vertices, n_el, curve_order, "
        "charge_host, charge_kind, charge_expo, ref_tet_pts, ref_tet_w, "
        "ref_tri_pts, ref_tri_w, curve_gl, curve_gw, ref_tet_pts_lo, "
        "ref_tet_w_lo, ref_tri_pts_lo, ref_tri_w_lo, ho_far_factor, "
        "image_masks, image_signs, reference_density)");
    const int n_el = PositiveInteger(prhs[5], "n_el");
    const int curve_order = PositiveInteger(prhs[6], "curve_order");
    if (curve_order != 2)
        BadArgument("curve_order must be 2 for the P2 curved tetrahedron contract");
    std::size_t cell_node_rows = 0, cell_node_cols = 0;
    auto cell_nodes = RealMatrix(prhs[1], cell_node_rows, cell_node_cols, "cell_nodes");
    if (cell_node_rows != 10u * static_cast<std::size_t>(n_el) || cell_node_cols != 3)
        BadArgument("cell_nodes must have shape (10*n_el)-by-3");
    std::size_t face_node_rows = 0, face_node_cols = 0;
    auto face_nodes = RealMatrix(prhs[2], face_node_rows, face_node_cols, "face_nodes");
    if (face_node_rows != 0 && (face_node_cols != 3 || face_node_rows % 6 != 0))
        BadArgument("face_nodes must have shape (6*n_face)-by-3");
    const int n_face = static_cast<int>(face_node_rows / 6);

    std::size_t cell_vertex_rows = 0, cell_vertex_cols = 0;
    auto cell_vertices = IntegerMatrix(
        prhs[3], cell_vertex_rows, cell_vertex_cols, "cell_vertices");
    if (cell_vertex_rows != static_cast<std::size_t>(n_el) || cell_vertex_cols != 4)
        BadArgument("cell_vertices must have shape n_el-by-4");
    std::size_t face_vertex_rows = 0, face_vertex_cols = 0;
    auto face_vertices = IntegerMatrix(
        prhs[4], face_vertex_rows, face_vertex_cols, "face_vertices");
    if (face_vertex_rows != static_cast<std::size_t>(n_face) ||
        (n_face != 0 && face_vertex_cols != 3))
        BadArgument("face_vertices must have shape n_face-by-3");
    for (int vertex : cell_vertices)
        if (vertex < 0)
            BadArgument("cell_vertices must use nonnegative zero-based vertex IDs");
    for (int vertex : face_vertices)
        if (vertex < 0)
            BadArgument("face_vertices must use nonnegative zero-based vertex IDs");

    auto charge_host = IntegerVector(prhs[7], "charge_host");
    auto charge_kind = IntegerVector(prhs[8], "charge_kind");
    std::size_t exponent_rows = 0, exponent_cols = 0;
    auto charge_expo = IntegerMatrix(
        prhs[9], exponent_rows, exponent_cols, "charge_expo");
    ValidateChargeDescriptors(charge_host, charge_kind, charge_expo,
                              exponent_rows, exponent_cols, n_el, n_face);

    std::vector<double> ref_tet_pts, ref_tet_w, ref_tri_pts, ref_tri_w;
    std::vector<double> ref_tet_pts_lo, ref_tet_w_lo, ref_tri_pts_lo, ref_tri_w_lo;
    ReadQuadrature(prhs[10], prhs[11], 3, true, "ref_tet_pts", "ref_tet_w",
                   ref_tet_pts, ref_tet_w);
    ReadQuadrature(prhs[12], prhs[13], 2, true, "ref_tri_pts", "ref_tri_w",
                   ref_tri_pts, ref_tri_w);
    auto curve_gl = RealVector(prhs[14], "curve_gl");
    auto curve_gw = RealVector(prhs[15], "curve_gw");
    if (curve_gl.empty() || curve_gl.size() != curve_gw.size())
        BadArgument("curve_gl and curve_gw must be nonempty vectors of equal length");
    for (double point : curve_gl)
        if (!std::isfinite(point) || point < 0.0 || point > 1.0)
            BadArgument("curve_gl points must be finite and lie in [0,1]");
    for (double weight : curve_gw)
        if (!std::isfinite(weight) || weight <= 0.0)
            BadArgument("curve_gw weights must be positive and finite");
    ReadQuadrature(prhs[16], prhs[17], 3, false, "ref_tet_pts_lo", "ref_tet_w_lo",
                   ref_tet_pts_lo, ref_tet_w_lo);
    ReadQuadrature(prhs[18], prhs[19], 2, false, "ref_tri_pts_lo", "ref_tri_w_lo",
                   ref_tri_pts_lo, ref_tri_w_lo);
    if (ref_tet_w_lo.empty() != ref_tri_w_lo.empty())
        BadArgument("low tet and triangle quadrature rules must be supplied together");
    const double ho_far_factor = Scalar(prhs[20], "ho_far_factor");
    if (ho_far_factor <= 0.0)
        BadArgument("ho_far_factor must be positive");
    auto image_masks = IntegerVector(prhs[21], "image_masks");
    auto image_signs = RealVector(prhs[22], "image_signs");
    if (image_masks.size() != image_signs.size())
        BadArgument("image_masks and image_signs must have equal lengths");
    for (double sign : image_signs)
        if (!std::isfinite(sign))
            BadArgument("image_signs must contain finite values");
    const bool reference_density = Boolean(prhs[23], "reference_density");

    auto holder = std::make_unique<ChargeGramHandle>();
    holder->n_dof = static_cast<int>(charge_host.size());
    holder->manager = std::make_unique<HACApKChargeGram>(
        std::move(cell_nodes), std::move(face_nodes),
        std::move(cell_vertices), std::move(face_vertices), n_el, curve_order,
        std::move(charge_host), std::move(charge_kind), std::move(charge_expo),
        std::move(ref_tet_pts), std::move(ref_tet_w),
        std::move(ref_tri_pts), std::move(ref_tri_w),
        std::move(curve_gl), std::move(curve_gw),
        std::move(ref_tet_pts_lo), std::move(ref_tet_w_lo),
        std::move(ref_tri_pts_lo), std::move(ref_tri_w_lo), ho_far_factor,
        std::move(image_masks), std::move(image_signs), reference_density);
    plhs[0] = Uint64Output(RegisterChargeGram(std::move(holder)));
}

void ChargeGramCreateHex(int nlhs, mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]) {
    CheckArity(nrhs, 24, nlhs, 1,
        "handle = radia_mex('hacapk.charge_gram.create_hex', hex_cell_nodes, "
        "quad_face_nodes, n_el, n_bf, charge_host, charge_kind, charge_expo, "
        "sym_tet_pts, sym_tet_w, sym_tri_pts, sym_tri_w, gl_out, gw_out, "
        "gl_in, gw_in, far_tet_pts, far_tet_w, far_tri_pts, far_tri_w, "
        "near_grade, far_inner_factor, image_masks, image_signs)");
    const int n_el = PositiveInteger(prhs[3], "n_el");
    const int n_bf = NonnegativeInteger(prhs[4], "n_bf");
    std::size_t cell_rows = 0, cell_cols = 0;
    auto cell_nodes = RealMatrix(prhs[1], cell_rows, cell_cols, "hex_cell_nodes");
    if (cell_rows != 27u * static_cast<std::size_t>(n_el) || cell_cols != 3)
        BadArgument("hex_cell_nodes must have shape (27*n_el)-by-3");
    std::size_t face_rows = 0, face_cols = 0;
    auto face_nodes = RealMatrix(prhs[2], face_rows, face_cols, "quad_face_nodes");
    if (face_rows != 9u * static_cast<std::size_t>(n_bf) ||
        (n_bf != 0 && face_cols != 3))
        BadArgument("quad_face_nodes must have shape (9*n_bf)-by-3");

    auto charge_host = IntegerVector(prhs[5], "charge_host");
    auto charge_kind = IntegerVector(prhs[6], "charge_kind");
    std::size_t exponent_rows = 0, exponent_cols = 0;
    auto charge_expo = IntegerMatrix(
        prhs[7], exponent_rows, exponent_cols, "charge_expo");
    ValidateChargeDescriptors(charge_host, charge_kind, charge_expo,
                              exponent_rows, exponent_cols, n_el, n_bf);

    std::vector<double> sym_tet_pts, sym_tet_w, sym_tri_pts, sym_tri_w;
    std::vector<double> far_tet_pts, far_tet_w, far_tri_pts, far_tri_w;
    std::vector<double> gl_out, gw_out, gl_in, gw_in;
    ReadQuadrature(prhs[8], prhs[9], 3, true, "sym_tet_pts", "sym_tet_w",
                   sym_tet_pts, sym_tet_w);
    ReadQuadrature(prhs[10], prhs[11], 2, true, "sym_tri_pts", "sym_tri_w",
                   sym_tri_pts, sym_tri_w);
    ReadRule1D(prhs[12], prhs[13], "gl_out", "gw_out", gl_out, gw_out);
    ReadRule1D(prhs[14], prhs[15], "gl_in", "gw_in", gl_in, gw_in);
    ReadQuadrature(prhs[16], prhs[17], 3, true, "far_tet_pts", "far_tet_w",
                   far_tet_pts, far_tet_w);
    ReadQuadrature(prhs[18], prhs[19], 2, true, "far_tri_pts", "far_tri_w",
                   far_tri_pts, far_tri_w);
    const double near_grade = Scalar(prhs[20], "near_grade");
    const double far_inner_factor = Scalar(prhs[21], "far_inner_factor");
    if (near_grade <= 0.0 || far_inner_factor <= 0.0)
        BadArgument("near_grade and far_inner_factor must be positive");
    auto image_masks = IntegerVector(prhs[22], "image_masks");
    auto image_signs = RealVector(prhs[23], "image_signs");
    if (image_masks.size() != image_signs.size())
        BadArgument("image_masks and image_signs must have equal lengths");
    for (double sign : image_signs)
        if (!std::isfinite(sign))
            BadArgument("image_signs must contain finite values");

    auto holder = std::make_unique<ChargeGramHandle>();
    holder->n_dof = static_cast<int>(charge_host.size());
    holder->manager = std::make_unique<HACApKChargeGram>(
        std::move(cell_nodes), std::move(face_nodes), n_el, n_bf,
        std::move(charge_host), std::move(charge_kind), std::move(charge_expo),
        std::move(sym_tet_pts), std::move(sym_tet_w),
        std::move(sym_tri_pts), std::move(sym_tri_w),
        std::move(gl_out), std::move(gw_out), std::move(gl_in), std::move(gw_in),
        std::move(far_tet_pts), std::move(far_tet_w),
        std::move(far_tri_pts), std::move(far_tri_w), near_grade,
        far_inner_factor, std::move(image_masks), std::move(image_signs));
    plhs[0] = Uint64Output(RegisterChargeGram(std::move(holder)));
}

void ChargeGramCreateWedge(int nlhs, mxArray* plhs[], int nrhs,
                           const mxArray* prhs[]) {
    CheckArity(nrhs, 27, nlhs, 1,
        "handle = radia_mex('hacapk.charge_gram.create_wedge', wedge_cell_nodes, "
        "face_nodes, face_type, n_el, n_bf, charge_host, charge_kind, charge_expo, "
        "sym_tet_pts, sym_tet_w, sym_tri_pts, sym_tri_w, field_tri_pts, "
        "field_tri_w, gl_out, gw_out, gl_in, gw_in, far_tet_pts, far_tet_w, "
        "far_tri_pts, far_tri_w, near_grade, far_inner_factor, image_masks, "
        "image_signs)");
    const int n_el = PositiveInteger(prhs[4], "n_el");
    const int n_bf = NonnegativeInteger(prhs[5], "n_bf");
    std::size_t cell_rows = 0, cell_cols = 0;
    auto cell_nodes = RealMatrix(prhs[1], cell_rows, cell_cols, "wedge_cell_nodes");
    if (cell_rows != 18u * static_cast<std::size_t>(n_el) || cell_cols != 3)
        BadArgument("wedge_cell_nodes must have shape (18*n_el)-by-3");
    std::size_t face_rows = 0, face_cols = 0;
    auto face_nodes = RealMatrix(prhs[2], face_rows, face_cols, "face_nodes");
    if (face_rows != 9u * static_cast<std::size_t>(n_bf) ||
        (n_bf != 0 && face_cols != 3))
        BadArgument("face_nodes must have shape (9*n_bf)-by-3");
    auto face_type = IntegerVector(prhs[3], "face_type");
    if (face_type.size() != static_cast<std::size_t>(n_bf))
        BadArgument("face_type must have n_bf entries");
    for (int type : face_type)
        if (type != 0 && type != 1)
            BadArgument("face_type entries must be 0 (triangle) or 1 (quadrilateral)");

    auto charge_host = IntegerVector(prhs[6], "charge_host");
    auto charge_kind = IntegerVector(prhs[7], "charge_kind");
    std::size_t exponent_rows = 0, exponent_cols = 0;
    auto charge_expo = IntegerMatrix(
        prhs[8], exponent_rows, exponent_cols, "charge_expo");
    ValidateChargeDescriptors(charge_host, charge_kind, charge_expo,
                              exponent_rows, exponent_cols, n_el, n_bf);
    for (int exponent : charge_expo)
        if (exponent > 2)
            BadArgument("wedge charge_expo entries must lie in {0,1,2}");

    std::vector<double> sym_tet_pts, sym_tet_w, sym_tri_pts, sym_tri_w;
    std::vector<double> field_tri_pts, field_tri_w;
    std::vector<double> far_tet_pts, far_tet_w, far_tri_pts, far_tri_w;
    std::vector<double> gl_out, gw_out, gl_in, gw_in;
    ReadQuadrature(prhs[9], prhs[10], 3, true, "sym_tet_pts", "sym_tet_w",
                   sym_tet_pts, sym_tet_w);
    ReadQuadrature(prhs[11], prhs[12], 2, true, "sym_tri_pts", "sym_tri_w",
                   sym_tri_pts, sym_tri_w);
    ReadQuadrature(prhs[13], prhs[14], 2, true, "field_tri_pts", "field_tri_w",
                   field_tri_pts, field_tri_w);
    ReadRule1D(prhs[15], prhs[16], "gl_out", "gw_out", gl_out, gw_out);
    ReadRule1D(prhs[17], prhs[18], "gl_in", "gw_in", gl_in, gw_in);
    ReadQuadrature(prhs[19], prhs[20], 3, true, "far_tet_pts", "far_tet_w",
                   far_tet_pts, far_tet_w);
    ReadQuadrature(prhs[21], prhs[22], 2, true, "far_tri_pts", "far_tri_w",
                   far_tri_pts, far_tri_w);
    const double near_grade = Scalar(prhs[23], "near_grade");
    const double far_inner_factor = Scalar(prhs[24], "far_inner_factor");
    if (near_grade <= 0.0 || far_inner_factor <= 0.0)
        BadArgument("near_grade and far_inner_factor must be positive");
    auto image_masks = IntegerVector(prhs[25], "image_masks");
    auto image_signs = RealVector(prhs[26], "image_signs");
    if (image_masks.size() != image_signs.size())
        BadArgument("image_masks and image_signs must have equal lengths");
    for (double sign : image_signs)
        if (!std::isfinite(sign))
            BadArgument("image_signs must contain finite values");

    auto holder = std::make_unique<ChargeGramHandle>();
    holder->n_dof = static_cast<int>(charge_host.size());
    holder->manager = std::make_unique<HACApKChargeGram>(
        std::move(cell_nodes), std::move(face_nodes), std::move(face_type),
        n_el, n_bf, std::move(charge_host), std::move(charge_kind),
        std::move(charge_expo), std::move(sym_tet_pts), std::move(sym_tet_w),
        std::move(sym_tri_pts), std::move(sym_tri_w),
        std::move(field_tri_pts), std::move(field_tri_w),
        std::move(gl_out), std::move(gw_out), std::move(gl_in), std::move(gw_in),
        std::move(far_tet_pts), std::move(far_tet_w),
        std::move(far_tri_pts), std::move(far_tri_w), near_grade,
        far_inner_factor, std::move(image_masks), std::move(image_signs));
    plhs[0] = Uint64Output(RegisterChargeGram(std::move(holder)));
}

void ChargeGramCreatePlanar2D(int nlhs, mxArray* plhs[], int nrhs,
                              const mxArray* prhs[]) {
    CheckArity(nrhs, 25, nlhs, 1,
        "handle = radia_mex('hacapk.charge_gram.create_planar_2d', dim2, "
        "geometry_order, cell_map, cell_type, edge_map, n_el, n_be, "
        "charge_host, charge_kind, charge_expo, sym_tri_pts, sym_tri_w, "
        "gl_quad, gw_quad, gl_edge, gw_edge, gl_in, gw_in, far_tri_pts, "
        "far_tri_w, near_grade, far_inner_factor, image_masks, image_signs)");
    const int dim2 = PositiveInteger(prhs[1], "dim2");
    if (dim2 != 2)
        BadArgument("dim2 must equal 2 for the planar logarithmic kernel");
    const int geometry_order = PositiveInteger(prhs[2], "geometry_order");
    if (geometry_order > 3)
        BadArgument("geometry_order must lie in {1,2,3}");
    const int n_el = PositiveInteger(prhs[6], "n_el");
    const int n_be = NonnegativeInteger(prhs[7], "n_be");

    std::size_t cell_rows = 0, cell_cols = 0;
    auto cell_map = RealMatrix(prhs[3], cell_rows, cell_cols, "cell_map");
    const std::size_t cell_coefficients =
        static_cast<std::size_t>(geometry_order + 1) *
        static_cast<std::size_t>(geometry_order + 1);
    if (cell_rows != static_cast<std::size_t>(n_el) * cell_coefficients ||
        cell_cols != 2)
        BadArgument("cell_map must have shape (n_el*(geometry_order+1)^2)-by-2");
    auto cell_type = IntegerVector(prhs[4], "cell_type");
    if (cell_type.size() != static_cast<std::size_t>(n_el))
        BadArgument("cell_type must have n_el entries");
    for (int type : cell_type)
        if (type != 0 && type != 1)
            BadArgument("cell_type entries must be 0 (triangle) or 1 (quadrilateral)");

    std::size_t edge_rows = 0, edge_cols = 0;
    auto edge_map = RealMatrix(prhs[5], edge_rows, edge_cols, "edge_map");
    if (edge_rows != static_cast<std::size_t>(n_be) *
                         static_cast<std::size_t>(geometry_order + 1) ||
        (n_be != 0 && edge_cols != 2))
        BadArgument("edge_map must have shape (n_be*(geometry_order+1))-by-2");
    if (!std::all_of(cell_map.begin(), cell_map.end(),
                     [](double value) { return std::isfinite(value); }) ||
        !std::all_of(edge_map.begin(), edge_map.end(),
                     [](double value) { return std::isfinite(value); }))
        BadArgument("cell_map and edge_map must contain finite values");

    auto charge_host = IntegerVector(prhs[8], "charge_host");
    auto charge_kind = IntegerVector(prhs[9], "charge_kind");
    std::size_t exponent_rows = 0, exponent_cols = 0;
    auto charge_expo = IntegerMatrix(
        prhs[10], exponent_rows, exponent_cols, "charge_expo");
    ValidateChargeDescriptors(charge_host, charge_kind, charge_expo,
                              exponent_rows, exponent_cols, n_el, n_be);
    for (std::size_t charge = 0; charge < charge_host.size(); ++charge) {
        const int ex = charge_expo[3 * charge];
        const int ey = charge_expo[3 * charge + 1];
        const int ez = charge_expo[3 * charge + 2];
        if (ex > 3 || ey > 3 || ez != 0)
            BadArgument("planar charge exponents must satisfy ex,ey <= 3 and ez = 0");
        if (charge_kind[charge] == 1 && ey != 0)
            BadArgument("planar edge charges support only the first reference exponent");
    }

    std::vector<double> sym_tri_pts, sym_tri_w, far_tri_pts, far_tri_w;
    std::vector<double> gl_quad, gw_quad, gl_edge, gw_edge, gl_in, gw_in;
    ReadQuadrature(prhs[11], prhs[12], 2, true, "sym_tri_pts", "sym_tri_w",
                   sym_tri_pts, sym_tri_w);
    ReadRule1D(prhs[13], prhs[14], "gl_quad", "gw_quad", gl_quad, gw_quad);
    ReadRule1D(prhs[15], prhs[16], "gl_edge", "gw_edge", gl_edge, gw_edge);
    ReadRule1D(prhs[17], prhs[18], "gl_in", "gw_in", gl_in, gw_in);
    ReadQuadrature(prhs[19], prhs[20], 2, true, "far_tri_pts", "far_tri_w",
                   far_tri_pts, far_tri_w);
    const double near_grade = Scalar(prhs[21], "near_grade");
    const double far_inner_factor = Scalar(prhs[22], "far_inner_factor");
    if (!std::isfinite(near_grade) || !std::isfinite(far_inner_factor) ||
        near_grade <= 0.0 || far_inner_factor <= 0.0)
        BadArgument("near_grade and far_inner_factor must be positive and finite");
    auto image_masks = IntegerVector(prhs[23], "image_masks");
    auto image_signs = RealVector(prhs[24], "image_signs");
    if (image_masks.size() != image_signs.size())
        BadArgument("image_masks and image_signs must have equal lengths");
    for (int mask : image_masks)
        if (mask < 1 || mask > 3)
            BadArgument("planar image_masks entries must use x/y bits 1..3");
    for (double sign : image_signs)
        if (!std::isfinite(sign))
            BadArgument("image_signs must contain finite values");

    auto holder = std::make_unique<ChargeGramHandle>();
    holder->n_dof = static_cast<int>(charge_host.size());
    holder->manager = std::make_unique<HACApKChargeGram>(
        dim2, geometry_order, std::move(cell_map), std::move(cell_type),
        std::move(edge_map), n_el, n_be, std::move(charge_host),
        std::move(charge_kind), std::move(charge_expo),
        std::move(sym_tri_pts), std::move(sym_tri_w),
        std::move(gl_quad), std::move(gw_quad),
        std::move(gl_edge), std::move(gw_edge),
        std::move(gl_in), std::move(gw_in),
        std::move(far_tri_pts), std::move(far_tri_w), near_grade,
        far_inner_factor, std::move(image_masks), std::move(image_signs));
    plhs[0] = Uint64Output(RegisterChargeGram(std::move(holder)));
}

void ChargeGramCreateCurvedPolytope(int nlhs, mxArray* plhs[], int nrhs,
                                    const mxArray* prhs[]) {
    CheckArity(nrhs, 16, nlhs, 1,
        "handle = radia_mex('hacapk.charge_gram.create_curved_polytope', "
        "cell_curved_nodes, cell_subtet_off, cell_cent, cell_meas, "
        "face_curved_nodes, face_subtri_off, face_cent, face_meas, "
        "ref_tet_pts, ref_tet_w, ref_tri_pts, ref_tri_w, curve_gl, "
        "curve_gw, n_el)");
    const int n_el = PositiveInteger(prhs[15], "n_el");
    std::size_t cell_node_rows = 0, cell_node_cols = 0;
    auto cell_nodes = RealMatrix(
        prhs[1], cell_node_rows, cell_node_cols, "cell_curved_nodes");
    if (cell_node_rows == 0 || cell_node_cols != 3 || cell_node_rows % 10 != 0)
        BadArgument("cell_curved_nodes must have shape (10*n_subtet)-by-3");
    const int n_subtet = static_cast<int>(cell_node_rows / 10);
    auto cell_offsets = IntegerVector(prhs[2], "cell_subtet_off");
    if (cell_offsets.size() != static_cast<std::size_t>(n_el + 1) ||
        cell_offsets.front() != 0 || cell_offsets.back() != n_subtet)
        BadArgument("cell_subtet_off must be a zero-based n_el+1 CSR vector ending at n_subtet");
    for (std::size_t i = 1; i < cell_offsets.size(); ++i)
        if (cell_offsets[i] <= cell_offsets[i - 1])
            BadArgument("every curved polytope cell must contain at least one subtet");

    std::size_t cell_cent_rows = 0, cell_cent_cols = 0;
    auto cell_cent = RealMatrix(prhs[3], cell_cent_rows, cell_cent_cols, "cell_cent");
    auto cell_meas = RealVector(prhs[4], "cell_meas");
    if (cell_cent_rows != static_cast<std::size_t>(n_el) || cell_cent_cols != 3 ||
        cell_meas.size() != static_cast<std::size_t>(n_el))
        BadArgument("cell_cent must be n_el-by-3 with n_el cell_meas entries");

    std::size_t face_node_rows = 0, face_node_cols = 0;
    auto face_nodes = RealMatrix(
        prhs[5], face_node_rows, face_node_cols, "face_curved_nodes");
    if (face_node_rows != 0 && (face_node_cols != 3 || face_node_rows % 6 != 0))
        BadArgument("face_curved_nodes must have shape (6*n_subtri)-by-3");
    const int n_subtri = static_cast<int>(face_node_rows / 6);
    auto face_offsets = IntegerVector(prhs[6], "face_subtri_off");
    if (face_offsets.empty())
        BadArgument("face_subtri_off must contain at least the initial zero");
    const int n_bf = static_cast<int>(face_offsets.size() - 1);
    if (face_offsets.front() != 0 || face_offsets.back() != n_subtri)
        BadArgument("face_subtri_off must be a zero-based CSR vector ending at n_subtri");
    for (std::size_t i = 1; i < face_offsets.size(); ++i)
        if (face_offsets[i] <= face_offsets[i - 1])
            BadArgument("every curved polytope face must contain at least one subtriangle");

    std::size_t face_cent_rows = 0, face_cent_cols = 0;
    auto face_cent = RealMatrix(prhs[7], face_cent_rows, face_cent_cols, "face_cent");
    auto face_meas = RealVector(prhs[8], "face_meas");
    if (face_cent_rows != static_cast<std::size_t>(n_bf) ||
        (n_bf != 0 && face_cent_cols != 3) ||
        face_meas.size() != static_cast<std::size_t>(n_bf))
        BadArgument("face_cent must be n_face-by-3 with n_face face_meas entries");

    auto finite_geometry = [](const std::vector<double>& values) {
        return std::all_of(values.begin(), values.end(),
                           [](double value) { return std::isfinite(value); });
    };
    if (!finite_geometry(cell_nodes) || !finite_geometry(cell_cent) ||
        !finite_geometry(face_nodes) || !finite_geometry(face_cent))
        BadArgument("curved polytope nodes and centroids must contain finite values");
    for (double measure : cell_meas)
        if (!std::isfinite(measure) || measure <= 0.0)
            BadArgument("cell_meas entries must be positive and finite");
    for (double measure : face_meas)
        if (!std::isfinite(measure) || measure <= 0.0)
            BadArgument("face_meas entries must be positive and finite");

    std::vector<double> ref_tet_pts, ref_tet_w, ref_tri_pts, ref_tri_w;
    std::vector<double> curve_gl, curve_gw;
    ReadQuadrature(prhs[9], prhs[10], 3, true, "ref_tet_pts", "ref_tet_w",
                   ref_tet_pts, ref_tet_w);
    ReadQuadrature(prhs[11], prhs[12], 2, true, "ref_tri_pts", "ref_tri_w",
                   ref_tri_pts, ref_tri_w);
    ReadRule1D(prhs[13], prhs[14], "curve_gl", "curve_gw", curve_gl, curve_gw);

    auto holder = std::make_unique<ChargeGramHandle>();
    holder->n_dof = n_el + n_bf;
    holder->manager = std::make_unique<HACApKChargeGram>(
        std::move(cell_nodes), std::move(cell_offsets), std::move(cell_cent),
        std::move(cell_meas), std::move(face_nodes), std::move(face_offsets),
        std::move(face_cent), std::move(face_meas), std::move(ref_tet_pts),
        std::move(ref_tet_w), std::move(ref_tri_pts), std::move(ref_tri_w),
        std::move(curve_gl), std::move(curve_gw), n_el);
    plhs[0] = Uint64Output(RegisterChargeGram(std::move(holder)));
}

void ChargeGramCreateLocalPolynomials(int nlhs, mxArray* plhs[], int nrhs,
                                      const mxArray* prhs[]) {
    CheckArity(nrhs, 8, nlhs, 1,
        "handle = radia_mex('hacapk.charge_gram.create_local_polynomials', "
        "cell_verts, n_el, charge_host, polynomial_coefficients, "
        "polynomial_exponents, ref_tet_pts, ref_tet_w)");
    const int n_el = PositiveInteger(prhs[2], "n_el");
    std::size_t cell_rows = 0, cell_cols = 0;
    auto cell_verts = RealMatrix(prhs[1], cell_rows, cell_cols, "cell_verts");
    if (cell_rows != 4u * static_cast<std::size_t>(n_el) || cell_cols != 3)
        BadArgument("cell_verts must have shape (4*n_el)-by-3");
    auto charge_host = IntegerVector(prhs[3], "charge_host");
    std::size_t coefficient_rows = 0, coefficient_cols = 0;
    auto coefficients = RealMatrix(
        prhs[4], coefficient_rows, coefficient_cols, "polynomial_coefficients");
    std::size_t exponent_rows = 0, exponent_cols = 0;
    auto exponent_values = RealMatrix(
        prhs[5], exponent_rows, exponent_cols, "polynomial_exponents");
    std::size_t ref_rows = 0, ref_cols = 0;
    auto ref_points = RealMatrix(prhs[6], ref_rows, ref_cols, "ref_tet_pts");
    auto ref_weights = RealVector(prhs[7], "ref_tet_w");
    if (charge_host.size() != coefficient_rows)
        BadArgument("charge_host and polynomial_coefficients must have equal row counts");
    if (exponent_cols != 3 || exponent_rows != coefficient_cols)
        BadArgument("polynomial_exponents must have shape (n_monomial,3)");
    if (ref_cols != 3 || ref_rows != ref_weights.size() || ref_rows == 0)
        BadArgument("ref_tet_pts must be N-by-3 with one ref_tet_w per row");
    for (int host : charge_host)
        if (host < 0 || host >= n_el)
            BadArgument("charge_host uses zero-based indices in [0,n_el)");
    std::vector<int> exponents(exponent_values.size());
    for (std::size_t i = 0; i < exponent_values.size(); ++i) {
        const double value = exponent_values[i];
        if (value < 0.0 || value != std::floor(value) ||
            value > static_cast<double>(std::numeric_limits<int>::max()))
            BadArgument("polynomial_exponents must contain nonnegative integers");
        exponents[i] = static_cast<int>(value);
    }
    auto holder = std::make_unique<ChargeGramHandle>();
    holder->n_dof = static_cast<int>(coefficient_rows);
    holder->manager = std::make_unique<HACApKChargeGram>(
        std::move(cell_verts), n_el, std::move(charge_host),
        std::move(coefficients), std::move(exponents),
        std::move(ref_points), std::move(ref_weights));
    plhs[0] = Uint64Output(RegisterChargeGram(std::move(holder)));
}

void ChargeGramBuild(int nlhs, mxArray* plhs[], int nrhs,
                     const mxArray* prhs[]) {
    CheckArity(nrhs, 7, nlhs, 1,
        "ok = radia_mex('hacapk.charge_gram.build', handle, aca_eps, "
        "leaf_size, eta, max_rank, print_level)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    RadHACApKParams params;
    const double aca_eps = Scalar(prhs[2], "aca_eps");
    const int leaf_size = IntegerScalar(prhs[3], "leaf_size");
    const double eta = Scalar(prhs[4], "eta");
    const int max_rank = IntegerScalar(prhs[5], "max_rank");
    params.print_level = IntegerScalar(prhs[6], "print_level");
    if (aca_eps > 0.0) params.aca_eps = aca_eps;
    if (leaf_size > 0) params.leaf_size = leaf_size;
    if (eta > 0.0) params.eta = eta;
    if (max_rank > 0) params.max_rank = max_rank;
    plhs[0] = mxCreateLogicalScalar(holder.manager->BuildHMatrix(params));
}

void ChargeGramMatVec(const std::string& command, int nlhs, mxArray* plhs[],
                      int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
        "y = radia_mex('hacapk.charge_gram.<matvec>', handle, x)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const auto x = RealVector(prhs[2], "x");
    const int n = holder.n_dof;
    if (!holder.manager->IsValid())
        BadArgument("charge-Gram matvec requires build() first");
    if (n <= 0 || x.size() != static_cast<std::size_t>(n))
        BadArgument("x must have one entry per charge degree of freedom after build");
    std::vector<double> y(static_cast<std::size_t>(n));
    if (command == "hacapk.charge_gram.matvec")
        holder.manager->MatVec(x, y);
    else if (command == "hacapk.charge_gram.matvec_transpose")
        holder.manager->MatVecTranspose(x, y);
    else
        holder.manager->MatVecSym(x, y);
    plhs[0] = RealColumn(y);
}

void ChargeGramEntry(int nlhs, mxArray* plhs[], int nrhs,
                     const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "value = radia_mex('hacapk.charge_gram.entry', handle, i, j)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const int n = holder.n_dof;
    const int i = PositiveInteger(prhs[2], "i");
    const int j = PositiveInteger(prhs[3], "j");
    if (i > n || j > n)
        BadArgument("charge entry indices are out of range");
    plhs[0] = mxCreateDoubleScalar(
        holder.manager->GetInteractionMatrixElement(i - 1, j - 1));
}

void ChargeGramHexVolumeSelfBlockDirectionalDerivative(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "block = radia_mex('hacapk.charge_gram."
        "hex_volume_self_block_directional_derivative', handle, host, node_velocity)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const int host = PositiveInteger(prhs[2], "host");
    std::size_t rows = 0, cols = 0;
    auto velocity = RealMatrix(prhs[3], rows, cols, "node_velocity");
    if (rows != 27 || cols != 3)
        BadArgument("node_velocity must have shape 27-by-3");
    if (!std::all_of(velocity.begin(), velocity.end(),
                     [](double value) { return std::isfinite(value); }))
        BadArgument("node_velocity must contain finite values");
    auto derivative = holder.manager->HexVolumeSelfBlockDirectionalDerivative(
        host - 1, velocity);
    const std::size_t block_size = static_cast<std::size_t>(
        std::llround(std::sqrt(static_cast<double>(derivative.size()))));
    if (block_size * block_size != derivative.size())
        throw std::runtime_error("HEX self-block derivative returned a nonsquare block");
    plhs[0] = RealMatrixOutput(derivative, block_size, block_size);
}

void ChargeGramHexFaceSelfBlockDirectionalDerivative(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "block = radia_mex('hacapk.charge_gram."
        "hex_face_self_block_directional_derivative', handle, host, node_velocity)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const int host = PositiveInteger(prhs[2], "host");
    std::size_t rows = 0, cols = 0;
    auto velocity = RealMatrix(prhs[3], rows, cols, "node_velocity");
    if (rows != 9 || cols != 3)
        BadArgument("node_velocity must have shape 9-by-3");
    if (!std::all_of(velocity.begin(), velocity.end(),
                     [](double value) { return std::isfinite(value); }))
        BadArgument("node_velocity must contain finite values");
    auto derivative = holder.manager->HexFaceSelfBlockDirectionalDerivative(
        host - 1, velocity);
    const std::size_t block_size = static_cast<std::size_t>(
        std::llround(std::sqrt(static_cast<double>(derivative.size()))));
    if (block_size * block_size != derivative.size())
        throw std::runtime_error(
            "HEX face self-block derivative returned a nonsquare block");
    plhs[0] = RealMatrixOutput(derivative, block_size, block_size);
}

void ChargeGramTetVolumeSelfBlockDirectionalDerivative(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "block = radia_mex('hacapk.charge_gram."
        "tet_volume_self_block_directional_derivative', handle, host, node_velocity)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const int host = PositiveInteger(prhs[2], "host");
    std::size_t rows = 0, cols = 0;
    auto velocity = RealMatrix(prhs[3], rows, cols, "node_velocity");
    if (rows != 4 || cols != 3)
        BadArgument("node_velocity must have shape 4-by-3");
    if (!std::all_of(velocity.begin(), velocity.end(),
                     [](double value) { return std::isfinite(value); }))
        BadArgument("node_velocity must contain finite values");
    auto derivative = holder.manager->TetVolumeSelfBlockDirectionalDerivative(
        host - 1, velocity);
    const std::size_t block_size = static_cast<std::size_t>(
        std::llround(std::sqrt(static_cast<double>(derivative.size()))));
    if (block_size * block_size != derivative.size())
        throw std::runtime_error(
            "TET volume self-block derivative returned a nonsquare block");
    plhs[0] = RealMatrixOutput(derivative, block_size, block_size);
}

void ChargeGramTetFaceSelfBlockDirectionalDerivative(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "block = radia_mex('hacapk.charge_gram."
        "tet_face_self_block_directional_derivative', handle, host, node_velocity)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const int host = PositiveInteger(prhs[2], "host");
    std::size_t rows = 0, cols = 0;
    auto velocity = RealMatrix(prhs[3], rows, cols, "node_velocity");
    if (rows != 3 || cols != 3)
        BadArgument("node_velocity must have shape 3-by-3");
    if (!std::all_of(velocity.begin(), velocity.end(),
                     [](double value) { return std::isfinite(value); }))
        BadArgument("node_velocity must contain finite values");
    auto derivative = holder.manager->TetFaceSelfBlockDirectionalDerivative(
        host - 1, velocity);
    const std::size_t block_size = static_cast<std::size_t>(
        std::llround(std::sqrt(static_cast<double>(derivative.size()))));
    if (block_size * block_size != derivative.size())
        throw std::runtime_error(
            "TET face self-block derivative returned a nonsquare block");
    plhs[0] = RealMatrixOutput(derivative, block_size, block_size);
}

void ChargeGramWedgeVolumeSelfBlockDirectionalDerivative(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "block = radia_mex('hacapk.charge_gram."
        "wedge_volume_self_block_directional_derivative', handle, host, node_velocity)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const int host = PositiveInteger(prhs[2], "host");
    std::size_t rows = 0, cols = 0;
    auto velocity = RealMatrix(prhs[3], rows, cols, "node_velocity");
    if (rows != 18 || cols != 3)
        BadArgument("node_velocity must have shape 18-by-3");
    if (!std::all_of(velocity.begin(), velocity.end(),
                     [](double value) { return std::isfinite(value); }))
        BadArgument("node_velocity must contain finite values");
    auto derivative = holder.manager->WedgeVolumeSelfBlockDirectionalDerivative(
        host - 1, velocity);
    const std::size_t block_size = static_cast<std::size_t>(
        std::llround(std::sqrt(static_cast<double>(derivative.size()))));
    if (block_size * block_size != derivative.size())
        throw std::runtime_error(
            "WEDGE volume self-block derivative returned a nonsquare block");
    plhs[0] = RealMatrixOutput(derivative, block_size, block_size);
}

void ChargeGramWedgeFaceSelfBlockDirectionalDerivative(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "block = radia_mex('hacapk.charge_gram."
        "wedge_face_self_block_directional_derivative', handle, host, node_velocity)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const int host = PositiveInteger(prhs[2], "host");
    std::size_t rows = 0, cols = 0;
    auto velocity = RealMatrix(prhs[3], rows, cols, "node_velocity");
    if ((rows != 6 && rows != 9) || cols != 3)
        BadArgument("node_velocity must have shape 6-by-3 or 9-by-3");
    if (!std::all_of(velocity.begin(), velocity.end(),
                     [](double value) { return std::isfinite(value); }))
        BadArgument("node_velocity must contain finite values");
    auto derivative = holder.manager->WedgeFaceSelfBlockDirectionalDerivative(
        host - 1, velocity);
    const std::size_t block_size = static_cast<std::size_t>(
        std::llround(std::sqrt(static_cast<double>(derivative.size()))));
    if (block_size * block_size != derivative.size())
        throw std::runtime_error(
            "WEDGE face self-block derivative returned a nonsquare block");
    plhs[0] = RealMatrixOutput(derivative, block_size, block_size);
}

void ChargeGramHexDirectionalDerivative(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "matrix = radia_mex('hacapk.charge_gram.hex_directional_derivative', "
        "handle, cell_node_velocity, face_node_velocity)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    std::size_t cell_count = 0, cell_nodes = 0, cell_components = 0;
    std::size_t face_count = 0, face_nodes = 0, face_components = 0;
    auto cell_velocity = RealTensor3(
        prhs[2], cell_count, cell_nodes, cell_components, "cell_node_velocity");
    auto face_velocity = RealTensor3(
        prhs[3], face_count, face_nodes, face_components, "face_node_velocity");
    if (cell_nodes != 27 || cell_components != 3)
        BadArgument("cell_node_velocity must have shape ncell-by-27-by-3");
    if (face_nodes != 9 || face_components != 3)
        BadArgument("face_node_velocity must have shape nface-by-9-by-3");
    const auto finite = [](double value) { return std::isfinite(value); };
    if (!std::all_of(cell_velocity.begin(), cell_velocity.end(), finite) ||
        !std::all_of(face_velocity.begin(), face_velocity.end(), finite))
        BadArgument("charge-Gram velocities must contain finite values");
    auto derivative = holder.manager->HexChargeGramDirectionalDerivative(
        cell_velocity, face_velocity);
    const std::size_t n = static_cast<std::size_t>(
        std::llround(std::sqrt(static_cast<double>(derivative.size()))));
    if (n * n != derivative.size() ||
        !std::all_of(derivative.begin(), derivative.end(), finite))
        throw std::runtime_error(
            "HEX charge-Gram derivative returned an invalid matrix");
    plhs[0] = RealMatrixOutput(derivative, n, n);
}

void ChargeGramTetDirectionalDerivative(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "matrix = radia_mex('hacapk.charge_gram.tet_directional_derivative', "
        "handle, cell_vertex_velocity, face_vertex_velocity)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    std::size_t cell_count = 0, cell_nodes = 0, cell_components = 0;
    std::size_t face_count = 0, face_nodes = 0, face_components = 0;
    auto cell_velocity = RealTensor3(
        prhs[2], cell_count, cell_nodes, cell_components, "cell_vertex_velocity");
    auto face_velocity = RealTensor3(
        prhs[3], face_count, face_nodes, face_components, "face_vertex_velocity");
    if (cell_nodes != 4 || cell_components != 3)
        BadArgument("cell_vertex_velocity must have shape ncell-by-4-by-3");
    if (face_nodes != 3 || face_components != 3)
        BadArgument("face_vertex_velocity must have shape nface-by-3-by-3");
    const auto finite = [](double value) { return std::isfinite(value); };
    if (!std::all_of(cell_velocity.begin(), cell_velocity.end(), finite) ||
        !std::all_of(face_velocity.begin(), face_velocity.end(), finite))
        BadArgument("charge-Gram velocities must contain finite values");
    auto derivative = holder.manager->TetChargeGramDirectionalDerivative(
        cell_velocity, face_velocity);
    const std::size_t n = static_cast<std::size_t>(
        std::llround(std::sqrt(static_cast<double>(derivative.size()))));
    if (n * n != derivative.size() ||
        !std::all_of(derivative.begin(), derivative.end(), finite))
        throw std::runtime_error(
            "TET charge-Gram derivative returned an invalid matrix");
    plhs[0] = RealMatrixOutput(derivative, n, n);
}

void ChargeGramTetChargeMapRowDirectionalRates(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "rates = radia_mex('hacapk.charge_gram."
        "tet_charge_map_row_directional_rates', handle, "
        "cell_vertex_velocity, face_vertex_velocity)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    std::size_t cell_count = 0, cell_nodes = 0, cell_components = 0;
    std::size_t face_count = 0, face_nodes = 0, face_components = 0;
    auto cell_velocity = RealTensor3(
        prhs[2], cell_count, cell_nodes, cell_components, "cell_vertex_velocity");
    auto face_velocity = RealTensor3(
        prhs[3], face_count, face_nodes, face_components, "face_vertex_velocity");
    if (cell_nodes != 4 || cell_components != 3)
        BadArgument("cell_vertex_velocity must have shape ncell-by-4-by-3");
    if (face_nodes != 3 || face_components != 3)
        BadArgument("face_vertex_velocity must have shape nface-by-3-by-3");
    const auto finite = [](double value) { return std::isfinite(value); };
    if (!std::all_of(cell_velocity.begin(), cell_velocity.end(), finite) ||
        !std::all_of(face_velocity.begin(), face_velocity.end(), finite))
        BadArgument("charge-map velocities must contain finite values");
    auto rates = holder.manager->TetChargeMapRowDirectionalRates(
        cell_velocity, face_velocity);
    if (!std::all_of(rates.begin(), rates.end(), finite))
        throw std::runtime_error("TET charge-map derivative returned invalid rates");
    plhs[0] = RealColumn(rates);
}

void ChargeGramWedgeDirectionalDerivative(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "matrix = radia_mex('hacapk.charge_gram.wedge_directional_derivative', "
        "handle, cell_node_velocity, face_node_velocity)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    std::size_t cell_count = 0, cell_nodes = 0, cell_components = 0;
    std::size_t face_count = 0, face_nodes = 0, face_components = 0;
    auto cell_velocity = RealTensor3(
        prhs[2], cell_count, cell_nodes, cell_components, "cell_node_velocity");
    auto face_velocity = RealTensor3(
        prhs[3], face_count, face_nodes, face_components, "face_node_velocity");
    if (cell_nodes != 18 || cell_components != 3)
        BadArgument("cell_node_velocity must have shape ncell-by-18-by-3");
    if (face_nodes != 9 || face_components != 3)
        BadArgument("face_node_velocity must have shape nface-by-9-by-3");
    const auto finite = [](double value) { return std::isfinite(value); };
    if (!std::all_of(cell_velocity.begin(), cell_velocity.end(), finite) ||
        !std::all_of(face_velocity.begin(), face_velocity.end(), finite))
        BadArgument("charge-Gram velocities must contain finite values");
    auto derivative = holder.manager->WedgeChargeGramDirectionalDerivative(
        cell_velocity, face_velocity);
    const std::size_t n = static_cast<std::size_t>(
        std::llround(std::sqrt(static_cast<double>(derivative.size()))));
    if (n * n != derivative.size() ||
        !std::all_of(derivative.begin(), derivative.end(), finite))
        throw std::runtime_error(
            "WEDGE charge-Gram derivative returned an invalid matrix");
    plhs[0] = RealMatrixOutput(derivative, n, n);
}

void ChargeGramDirectionalDerivativeOperatorCreate(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 8, nlhs, 1,
        "handle = radia_mex('hacapk.charge_gram.directional_derivative_operator', "
        "parent, family, cell_velocity, face_velocity, aca_eps, leaf_size, eta)");
    ChargeGramHandle& parent = ChargeGram(Handle(prhs[1]));
    const std::string family_name = Text(prhs[2], "family");
    ChargeDerivativeFamily family;
    std::size_t expected_cell_nodes = 0;
    std::size_t expected_face_nodes = 0;
    if (family_name == "hex") {
        family = ChargeDerivativeFamily::Hex;
        expected_cell_nodes = 27;
        expected_face_nodes = 9;
    } else if (family_name == "tet") {
        family = ChargeDerivativeFamily::Tet;
        expected_cell_nodes = 4;
        expected_face_nodes = 3;
    } else if (family_name == "wedge") {
        family = ChargeDerivativeFamily::Wedge;
        expected_cell_nodes = 18;
        expected_face_nodes = 9;
    } else {
        BadArgument("family must be 'hex', 'tet', or 'wedge'");
    }
    std::size_t cell_count = 0, cell_nodes = 0, cell_components = 0;
    std::size_t face_count = 0, face_nodes = 0, face_components = 0;
    auto cell_velocity = RealTensor3(
        prhs[3], cell_count, cell_nodes, cell_components, "cell_velocity");
    auto face_velocity = RealTensor3(
        prhs[4], face_count, face_nodes, face_components, "face_velocity");
    if (cell_nodes != expected_cell_nodes || cell_components != 3)
        BadArgument("cell_velocity has the wrong node count for the selected family");
    if (face_nodes != expected_face_nodes || face_components != 3)
        BadArgument("face_velocity has the wrong node count for the selected family");
    const auto finite = [](double value) { return std::isfinite(value); };
    if (!std::all_of(cell_velocity.begin(), cell_velocity.end(), finite) ||
        !std::all_of(face_velocity.begin(), face_velocity.end(), finite))
        BadArgument("charge-Gram velocities must contain finite values");
    RadHACApKParams params;
    const double aca_eps = Scalar(prhs[5], "aca_eps");
    const int leaf_size = IntegerScalar(prhs[6], "leaf_size");
    const double eta = Scalar(prhs[7], "eta");
    if (!(aca_eps > 0.0) || leaf_size <= 0 || !(eta > 0.0))
        BadArgument("aca_eps, leaf_size, and eta must be positive");
    params.aca_eps = aca_eps;
    params.leaf_size = leaf_size;
    params.eta = eta;
    auto holder = std::make_unique<ChargeGramDerivativeHandle>();
    holder->parent = parent.manager;
    holder->manager = parent.manager->BuildDirectionalDerivativeOperator(
        family, cell_velocity, face_velocity, params);
    holder->n_dof = holder->manager->GetNDOF();
    plhs[0] = Uint64Output(RegisterChargeGramDerivative(std::move(holder)));
}

void ChargeGramDirectionalDerivativeContractions(
    int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 7, nlhs, 1,
        "values = radia_mex('hacapk.charge_gram."
        "directional_derivative_contractions', parent, family, cell_velocity, "
        "face_velocity, left, right)");
    ChargeGramHandle& parent = ChargeGram(Handle(prhs[1]));
    const std::string family_name = Text(prhs[2], "family");
    ChargeDerivativeFamily family;
    std::size_t expected_cell_nodes = 0;
    std::size_t expected_face_nodes = 0;
    if (family_name == "hex") {
        family = ChargeDerivativeFamily::Hex;
        expected_cell_nodes = 27;
        expected_face_nodes = 9;
    } else if (family_name == "tet") {
        family = ChargeDerivativeFamily::Tet;
        expected_cell_nodes = 4;
        expected_face_nodes = 3;
    } else if (family_name == "wedge") {
        family = ChargeDerivativeFamily::Wedge;
        expected_cell_nodes = 18;
        expected_face_nodes = 9;
    } else {
        BadArgument("family must be 'hex', 'tet', or 'wedge'");
    }

    std::size_t n_mode = 0, cell_count = 0, cell_nodes = 0;
    std::size_t cell_components = 0, face_modes = 0, face_count = 0;
    std::size_t face_nodes = 0, face_components = 0;
    auto cell_velocity = RealTensor4(prhs[3], n_mode, cell_count, cell_nodes,
                                     cell_components, "cell_velocity");
    auto face_velocity = RealTensor4(prhs[4], face_modes, face_count, face_nodes,
                                     face_components, "face_velocity");
    if (n_mode == 0 || n_mode != face_modes)
        BadArgument("cell_velocity and face_velocity must have the same positive mode count");
    if (cell_nodes != expected_cell_nodes || cell_components != 3)
        BadArgument("cell_velocity has the wrong node count for the selected family");
    if (face_nodes != expected_face_nodes || face_components != 3)
        BadArgument("face_velocity has the wrong node count for the selected family");
    if (n_mode > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("the derivative mode count exceeds the native integer range");
    auto left = RealVector(prhs[5], "left");
    auto right = RealVector(prhs[6], "right");
    if (left.size() != static_cast<std::size_t>(parent.n_dof) ||
        right.size() != static_cast<std::size_t>(parent.n_dof))
        BadArgument("left and right must contain one entry per charge degree of freedom");
    const auto finite = [](double value) { return std::isfinite(value); };
    if (!std::all_of(cell_velocity.begin(), cell_velocity.end(), finite) ||
        !std::all_of(face_velocity.begin(), face_velocity.end(), finite) ||
        !std::all_of(left.begin(), left.end(), finite) ||
        !std::all_of(right.begin(), right.end(), finite))
        BadArgument("derivative velocities and contraction vectors must be finite");
    auto values = parent.manager->DirectionalDerivativeContractions(
        family, static_cast<int>(n_mode), cell_velocity, face_velocity,
        left, right);
    if (!std::all_of(values.begin(), values.end(), finite))
        throw std::runtime_error("charge-Gram derivative contractions returned invalid values");
    plhs[0] = RealColumn(values);
}

void ChargeGramDerivativeInfo(int nlhs, mxArray* plhs[], int nrhs,
                              const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
        "info = radia_mex('hacapk.charge_gram_derivative.info', handle)");
    ChargeGramDerivativeHandle& holder = ChargeGramDerivative(Handle(prhs[1]));
    const auto& stats = holder.manager->GetStats();
    const char* fields[] = {
        "n_dof", "valid", "n_lowrank", "n_dense", "max_rank", "n_leaves",
        "compression", "build_time", "memory_mb", "dense_memory_mb"};
    plhs[0] = mxCreateStructMatrix(1, 1, 10, fields);
    mxSetField(plhs[0], 0, "n_dof", mxCreateDoubleScalar(holder.n_dof));
    mxSetField(plhs[0], 0, "valid",
               mxCreateLogicalScalar(holder.manager->IsValid()));
    mxSetField(plhs[0], 0, "n_lowrank", mxCreateDoubleScalar(stats.n_lowrank));
    mxSetField(plhs[0], 0, "n_dense", mxCreateDoubleScalar(stats.n_dense));
    mxSetField(plhs[0], 0, "max_rank", mxCreateDoubleScalar(stats.max_rank));
    mxSetField(plhs[0], 0, "n_leaves", mxCreateDoubleScalar(stats.n_leaves));
    mxSetField(plhs[0], 0, "compression", mxCreateDoubleScalar(stats.compression));
    mxSetField(plhs[0], 0, "build_time", mxCreateDoubleScalar(stats.build_time));
    mxSetField(plhs[0], 0, "memory_mb", mxCreateDoubleScalar(stats.memory_mb));
    mxSetField(plhs[0], 0, "dense_memory_mb",
               mxCreateDoubleScalar(stats.dense_memory_mb));
}

void ChargeGramDerivativeEntry(int nlhs, mxArray* plhs[], int nrhs,
                               const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "value = radia_mex('hacapk.charge_gram_derivative.entry', handle, i, j)");
    ChargeGramDerivativeHandle& holder = ChargeGramDerivative(Handle(prhs[1]));
    const int i = PositiveInteger(prhs[2], "i");
    const int j = PositiveInteger(prhs[3], "j");
    if (i > holder.n_dof || j > holder.n_dof)
        BadArgument("charge-Gram derivative entry indices are out of range");
    plhs[0] = mxCreateDoubleScalar(
        holder.manager->GetInteractionMatrixElement(i - 1, j - 1));
}

void ChargeGramDerivativeMatVecSym(int nlhs, mxArray* plhs[], int nrhs,
                                   const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
        "y = radia_mex('hacapk.charge_gram_derivative.matvec_sym', handle, x)");
    ChargeGramDerivativeHandle& holder = ChargeGramDerivative(Handle(prhs[1]));
    const auto x = RealVector(prhs[2], "x");
    if (x.size() != static_cast<std::size_t>(holder.n_dof))
        BadArgument("x must have one entry per charge degree of freedom");
    std::vector<double> y(static_cast<std::size_t>(holder.n_dof));
    holder.manager->MatVecSym(x, y);
    plhs[0] = RealColumn(y);
}

void ChargeGramInfo(int nlhs, mxArray* plhs[], int nrhs,
                    const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
        "info = radia_mex('hacapk.charge_gram.info', handle)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const auto& stats = holder.manager->GetStats();
    const char* fields[] = {
        "n_dof", "valid", "n_lowrank", "n_dense", "max_rank", "n_leaves",
        "compression", "build_time", "memory_mb", "dense_memory_mb"};
    plhs[0] = mxCreateStructMatrix(1, 1, 10, fields);
    mxSetField(plhs[0], 0, "n_dof",
               mxCreateDoubleScalar(holder.n_dof));
    mxSetField(plhs[0], 0, "valid",
               mxCreateLogicalScalar(holder.manager->IsValid()));
    mxSetField(plhs[0], 0, "n_lowrank", mxCreateDoubleScalar(stats.n_lowrank));
    mxSetField(plhs[0], 0, "n_dense", mxCreateDoubleScalar(stats.n_dense));
    mxSetField(plhs[0], 0, "max_rank", mxCreateDoubleScalar(stats.max_rank));
    mxSetField(plhs[0], 0, "n_leaves", mxCreateDoubleScalar(stats.n_leaves));
    mxSetField(plhs[0], 0, "compression", mxCreateDoubleScalar(stats.compression));
    mxSetField(plhs[0], 0, "build_time", mxCreateDoubleScalar(stats.build_time));
    mxSetField(plhs[0], 0, "memory_mb", mxCreateDoubleScalar(stats.memory_mb));
    mxSetField(plhs[0], 0, "dense_memory_mb",
               mxCreateDoubleScalar(stats.dense_memory_mb));
}


void ChargeGramHexStateCheck(int nlhs, mxArray* plhs[], int nrhs,
                             const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "check = radia_mex('hacapk.charge_gram.hex_state_check', handle)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const char* fields[] = {"ctor", "now"};
    plhs[0] = mxCreateStructMatrix(1, 1, 2, fields);
    mxSetField(plhs[0], 0, "ctor", mxCreateDoubleScalar(holder.manager->HexStateCtorChecksum()));
    mxSetField(plhs[0], 0, "now", mxCreateDoubleScalar(holder.manager->HexStateChecksum()));
}

void ChargeGramHexStoredNodes(int nlhs, mxArray* plhs[], int nrhs,
                              const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "nodes = radia_mex('hacapk.charge_gram.hex_stored_nodes', handle)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const char* fields[] = {"cell_nodes", "face_nodes"};
    plhs[0] = mxCreateStructMatrix(1, 1, 2, fields);
    mxSetField(plhs[0], 0, "cell_nodes", RealColumn(holder.manager->HexStoredCellNodes()));
    mxSetField(plhs[0], 0, "face_nodes", RealColumn(holder.manager->HexStoredFaceNodes()));
}

void ChargeGramHexStateBreakdown(int nlhs, mxArray* plhs[], int nrhs,
                                 const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "breakdown = radia_mex('hacapk.charge_gram.hex_state_breakdown', handle)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    plhs[0] = PairStructOutput(holder.manager->HexStateBreakdown());
}

void ChargeGramConfigureChargeMap(int nlhs, mxArray* plhs[], int nrhs,
                                  const mxArray* prhs[]) {
    CheckArity(nrhs, 6, nlhs, 0,
               "radia_mex('hacapk.charge_gram.configure_charge_map', handle, B_indptr, B_indices, B_data, n_face)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    holder.manager->ConfigureChargeMap(
        IntegerVector(prhs[2], "B_indptr"), IntegerVector(prhs[3], "B_indices"),
        RealVector(prhs[4], "B_data"), NonnegativeInteger(prhs[5], "n_face"));
}

void ChargeGramConfigureVectorChargeMap(int nlhs, mxArray* plhs[], int nrhs,
                                        const mxArray* prhs[]) {
    CheckArity(nrhs, 7, nlhs, 0,
               "radia_mex('hacapk.charge_gram.configure_vector_charge_map', handle, B_indptr, B_indices, B_data, n_face, n_components)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    holder.manager->ConfigureVectorChargeMap(
        IntegerVector(prhs[2], "B_indptr"), IntegerVector(prhs[3], "B_indices"),
        RealVector(prhs[4], "B_data"), NonnegativeInteger(prhs[5], "n_face"),
        PositiveInteger(prhs[6], "n_components"));
}

void ChargeGramConfigureMass(const std::string& command, int nlhs,
                             mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 6, nlhs, 0,
               "radia_mex('hacapk.charge_gram.configure_<mass>', handle, rows, cols, values, n_face)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    auto rows = IntegerVector(prhs[2], "mI");
    auto cols = IntegerVector(prhs[3], "mJ");
    auto values = RealVector(prhs[4], "mV");
    if (rows.size() != cols.size() || rows.size() != values.size())
        BadArgument("COO rows, cols, and values must have the same length");
    const int n_face = NonnegativeInteger(prhs[5], "n_face");
    if (command == "hacapk.charge_gram.configure_mass_matrix")
        holder.manager->ConfigureMassMatrix(std::move(rows), std::move(cols),
                                            std::move(values), n_face);
    else
        holder.manager->ConfigureGeometryMassMatrix(std::move(rows), std::move(cols),
                                                    std::move(values), n_face);
}

void ChargeGramConfigureMassNGSolve(const std::string& command, int nlhs,
                                    mxArray* plhs[], int nrhs,
                                    const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 0,
               "radia_mex('hacapk.charge_gram.configure_<mass>_ngsolve', charge_handle, matrix_handle)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const auto& matrix = Matrix(Handle(prhs[2]));
    auto mass = ExtractNGSolveScalarSparse(matrix.matrix, command.c_str());
    if (command == "hacapk.charge_gram.configure_mass_matrix_ngsolve")
        holder.manager->ConfigureMassMatrix(
            std::move(mass.rows), std::move(mass.cols),
            std::move(mass.values), mass.size);
    else
        holder.manager->ConfigureGeometryMassMatrix(
            std::move(mass.rows), std::move(mass.cols),
            std::move(mass.values), mass.size);
}

void ChargeGramRestoreGeometryMassMatrix(int nlhs, mxArray* plhs[], int nrhs,
                                         const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "changed = radia_mex('hacapk.charge_gram.restore_geometry_mass_matrix', handle)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    plhs[0] = mxCreateLogicalScalar(holder.manager->RestoreGeometryMassMatrix());
}

void ChargeGramOperatorInfo(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "info = radia_mex('hacapk.charge_gram.operator_info', handle)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const char* fields[] = {"operator_configured", "charge_map_configured",
                            "mass_matrix_configured", "geometry_mass_configured",
                            "n_face", "constraint_count"};
    plhs[0] = mxCreateStructMatrix(1, 1, 6, fields);
    mxSetField(plhs[0], 0, "operator_configured",
               mxCreateLogicalScalar(holder.manager->HasConfiguredChargeMap() &&
                                     holder.manager->HasConfiguredMassMatrix() &&
                                     holder.manager->HasConfiguredGeometryMassMatrix()));
    mxSetField(plhs[0], 0, "charge_map_configured",
               mxCreateLogicalScalar(holder.manager->HasConfiguredChargeMap()));
    mxSetField(plhs[0], 0, "mass_matrix_configured",
               mxCreateLogicalScalar(holder.manager->HasConfiguredMassMatrix()));
    mxSetField(plhs[0], 0, "geometry_mass_configured",
               mxCreateLogicalScalar(holder.manager->HasConfiguredGeometryMassMatrix()));
    mxSetField(plhs[0], 0, "n_face",
               mxCreateDoubleScalar(holder.manager->ConfiguredNFace()));
    mxSetField(plhs[0], 0, "constraint_count",
               mxCreateDoubleScalar(holder.manager->ConfiguredConstraintCount()));
}

void ChargeGramDemagMatrix(int nlhs, mxArray* plhs[], int nrhs,
                           const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
        "matrix = radia_mex('hacapk.charge_gram.demag_matrix', handle)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    auto matrix = std::make_shared<radia::ngsolve_bridge::HDivDemagMatrix>(
        holder.manager);
    plhs[0] = Uint64Output(RegisterMatrix(MakeNGSolveMatrixHandle(
        std::move(matrix), nullptr, "hdiv_demag_matrix")));
}

void ChargeGramConfiguredApply(const std::string& command, int nlhs,
                               mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (command == "hacapk.charge_gram.demag_apply") {
        if ((nrhs != 3 && nrhs != 4) || nlhs != 1)
            BadArgument("usage: y = radia_mex('hacapk.charge_gram.demag_apply', handle, x, symmetric)");
    } else {
        CheckArity(nrhs, 3, nlhs, 1,
                   "y = radia_mex('hacapk.charge_gram.<configured_apply>', handle, x)");
    }
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const int n = holder.manager->ConfiguredNFace();
    auto x = RealVector(prhs[2], "x");
    if (n < 0 || x.size() != static_cast<std::size_t>(n))
        BadArgument("x must have one entry per configured face degree of freedom");
    std::vector<double> y(static_cast<std::size_t>(n), 0.0);
    if (command == "hacapk.charge_gram.demag_apply") {
        const bool symmetric = nrhs == 4 ? Boolean(prhs[3], "symmetric") : true;
        holder.manager->ApplyConfiguredDemag(x.data(), y.data(), symmetric);
    } else if (command == "hacapk.charge_gram.geometry_mass_apply") {
        holder.manager->ApplyConfiguredGeometryMass(x.data(), y.data());
    } else {
        y = holder.manager->ApplyConfiguredMassRiesz(x);
    }
    plhs[0] = RealColumn(y);
}

void ChargeGramSolveConfigured(const std::string& command, int nlhs,
                               mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nlhs != 1 || nrhs < 4 || nrhs > 8)
        BadArgument("usage: result = radia_mex('hacapk.charge_gram.solve_configured_linear_material[_auto_prec]', handle, inv_chi, rhs, tol, maxit, symmetric, x0)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const double inv_chi = Scalar(prhs[2], "inv_chi");
    auto rhs = RealVector(prhs[3], "rhs");
    const double tol = nrhs >= 5 ? Scalar(prhs[4], "tol") : 1.0e-8;
    const int maxit = nrhs >= 6 ? PositiveInteger(prhs[5], "maxit") : 5000;
    const bool symmetric = nrhs >= 7 ? Boolean(prhs[6], "symmetric") : true;
    std::vector<double> x0;
    const std::vector<double>* x0_ptr = nullptr;
    if (nrhs >= 8 && mxGetNumberOfElements(prhs[7]) != 0) {
        x0 = RealVector(prhs[7], "x0");
        x0_ptr = &x0;
    }
    int iters = 0;
    double prec_min = 0.0, prec_max = 0.0;
    std::vector<double> result;
    if (command == "hacapk.charge_gram.solve_configured_linear_material")
        result = holder.manager->SolveConfiguredLinearMaterial(
            inv_chi, rhs, tol, maxit, iters, true, symmetric, x0_ptr);
    else
        result = holder.manager->SolveConfiguredLinearMaterialAutoPrec(
            inv_chi, rhs, tol, maxit, iters, prec_min, prec_max, x0_ptr);
    const char* fields[] = {"m", "iters", "prec_min", "prec_max", "timings"};
    plhs[0] = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(plhs[0], 0, "m", RealColumn(result));
    mxSetField(plhs[0], 0, "iters", mxCreateDoubleScalar(iters));
    mxSetField(plhs[0], 0, "prec_min", mxCreateDoubleScalar(prec_min));
    mxSetField(plhs[0], 0, "prec_max", mxCreateDoubleScalar(prec_max));
    mxSetField(plhs[0], 0, "timings", PairStructOutput(holder.manager->LastSolveTimings()));
}

void ChargeGramCreateFieldEvaluator(int nlhs, mxArray* plhs[], int nrhs,
                                    const mxArray* prhs[]) {
    if (nlhs != 1 || (nrhs != 3 && nrhs != 9))
        BadArgument("usage: handle = radia_mex('hacapk.charge_gram.create_field_evaluator', handle, magnetization [, leaf_size, theta, tree_min_sources, auto_min_work, tree_relative_tolerance, probe_count])");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    auto magnetization = RealVector(prhs[2], "magnetization");
    const auto options = FieldOptions(nrhs, prhs, 3,
        "handle = radia_mex('hacapk.charge_gram.create_field_evaluator', handle, magnetization [, options])");
    plhs[0] = Uint64Output(RegisterField(
        holder.manager->CreateConfiguredFieldEvaluator(magnetization, options)));
}

void ChargeGramCreatePlanarFieldEvaluator(int nlhs, mxArray* plhs[], int nrhs,
                                          const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
               "handle = radia_mex('hacapk.charge_gram.create_planar_field_evaluator', handle, magnetization)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    auto magnetization = RealVector(prhs[2], "magnetization");
    plhs[0] = Uint64Output(RegisterPlanar(
        holder.manager->CreateConfiguredPlanarFieldEvaluator(magnetization)));
}

void ChargeGramStats(int nlhs, mxArray* plhs[], int nrhs,
                     const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "stats = radia_mex('hacapk.charge_gram.stats', handle)");
    ChargeGramHandle& holder = ChargeGram(Handle(prhs[1]));
    const auto& stats = holder.manager->GetStats();
    std::vector<std::pair<std::string, double>> values = {
        {"n_dof", static_cast<double>(stats.n_dof)},
        {"n_leaves", static_cast<double>(stats.n_leaves)},
        {"n_lowrank", static_cast<double>(stats.n_lowrank)},
        {"n_dense", static_cast<double>(stats.n_dense)},
        {"max_rank", static_cast<double>(stats.max_rank)},
        {"compression", stats.compression}, {"build_time", stats.build_time},
        {"memory_mb", stats.memory_mb}, {"dense_memory_mb", stats.dense_memory_mb}};
    const auto cache = holder.manager->HexCacheStats();
    values.insert(values.end(), cache.begin(), cache.end());
    plhs[0] = PairStructOutput(values);
}

void HDivFieldFromTet(int nlhs, mxArray* plhs[], int nrhs,
                      const mxArray* prhs[]) {
    if (nlhs != 1 && nlhs != 2)
        BadArgument("usage: handle = radia_mex('hdiv.field_evaluator.from_tet', volume, surface, image_masks, image_signs [, options])");
    if (nrhs != 5 && nrhs != 11)
        BadArgument("usage: handle = radia_mex('hdiv.field_evaluator.from_tet', volume, surface, image_masks, image_signs [, leaf_size, theta, tree_min_sources, auto_min_work, tree_relative_tolerance, probe_count])");
    std::size_t volume_rows = 0, volume_cols = 0;
    std::size_t surface_rows = 0, surface_cols = 0;
    auto volume = RealMatrix(prhs[1], volume_rows, volume_cols, "volume");
    auto surface = RealMatrix(prhs[2], surface_rows, surface_cols, "surface");
    if (volume_cols != 16 || surface_cols != 22)
        BadArgument("volume must be N-by-16 and surface must be M-by-22");
    auto image_masks = IntegerVector(prhs[3], "image_masks");
    auto image_signs = RealVector(prhs[4], "image_signs");
    const auto options = FieldOptions(nrhs, prhs, 5,
        "handle = radia_mex('hdiv.field_evaluator.from_tet', volume, surface, image_masks, image_signs [, options])");
    auto evaluator = HDivFieldEvaluator::FromTet(
        std::move(volume), std::move(surface), std::move(image_masks),
        std::move(image_signs), options);
    plhs[0] = Uint64Output(RegisterField(std::move(evaluator)));
    if (nlhs == 2) {
        const char* fields[] = {"volume_count", "surface_count"};
        plhs[1] = mxCreateStructMatrix(1, 1, 2, fields);
        mxSetField(plhs[1], 0, "volume_count", mxCreateDoubleScalar(volume_rows));
        mxSetField(plhs[1], 0, "surface_count", mxCreateDoubleScalar(surface_rows));
    }
}

void HDivFieldFromCloud(int nlhs, mxArray* plhs[], int nrhs,
                        const mxArray* prhs[]) {
    if (nlhs != 1 || (nrhs != 5 && nrhs != 11))
        BadArgument("usage: handle = radia_mex('hdiv.field_evaluator.from_cloud', xyz, strength, image_masks, image_signs [, leaf_size, theta, tree_min_sources, auto_min_work, tree_relative_tolerance, probe_count])");
    std::size_t rows = 0, cols = 0;
    auto xyz = RealMatrix(prhs[1], rows, cols, "xyz");
    if (cols != 3)
        BadArgument("xyz must be N-by-3");
    auto strength = RealVector(prhs[2], "strength");
    if (strength.size() != rows)
        BadArgument("strength must have one entry per xyz row");
    const auto options = FieldOptions(nrhs, prhs, 5,
        "handle = radia_mex('hdiv.field_evaluator.from_cloud', xyz, strength, image_masks, image_signs [, options])");
    plhs[0] = Uint64Output(RegisterField(HDivFieldEvaluator::FromCloud(
        std::move(xyz), std::move(strength), IntegerVector(prhs[3], "image_masks"),
        RealVector(prhs[4], "image_signs"), options)));
}

void HDivFieldFromCurvedTet(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    if (nlhs != 1 || (nrhs != 7 && nrhs != 13))
        BadArgument("usage: handle = radia_mex('hdiv.field_evaluator.from_curved_tet', volume, surface, gauss_points, gauss_weights, image_masks, image_signs [, options])");
    std::size_t volume_rows = 0, volume_cols = 0;
    std::size_t surface_rows = 0, surface_cols = 0;
    auto volume = RealMatrix(prhs[1], volume_rows, volume_cols, "volume");
    auto surface = RealMatrix(prhs[2], surface_rows, surface_cols, "surface");
    if (volume_cols != 34 || surface_cols != 24)
        BadArgument("curved volume must be N-by-34 and curved surface must be M-by-24");
    auto gauss_points = RealVector(prhs[3], "gauss_points");
    auto gauss_weights = RealVector(prhs[4], "gauss_weights");
    if (gauss_points.empty() || gauss_points.size() != gauss_weights.size())
        BadArgument("gauss_points and gauss_weights must be nonempty and have equal lengths");
    auto options = FieldOptions(nrhs, prhs, 7,
        "handle = radia_mex('hdiv.field_evaluator.from_curved_tet', volume, surface, gauss_points, gauss_weights, image_masks, image_signs [, options])");
    plhs[0] = Uint64Output(RegisterField(HDivFieldEvaluator::FromCurvedTet(
        std::move(volume), std::move(surface), std::move(gauss_points),
        std::move(gauss_weights), IntegerVector(prhs[5], "image_masks"),
        RealVector(prhs[6], "image_signs"), options)));
}

void HDivFieldDestroy(int nlhs, mxArray* plhs[], int nrhs,
                      const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 0,
               "radia_mex('hdiv.field_evaluator.destroy', handle)");
    DestroyField(Handle(prhs[1]));
}

void HDivFieldField(int nlhs, mxArray* plhs[], int nrhs,
                    const mxArray* prhs[]) {
    if ((nrhs != 3 && nrhs != 4) || nlhs != 1)
        BadArgument("usage: value = radia_mex('hdiv.field_evaluator.field', handle, observations [, algorithm])");
    const auto evaluator = Field(Handle(prhs[1]));
    std::size_t rows = 0, cols = 0;
    auto observations = RealMatrix(prhs[2], rows, cols, "observations");
    if (cols != 3) BadArgument("observations must be N-by-3");
    const auto algorithm = nrhs == 4
        ? HDivFieldEvaluator::ParseAlgorithm(Text(prhs[3], "algorithm"))
        : HDivFieldEvaluator::Algorithm::Auto;
    std::vector<double> output(rows * 3, 0.0);
    evaluator->Evaluate(observations.data(), rows, output.data(), algorithm);
    plhs[0] = RealMatrixOutput(output, rows, 3);
}

void HDivFieldCandidateAlgorithm(int nlhs, mxArray* plhs[], int nrhs,
                                 const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
               "algorithm = radia_mex('hdiv.field_evaluator.candidate_algorithm', handle, n_observations)");
    const auto evaluator = Field(Handle(prhs[1]));
    plhs[0] = TextOutput(HDivFieldEvaluator::AlgorithmName(
        evaluator->AlgorithmFor(static_cast<std::size_t>(NonnegativeLong(
            prhs[2], "n_observations")))));
}

void HDivFieldLastAlgorithm(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "algorithm = radia_mex('hdiv.field_evaluator.last_algorithm', handle)");
    const auto evaluator = Field(Handle(prhs[1]));
    plhs[0] = TextOutput(HDivFieldEvaluator::AlgorithmName(evaluator->LastAlgorithm()));
}

void HDivFieldStats(int nlhs, mxArray* plhs[], int nrhs,
                    const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "stats = radia_mex('hdiv.field_evaluator.stats', handle)");
    const auto evaluator = Field(Handle(prhs[1]));
    double lower[3] = {0.0, 0.0, 0.0}, upper[3] = {0.0, 0.0, 0.0};
    evaluator->Bounds(lower, upper);
    const char* fields[] = {"source_count", "image_count", "tree_nodes", "leaf_size",
                            "theta", "tree_min_sources", "auto_min_work",
                            "tree_relative_tolerance", "probe_count",
                            "source_representation", "last_algorithm", "bounds_min", "bounds_max"};
    plhs[0] = mxCreateStructMatrix(1, 1, 13, fields);
    mxSetField(plhs[0], 0, "source_count", mxCreateDoubleScalar(evaluator->SourceCount()));
    mxSetField(plhs[0], 0, "image_count", mxCreateDoubleScalar(evaluator->ImageCount()));
    mxSetField(plhs[0], 0, "tree_nodes", mxCreateDoubleScalar(evaluator->TreeNodeCount()));
    mxSetField(plhs[0], 0, "leaf_size", mxCreateDoubleScalar(evaluator->LeafSize()));
    mxSetField(plhs[0], 0, "theta", mxCreateDoubleScalar(evaluator->Theta()));
    mxSetField(plhs[0], 0, "tree_min_sources", mxCreateDoubleScalar(evaluator->TreeMinSources()));
    mxSetField(plhs[0], 0, "auto_min_work", mxCreateDoubleScalar(evaluator->AutoMinWork()));
    mxSetField(plhs[0], 0, "tree_relative_tolerance",
               mxCreateDoubleScalar(evaluator->TreeRelativeTolerance()));
    mxSetField(plhs[0], 0, "probe_count", mxCreateDoubleScalar(evaluator->ProbeCount()));
    mxSetField(plhs[0], 0, "source_representation", TextOutput(evaluator->SourceRepresentation()));
    mxSetField(plhs[0], 0, "last_algorithm",
               TextOutput(HDivFieldEvaluator::AlgorithmName(evaluator->LastAlgorithm())));
    mxSetField(plhs[0], 0, "bounds_min", RealRow({lower[0], lower[1], lower[2]}));
    mxSetField(plhs[0], 0, "bounds_max", RealRow({upper[0], upper[1], upper[2]}));
}

void HDivFieldAsCoefficient(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    if ((nrhs != 2 && nrhs != 3) || nlhs != 1)
        BadArgument(
            "handle = radia_mex('hdiv.field_evaluator.as_coefficient', "
            "evaluator_handle [, algorithm])");
    const auto evaluator = Field(Handle(prhs[1]));
    const std::string algorithm =
        nrhs == 3 ? Text(prhs[2], "algorithm") : "direct";
    auto coefficient =
        std::make_shared<radia::ngsolve_bridge::HDivFieldCoefficient>(
            evaluator, algorithm);
    plhs[0] = Uint64Output(RegisterCoefficient(std::move(coefficient)));
}

void PlanarFieldCreate(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    CheckArity(nrhs, 5, nlhs, 1,
               "handle = radia_mex('hdiv.planar_evaluator.create', positions, strengths, image_masks, image_signs)");
    std::size_t rows = 0, cols = 0;
    auto positions = RealMatrix(prhs[1], rows, cols, "positions");
    auto strengths = RealVector(prhs[2], "strengths");
    if (cols != 2 || strengths.size() != rows)
        BadArgument("positions must be N-by-2 with one strength per row");
    auto image_masks = IntegerVector(prhs[3], "image_masks");
    auto image_signs = RealVector(prhs[4], "image_signs");
    if (image_masks.size() != image_signs.size())
        BadArgument("image_masks and image_signs must have equal lengths");
    plhs[0] = Uint64Output(RegisterPlanar(std::make_shared<PlanarFieldEvaluator>(
        std::move(positions), std::move(strengths), std::move(image_masks),
        std::move(image_signs))));
}

void PlanarFieldDestroy(int nlhs, mxArray* plhs[], int nrhs,
                        const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 0,
               "radia_mex('hdiv.planar_evaluator.destroy', handle)");
    DestroyPlanar(Handle(prhs[1]));
}

void PlanarFieldEvaluate(const std::string& command, int nlhs, mxArray* plhs[],
                         int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
               "value = radia_mex('hdiv.planar_evaluator.<field|az>', handle, points)");
    const auto evaluator = Planar(Handle(prhs[1]));
    std::size_t rows = 0, cols = 0;
    auto points = RealMatrix(prhs[2], rows, cols, "points");
    if (cols != 2) BadArgument("points must be N-by-2");
    if (command == "hdiv.planar_evaluator.field") {
        std::vector<double> output(rows * 2, 0.0);
        evaluator->EvaluateField(points.data(), rows, output.data());
        plhs[0] = RealMatrixOutput(output, rows, 2);
    } else {
        std::vector<double> output(rows, 0.0);
        evaluator->EvaluateAz(points.data(), rows, output.data());
        plhs[0] = RealColumn(output);
    }
}

void PlanarFieldStats(int nlhs, mxArray* plhs[], int nrhs,
                      const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
               "stats = radia_mex('hdiv.planar_evaluator.stats', handle)");
    const auto evaluator = Planar(Handle(prhs[1]));
    const char* fields[] = {"base_source_count", "source_count", "image_count"};
    plhs[0] = mxCreateStructMatrix(1, 1, 3, fields);
    mxSetField(plhs[0], 0, "base_source_count", mxCreateDoubleScalar(evaluator->BaseSourceCount()));
    mxSetField(plhs[0], 0, "source_count", mxCreateDoubleScalar(evaluator->SourceCount()));
    mxSetField(plhs[0], 0, "image_count", mxCreateDoubleScalar(evaluator->ImageCount()));
}

void PlanarFieldAsCoefficient(int nlhs, mxArray* plhs[], int nrhs,
                              const mxArray* prhs[]) {
    CheckArity(nrhs, 6, nlhs, 1,
        "handle = radia_mex('hdiv.planar_evaluator.as_coefficient', "
        "evaluator_handle, source_angle, target_angle, center_x, center_y)");
    const auto evaluator = Planar(Handle(prhs[1]));
    auto coefficient =
        std::make_shared<radia::ngsolve_bridge::PlanarHDivFieldCoefficient>(
            evaluator, Scalar(prhs[2], "source_angle"),
            Scalar(prhs[3], "target_angle"),
            Scalar(prhs[4], "center_x"), Scalar(prhs[5], "center_y"));
    plhs[0] = Uint64Output(RegisterCoefficient(std::move(coefficient)));
}

int CreatePolyhedron(const mxArray* vertices_value, const mxArray* magnetization_value,
                     int expected_vertices, const std::vector<int>& faces,
                     const std::vector<int>& face_lengths) {
    std::size_t rows = 0, cols = 0;
    auto vertices = RealMatrix(vertices_value, rows, cols, "vertices");
    if (rows != static_cast<std::size_t>(expected_vertices) || cols != 3)
        BadArgument("vertices have the wrong shape for this element type");
    auto magnetization = FixedRealVector(magnetization_value, 3, "magnetization");
    std::vector<int> mutable_faces = faces;
    std::vector<int> mutable_lengths = face_lengths;
    double linear_m[9] = {0};
    double current[3] = {0};
    double linear_current[9] = {0};
    int handle = 0;
    CheckRadia(RadObjPolyhdr(
        &handle, vertices.data(), expected_vertices, mutable_faces.data(),
        mutable_lengths.data(), static_cast<int>(mutable_lengths.size()),
        magnetization.data(), linear_m, current, linear_current));
    return handle;
}

void RadiaPolyhedron(const std::string& command, int nlhs, mxArray* plhs[],
                     int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
        "handle = radia_mex('radia.Obj<Element>', vertices, magnetization)");
    int handle = 0;
    if (command == "radia.ObjHexahedron") {
        handle = CreatePolyhedron(prhs[1], prhs[2], 8,
            {1,2,3,4, 2,6,7,3, 1,5,6,2, 1,4,8,5, 3,7,8,4, 5,8,7,6},
            {4,4,4,4,4,4});
    } else if (command == "radia.ObjTetrahedron") {
        handle = CreatePolyhedron(prhs[1], prhs[2], 4,
            {1,3,2, 1,2,4, 2,3,4, 3,1,4}, {3,3,3,3});
    } else if (command == "radia.ObjWedge") {
        handle = CreatePolyhedron(prhs[1], prhs[2], 6,
            {1,3,2, 4,5,6, 1,2,5,4, 2,3,6,5, 3,1,4,6}, {3,3,4,4,4});
    } else {
        handle = CreatePolyhedron(prhs[1], prhs[2], 5,
            {1,4,3,2, 1,2,5, 2,3,5, 3,4,5, 4,1,5}, {4,3,3,3,3});
    }
    plhs[0] = mxCreateDoubleScalar(handle);
}

std::vector<char> MutableText(const mxArray* value, const char* name,
                              std::size_t max_length = 0) {
    const std::string text = Text(value, name);
    if (max_length != 0 && text.size() > max_length)
        BadArgument(std::string(name) + " is too long");
    return std::vector<char>(text.begin(), text.end());
}

char SingleCharacter(const mxArray* value, const char* name) {
    const std::string text = Text(value, name);
    if (text.size() != 1)
        BadArgument(std::string(name) + " must contain one character");
    return text[0];
}

void RadiaExtendedObject(const std::string& command, int nlhs, mxArray* plhs[],
                         int nrhs, const mxArray* prhs[]) {
    int handle = 0;
    if (command == "radia.ObjBckg") {
        CheckArity(nrhs, 2, nlhs, 1,
            "handle = radia_mex('radia.ObjBckg', field_vector)");
        auto field = FixedRealVector(prhs[1], 3, "field_vector");
        CheckRadia(RadObjBckg(&handle, field.data()));
    } else if (command == "radia.ObjArcPgnMag") {
        if ((nrhs < 6 || nrhs > 8) || nlhs != 1)
            BadArgument("usage: handle = radia_mex('radia.ObjArcPgnMag', center, axis, vertices, phi_range, nseg [, sym_nosym [, magnetization]])");
        auto center = FixedRealVector(prhs[1], 2, "center");
        auto axis = SingleCharacter(prhs[2], "axis");
        std::size_t rows = 0, cols = 0;
        auto vertices = RealMatrix(prhs[3], rows, cols, "vertices");
        if (rows < 3 || cols != 2)
            BadArgument("vertices must have shape count-by-2 with at least 3 points");
        auto phi = FixedRealVector(prhs[4], 2, "phi_range");
        const int nseg = PositiveInteger(prhs[5], "nseg");
        const char sym_no = (nrhs >= 7 && Text(prhs[6], "sym_nosym") == "sym") ? 's' : 'n';
        double magnetization[3] = {0.0, 0.0, 0.0};
        if (nrhs == 8 && mxGetNumberOfElements(prhs[7]) != 0) {
            auto value = FixedRealVector(prhs[7], 3, "magnetization");
            std::copy(value.begin(), value.end(), magnetization);
        }
        CheckRadia(RadObjArcPgnMag(&handle, center.data(), axis, vertices.data(),
                                   MatrixDimension(rows, "vertices"), phi.data(),
                                   nseg, sym_no, magnetization));
    } else if (command == "radia.ObjThckPgn") {
        CheckArity(nrhs, 6, nlhs, 1,
            "handle = radia_mex('radia.ObjThckPgn', xc, lx, polygon, axis, magnetization)");
        std::size_t rows = 0, cols = 0;
        auto polygon = RealMatrix(prhs[3], rows, cols, "polygon");
        if (rows < 3 || cols != 2)
            BadArgument("polygon must have shape count-by-2 with at least 3 points");
        auto axis = SingleCharacter(prhs[4], "axis");
        auto magnetization = FixedRealVector(prhs[5], 3, "magnetization");
        CheckRadia(RadObjThckPgn(&handle, Scalar(prhs[1], "xc"),
                                  Scalar(prhs[2], "lx"), polygon.data(),
                                  MatrixDimension(rows, "polygon"), axis,
                                  magnetization.data()));
    } else if (command == "radia.ObjCylMag") {
        CheckArity(nrhs, 7, nlhs, 1,
            "handle = radia_mex('radia.ObjCylMag', center, radius, height, nseg, axis, magnetization)");
        auto center = FixedRealVector(prhs[1], 3, "center");
        auto axis = SingleCharacter(prhs[5], "axis");
        auto magnetization = FixedRealVector(prhs[6], 3, "magnetization");
        CheckRadia(RadObjCylMag(&handle, center.data(), Scalar(prhs[2], "radius"),
                                 Scalar(prhs[3], "height"),
                                 PositiveInteger(prhs[4], "nseg"), axis,
                                 magnetization.data()));
    } else if (command == "radia.ObjRecCur") {
        CheckArity(nrhs, 4, nlhs, 1,
            "handle = radia_mex('radia.ObjRecCur', center, dimensions, current_density)");
        auto center = FixedRealVector(prhs[1], 3, "center");
        auto dimensions = FixedRealVector(prhs[2], 3, "dimensions");
        auto current = FixedRealVector(prhs[3], 3, "current_density");
        CheckRadia(RadObjRecCur(&handle, center.data(), dimensions.data(), current.data()));
    } else if (command == "radia.ObjArcCur") {
        CheckArity(nrhs, 9, nlhs, 1,
            "handle = radia_mex('radia.ObjArcCur', center, radii, phi, height, nseg, man_auto, axis, current_density)");
        auto center = FixedRealVector(prhs[1], 3, "center");
        auto radii = FixedRealVector(prhs[2], 2, "radii");
        auto phi = FixedRealVector(prhs[3], 2, "phi");
        auto man_auto = SingleCharacter(prhs[6], "man_auto");
        auto axis = SingleCharacter(prhs[7], "axis");
        CheckRadia(RadObjArcCur(&handle, center.data(), radii.data(), phi.data(),
                                Scalar(prhs[4], "height"), PositiveInteger(prhs[5], "nseg"),
                                man_auto, axis, Scalar(prhs[8], "current_density")));
    } else if (command == "radia.ObjRaceTrk") {
        CheckArity(nrhs, 9, nlhs, 1,
            "handle = radia_mex('radia.ObjRaceTrk', center, radii, lengths, height, nseg, man_auto, axis, current_density)");
        auto center = FixedRealVector(prhs[1], 3, "center");
        auto radii = FixedRealVector(prhs[2], 2, "radii");
        auto lengths = FixedRealVector(prhs[3], 2, "lengths");
        auto man_auto = SingleCharacter(prhs[6], "man_auto");
        auto axis = SingleCharacter(prhs[7], "axis");
        CheckRadia(RadObjRaceTrk(&handle, center.data(), radii.data(), lengths.data(),
                                 Scalar(prhs[4], "height"), PositiveInteger(prhs[5], "nseg"),
                                 man_auto, axis, Scalar(prhs[8], "current_density")));
    } else {
        CheckArity(nrhs, 3, nlhs, 1,
            "handle = radia_mex('radia.ObjFlmCur', points, current)");
        std::size_t rows = 0, cols = 0;
        auto points = RealMatrix(prhs[1], rows, cols, "points");
        if (rows < 2 || cols != 3)
            BadArgument("points must have shape count-by-3 with at least 2 points");
        CheckRadia(RadObjFlmCur(&handle, points.data(), MatrixDimension(rows, "points"),
                                Scalar(prhs[2], "current")));
    }
    plhs[0] = mxCreateDoubleScalar(handle);
}

void RadiaExtendedField(const std::string& command, int nlhs, mxArray* plhs[],
                        int nrhs, const mxArray* prhs[]) {
    if (command == "radia.FldFrcShpRtg") {
        CheckArity(nrhs, 3, nlhs, 1,
            "shape = radia_mex('radia.FldFrcShpRtg', center, dimensions)");
        auto center = FixedRealVector(prhs[1], 3, "center");
        auto dimensions = FixedRealVector(prhs[2], 2, "dimensions");
        int shape = 0;
        CheckRadia(RadFldFrcShpRtg(&shape, center.data(), dimensions.data()));
        plhs[0] = mxCreateDoubleScalar(shape);
    } else if (command == "radia.FldFrc") {
        CheckArity(nrhs, 3, nlhs, 1,
            "force = radia_mex('radia.FldFrc', object, shape)");
        double force[6] = {0.0};
        int count = 6;
        CheckRadia(RadFldFrc(force, &count, PositiveInteger(prhs[1], "object"),
                             PositiveInteger(prhs[2], "shape")));
        if (count < 0 || count > 6)
            throw std::runtime_error("Radia returned an invalid force size");
        plhs[0] = RealRow(std::vector<double>(force, force + count));
    } else if (command == "radia.FldLst") {
        if ((nrhs < 6 || nrhs > 8) || nlhs != 1)
            BadArgument("usage: values = radia_mex('radia.FldLst', object, field_id, p1, p2, np [, arg_opt [, start]])");
        const int object = PositiveInteger(prhs[1], "object");
        const std::string field_id = Text(prhs[2], "field_id");
        auto p1 = FixedRealVector(prhs[3], 3, "p1");
        auto p2 = FixedRealVector(prhs[4], 3, "p2");
        const int np = PositiveInteger(prhs[5], "np");
        const std::string arg_opt = nrhs >= 7 ? Text(prhs[6], "arg_opt") : "noarg";
        const double start = nrhs >= 8 ? Scalar(prhs[7], "start") : 0.0;
        std::string field_lower = field_id;
        std::transform(field_lower.begin(), field_lower.end(), field_lower.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        std::string arg_lower = arg_opt;
        std::transform(arg_lower.begin(), arg_lower.end(), arg_lower.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (arg_lower != "arg" && arg_lower != "noarg")
            BadArgument("arg_opt must be 'arg' or 'noarg'");
        const bool with_arg = arg_lower == "arg";
        const bool is_vector = field_lower.empty() || field_lower == "b" ||
                               field_lower == "h" || field_lower == "a" ||
                               field_lower == "m";
        const std::size_t values_per_point = is_vector ? 3 : 1;
        const std::size_t stride = values_per_point + (with_arg ? 1 : 0);
        std::vector<double> result(static_cast<std::size_t>(np) * stride, 0.0);
        std::vector<char> id_buffer(field_id.begin(), field_id.end());
        id_buffer.push_back('\0');
        std::vector<char> arg_buffer(arg_lower.begin(), arg_lower.end());
        arg_buffer.push_back('\0');
        int count = 0;
        CheckRadia(RadFldLst(result.data(), &count, object, id_buffer.data(),
                             p1.data(), p2.data(), np, arg_buffer.data(), start));
        plhs[0] = RealMatrixOutput(result, static_cast<std::size_t>(np), stride);
    } else if (command == "radia.FldInt") {
        CheckArity(nrhs, 6, nlhs, 1,
            "integral = radia_mex('radia.FldInt', object, inf_fin, field_id, p1, p2)");
        auto p1 = FixedRealVector(prhs[4], 3, "p1");
        auto p2 = FixedRealVector(prhs[5], 3, "p2");
        auto inf_fin = MutableText(prhs[2], "inf_fin", 15);
        auto field_id = MutableText(prhs[3], "field_id", 15);
        inf_fin.push_back('\0');
        field_id.push_back('\0');
        double result[3] = {0.0, 0.0, 0.0};
        int count = 0;
        CheckRadia(RadFldInt(result, &count, PositiveInteger(prhs[1], "object"),
                             inf_fin.data(), field_id.data(), p1.data(), p2.data()));
        if (count < 0 || count > 3)
            throw std::runtime_error("Radia returned an invalid field-integral size");
        if (count == 1)
            plhs[0] = mxCreateDoubleScalar(result[0]);
        else
            plhs[0] = RealRow(std::vector<double>(result, result + count));
    } else if (command == "radia.ObjCenFld") {
        if ((nrhs != 2 && nrhs != 3) || nlhs != 1)
            BadArgument("usage: result = radia_mex('radia.ObjCenFld', object [, field_type])");
        const std::string field_type = nrhs == 3 ? Text(prhs[2], "field_type") : "B";
        const char type = field_type.empty() ? 'B' : field_type[0];
        double result[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        int mesh[3] = {1, 1, 1};
        CheckRadia(RadObjCenFld(result, mesh, PositiveInteger(prhs[1], "object"), type));
        const char* fields[] = {"center", "field"};
        plhs[0] = mxCreateStructMatrix(1, 1, 2, fields);
        mxSetField(plhs[0], 0, "center", RealRow({result[0], result[1], result[2]}));
        mxSetField(plhs[0], 0, "field", RealRow({result[3], result[4], result[5]}));
    } else if (command == "radia.FldCmpCrt") {
        CheckArity(nrhs, 7, nlhs, 0,
            "radia_mex('radia.FldCmpCrt', prcB, prcA, prcBInt, prcFrc, prcTrjCrd, prcTrjAng)");
        int result = 0;
        CheckRadia(RadFldCmpCrt(&result, Scalar(prhs[1], "prcB"),
                                 Scalar(prhs[2], "prcA"), Scalar(prhs[3], "prcBInt"),
                                 Scalar(prhs[4], "prcFrc"), Scalar(prhs[5], "prcTrjCrd"),
                                 Scalar(prhs[6], "prcTrjAng")));
    } else if (command == "radia.FldCmpPrc") {
        CheckArity(nrhs, 2, nlhs, 0,
            "radia_mex('radia.FldCmpPrc', options)");
        auto options = MutableText(prhs[1], "options");
        options.push_back('\0');
        int result = 0;
        CheckRadia(RadFldCmpPrc(&result, options.data()));
    } else if (command == "radia.FldLenRndSw") {
        CheckArity(nrhs, 2, nlhs, 0,
            "radia_mex('radia.FldLenRndSw', on_off)");
        auto option = MutableText(prhs[1], "on_off", 15);
        option.push_back('\0');
        int result = 0;
        CheckRadia(RadFldLenRndSw(&result, option.data()));
    } else {
        if ((nrhs != 3 && nrhs != 4) || nlhs != 0)
            BadArgument("usage: radia_mex('radia.FldLenTol', abs_val, rel_val [, zero_val])");
        int result = 0;
        CheckRadia(RadFldLenTol(&result, Scalar(prhs[1], "abs_val"),
                                Scalar(prhs[2], "rel_val"),
                                nrhs == 4 ? Scalar(prhs[3], "zero_val") : 0.0));
    }
}

void RadiaObjCnt(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1, "handle = radia_mex('radia.ObjCnt', handles)");
    auto handles = IntegerVector(prhs[1], "handles");
    int result = 0;
    CheckRadia(RadObjCnt(&result, handles.empty() ? nullptr : handles.data(),
                         static_cast<int>(handles.size())));
    plhs[0] = mxCreateDoubleScalar(result);
}

void RadiaMatLin(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if ((nrhs != 2 && nrhs != 3) || nlhs != 1)
        BadArgument("usage: mat = radia_mex('radia.MatLin', mu_r [, easy_axis])");
    auto mu = RealVector(prhs[1], "mu_r");
    int material = 0;
    if (mu.size() == 1) {
        if (mu[0] < 1.0)
            BadArgument("mu_r must be at least 1");
        CheckRadia(RadMatLinIso(&material, mu[0] - 1.0));
    } else if (mu.size() == 2) {
        if (nrhs != 3)
            BadArgument("anisotropic mu_r requires easy_axis");
        if (mu[0] < 1.0 || mu[1] < 1.0)
            BadArgument("mu_r entries must be at least 1");
        auto axis = FixedRealVector(prhs[2], 3, "easy_axis");
        double susceptibility[2] = {mu[0] - 1.0, mu[1] - 1.0};
        CheckRadia(RadMatLinAniso(&material, susceptibility, axis.data()));
    } else {
        BadArgument("mu_r must be scalar or have two entries");
    }
    plhs[0] = mxCreateDoubleScalar(material);
}

void RadiaMatSatIsoTab(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 1,
        "mat = radia_mex('radia.MatSatIsoTab', bh_data)");
    std::size_t rows = 0, cols = 0;
    auto table = RealMatrix(prhs[1], rows, cols, "bh_data");
    if (rows == 0 || cols != 2)
        BadArgument("bh_data must have shape count-by-2 with [H,B] rows");
    int material = 0;
    CheckRadia(RadMatSatIsoTab(&material, table.data(), static_cast<int>(rows)));
    plhs[0] = mxCreateDoubleScalar(material);
}

void RadiaExtendedMaterial(const std::string& command, int nlhs, mxArray* plhs[],
                           int nrhs, const mxArray* prhs[]) {
    if (command == "radia.MatSatIsoFrm") {
        CheckArity(nrhs, 2, nlhs, 1,
            "mat = radia_mex('radia.MatSatIsoFrm', params)");
        std::size_t rows = 0, cols = 0;
        auto params = RealMatrix(prhs[1], rows, cols, "params");
        if (rows > 3 || cols != 2)
            BadArgument("params must have at most 3 rows and 2 columns");
        double ksi_ms[3][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
        for (std::size_t i = 0; i < rows; ++i) {
            ksi_ms[i][0] = params[i * 2];
            ksi_ms[i][1] = params[i * 2 + 1];
        }
        int material = 0;
        CheckRadia(RadMatSatIsoFrm(&material, ksi_ms[0], ksi_ms[1], ksi_ms[2]));
        plhs[0] = mxCreateDoubleScalar(material);
    } else if (command == "radia.MatSatAniso") {
        CheckArity(nrhs, 3, nlhs, 1,
            "mat = radia_mex('radia.MatSatAniso', data_par, data_per)");
        auto data_par = RealVector(prhs[1], "data_par");
        auto data_per = RealVector(prhs[2], "data_per");
        if (data_par.empty() || data_per.empty())
            BadArgument("data_par and data_per must be non-empty coefficient vectors");
        int material = 0;
        CheckRadia(RadMatSatAniso(&material, data_par.data(),
                                  MatrixDimension(data_par.size(), "data_par"),
                                  data_per.data(), MatrixDimension(data_per.size(), "data_per")));
        plhs[0] = mxCreateDoubleScalar(material);
    } else if (command == "radia.MatSatLamTab") {
        CheckArity(nrhs, 4, nlhs, 1,
            "mat = radia_mex('radia.MatSatLamTab', mh_data, packing_factor, normal)");
        std::size_t rows = 0, cols = 0;
        auto mh_data = RealMatrix(prhs[1], rows, cols, "mh_data");
        if (rows == 0 || cols != 2)
            BadArgument("mh_data must be a non-empty count-by-2 matrix");
        auto normal = FixedRealVector(prhs[3], 3, "normal");
        int material = 0;
        CheckRadia(RadMatSatLamTab(&material, mh_data.data(),
                                   MatrixDimension(rows, "mh_data"),
                                   Scalar(prhs[2], "packing_factor"), normal.data()));
        plhs[0] = mxCreateDoubleScalar(material);
    } else if (command == "radia.MatSatLamFrm") {
        CheckArity(nrhs, 6, nlhs, 1,
            "mat = radia_mex('radia.MatSatLamFrm', ksi_ms1, ksi_ms2, ksi_ms3, packing, normal)");
        auto ksi_ms1 = FixedRealVector(prhs[1], 2, "ksi_ms1");
        auto ksi_ms2 = FixedRealVector(prhs[2], 2, "ksi_ms2");
        auto ksi_ms3 = FixedRealVector(prhs[3], 2, "ksi_ms3");
        auto normal = FixedRealVector(prhs[5], 3, "normal");
        int material = 0;
        CheckRadia(RadMatSatLamFrm(&material, ksi_ms1.data(), ksi_ms2.data(),
                                   ksi_ms3.data(), Scalar(prhs[4], "packing"),
                                   normal.data()));
        plhs[0] = mxCreateDoubleScalar(material);
    } else {
        CheckArity(nrhs, 4, nlhs, 1,
            "magnetization = radia_mex('radia.MatMvsH', material, component, h_field)");
        auto component = MutableText(prhs[2], "component", 15);
        component.push_back('\0');
        auto h_field = FixedRealVector(prhs[3], 3, "h_field");
        double magnetization[3] = {0.0, 0.0, 0.0};
        int count = 3;
        CheckRadia(RadMatMvsH(magnetization, &count,
                              PositiveInteger(prhs[1], "material"), component.data(),
                              h_field.data()));
        if (count < 0 || count > 3)
            throw std::runtime_error("Radia returned an invalid magnetization size");
        plhs[0] = RealRow(std::vector<double>(magnetization, magnetization + count));
    }
}

void RadiaHysteresis(const std::string& command, int nlhs, mxArray* plhs[],
                     int nrhs, const mxArray* prhs[]) {
    if (command == "radia.MatEnergyHysteresis" || command == "radia.MatPlayHysteresis") {
        const bool energy = command == "radia.MatEnergyHysteresis";
        if ((energy && (nrhs != 6 && nrhs != 7)) || (!energy && nrhs != 6) || nlhs != 1)
            BadArgument(energy
                ? "usage: mat = radia_mex('radia.MatEnergyHysteresis', K, chi, r_flat, f_flat, table_sizes [, eps])"
                : "usage: mat = radia_mex('radia.MatPlayHysteresis', K, eta, r_flat, f_flat, table_sizes)");
        const int K = PositiveInteger(prhs[1], "K");
        auto thresholds = RealVector(prhs[2], energy ? "chi" : "eta");
        auto r_flat = RealVector(prhs[3], "r_flat");
        auto f_flat = RealVector(prhs[4], "f_flat");
        auto table_sizes_raw = IntegerVector(prhs[5], "table_sizes");
        if (thresholds.size() != static_cast<std::size_t>(K))
            BadArgument("threshold vector must have length K");
        if (r_flat.size() != f_flat.size())
            BadArgument("r_flat and f_flat must have the same length");
        if (table_sizes_raw.size() != static_cast<std::size_t>(K))
            BadArgument("table_sizes must have length K");
        std::size_t total = 0;
        for (int size : table_sizes_raw) {
            if (size <= 0)
                BadArgument("table_sizes entries must be positive");
            total += static_cast<std::size_t>(size);
        }
        if (total != r_flat.size())
            BadArgument("table_sizes must sum to the flattened table length");
        int material = 0;
        if (energy) {
            const double eps = nrhs == 7 ? Scalar(prhs[6], "eps") : 1e-8;
            CheckRadia(RadMatEnergyHysteresis(&material, K, thresholds.data(),
                                               r_flat.data(), f_flat.data(),
                                               table_sizes_raw.data(), eps));
        } else {
            CheckRadia(RadMatPlayHysteresis(&material, K, thresholds.data(),
                                             r_flat.data(), f_flat.data(),
                                             table_sizes_raw.data()));
        }
        plhs[0] = mxCreateDoubleScalar(material);
        return;
    }

    if (command == "radia.MatHysSaveState") {
        CheckArity(nrhs, 2, nlhs, 1,
            "state = radia_mex('radia.MatHysSaveState', material)");
        int length = 0;
        CheckRadia(RadMatHysSaveState(PositiveInteger(prhs[1], "material"),
                                      nullptr, &length));
        if (length < 0)
            throw std::runtime_error("Radia returned an invalid hysteresis state size");
        std::vector<double> state(static_cast<std::size_t>(length));
        CheckRadia(RadMatHysSaveState(PositiveInteger(prhs[1], "material"),
                                      state.empty() ? nullptr : state.data(), &length));
        plhs[0] = RealColumn(state);
    } else if (command == "radia.MatHysRestoreState") {
        CheckArity(nrhs, 3, nlhs, 0,
            "radia_mex('radia.MatHysRestoreState', material, state)");
        auto state = RealVector(prhs[2], "state");
        CheckRadia(RadMatHysRestoreState(PositiveInteger(prhs[1], "material"),
                                          state.empty() ? nullptr : state.data(),
                                          MatrixDimension(state.size(), "state")));
    } else if (command == "radia.MatHysCommitState") {
        CheckArity(nrhs, 2, nlhs, 0,
            "radia_mex('radia.MatHysCommitState', material)");
        CheckRadia(RadMatHysCommitState(PositiveInteger(prhs[1], "material")));
    } else if (command == "radia.MatHysGetNuRev") {
        CheckArity(nrhs, 2, nlhs, 1,
            "nu_rev = radia_mex('radia.MatHysGetNuRev', material)");
        double nu_rev = 0.0;
        CheckRadia(RadMatHysGetNuRev(PositiveInteger(prhs[1], "material"), &nu_rev));
        plhs[0] = mxCreateDoubleScalar(nu_rev);
    } else {
        CheckArity(nrhs, 3, nlhs, 1,
            "h_irr = radia_mex('radia.MatHysIrreversible', material, B)");
        auto B = FixedRealVector(prhs[2], 3, "B");
        double h_irr[3] = {0.0, 0.0, 0.0};
        CheckRadia(RadMatHysIrreversible(PositiveInteger(prhs[1], "material"),
                                          B.data(), h_irr));
        plhs[0] = RealRow({h_irr[0], h_irr[1], h_irr[2]});
    }
}

void RadiaHysteresisBatch(const std::string& command, int nlhs, mxArray* plhs[],
                          int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "output = radia_mex('radia.MatHysForwardBatch|radia.MatHysCommitBatch', material, B, states)");
    const int material = PositiveInteger(prhs[1], "material");
    std::size_t b_rows = 0, b_cols = 0, state_rows = 0, state_cols = 0;
    auto B = RealMatrix(prhs[2], b_rows, b_cols, "B");
    auto states = RealMatrix(prhs[3], state_rows, state_cols, "states");
    if (b_cols != 3)
        BadArgument("B must have shape count-by-3");
    if (state_rows != b_rows)
        BadArgument("states must have one row per B row");
    const int count = MatrixDimension(b_rows, "B");
    const int state_length = MatrixDimension(state_cols, "states");

    if (command == "radia.MatHysForwardBatch") {
        double nu_rev = 0.0;
        CheckRadia(RadMatHysGetNuRev(material, &nu_rev));
        std::vector<double> H(b_rows * 3, 0.0);
        for (int i = 0; i < count; ++i) {
            CheckRadia(RadMatHysRestoreState(
                material, states.data() + static_cast<std::size_t>(i) * state_cols,
                state_length));
            double pB[3] = {
                B[static_cast<std::size_t>(i) * 3 + 0],
                B[static_cast<std::size_t>(i) * 3 + 1],
                B[static_cast<std::size_t>(i) * 3 + 2]};
            double h_irr[3] = {0.0, 0.0, 0.0};
            CheckRadia(RadMatHysIrreversible(material, pB, h_irr));
            H[static_cast<std::size_t>(i) * 3 + 0] = nu_rev * pB[0] + h_irr[0];
            H[static_cast<std::size_t>(i) * 3 + 1] = nu_rev * pB[1] + h_irr[1];
            H[static_cast<std::size_t>(i) * 3 + 2] = nu_rev * pB[2] + h_irr[2];
        }
        plhs[0] = RealMatrixOutput(H, b_rows, 3);
        return;
    }

    int saved_length = 0;
    CheckRadia(RadMatHysSaveState(material, nullptr, &saved_length));
    if (saved_length != state_length)
        BadArgument("states column count does not match the material state size");
    std::vector<double> output(b_rows * state_cols, 0.0);
    for (int i = 0; i < count; ++i) {
        const double* input_state = states.data() + static_cast<std::size_t>(i) * state_cols;
        double* output_state = output.data() + static_cast<std::size_t>(i) * state_cols;
        CheckRadia(RadMatHysRestoreState(material, input_state, state_length));
        double pB[3] = {
            B[static_cast<std::size_t>(i) * 3 + 0],
            B[static_cast<std::size_t>(i) * 3 + 1],
            B[static_cast<std::size_t>(i) * 3 + 2]};
        double h_irr[3] = {0.0, 0.0, 0.0};
        CheckRadia(RadMatHysIrreversible(material, pB, h_irr));
        CheckRadia(RadMatHysCommitState(material));
        int output_length = state_length;
        CheckRadia(RadMatHysSaveState(material, output_state, &output_length));
        if (output_length != state_length)
            throw std::runtime_error("Radia returned an inconsistent hysteresis state size");
    }
    plhs[0] = RealMatrixOutput(output, b_rows, state_cols);
}

void RadiaMatApl(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
        "obj = radia_mex('radia.MatApl', object, material)");
    int result = 0;
    CheckRadia(RadMatApl(&result, PositiveInteger(prhs[1], "object"),
                         PositiveInteger(prhs[2], "material")));
    plhs[0] = mxCreateDoubleScalar(result);
}

void RadiaSolve(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 6, nlhs, 1,
        "result = radia_mex('radia.Solve', object, precision, max_iter, method, image)");
    const int object = PositiveInteger(prhs[1], "object");
    const double precision = Scalar(prhs[2], "precision");
    const int max_iter = PositiveInteger(prhs[3], "max_iter");
    const int method = IntegerScalar(prhs[4], "method");
    const std::string image = Text(prhs[5], "image");
    double result[4] = {0};
    int count = 0;
    CheckRadia(RadSolve(result, &count, object, precision, max_iter, method,
                        image.empty() ? nullptr : image.c_str()));
    plhs[0] = RealRow(std::vector<double>(result, result + std::min(4, count)));
}

void RadiaSolveNonlinear(int nlhs, mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]) {
    if ((nrhs != 6 && nrhs != 7) || nlhs != 1)
        BadArgument("usage: result = radia_mex('radia.SolveNonl', object, precision, max_iter, method, nonl_method [, image])");
    const int object = PositiveInteger(prhs[1], "object");
    const double precision = Scalar(prhs[2], "precision");
    const int max_iter = PositiveInteger(prhs[3], "max_iter");
    const int method = IntegerScalar(prhs[4], "method");
    const int nonl_method = IntegerScalar(prhs[5], "nonl_method");
    const std::string image = nrhs == 7 ? Text(prhs[6], "image") : "";
    double result[4] = {0.0, 0.0, 0.0, 0.0};
    int count = 4;
    CheckRadia(RadSolveNonl(result, &count, object, precision, max_iter,
                             method, nonl_method,
                             image.empty() ? nullptr : image.c_str()));
    if (count < 0 || count > 4)
        throw std::runtime_error("Radia returned an invalid nonlinear-solve result size");
    plhs[0] = RealRow(std::vector<double>(result, result + count));
}

void RadiaSolveStats(int nlhs, mxArray* plhs[], int nrhs) {
    CheckArity(nrhs, 1, nlhs, 1,
        "stats = radia_mex('radia.GetSolveStats')");
    double stats[20] = {0.0};
    int count = 0;
    CheckRadia(RadGetSolveStats(stats, &count));
    if (count < 0 || count > 20)
        throw std::runtime_error("Radia returned an invalid solve-statistics size");
    const char* fields[] = {
        "t_matrix_build", "t_linear_solve", "linear_iterations", "nonl_iterations",
        "taskmanager_enabled", "num_threads", "t_lu_decomp", "t_hmatrix_build",
        "t_moment_fieldgrad", "t_moment_system_build"};
    plhs[0] = mxCreateStructMatrix(1, 1, 10, fields);
    for (int i = 0; i < 10; ++i)
        mxSetField(plhs[0], 0, fields[i], mxCreateDoubleMatrix(0, 0, mxREAL));
    if (count >= 4) {
        mxSetField(plhs[0], 0, "t_matrix_build", mxCreateDoubleScalar(stats[0]));
        mxSetField(plhs[0], 0, "t_linear_solve", mxCreateDoubleScalar(stats[1]));
        mxSetField(plhs[0], 0, "linear_iterations", mxCreateDoubleScalar(stats[2]));
        mxSetField(plhs[0], 0, "nonl_iterations", mxCreateDoubleScalar(stats[3]));
    }
    if (count >= 6) {
        mxSetField(plhs[0], 0, "taskmanager_enabled",
                   mxCreateLogicalScalar(stats[4] > 0.5));
        mxSetField(plhs[0], 0, "num_threads", mxCreateDoubleScalar(stats[5]));
    }
    if (count >= 7)
        mxSetField(plhs[0], 0, "t_lu_decomp", mxCreateDoubleScalar(stats[6]));
    if (count >= 8)
        mxSetField(plhs[0], 0, "t_hmatrix_build", mxCreateDoubleScalar(stats[7]));
    if (count >= 11)
        mxSetField(plhs[0], 0, "t_moment_fieldgrad", mxCreateDoubleScalar(stats[10]));
    if (count >= 12)
        mxSetField(plhs[0], 0, "t_moment_system_build", mxCreateDoubleScalar(stats[11]));
}

void RadiaInteraction(const std::string& command, int nlhs, mxArray* plhs[],
                      int nrhs, const mxArray* prhs[]) {
    if (command == "radia.BuildMatrix") {
        if ((nrhs != 2 && nrhs != 3) || nlhs != 1)
            BadArgument("usage: handle = radia_mex('radia.BuildMatrix', object [, image])");
        const std::string image = nrhs == 3 ? Text(prhs[2], "image") : "";
        int handle = 0;
        CheckRadia(RadBuildMatrix(&handle, PositiveInteger(prhs[1], "object"),
                                  image.empty() ? nullptr : image.c_str()));
        plhs[0] = mxCreateDoubleScalar(handle);
    } else if (command == "radia.GetInteractMatrix") {
        if ((nrhs != 2) || (nlhs != 1 && nlhs != 2))
            BadArgument("usage: [matrix, dof] = radia_mex('radia.GetInteractMatrix', handle)");
        int dof = 0;
        CheckRadia(RadGetInteractMatrix(nullptr, &dof,
                                        PositiveInteger(prhs[1], "handle")));
        if (dof <= 0)
            throw std::runtime_error("Radia returned an empty interaction matrix");
        const std::size_t size = static_cast<std::size_t>(dof) * static_cast<std::size_t>(dof);
        std::vector<double> matrix(size, 0.0);
        int returned_dof = dof;
        CheckRadia(RadGetInteractMatrix(matrix.data(), &returned_dof,
                                        PositiveInteger(prhs[1], "handle")));
        if (returned_dof != dof)
            throw std::runtime_error("Radia interaction-matrix size changed between queries");
        plhs[0] = RealMatrixOutput(matrix, static_cast<std::size_t>(dof),
                                   static_cast<std::size_t>(dof));
        if (nlhs == 2)
            plhs[1] = mxCreateDoubleScalar(dof);
    } else {
        if ((nrhs != 2) || (nlhs != 1 && nlhs != 2))
            BadArgument("usage: [geometry, dof] = radia_mex('radia.GetFaceGeom', handle)");
        int dof = 0;
        CheckRadia(RadGetFaceGeom(nullptr, &dof,
                                  PositiveInteger(prhs[1], "handle")));
        if (dof <= 0)
            throw std::runtime_error("Radia returned empty face geometry");
        std::vector<double> geometry(static_cast<std::size_t>(dof) * 11, 0.0);
        int returned_dof = dof;
        CheckRadia(RadGetFaceGeom(geometry.data(), &returned_dof,
                                  PositiveInteger(prhs[1], "handle")));
        if (returned_dof != dof)
            throw std::runtime_error("Radia face-geometry size changed between queries");
        plhs[0] = RealMatrixOutput(geometry, static_cast<std::size_t>(dof), 11);
        if (nlhs == 2)
            plhs[1] = mxCreateDoubleScalar(dof);
    }
}

void RadiaPlanar(const std::string& command, int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[]) {
    if (command == "radia.PlanarChargeField" || command == "radia.PlanarChargeAz") {
        CheckArity(nrhs, 4, nlhs, 1,
            "value = radia_mex('radia.PlanarChargeField|radia.PlanarChargeAz', Xq, Q, P)");
        std::size_t nq_rows = 0, nq_cols = 0, np_rows = 0, np_cols = 0;
        auto xq = RealMatrix(prhs[1], nq_rows, nq_cols, "Xq");
        auto q = RealVector(prhs[2], "Q");
        auto p = RealMatrix(prhs[3], np_rows, np_cols, "P");
        if (nq_cols != 2)
            BadArgument("Xq must have shape count-by-2");
        if (q.size() != nq_rows)
            BadArgument("Q must have one entry per Xq row");
        if (np_cols != 2)
            BadArgument("P must have shape count-by-2");
        const int nq = MatrixDimension(nq_rows, "Xq");
        const int np = MatrixDimension(np_rows, "P");
        if (command == "radia.PlanarChargeField") {
            std::vector<double> result(np_rows * 2, 0.0);
            rad_planar_charges::Field(nq, xq.data(), q.data(), np, p.data(), result.data());
            plhs[0] = RealMatrixOutput(result, np_rows, 2);
        } else {
            std::vector<double> result(np_rows, 0.0);
            rad_planar_charges::FieldAz(nq, xq.data(), q.data(), np, p.data(), result.data());
            plhs[0] = RealColumn(result);
        }
        return;
    }

    if (command == "radia.PlanarMaxwellTorqueCircle") {
        CheckArity(nrhs, 9, nlhs, 1,
            "torque = radia_mex('radia.PlanarMaxwellTorqueCircle', Xq, Q, Rc, cx, cy, n, hextx, hexty)");
        std::size_t rows = 0, cols = 0;
        auto xq = RealMatrix(prhs[1], rows, cols, "Xq");
        auto q = RealVector(prhs[2], "Q");
        if (cols != 2 || q.size() != rows)
            BadArgument("Xq must be count-by-2 and Q must match its row count");
        const double rc = Scalar(prhs[3], "Rc");
        if (rc <= 0.0)
            BadArgument("Rc must be positive");
        const int n = PositiveInteger(prhs[6], "n");
        plhs[0] = mxCreateDoubleScalar(rad_planar_charges::MaxwellTorqueCircle(
            MatrixDimension(rows, "Xq"), xq.data(), q.data(), rc,
            Scalar(prhs[4], "cx"), Scalar(prhs[5], "cy"), n,
            Scalar(prhs[7], "hextx"), Scalar(prhs[8], "hexty")));
        return;
    }

    CheckArity(nrhs, 9, nlhs, 1,
        "force = radia_mex('radia.PlanarMaxwellForceCircle', Xq, Q, Rc, cx, cy, n, hextx, hexty)");
    std::size_t rows = 0, cols = 0;
    auto xq = RealMatrix(prhs[1], rows, cols, "Xq");
    auto q = RealVector(prhs[2], "Q");
    if (cols != 2 || q.size() != rows)
        BadArgument("Xq must be count-by-2 and Q must match its row count");
    const double rc = Scalar(prhs[3], "Rc");
    if (rc <= 0.0)
        BadArgument("Rc must be positive");
    const int n = PositiveInteger(prhs[6], "n");
    std::vector<double> force(2, 0.0);
    rad_planar_charges::MaxwellForceCircle(
        MatrixDimension(rows, "Xq"), xq.data(), q.data(), rc,
        Scalar(prhs[4], "cx"), Scalar(prhs[5], "cy"), n,
        Scalar(prhs[7], "hextx"), Scalar(prhs[8], "hexty"), force.data());
    plhs[0] = RealRow(force);
}

void RadiaAverageField(const std::string& command, int nlhs, mxArray* plhs[],
                       int nrhs, const mxArray* prhs[]) {
    if (command == "radia.AverageBInBox") {
        CheckArity(nrhs, 6, nlhs, 1,
            "B = radia_mex('radia.AverageBInBox', M, src_min, src_max, tgt_min, tgt_max)");
        auto M = FixedRealVector(prhs[1], 3, "M");
        auto src_min = FixedRealVector(prhs[2], 3, "src_min");
        auto src_max = FixedRealVector(prhs[3], 3, "src_max");
        auto tgt_min = FixedRealVector(prhs[4], 3, "tgt_min");
        auto tgt_max = FixedRealVector(prhs[5], 3, "tgt_max");
        double B[3] = {0.0, 0.0, 0.0};
        radia::average_field::AverageBInBox(
            M.data(), src_min.data(), src_max.data(),
            tgt_min.data(), tgt_max.data(), B);
        plhs[0] = RealRow({B[0], B[1], B[2]});
        return;
    }

    CheckArity(nrhs, 5, nlhs, 1,
        "A = radia_mex('radia.AverageDemagTensor', src_min, src_max, tgt_min, tgt_max)");
    auto src_min = FixedRealVector(prhs[1], 3, "src_min");
    auto src_max = FixedRealVector(prhs[2], 3, "src_max");
    auto tgt_min = FixedRealVector(prhs[3], 3, "tgt_min");
    auto tgt_max = FixedRealVector(prhs[4], 3, "tgt_max");
    double A[9] = {0.0};
    radia::average_field::AverageDemagTensor(
        src_min.data(), src_max.data(), tgt_min.data(), tgt_max.data(), A);
    plhs[0] = RealMatrixOutput(std::vector<double>(A, A + 9), 3, 3);
}

void RadiaEquivalenceSource(const std::string& command, int nlhs, mxArray* plhs[],
                            int nrhs, const mxArray* prhs[]) {
    if (command == "equivalence.static_h") {
        if ((nrhs != 6 && nrhs != 7) || nlhs != 1)
            BadArgument("usage: H = radia_mex('equivalence.static_h', centroids, normals, areas, H_surf, obs [, n_threads])");
        std::size_t face_rows = 0, face_cols = 0, normal_rows = 0, normal_cols = 0;
        std::size_t field_rows = 0, field_cols = 0, obs_rows = 0, obs_cols = 0;
        auto centroids = RealMatrix(prhs[1], face_rows, face_cols, "centroids");
        auto normals = RealMatrix(prhs[2], normal_rows, normal_cols, "normals");
        auto areas = RealVector(prhs[3], "areas");
        auto H_surf = RealMatrix(prhs[4], field_rows, field_cols, "H_surf");
        auto obs = RealMatrix(prhs[5], obs_rows, obs_cols, "obs");
        if (face_cols != 3 || normal_cols != 3 || field_cols != 3 || obs_cols != 3 ||
            normal_rows != face_rows || areas.size() != face_rows || field_rows != face_rows)
            BadArgument("centroids, normals, H_surf, and obs must be count-by-3; face arrays must match");
        const int n_threads = nrhs == 7 ? NonnegativeInteger(prhs[6], "n_threads") : 0;
        std::vector<double> output(obs_rows * 3, 0.0);
        radia::eqsrc::EvaluateStaticH(
            centroids.data(), normals.data(), areas.data(), H_surf.data(),
            MatrixDimension(face_rows, "centroids"), obs.data(),
            MatrixDimension(obs_rows, "obs"), output.data(), n_threads);
        plhs[0] = RealMatrixOutput(output, obs_rows, 3);
        return;
    }

    if ((nrhs != 10 && nrhs != 11) || nlhs != 4)
        BadArgument("usage: [Ere,Eim,Hre,Him] = radia_mex('equivalence.harmonic', centroids, normals, areas, Ere, Eim, Hre, Him, obs, omega [, n_threads])");
    std::size_t face_rows = 0, face_cols = 0, normal_rows = 0, normal_cols = 0;
    std::size_t er_rows = 0, er_cols = 0, ei_rows = 0, ei_cols = 0;
    std::size_t hr_rows = 0, hr_cols = 0, hi_rows = 0, hi_cols = 0;
    std::size_t obs_rows = 0, obs_cols = 0;
    auto centroids = RealMatrix(prhs[1], face_rows, face_cols, "centroids");
    auto normals = RealMatrix(prhs[2], normal_rows, normal_cols, "normals");
    auto areas = RealVector(prhs[3], "areas");
    auto E_re = RealMatrix(prhs[4], er_rows, er_cols, "E_re");
    auto E_im = RealMatrix(prhs[5], ei_rows, ei_cols, "E_im");
    auto H_re = RealMatrix(prhs[6], hr_rows, hr_cols, "H_re");
    auto H_im = RealMatrix(prhs[7], hi_rows, hi_cols, "H_im");
    auto obs = RealMatrix(prhs[8], obs_rows, obs_cols, "obs");
    const double omega = Scalar(prhs[9], "omega");
    if (omega <= 0.0)
        BadArgument("omega must be positive for equivalence.harmonic");
    if (face_cols != 3 || normal_cols != 3 || er_cols != 3 || ei_cols != 3 ||
        hr_cols != 3 || hi_cols != 3 || obs_cols != 3 ||
        normal_rows != face_rows || areas.size() != face_rows ||
        er_rows != face_rows || ei_rows != face_rows ||
        hr_rows != face_rows || hi_rows != face_rows)
        BadArgument("all equivalence face arrays must have matching count-by-3 shapes");
    const int n_threads = nrhs == 11 ? NonnegativeInteger(prhs[10], "n_threads") : 0;
    std::vector<double> E_re_out(obs_rows * 3, 0.0), E_im_out(obs_rows * 3, 0.0);
    std::vector<double> H_re_out(obs_rows * 3, 0.0), H_im_out(obs_rows * 3, 0.0);
    radia::eqsrc::EvaluateHarmonic(
        centroids.data(), normals.data(), areas.data(), E_re.data(), E_im.data(),
        H_re.data(), H_im.data(), MatrixDimension(face_rows, "centroids"),
        obs.data(), MatrixDimension(obs_rows, "obs"), omega,
        E_re_out.data(), E_im_out.data(), H_re_out.data(), H_im_out.data(), n_threads);
    plhs[0] = RealMatrixOutput(E_re_out, obs_rows, 3);
    plhs[1] = RealMatrixOutput(E_im_out, obs_rows, 3);
    plhs[2] = RealMatrixOutput(H_re_out, obs_rows, 3);
    plhs[3] = RealMatrixOutput(H_im_out, obs_rows, 3);
}

mxArray* HLUKindBreakdown(const long* values) {
    static const char* fields[] = {
        "internalxinternal", "internalxrk", "internalxdense",
        "rkxinternal", "rkxrk", "rkxdense",
        "densexinternal", "densexrk", "densexdense"};
    mxArray* result = mxCreateStructMatrix(1, 1, 9, fields);
    for (int i = 0; i < 9; ++i)
        mxSetField(result, 0, fields[i], mxCreateDoubleScalar(
            static_cast<double>(values[i])));
    return result;
}

void RadiaHLU(const std::string& command, int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[]) {
    if (command == "hlu.set_trunc_tol") {
        CheckArity(nrhs, 2, nlhs, 0,
            "radia_mex('hlu.set_trunc_tol', tol)");
        cHACApK_hlu_set_trunc_tol(Scalar(prhs[1], "tol"));
    } else if (command == "hlu.get_trunc_tol") {
        CheckArity(nrhs, 1, nlhs, 1,
            "tol = radia_mex('hlu.get_trunc_tol')");
        plhs[0] = mxCreateDoubleScalar(cHACApK_hlu_get_trunc_tol());
    } else if (command == "hlu.last_timings") {
        CheckArity(nrhs, 1, nlhs, 1,
            "stats = radia_mex('hlu.last_timings')");
        double t_decomp = 0.0, t_solve = 0.0;
        long n_lu = 0, n_gemm = 0;
        cHACApK_hlu_get_timings(&t_decomp, &t_solve, &n_lu, &n_gemm);
        const char* fields[] = {"t_decomp_sec", "t_solve_sec", "n_dense_lu", "n_dense_gemm"};
        plhs[0] = mxCreateStructMatrix(1, 1, 4, fields);
        mxSetField(plhs[0], 0, "t_decomp_sec", mxCreateDoubleScalar(t_decomp));
        mxSetField(plhs[0], 0, "t_solve_sec", mxCreateDoubleScalar(t_solve));
        mxSetField(plhs[0], 0, "n_dense_lu", mxCreateDoubleScalar(static_cast<double>(n_lu)));
        mxSetField(plhs[0], 0, "n_dense_gemm", mxCreateDoubleScalar(static_cast<double>(n_gemm)));
    } else if (command == "hlu.materialize_stats") {
        CheckArity(nrhs, 1, nlhs, 1,
            "stats = radia_mex('hlu.materialize_stats')");
        long n_calls = 0, n_elems = 0, n_internal = 0, n_leaf = 0;
        cHACApK_hlu_get_materialize_stats(&n_calls, &n_elems);
        cHACApK_hlu_get_materialize_split(&n_internal, &n_leaf);
        const char* fields[] = {"n_calls", "n_elems", "n_internal", "n_leaf"};
        plhs[0] = mxCreateStructMatrix(1, 1, 4, fields);
        mxSetField(plhs[0], 0, "n_calls", mxCreateDoubleScalar(static_cast<double>(n_calls)));
        mxSetField(plhs[0], 0, "n_elems", mxCreateDoubleScalar(static_cast<double>(n_elems)));
        mxSetField(plhs[0], 0, "n_internal", mxCreateDoubleScalar(static_cast<double>(n_internal)));
        mxSetField(plhs[0], 0, "n_leaf", mxCreateDoubleScalar(static_cast<double>(n_leaf)));
    } else if (command == "hlu.set_parallel") {
        CheckArity(nrhs, 2, nlhs, 0,
            "radia_mex('hlu.set_parallel', on)");
        cHACApK_hlu_set_parallel(Boolean(prhs[1], "on") ? 1 : 0);
    } else if (command == "hlu.get_parallel") {
        CheckArity(nrhs, 1, nlhs, 1,
            "on = radia_mex('hlu.get_parallel')");
        plhs[0] = mxCreateLogicalScalar(cHACApK_hlu_get_parallel() != 0);
    } else if (command == "hlu.set_par_cutoff") {
        CheckArity(nrhs, 2, nlhs, 0,
            "radia_mex('hlu.set_par_cutoff', cutoff)");
        cHACApK_hlu_set_par_cutoff(NonnegativeLong(prhs[1], "cutoff"));
    } else if (command == "hlu.max_threads") {
        CheckArity(nrhs, 1, nlhs, 1,
            "n = radia_mex('hlu.max_threads')");
        plhs[0] = mxCreateDoubleScalar(static_cast<double>(chacapk_max_threads()));
    } else if (command == "hlu.set_accum_cap") {
        CheckArity(nrhs, 2, nlhs, 0,
            "radia_mex('hlu.set_accum_cap', cap)");
        cHACApK_hlu_set_accum_cap(NonnegativeInteger(prhs[1], "cap"));
    } else if (command == "hlu.get_accum_cap") {
        CheckArity(nrhs, 1, nlhs, 1,
            "cap = radia_mex('hlu.get_accum_cap')");
        plhs[0] = mxCreateDoubleScalar(static_cast<double>(cHACApK_hlu_get_accum_cap()));
    } else if (command == "hlu.mixed_breakdown") {
        CheckArity(nrhs, 1, nlhs, 1,
            "breakdown = radia_mex('hlu.mixed_breakdown')");
        long addmul[9] = {0}, lln[9] = {0}, run[9] = {0};
        cHACApK_hlu_get_mixed_breakdown(addmul, lln, run);
        const char* fields[] = {"addmul", "lln", "run"};
        plhs[0] = mxCreateStructMatrix(1, 1, 3, fields);
        mxSetField(plhs[0], 0, "addmul", HLUKindBreakdown(addmul));
        mxSetField(plhs[0], 0, "lln", HLUKindBreakdown(lln));
        mxSetField(plhs[0], 0, "run", HLUKindBreakdown(run));
    } else if (command == "hlu.cluster_strategy" ||
               command == "radia.GetClusterStrategy") {
        CheckArity(nrhs, 1, nlhs, 1,
            "strategy = radia_mex('hlu.cluster_strategy')");
        plhs[0] = mxCreateDoubleScalar(static_cast<double>(cHACApK_get_cluster_strategy()));
    } else if (command == "hlu.self_test") {
        if ((nrhs != 1 && nrhs != 3) || nlhs != 1)
            BadArgument("usage: error = radia_mex('hlu.self_test' [, depth, n_per_block])");
        const int depth = nrhs == 3 ? PositiveInteger(prhs[1], "depth") : 1;
        const int n_per_block = nrhs == 3 ? PositiveInteger(prhs[2], "n_per_block") : 100;
        plhs[0] = mxCreateDoubleScalar(cHACApK_harith_self_test(depth, n_per_block));
    } else if (command == "hlu.self_test_rk") {
        if ((nrhs != 1 && nrhs != 3) || nlhs != 1)
            BadArgument("usage: error = radia_mex('hlu.self_test_rk' [, n_per_block, rk_rank])");
        const int n_per_block = nrhs == 3 ? PositiveInteger(prhs[1], "n_per_block") : 100;
        const int rk_rank = nrhs == 3 ? PositiveInteger(prhs[2], "rk_rank") : 5;
        plhs[0] = mxCreateDoubleScalar(cHACApK_harith_self_test_rk(n_per_block, rk_rank));
    } else if (command == "hlu.self_test_addmul_rkrk") {
        if ((nrhs != 1 && nrhs != 7) || nlhs != 1)
            BadArgument("usage: error = radia_mex('hlu.self_test_addmul_rkrk' [, m,n,inner,kA,kB,kC])");
        const int m = nrhs == 7 ? PositiveInteger(prhs[1], "m") : 64;
        const int n = nrhs == 7 ? PositiveInteger(prhs[2], "n") : 64;
        const int inner = nrhs == 7 ? PositiveInteger(prhs[3], "inner") : 64;
        const int kA = nrhs == 7 ? PositiveInteger(prhs[4], "kA") : 5;
        const int kB = nrhs == 7 ? PositiveInteger(prhs[5], "kB") : 5;
        const int kC = nrhs == 7 ? PositiveInteger(prhs[6], "kC") : 5;
        plhs[0] = mxCreateDoubleScalar(cHACApK_harith_self_test_addmul_rkrk(
            m, n, inner, kA, kB, kC));
    } else if (command == "hlu.self_test_radia_exact_with_matrix") {
        CheckArity(nrhs, 3, nlhs, 1,
            "error = radia_mex('hlu.self_test_radia_exact_with_matrix', A, b)");
        std::size_t rows = 0, cols = 0;
        auto a = RealMatrix(prhs[1], rows, cols, "A");
        auto b = RealVector(prhs[2], "b");
        if (rows != 162 || cols != 162 || b.size() != 162)
            BadArgument("A must be 162-by-162 and b must have 162 entries");
        std::vector<double> a_colmajor(162 * 162, 0.0);
        for (int i = 0; i < 162; ++i)
            for (int j = 0; j < 162; ++j)
                a_colmajor[static_cast<std::size_t>(i) + static_cast<std::size_t>(j) * 162] =
                    a[static_cast<std::size_t>(i) * 162 + j];
        plhs[0] = mxCreateDoubleScalar(cHACApK_harith_self_test_radia_exact_with_matrix(
            a_colmajor.data(), b.data()));
    } else if (command == "hlu.self_test_radia_exact_diag") {
        if ((nrhs != 1 && nrhs != 2) || nlhs != 1)
            BadArgument("usage: error = radia_mex('hlu.self_test_radia_exact_diag' [, diag_boost])");
        plhs[0] = mxCreateDoubleScalar(cHACApK_harith_self_test_radia_exact_diag(
            nrhs == 2 ? Scalar(prhs[1], "diag_boost") : 2.0));
    } else if (command == "hlu.self_test_radia_exact") {
        CheckArity(nrhs, 1, nlhs, 1,
            "error = radia_mex('hlu.self_test_radia_exact')");
        plhs[0] = mxCreateDoubleScalar(cHACApK_harith_self_test_radia_exact());
    } else if (command == "hlu.self_test_depth3_asymmetric") {
        if ((nrhs != 1 && nrhs != 2) || nlhs != 1)
            BadArgument("usage: error = radia_mex('hlu.self_test_depth3_asymmetric' [, nb_tiny])");
        plhs[0] = mxCreateDoubleScalar(cHACApK_harith_self_test_depth3_asymmetric(
            nrhs == 2 ? PositiveInteger(prhs[1], "nb_tiny") : 3));
    } else if (command == "hlu.self_test_mixed_sibling_via_conversion") {
        if ((nrhs != 1 && nrhs != 2) || nlhs != 1)
            BadArgument("usage: error = radia_mex('hlu.self_test_mixed_sibling_via_conversion' [, nb_small])");
        plhs[0] = mxCreateDoubleScalar(cHACApK_harith_self_test_mixed_sibling_via_conversion(
            nrhs == 2 ? PositiveInteger(prhs[1], "nb_small") : 5));
    } else if (command == "hlu.self_test_mixed_sibling_nonuniform") {
        if ((nrhs != 1 && nrhs != 5) || nlhs != 1)
            BadArgument("usage: error = radia_mex('hlu.self_test_mixed_sibling_nonuniform' [, n1,n2,m1,m3])");
        const int n1 = nrhs == 5 ? PositiveInteger(prhs[1], "n1") : 5;
        const int n2 = nrhs == 5 ? PositiveInteger(prhs[2], "n2") : 7;
        const int m1 = nrhs == 5 ? PositiveInteger(prhs[3], "m1") : 2;
        const int m3 = nrhs == 5 ? PositiveInteger(prhs[4], "m3") : 3;
        plhs[0] = mxCreateDoubleScalar(cHACApK_harith_self_test_mixed_sibling_nonuniform(
            n1, n2, m1, m3));
    } else if (command == "hlu.self_test_mixed_sibling") {
        if ((nrhs != 1 && nrhs != 2) || nlhs != 1)
            BadArgument("usage: error = radia_mex('hlu.self_test_mixed_sibling' [, nb_small])");
        plhs[0] = mxCreateDoubleScalar(cHACApK_harith_self_test_mixed_sibling(
            nrhs == 2 ? PositiveInteger(prhs[1], "nb_small") : 5));
    } else {
        if ((nrhs != 1 && nrhs != 3) || nlhs != 1)
            BadArgument("usage: error = radia_mex('hlu.self_test_rk_deep' [, n_per_block, rk_rank])");
        const int n_per_block = nrhs == 3 ? PositiveInteger(prhs[1], "n_per_block") : 100;
        const int rk_rank = nrhs == 3 ? PositiveInteger(prhs[2], "rk_rank") : 5;
        plhs[0] = mxCreateDoubleScalar(cHACApK_harith_self_test_rk_deep(
            n_per_block, rk_rank));
    }
}

const mxArray* SolverConfigField(const mxArray* config, const char* name) {
    return mxGetField(config, 0, name);
}

void RadiaSolverConfig(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 0,
        "radia_mex('radia.SolverConfig', options_struct)");
    const mxArray* config = prhs[1];
    if (!mxIsStruct(config) || mxGetNumberOfElements(config) != 1)
        BadArgument("options_struct must be a scalar MATLAB struct");
    static const char* allowed[] = {
        "relax_param", "newton_method", "newton_damping",
        "newton_damping_max_iter", "newton_damping_min_omega",
        "b_input_newton", "b_input_hantila", "hantila_alpha",
        "hantila_relax", "keep_magnetization"};
    const int allowed_count = static_cast<int>(sizeof(allowed) / sizeof(allowed[0]));
    for (int i = 0; i < mxGetNumberOfFields(config); ++i) {
        const char* field = mxGetFieldNameByNumber(config, i);
        bool known = false;
        for (int j = 0; j < allowed_count; ++j)
            known = known || std::string(field) == allowed[j];
        if (!known)
            BadArgument(std::string("unknown SolverConfig option: ") + field);
    }
    int result = 0;
    if (const mxArray* value = SolverConfigField(config, "relax_param"))
        CheckRadia(RadSetRelaxParam(&result, Scalar(value, "relax_param")));
    if (const mxArray* value = SolverConfigField(config, "newton_method"))
        CheckRadia(RadSetNewtonMethod(&result, Boolean(value, "newton_method") ? 1 : 0));
    if (SolverConfigField(config, "newton_damping") ||
        SolverConfigField(config, "newton_damping_max_iter") ||
        SolverConfigField(config, "newton_damping_min_omega")) {
        const bool enabled = SolverConfigField(config, "newton_damping")
            ? Boolean(SolverConfigField(config, "newton_damping"), "newton_damping") : true;
        const int max_iter = SolverConfigField(config, "newton_damping_max_iter")
            ? PositiveInteger(SolverConfigField(config, "newton_damping_max_iter"),
                              "newton_damping_max_iter") : 5;
        const double min_omega = SolverConfigField(config, "newton_damping_min_omega")
            ? Scalar(SolverConfigField(config, "newton_damping_min_omega"),
                     "newton_damping_min_omega") : 0.01;
        CheckRadia(RadSetNewtonDamping(&result, enabled ? 1 : 0, max_iter, min_omega));
    }
    if (const mxArray* value = SolverConfigField(config, "b_input_newton"))
        CheckRadia(RadSetBInputNewton(&result, Boolean(value, "b_input_newton") ? 1 : 0));
    if (const mxArray* value = SolverConfigField(config, "b_input_hantila"))
        CheckRadia(RadSetBInputHantila(&result, Boolean(value, "b_input_hantila") ? 1 : 0));
    if (const mxArray* value = SolverConfigField(config, "hantila_alpha"))
        CheckRadia(RadSetHantilaAlpha(&result, Scalar(value, "hantila_alpha")));
    if (const mxArray* value = SolverConfigField(config, "hantila_relax"))
        CheckRadia(RadSetHantilaRelax(&result, Scalar(value, "hantila_relax")));
    if (const mxArray* value = SolverConfigField(config, "keep_magnetization"))
        CheckRadia(RadSetKeepMagnetization(&result,
                                           Boolean(value, "keep_magnetization") ? 1 : 0));
}

void RadiaGetSolverConfig(int nlhs, mxArray* plhs[], int nrhs) {
    CheckArity(nrhs, 1, nlhs, 1,
        "config = radia_mex('radia.GetSolverConfig')");
    const char* fields[] = {
        "relax_param", "keep_magnetization", "newton_method", "newton_damping",
        "newton_damping_max_iter", "newton_damping_min_omega", "b_input_newton",
        "b_input_hantila", "hantila_alpha", "hantila_relax"};
    plhs[0] = mxCreateStructMatrix(1, 1, 10, fields);
    double relax = 0.0;
    int keep = 0, newton = 0, damping = 0, max_iter = 0;
    double min_omega = 0.0, alpha = 0.0, hantila_relax = 0.0;
    int b_input_newton = 0, b_input_hantila = 0;
    CheckRadia(RadGetRelaxParam(&relax));
    CheckRadia(RadGetKeepMagnetization(&keep));
    CheckRadia(RadGetNewtonMethod(&newton));
    CheckRadia(RadGetNewtonDampingStats(&damping, &max_iter, &min_omega));
    CheckRadia(RadGetBInputNewton(&b_input_newton));
    CheckRadia(RadGetBInputHantila(&b_input_hantila));
    CheckRadia(RadGetHantilaAlpha(&alpha));
    CheckRadia(RadGetHantilaRelax(&hantila_relax));
    mxSetField(plhs[0], 0, "relax_param", mxCreateDoubleScalar(relax));
    mxSetField(plhs[0], 0, "keep_magnetization", mxCreateLogicalScalar(keep != 0));
    mxSetField(plhs[0], 0, "newton_method", mxCreateLogicalScalar(newton != 0));
    mxSetField(plhs[0], 0, "newton_damping", mxCreateLogicalScalar(damping != 0));
    mxSetField(plhs[0], 0, "newton_damping_max_iter", mxCreateDoubleScalar(max_iter));
    mxSetField(plhs[0], 0, "newton_damping_min_omega", mxCreateDoubleScalar(min_omega));
    mxSetField(plhs[0], 0, "b_input_newton", mxCreateLogicalScalar(b_input_newton != 0));
    mxSetField(plhs[0], 0, "b_input_hantila", mxCreateLogicalScalar(b_input_hantila != 0));
    mxSetField(plhs[0], 0, "hantila_alpha", mxCreateDoubleScalar(alpha));
    mxSetField(plhs[0], 0, "hantila_relax", mxCreateDoubleScalar(hantila_relax));
}

void RadiaFld(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "field = radia_mex('radia.Fld', object, field_type, points)");
    const int object = PositiveInteger(prhs[1], "object");
    std::string field_type = Text(prhs[2], "field_type");
    std::transform(field_type.begin(), field_type.end(), field_type.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    std::size_t rows = 0, cols = 0;
    auto points = RealMatrix(prhs[3], rows, cols, "points");
    if (rows == 0 || cols != 3)
        BadArgument("points must have shape count-by-3");
    if (rows > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        BadArgument("too many field points");
    const int count = static_cast<int>(rows);

    if (field_type == "b" || field_type == "h") {
        std::vector<double> b(rows * 3), h(rows * 3);
        CheckRadia(RadFldBatch(b.data(), h.data(), count, points.data(), object));
        plhs[0] = RealMatrixOutput(field_type == "h" ? h : b, rows, 3);
        return;
    }
    if (field_type == "a") {
        std::vector<double> result(rows * 3);
        CheckRadia(RadFldA(result.data(), count, points.data(), object));
        plhs[0] = RealMatrixOutput(result, rows, 3);
        return;
    }
    if (field_type == "phi") {
        std::vector<double> result(rows);
        CheckRadia(RadFldPhi(result.data(), count, points.data(), object));
        plhs[0] = RealMatrixOutput(result, rows, 1);
        return;
    }

    std::vector<double> first(6, 0.0);
    int result_size = 0;
    CheckRadia(RadFld(first.data(), &result_size, object, field_type.data(),
                      points.data(), 1));
    if (result_size <= 0 || result_size > 6)
        throw std::runtime_error("Radia returned an invalid field result size");
    std::vector<double> result(rows * static_cast<std::size_t>(result_size));
    std::copy(first.begin(), first.begin() + result_size, result.begin());
    for (int i = 1; i < count; ++i) {
        int current_size = 0;
        CheckRadia(RadFld(result.data() + static_cast<std::size_t>(i) * result_size,
                          &current_size, object, field_type.data(),
                          points.data() + static_cast<std::size_t>(i) * 3, 1));
        if (current_size != result_size)
            throw std::runtime_error("field result size changed between points");
    }
    plhs[0] = RealMatrixOutput(result, rows, static_cast<std::size_t>(result_size));
}

void RadiaUtility(const std::string& command, int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[]) {
    if (command == "radia.UtiDelAll") {
        CheckArity(nrhs, 1, nlhs, 0, "radia_mex('radia.UtiDelAll')");
        int result = 0;
        CheckRadia(RadUtiDelAll(&result));
    } else if (command == "radia.UtiDel") {
        CheckArity(nrhs, 2, nlhs, 0, "radia_mex('radia.UtiDel', object)");
        int result = 0;
        CheckRadia(RadUtiDel(&result, PositiveInteger(prhs[1], "object")));
    } else if (command == "radia.ObjGeoVol") {
        CheckArity(nrhs, 2, nlhs, 1, "volume = radia_mex('radia.ObjGeoVol', object)");
        double result = 0.0;
        CheckRadia(RadObjGeoVol(&result, PositiveInteger(prhs[1], "object")));
        plhs[0] = mxCreateDoubleScalar(result);
    } else {
        CheckArity(nrhs, 2, nlhs, 1, "ndof = radia_mex('radia.ObjDegFre', object)");
        int result = 0;
        CheckRadia(RadObjDegFre(&result, PositiveInteger(prhs[1], "object")));
        plhs[0] = mxCreateDoubleScalar(result);
    }
}

void RadiaContainer(const std::string& command, int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[]) {
    if (command == "radia.ObjAddToCnt") {
        CheckArity(nrhs, 3, nlhs, 0,
            "radia_mex('radia.ObjAddToCnt', container, objects)");
        auto objects = IntegerVector(prhs[2], "objects");
        CheckRadia(RadObjAddToCnt(PositiveInteger(prhs[1], "container"),
            objects.empty() ? nullptr : objects.data(),
            static_cast<int>(objects.size())));
    } else if (command == "radia.ObjCntSize") {
        CheckArity(nrhs, 2, nlhs, 1,
            "count = radia_mex('radia.ObjCntSize', container)");
        int count = 0;
        CheckRadia(RadObjCntSize(&count, PositiveInteger(prhs[1], "container")));
        plhs[0] = mxCreateDoubleScalar(count);
    } else {
        CheckArity(nrhs, 2, nlhs, 1,
            "objects = radia_mex('radia.ObjCntStuf', container)");
        const int container = PositiveInteger(prhs[1], "container");
        int count = 0;
        CheckRadia(RadObjCntSize(&count, container));
        if (count < 0)
            throw std::runtime_error("Radia returned a negative container size");
        std::vector<int> objects(static_cast<std::size_t>(count));
        CheckRadia(RadObjCntStuf(objects.empty() ? nullptr : objects.data(), container));
        std::vector<double> result(objects.begin(), objects.end());
        plhs[0] = RealRow(result);
    }
}

void RadiaObjectState(const std::string& command, int nlhs, mxArray* plhs[],
                      int nrhs, const mxArray* prhs[]) {
    if (command == "radia.ObjDpl") {
        if ((nrhs != 2 && nrhs != 3) || nlhs != 1)
            BadArgument("usage: copy = radia_mex('radia.ObjDpl', object [, option])");
        const std::string option = nrhs == 3 ? Text(prhs[2], "option") : "";
        std::vector<char> option_buffer(option.begin(), option.end());
        option_buffer.push_back('\0');
        int result = 0;
        CheckRadia(RadObjDpl(&result, PositiveInteger(prhs[1], "object"),
                             option_buffer.data()));
        plhs[0] = mxCreateDoubleScalar(result);
    } else if (command == "radia.ObjSetM") {
        CheckArity(nrhs, 3, nlhs, 0,
            "radia_mex('radia.ObjSetM', object, magnetization)");
        auto magnetization = FixedRealVector(prhs[2], 3, "magnetization");
        CheckRadia(RadObjSetM(PositiveInteger(prhs[1], "object"),
                              magnetization.data()));
    } else if (command == "radia.ObjScaleCur") {
        CheckArity(nrhs, 3, nlhs, 0,
            "radia_mex('radia.ObjScaleCur', object, scale)");
        CheckRadia(RadObjScaleCur(PositiveInteger(prhs[1], "object"),
                                  Scalar(prhs[2], "scale")));
    } else {
        CheckArity(nrhs, 2, nlhs, 1,
            "state = radia_mex('radia.ObjM', object)");
        const int object = PositiveInteger(prhs[1], "object");
        int ndof = 0;
        CheckRadia(RadObjDegFre(&ndof, object));
        if (ndof < 0)
            throw std::runtime_error("Radia returned a negative degree count");
        std::vector<double> raw(static_cast<std::size_t>(ndof) * 2 + 6, 0.0);
        int shape[20] = {0};
        CheckRadia(RadObjM(raw.data(), shape, object));
        const int count = shape[0] >= 3 ? shape[3] : 1;
        if (count < 0 || static_cast<std::size_t>(count) * 6 > raw.size())
            throw std::runtime_error("Radia returned an invalid magnetization shape");
        std::vector<double> centers(static_cast<std::size_t>(count) * 3);
        std::vector<double> magnetizations(static_cast<std::size_t>(count) * 3);
        for (int i = 0; i < count; ++i) {
            std::copy_n(raw.data() + static_cast<std::size_t>(i) * 6, 3,
                        centers.data() + static_cast<std::size_t>(i) * 3);
            std::copy_n(raw.data() + static_cast<std::size_t>(i) * 6 + 3, 3,
                        magnetizations.data() + static_cast<std::size_t>(i) * 3);
        }
        const char* fields[] = {"center", "magnetization"};
        plhs[0] = mxCreateStructMatrix(1, 1, 2, fields);
        mxSetField(plhs[0], 0, "center", RealMatrixOutput(centers, count, 3));
        mxSetField(plhs[0], 0, "magnetization",
                   RealMatrixOutput(magnetizations, count, 3));
    }
}

void RadiaTransform(const std::string& command, int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[]) {
    int result = 0;
    if (command == "radia.TrfTrsl") {
        CheckArity(nrhs, 2, nlhs, 1,
            "transform = radia_mex('radia.TrfTrsl', vector)");
        auto vector = FixedRealVector(prhs[1], 3, "vector");
        CheckRadia(RadTrfTrsl(&result, vector.data()));
    } else if (command == "radia.TrfRot") {
        CheckArity(nrhs, 4, nlhs, 1,
            "transform = radia_mex('radia.TrfRot', point, axis, angle)");
        auto point = FixedRealVector(prhs[1], 3, "point");
        auto axis = FixedRealVector(prhs[2], 3, "axis");
        CheckRadia(RadTrfRot(&result, point.data(), axis.data(),
                             Scalar(prhs[3], "angle")));
    } else if (command == "radia.TrfInv") {
        CheckArity(nrhs, 1, nlhs, 1, "transform = radia_mex('radia.TrfInv')");
        CheckRadia(RadTrfInv(&result));
    } else if (command == "radia.TrfCmbL" || command == "radia.TrfCmbR") {
        CheckArity(nrhs, 3, nlhs, 1,
            "transform = radia_mex('radia.TrfCmb<L|R>', original, transform)");
        const int original = PositiveInteger(prhs[1], "original_transform");
        const int transform = PositiveInteger(prhs[2], "transform");
        if (command == "radia.TrfCmbL")
            CheckRadia(RadTrfCmbL(&result, original, transform));
        else
            CheckRadia(RadTrfCmbR(&result, original, transform));
    } else {
        CheckArity(nrhs, 3, nlhs, 1,
            "object = radia_mex('radia.TrfOrnt', object, transform)");
        CheckRadia(RadTrfOrnt(&result, PositiveInteger(prhs[1], "object"),
                              PositiveInteger(prhs[2], "transform")));
    }
    plhs[0] = mxCreateDoubleScalar(result);
}

void RadiaMatPM(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "material = radia_mex('radia.MatPM', Br, Hc, easy_axis)");
    auto axis = FixedRealVector(prhs[3], 3, "easy_axis");
    int material = 0;
    CheckRadia(::RadMatPM(&material, Scalar(prhs[1], "Br"),
                          Scalar(prhs[2], "Hc"), axis.data()));
    plhs[0] = mxCreateDoubleScalar(material);
}

void RadiaVersion(int nlhs, mxArray* plhs[], int nrhs) {
    CheckArity(nrhs, 1, nlhs, 1, "version = radia_mex('radia.UtiVer')");
    double version = 0.0;
    CheckRadia(RadUtiVer(&version));
    plhs[0] = mxCreateDoubleScalar(version);
}

bool DispatchChargeGramCreate(const std::string& command, int nlhs,
                              mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (command == "hacapk.charge_gram.create_monopole") {
        ChargeGramCreateMonopole(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "hacapk.charge_gram.create_sampled_laplace" ||
        command == "hacapk.charge_gram.create_sampled_planar_log") {
        ChargeGramCreateSampled(command, nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "hacapk.charge_gram.create_local_polynomials") {
        ChargeGramCreateLocalPolynomials(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "hacapk.charge_gram.create_analytic_tet") {
        ChargeGramCreateAnalyticTet(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "hacapk.charge_gram.create_analytic_polytope") {
        ChargeGramCreateAnalyticPolytope(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "hacapk.charge_gram.create_high_order_tet") {
        ChargeGramCreateHighOrderTet(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "hacapk.charge_gram.create_curved_high_order_tet") {
        ChargeGramCreateCurvedHighOrderTet(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "hacapk.charge_gram.create_hex") {
        ChargeGramCreateHex(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "hacapk.charge_gram.create_wedge") {
        ChargeGramCreateWedge(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "hacapk.charge_gram.create_planar_2d") {
        ChargeGramCreatePlanar2D(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "hacapk.charge_gram.create_curved_polytope") {
        ChargeGramCreateCurvedPolytope(nlhs, plhs, nrhs, prhs);
        return true;
    }
    return false;
}

void Dispatch(const std::string& command, int nlhs, mxArray* plhs[], int nrhs,
              const mxArray* prhs[]) {
    if (command == "ngsolve.mesh.create") {
        NGSolveMeshCreate(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.mesh.info") {
        NGSolveMeshInfo(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.mesh.set_deformation") {
        NGSolveMeshSetDeformation(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.mesh.unset_deformation") {
        NGSolveMeshUnsetDeformation(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.mesh.trafo_quality") {
        NGSolveMeshTrafoQuality(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.mesh.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('ngsolve.mesh.destroy', handle)");
        DestroyMesh(Handle(prhs[1]));
        return;
    }
    if (command == "ngsolve.grid_function.from_fespace") {
        NGSolveGridFunctionFromFESpace(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.linear_form.create") {
        NGSolveLinearFormCreate(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.linear_form.create_from_coefficient") {
        NGSolveLinearFormCreateFromCoefficient(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.linear_form.create_boundary_from_coefficient") {
        NGSolveLinearFormCreateBoundaryFromCoefficient(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.linear_form.info") {
        NGSolveLinearFormInfo(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.linear_form.vector") {
        NGSolveLinearFormVector(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.linear_form.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('ngsolve.linear_form.destroy', handle)");
        DestroyLinearForm(Handle(prhs[1]));
        return;
    }
    if (command == "ngsolve.fespace.create") {
        NGSolveFESpaceCreate(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.fespace.info") {
        NGSolveFESpaceInfo(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.fespace.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('ngsolve.fespace.destroy', handle)");
        DestroyFESpace(Handle(prhs[1]));
        return;
    }
    if (command == "ngsolve.bilinear_form.create") {
        NGSolveBilinearFormCreate(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.bilinear_form.create_from_coefficient") {
        NGSolveBilinearFormCreateFromCoefficient(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.bilinear_form.create_boundary_from_coefficient") {
        NGSolveBilinearFormCreateBoundaryFromCoefficient(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.bilinear_form.info") {
        NGSolveBilinearFormInfo(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.bilinear_form.matrix") {
        NGSolveBilinearFormMatrix(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.bilinear_form.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('ngsolve.bilinear_form.destroy', handle)");
        DestroyBilinearForm(Handle(prhs[1]));
        return;
    }
    if (command == "ngsolve.matrix.info") {
        NGSolveMatrixInfo(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.matrix.values") {
        NGSolveMatrixValues(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.matrix.vector") {
        NGSolveMatrixVector(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.matrix.matvec") {
        NGSolveMatrixMatVec(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.matrix.inverse") {
        NGSolveMatrixInverse(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.matrix.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('ngsolve.matrix.destroy', handle)");
        DestroyMatrix(Handle(prhs[1]));
        return;
    }
    if (command == "ngsolve.solver.create") {
        NGSolveSolverCreate(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.solver.info") {
        NGSolveSolverInfo(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.solver.solve") {
        NGSolveSolverSolve(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.solver.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('ngsolve.solver.destroy', handle)");
        DestroySolver(Handle(prhs[1]));
        return;
    }
    // Keep the Simulink state-space commands outside the legacy dispatch
    // chain.  The latter is intentionally broad for API compatibility and
    // exceeds MSVC's nested-block limit when new commands are appended.
    if (command == "simulink.state_space.create") {
        SimulinkStateSpaceCreate(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "simulink.state_space.info") {
        SimulinkStateSpaceInfo(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "simulink.state_space.step") {
        SimulinkStateSpaceStep(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "simulink.state_space.reset") {
        SimulinkStateSpaceReset(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "simulink.state_space.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('simulink.state_space.destroy', handle)");
        DestroyStateSpace(Handle(prhs[1]));
        return;
    }
    // New parity commands use flat early-return dispatch.  Keep additions out
    // of the legacy compatibility chain, which is already at MSVC's nested
    // statement limit.
    if (command == "hdiv.affine_cell_self_energy_shape_derivative") {
        AffineCellSelfEnergyShapeDerivative(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "stream.aca_tsvd") {
        StreamACATSVD(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "hdiv.field_evaluator.as_coefficient") {
        HDivFieldAsCoefficient(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "hdiv.planar_evaluator.as_coefficient") {
        PlanarFieldAsCoefficient(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.radia_field.create") {
        NGSolveRadiaFieldCreate(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.radia_field.info") {
        NGSolveRadiaFieldInfo(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.radia_field.prepare_cache") {
        NGSolveRadiaFieldPrepareCache(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.radia_field.clear_cache") {
        NGSolveRadiaFieldClearCache(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.radia_field.cache_stats") {
        NGSolveRadiaFieldCacheStats(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.radia_field.as_voxel_coefficient") {
        NGSolveRadiaFieldAsVoxelCoefficient(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.matrix.projected_create") {
        NGSolveProjectedMatrixCreate(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.matrix.reduced_block_create") {
        NGSolveReducedBlockMatrixCreate(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.matrix.diagonal_preconditioner") {
        NGSolveMatrixDiagonalPreconditioner(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "ngsolve.matrix.term_count") {
        NGSolveMatrixTermCount(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "hacapk.charge_gram.demag_matrix") {
        ChargeGramDemagMatrix(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command ==
        "hacapk.charge_gram.hex_face_self_block_directional_derivative") {
        ChargeGramHexFaceSelfBlockDirectionalDerivative(
            nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "hacapk.charge_gram.hex_directional_derivative") {
        ChargeGramHexDirectionalDerivative(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command ==
        "hacapk.charge_gram.tet_volume_self_block_directional_derivative") {
        ChargeGramTetVolumeSelfBlockDirectionalDerivative(
            nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command ==
        "hacapk.charge_gram.tet_face_self_block_directional_derivative") {
        ChargeGramTetFaceSelfBlockDirectionalDerivative(
            nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "hacapk.charge_gram.tet_directional_derivative") {
        ChargeGramTetDirectionalDerivative(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command ==
        "hacapk.charge_gram.tet_charge_map_row_directional_rates") {
        ChargeGramTetChargeMapRowDirectionalRates(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command ==
        "hacapk.charge_gram.wedge_volume_self_block_directional_derivative") {
        ChargeGramWedgeVolumeSelfBlockDirectionalDerivative(
            nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command ==
        "hacapk.charge_gram.wedge_face_self_block_directional_derivative") {
        ChargeGramWedgeFaceSelfBlockDirectionalDerivative(
            nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "hacapk.charge_gram.wedge_directional_derivative") {
        ChargeGramWedgeDirectionalDerivative(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "hacapk.charge_gram.directional_derivative_operator") {
        ChargeGramDirectionalDerivativeOperatorCreate(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "hacapk.charge_gram.directional_derivative_contractions") {
        ChargeGramDirectionalDerivativeContractions(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "hacapk.charge_gram_derivative.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
            "radia_mex('hacapk.charge_gram_derivative.destroy', handle)");
        DestroyChargeGramDerivative(Handle(prhs[1]));
        return;
    }
    if (command == "hacapk.charge_gram_derivative.info") {
        ChargeGramDerivativeInfo(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "hacapk.charge_gram_derivative.entry") {
        ChargeGramDerivativeEntry(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "hacapk.charge_gram_derivative.matvec_sym") {
        ChargeGramDerivativeMatVecSym(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "api.info")
        ApiInfo(nlhs, plhs, nrhs);
    else if (command == "api.commands") {
        CheckArity(nrhs, 1, nlhs, 1, "commands = radia_mex('api.commands')");
        plhs[0] = Commands();
    } else if (command == "taskmanager.probe")
        TaskProbe(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.space_info")
        SpaceInfo(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.matrix_dump")
        NGSolveMatrixDump(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.coefficient_function.constant_create")
        NGSolveCoefficientConstantCreate(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.coefficient_function.add" ||
             command == "ngsolve.coefficient_function.subtract" ||
             command == "ngsolve.coefficient_function.multiply")
        NGSolveCoefficientBinary(command, nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.coefficient_function.scale")
        NGSolveCoefficientScale(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.coefficient_function.info")
        NGSolveCoefficientInfo(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.coefficient_function.evaluate")
        NGSolveCoefficientEvaluate(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.coefficient_function.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('ngsolve.coefficient_function.destroy', handle)");
        DestroyCoefficient(Handle(prhs[1]));
    } else if (command == "ngsolve.grid_function.create")
        NGSolveGridFunctionCreate(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.grid_function.info")
        NGSolveGridFunctionInfo(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.grid_function.vector")
        NGSolveGridFunctionVector(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.grid_function.set_vector")
        NGSolveGridFunctionSetVector(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.grid_function.interpolate")
        NGSolveGridFunctionInterpolate(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.grid_function.as_coefficient")
        NGSolveGridFunctionAsCoefficient(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.grid_function.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('ngsolve.grid_function.destroy', handle)");
        DestroyGridFunction(Handle(prhs[1]));
    }
    else if (command == "ngsolve.grid_function.vector_handle")
        NGSolveGridFunctionVectorHandle(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.vector.info")
        NGSolveVectorInfo(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.vector.copy")
        NGSolveVectorCopy(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.vector.set_zero")
        NGSolveVectorSetZero(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.vector.scale")
        NGSolveVectorScale(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.vector.axpy")
        NGSolveVectorAxpy(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.vector.dot")
        NGSolveVectorDot(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.vector.norm")
        NGSolveVectorNorm(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.vector.values")
        NGSolveVectorValues(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.vector.set_values")
        NGSolveVectorSetValues(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.vector.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('ngsolve.vector.destroy', handle)");
        DestroyVector(Handle(prhs[1]));
    }
    else if (command == "hcurl.eddy_cln.native_basis")
        HCurlEddyCLNNativeBasis(nlhs, plhs, nrhs, prhs);
    else if (command == "energy_stop.create")
        EnergyCreate(nlhs, plhs, nrhs, prhs);
    else if (command == "energy_stop.destroy") {
        CheckArity(nrhs, 2, nlhs, 0, "radia_mex('energy_stop.destroy', handle)");
        Destroy(Handle(prhs[1]));
    } else if (command == "energy_stop.info")
        EnergyInfo(nlhs, plhs, nrhs, prhs);
    else if (command == "energy_stop.state0")
        EnergyState0(nlhs, plhs, nrhs, prhs);
    else if (command == "energy_stop.forward")
        EnergyBatch(BatchOperation::Forward, nlhs, plhs, nrhs, prhs);
    else if (command == "energy_stop.commit")
        EnergyBatch(BatchOperation::Commit, nlhs, plhs, nrhs, prhs);
    else if (command == "energy_stop.stored_energy")
        EnergyBatch(BatchOperation::StoredEnergy, nlhs, plhs, nrhs, prhs);
    else if (command == "hybrid_vim.solve")
        HybridSolve(nlhs, plhs, nrhs, prhs);
    else if (command == "hybrid_vim.schur")
        HybridSchur(nlhs, plhs, nrhs, prhs);
    else if (command == "hybrid_vim.skin_impedance" ||
             command == "hybrid_vim.sibc_admittance_tail" ||
             command == "hybrid_vim.sibc_termination_impedance" ||
             command == "hybrid_vim.sibc_termination_admittance")
        HybridScalar(command, nlhs, plhs, nrhs, prhs);
    else if (command == "cln.lanczos")
        CLNLanczos(nlhs, plhs, nrhs, prhs);
    else if (command == "cln.build_tridiagonal")
        CLNBuildTridiagonal(nlhs, plhs, nrhs, prhs);
    else if (command == "cln.impedance")
        CLNImpedance(nlhs, plhs, nrhs, prhs);
    else if (command == "cln.impedance_sweep")
        CLNImpedanceSweep(nlhs, plhs, nrhs, prhs);
    else if (command == "cln.transform_coupling")
        CLNTransformCoupling(nlhs, plhs, nrhs, prhs);
    else if (command == "cln.transform_port")
        CLNTransformPort(nlhs, plhs, nrhs, prhs);
    else if (command == "cln.aca_compress")
        CLNACACompress(nlhs, plhs, nrhs, prhs);
    else if (command == "evrs.tmethod")
        EVRSTMethodAlgebra(nlhs, plhs, nrhs, prhs);
    else if (command == "hcurl.tet_reduced_gram")
        TetHCurlReducedGram(nlhs, plhs, nrhs, prhs);
    else if (command == "biot_savart.h_segments_complex" ||
             command == "biot_savart.a_segments_complex")
        BiotSavartSegments(command, nlhs, plhs, nrhs, prhs);
    else if (command == "biot_savart.a_triangles_complex" ||
             command == "biot_savart.b_triangles_complex")
        BiotSavartTriangles(command, nlhs, plhs, nrhs, prhs);
    else if (command == "bem.assemble_sldl" ||
             command == "bem.assemble_sldl_p2")
        BemGalerkin(command, nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.bem.create")
        BEMCreate(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.bem.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('hacapk.bem.destroy', handle)");
        DestroyBEM(Handle(prhs[1]));
    } else if (command == "hacapk.bem.build")
        BEMBuild(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.bem.matvec")
        BEMMatVec(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.bem.info")
        BEMInfo(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.peec.create")
        PEECCreate(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.peec.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('hacapk.peec.destroy', handle)");
        DestroyPEEC(Handle(prhs[1]));
    } else if (command == "hacapk.peec.build")
        PEECBuild(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.peec.matvec")
        PEECMatVec(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.peec.info")
        PEECInfo(nlhs, plhs, nrhs, prhs);
    else if (DispatchChargeGramCreate(command, nlhs, plhs, nrhs, prhs)) {
    }
    else if (command == "hacapk.charge_gram.destroy") {
        CheckArity(nrhs, 2, nlhs, 0,
                   "radia_mex('hacapk.charge_gram.destroy', handle)");
        DestroyChargeGram(Handle(prhs[1]));
    } else if (command == "hacapk.charge_gram.build")
        ChargeGramBuild(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.matvec" ||
             command == "hacapk.charge_gram.matvec_transpose" ||
             command == "hacapk.charge_gram.matvec_sym")
        ChargeGramMatVec(command, nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.entry")
        ChargeGramEntry(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.hex_volume_self_block_directional_derivative")
        ChargeGramHexVolumeSelfBlockDirectionalDerivative(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.info")
        ChargeGramInfo(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.hex_state_check")
        ChargeGramHexStateCheck(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.hex_stored_nodes")
        ChargeGramHexStoredNodes(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.hex_state_breakdown")
        ChargeGramHexStateBreakdown(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.configure_charge_map")
        ChargeGramConfigureChargeMap(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.configure_vector_charge_map")
        ChargeGramConfigureVectorChargeMap(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.configure_mass_matrix" ||
             command == "hacapk.charge_gram.configure_geometry_mass_matrix")
        ChargeGramConfigureMass(command, nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.configure_mass_matrix_ngsolve" ||
             command == "hacapk.charge_gram.configure_geometry_mass_matrix_ngsolve")
        ChargeGramConfigureMassNGSolve(command, nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.restore_geometry_mass_matrix")
        ChargeGramRestoreGeometryMassMatrix(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.operator_info")
        ChargeGramOperatorInfo(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.demag_apply" ||
             command == "hacapk.charge_gram.geometry_mass_apply" ||
             command == "hacapk.charge_gram.mass_riesz")
        ChargeGramConfiguredApply(command, nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.solve_configured_linear_material" ||
             command == "hacapk.charge_gram.solve_configured_linear_material_auto_prec")
        ChargeGramSolveConfigured(command, nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.create_field_evaluator")
        ChargeGramCreateFieldEvaluator(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.create_planar_field_evaluator")
        ChargeGramCreatePlanarFieldEvaluator(nlhs, plhs, nrhs, prhs);
    else if (command == "hacapk.charge_gram.stats")
        ChargeGramStats(nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.field_evaluator.from_tet")
        HDivFieldFromTet(nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.field_evaluator.from_cloud")
        HDivFieldFromCloud(nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.field_evaluator.from_curved_tet")
        HDivFieldFromCurvedTet(nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.field_evaluator.destroy")
        HDivFieldDestroy(nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.field_evaluator.field")
        HDivFieldField(nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.field_evaluator.candidate_algorithm")
        HDivFieldCandidateAlgorithm(nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.field_evaluator.last_algorithm")
        HDivFieldLastAlgorithm(nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.field_evaluator.stats")
        HDivFieldStats(nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.planar_evaluator.create")
        PlanarFieldCreate(nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.planar_evaluator.destroy")
        PlanarFieldDestroy(nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.planar_evaluator.field" ||
             command == "hdiv.planar_evaluator.az")
        PlanarFieldEvaluate(command, nlhs, plhs, nrhs, prhs);
    else if (command == "hdiv.planar_evaluator.stats")
        PlanarFieldStats(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.ObjHexahedron" ||
             command == "radia.ObjTetrahedron" ||
             command == "radia.ObjWedge" || command == "radia.ObjPyramid")
        RadiaPolyhedron(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.ObjThckPgn" || command == "radia.ObjCylMag" ||
             command == "radia.ObjRecCur" || command == "radia.ObjArcCur" ||
             command == "radia.ObjRaceTrk" || command == "radia.ObjFlmCur" ||
             command == "radia.ObjArcPgnMag" || command == "radia.ObjBckg")
        RadiaExtendedObject(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.ObjCnt")
        RadiaObjCnt(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.MatLin")
        RadiaMatLin(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.MatSatIsoTab")
        RadiaMatSatIsoTab(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.MatSatIsoFrm" || command == "radia.MatSatAniso" ||
             command == "radia.MatSatLamTab" || command == "radia.MatSatLamFrm" ||
             command == "radia.MatMvsH")
        RadiaExtendedMaterial(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.MatEnergyHysteresis" ||
             command == "radia.MatPlayHysteresis" ||
             command == "radia.MatHysSaveState" ||
             command == "radia.MatHysRestoreState" ||
             command == "radia.MatHysCommitState" ||
             command == "radia.MatHysGetNuRev" ||
             command == "radia.MatHysIrreversible")
        RadiaHysteresis(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.MatHysForwardBatch" ||
             command == "radia.MatHysCommitBatch")
        RadiaHysteresisBatch(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.MatApl")
        RadiaMatApl(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.Solve")
        RadiaSolve(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.SolveNonl")
        RadiaSolveNonlinear(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.GetSolveStats")
        RadiaSolveStats(nlhs, plhs, nrhs);
    else if (command == "radia.SolverConfig")
        RadiaSolverConfig(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.GetSolverConfig")
        RadiaGetSolverConfig(nlhs, plhs, nrhs);
    else if (command == "radia.BuildMatrix" || command == "radia.GetInteractMatrix" ||
             command == "radia.GetFaceGeom")
        RadiaInteraction(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.PlanarChargeField" || command == "radia.PlanarChargeAz" ||
             command == "radia.PlanarMaxwellTorqueCircle" ||
             command == "radia.PlanarMaxwellForceCircle")
        RadiaPlanar(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.AverageBInBox" ||
             command == "radia.AverageDemagTensor")
        RadiaAverageField(command, nlhs, plhs, nrhs, prhs);
    else if (command == "equivalence.static_h" ||
             command == "equivalence.harmonic")
        RadiaEquivalenceSource(command, nlhs, plhs, nrhs, prhs);
    else if (command.rfind("hlu.", 0) == 0 ||
             command == "radia.GetClusterStrategy")
        RadiaHLU(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.Fld")
        RadiaFld(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.FldFrcShpRtg" || command == "radia.FldFrc" ||
             command == "radia.FldLst" || command == "radia.FldInt" ||
             command == "radia.ObjCenFld" || command == "radia.FldCmpCrt" ||
             command == "radia.FldCmpPrc" || command == "radia.FldLenRndSw" ||
             command == "radia.FldLenTol")
        RadiaExtendedField(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.ObjGeoVol" || command == "radia.ObjDegFre" ||
             command == "radia.UtiDel" || command == "radia.UtiDelAll")
        RadiaUtility(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.ObjAddToCnt" || command == "radia.ObjCntSize" ||
             command == "radia.ObjCntStuf")
        RadiaContainer(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.ObjDpl" || command == "radia.ObjM" ||
             command == "radia.ObjSetM" || command == "radia.ObjScaleCur")
        RadiaObjectState(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.TrfTrsl" || command == "radia.TrfRot" ||
             command == "radia.TrfInv" || command == "radia.TrfCmbL" ||
             command == "radia.TrfCmbR" || command == "radia.TrfOrnt")
        RadiaTransform(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.MatPM")
        RadiaMatPM(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.UtiVer")
        RadiaVersion(nlhs, plhs, nrhs);
    else
        BadArgument("unknown radia_mex command: " + command);
}

} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    try {
        if (nrhs < 1)
            BadArgument("first input must be a command character vector");
        EnsureExitHandler();
        Dispatch(Text(prhs[0], "command"), nlhs, plhs, nrhs, prhs);
    } catch (const std::exception& exception) {
        mexErrMsgIdAndTxt("radia:mex:Exception", "%s", exception.what());
    } catch (...) {
        mexErrMsgIdAndTxt("radia:mex:UnknownException", "unknown C++ exception");
    }
}
