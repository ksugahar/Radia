#include "mex.h"
#include "optuna_cuda.hpp"
#include "optuna_runtime.hpp"
#ifndef _OPENMP
#error "Optuna shared kernels require the OpenMP-enabled host compiler"
#endif


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
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using Complex = std::complex<double>;

struct OptunaTPEHistoryTrial {
    int trial_number = -1;
    double objective = std::numeric_limits<double>::quiet_NaN();
    std::unordered_map<int, double> parameters;
};

struct OptunaRandomStateHandle {
    std::array<std::uint32_t, 624> mt{};
    std::size_t index = 624;
    bool has_gauss = false;
    double gauss = 0.0;
    std::vector<OptunaTPEHistoryTrial> tpe_history;
    // Computation cache only; MATLAB tables own the durable trial history.
    std::unordered_map<int, std::size_t> tpe_history_rows;
    std::vector<std::size_t> tpe_chronological;
};

// A failed CUDA operation must not consume the caller's accepted RNG state.
struct RandomStateCheckpoint {
    OptunaRandomStateHandle& state;
    std::array<std::uint32_t, 624> mt;
    std::size_t index;
    bool has_gauss;
    double gauss;
    bool accepted = false;
    explicit RandomStateCheckpoint(OptunaRandomStateHandle& s)
        : state(s), mt(s.mt), index(s.index), has_gauss(s.has_gauss), gauss(s.gauss) {}
    ~RandomStateCheckpoint() {
        if (!accepted) { state.mt=mt; state.index=index; state.has_gauss=has_gauss; state.gauss=gauss; }
    }
};

std::mutex registry_mutex;
std::unordered_map<std::uint64_t, std::unique_ptr<OptunaRandomStateHandle>>
    optuna_random_state_registry;
std::uint64_t next_handle = 1;
std::size_t lock_count = 0;
bool exit_handler_registered = false;


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

std::vector<Complex> ComplexTensor3(const mxArray* value, std::size_t& dim0,
                                    std::size_t& dim1, std::size_t& dim2,
                                    const char* name) {
    if (!mxIsDouble(value) || mxGetNumberOfDimensions(value) != 3)
        BadArgument(std::string(name) + " must be a double three-dimensional array");
    const mwSize* dims = mxGetDimensions(value);
    dim0 = static_cast<std::size_t>(dims[0]);
    dim1 = static_cast<std::size_t>(dims[1]);
    dim2 = static_cast<std::size_t>(dims[2]);
    std::vector<Complex> result(dim0 * dim1 * dim2);
    if (mxIsComplex(value)) {
        const auto* data = mxGetComplexDoubles(value);
        for (std::size_t i = 0; i < dim0; ++i)
            for (std::size_t j = 0; j < dim1; ++j)
                for (std::size_t k = 0; k < dim2; ++k) {
                    const auto& item = data[i + dim0 * (j + dim1 * k)];
                    result[(i * dim1 + j) * dim2 + k] =
                        Complex(item.real, item.imag);
                }
    } else {
        const double* data = mxGetDoubles(value);
        for (std::size_t i = 0; i < dim0; ++i)
            for (std::size_t j = 0; j < dim1; ++j)
                for (std::size_t k = 0; k < dim2; ++k)
                    result[(i * dim1 + j) * dim2 + k] =
                        Complex(data[i + dim0 * (j + dim1 * k)], 0.0);
    }
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

mxArray* Int32Column(const std::vector<int>& values) {
    static_assert(sizeof(int) == sizeof(std::int32_t));
    mxArray* result = mxCreateNumericMatrix(
        values.size(), 1, mxINT32_CLASS, mxREAL);
    auto* data = static_cast<std::int32_t*>(mxGetData(result));
    std::transform(values.begin(), values.end(), data,
                   [](int value) { return static_cast<std::int32_t>(value); });
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

mxArray* RealTensor3Output(const std::vector<double>& values,
                           std::size_t dim0, std::size_t dim1,
                           std::size_t dim2) {
    const mwSize dims[] = {static_cast<mwSize>(dim0),
                           static_cast<mwSize>(dim1),
                           static_cast<mwSize>(dim2)};
    mxArray* result = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double* data = mxGetDoubles(result);
    for (std::size_t i = 0; i < dim0; ++i)
        for (std::size_t j = 0; j < dim1; ++j)
            for (std::size_t k = 0; k < dim2; ++k)
                data[i + dim0 * (j + dim1 * k)] =
                    values[(i * dim1 + j) * dim2 + k];
    return result;
}

mxArray* RealTensor4Output(const std::vector<double>& values,
                           std::size_t dim0, std::size_t dim1,
                           std::size_t dim2, std::size_t dim3) {
    const mwSize dims[] = {static_cast<mwSize>(dim0),
                           static_cast<mwSize>(dim1),
                           static_cast<mwSize>(dim2),
                           static_cast<mwSize>(dim3)};
    mxArray* result = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    double* data = mxGetDoubles(result);
    for (std::size_t i = 0; i < dim0; ++i)
        for (std::size_t j = 0; j < dim1; ++j)
            for (std::size_t k = 0; k < dim2; ++k)
                for (std::size_t l = 0; l < dim3; ++l)
                    data[i + dim0 * (j + dim1 * (k + dim2 * l))] =
                        values[((i * dim1 + j) * dim2 + k) * dim3 + l];
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

template <class Registry>
void DestroyRegistered(Registry& registry, std::uint64_t handle,
                       const char* stale_handle_message) {
    bool erased = false;
    {
        std::lock_guard<std::mutex> guard(registry_mutex);
        erased = registry.erase(handle) != 0;
        if (erased)
            --lock_count;
    }
    if (!erased)
        BadArgument(stale_handle_message);

    // Keep MATLAB API calls outside the registry critical section; cleanup
    // timing must not be part of this mutex's locking contract.
    mexUnlock();
}


void Cleanup() {
    try {
        optuna_cuda::shutdown_runtime();
    } catch (const std::exception& error) {
        mexErrMsgIdAndTxt("radia:optuna:mex:LiveBuffers", "%s", error.what());
    }
    std::size_t unlock_count = 0;
    {
        std::lock_guard<std::mutex> guard(registry_mutex);
        optuna_random_state_registry.clear();
        unlock_count = lock_count;
        lock_count = 0;
    }
    while (unlock_count > 0) {
        mexUnlock();
        --unlock_count;
    }
}

bool HandleInUse(std::uint64_t handle) {
    return handle == 0 || optuna_random_state_registry.count(handle) != 0;
}

void EnsureExitHandler() {
    optuna_cuda::install_matlab_runtime();
    if (!exit_handler_registered) {
        mexAtExit(Cleanup);
        exit_handler_registered = true;
    }
}


std::uint64_t RegisterOptunaRandomState(
    std::unique_ptr<OptunaRandomStateHandle> random_state) {
    if (!random_state)
        BadArgument("Optuna random-state construction returned null");
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (HandleInUse(next_handle))
        ++next_handle;
    const std::uint64_t handle = next_handle++;
    optuna_random_state_registry.emplace(handle, std::move(random_state));
    mexLock();
    ++lock_count;
    return handle;
}

OptunaRandomStateHandle& OptunaRandomState(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = optuna_random_state_registry.find(handle);
    if (found == optuna_random_state_registry.end())
        BadArgument("invalid or stale Optuna random-state handle");
    return *found->second;
}

void DestroyOptunaRandomState(std::uint64_t handle) {
    DestroyRegistered(optuna_random_state_registry, handle,
                      "invalid or stale Optuna random-state handle");
}

mxArray* Commands() {
    static const char* names[] = {
        "api.info", "api.commands", "api.execution",
        "optuna.pareto.rank_crowding",
        "optuna.parzen.log_pdf_numerical",
        "optuna.parzen.log_pdf_categorical",
        "optuna.tpe.best_numerical",
        "optuna.tpe.best_joint",
        "optuna.tpe.best_joint_observations",
        "optuna.tpe.best_numerical_observations",
        "optuna.tpe.history.reset",
        "optuna.tpe.history.append_complete",
        "optuna.tpe.best_grouped_history",
        "optuna.sobol.points",
        "optuna.random_state.create",
        "optuna.random_state.rand",
        "optuna.random_state.randn",
        "optuna.random_state.randi",
        "optuna.random_state.randperm",
        "optuna.random_state.snapshot",
        "optuna.random_state.restore",
        "optuna.random_state.destroy"
    };
    constexpr std::size_t count = sizeof(names) / sizeof(names[0]);
    mxArray* result = mxCreateCellMatrix(1, count);
    for (std::size_t i = 0; i < count; ++i)
        mxSetCell(result, i, mxCreateString(names[i]));
    return result;
}


void OptunaParetoRankCrowding(int nlhs, mxArray* plhs[], int nrhs,
                              const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 2,
        "[rank, crowding] = radia_mex('optuna.pareto.rank_crowding', values, signs)");

    std::size_t count = 0, objective_count = 0;
    const auto values = RealMatrix(prhs[1], count, objective_count, "values");
    const auto signs = RealVector(prhs[2], "signs");
    if (objective_count == 0 || signs.size() != objective_count)
        BadArgument("signs must contain one entry per objective");
    for (double sign : signs)
        if (!std::isfinite(sign) || (sign != 1.0 && sign != -1.0))
            BadArgument("signs entries must be +1 (minimize) or -1 (maximize)");

    std::vector<double> normalized(values.size());
    for (std::size_t i = 0; i < count; ++i)
        for (std::size_t j = 0; j < objective_count; ++j) {
            const double value = values[i * objective_count + j];
            if (!std::isfinite(value))
                BadArgument("values must contain finite entries");
            normalized[i * objective_count + j] = value * signs[j];
        }

    auto dominates = [&](std::size_t left, std::size_t right) {
        bool strictly_better = false;
        for (std::size_t objective = 0; objective < objective_count; ++objective) {
            const double lhs = normalized[left * objective_count + objective];
            const double rhs = normalized[right * objective_count + objective];
            if (lhs > rhs)
                return false;
            strictly_better = strictly_better || lhs < rhs;
        }
        return strictly_better;
    };

    std::vector<std::vector<std::size_t>> dominated(count);
    std::vector<std::size_t> domination_count(count, 0);
    for (std::size_t left = 0; left < count; ++left)
        for (std::size_t right = left + 1; right < count; ++right) {
            if (dominates(left, right)) {
                dominated[left].push_back(right);
                ++domination_count[right];
            } else if (dominates(right, left)) {
                dominated[right].push_back(left);
                ++domination_count[left];
            }
        }

    std::vector<double> rank(count, 0.0);
    std::vector<std::vector<std::size_t>> fronts;
    fronts.emplace_back();
    for (std::size_t i = 0; i < count; ++i)
        if (domination_count[i] == 0) {
            rank[i] = 1.0;
            fronts.front().push_back(i);
        }
    for (std::size_t level = 0; level < fronts.size(); ++level) {
        std::vector<std::size_t> next;
        for (std::size_t left : fronts[level])
            for (std::size_t right : dominated[left]) {
                if (--domination_count[right] == 0) {
                    rank[right] = static_cast<double>(level + 2);
                    next.push_back(right);
                }
            }
        if (!next.empty()) {
            std::sort(next.begin(), next.end());
            fronts.push_back(std::move(next));
        }
    }

    std::vector<double> crowding(count, 0.0);
    const double infinity = std::numeric_limits<double>::infinity();
    for (const auto& front : fronts) {
        if (front.empty())
            continue;
        if (front.size() <= 2) {
            for (std::size_t index : front)
                crowding[index] = infinity;
            continue;
        }
        for (std::size_t objective = 0; objective < objective_count; ++objective) {
            auto order = front;
            std::stable_sort(order.begin(), order.end(),
                [&](std::size_t left, std::size_t right) {
                    return normalized[left * objective_count + objective] <
                           normalized[right * objective_count + objective];
                });
            crowding[order.front()] = infinity;
            crowding[order.back()] = infinity;
            const double minimum = normalized[order.front() * objective_count + objective];
            const double maximum = normalized[order.back() * objective_count + objective];
            const double span = maximum - minimum;
            if (span <= 0.0)
                continue;
            for (std::size_t position = 1; position + 1 < order.size(); ++position) {
                const std::size_t index = order[position];
                if (!std::isfinite(crowding[index]))
                    continue;
                const double previous = normalized[
                    order[position - 1] * objective_count + objective];
                const double next = normalized[
                    order[position + 1] * objective_count + objective];
                crowding[index] += (next - previous) / span;
            }
        }
    }

    plhs[0] = RealColumn(rank);
    plhs[1] = RealColumn(crowding);
}

double OptunaLogNormalMass(double lower, double upper) {
    constexpr double inverse_sqrt_two = 0.70710678118654752440;
    double mass = 0.0;
    if (upper <= 0.0) {
        mass = 0.5 * (std::erfc(-upper * inverse_sqrt_two) -
                      std::erfc(-lower * inverse_sqrt_two));
    } else if (lower >= 0.0) {
        mass = 0.5 * (std::erfc(lower * inverse_sqrt_two) -
                      std::erfc(upper * inverse_sqrt_two));
    } else {
        mass = 0.5 * (std::erf(upper * inverse_sqrt_two) -
                      std::erf(lower * inverse_sqrt_two));
    }
    return std::log(std::max(mass, std::numeric_limits<double>::min()));
}

double OptunaLogSumExp(const std::vector<double>& values) {
    const double maximum = *std::max_element(values.begin(), values.end());
    double total = 0.0;
    for (double value : values)
        total += std::exp(value - maximum);
    return maximum + std::log(total);
}

void OptunaParzenLogPdfNumerical(int nlhs, mxArray* plhs[], int nrhs,
                                 const mxArray* prhs[]) {
    CheckArity(nrhs, 9, nlhs, 1,
        "value = radia_mex('optuna.parzen.log_pdf_numerical', samples, weights, mu, sigma, internal_low, internal_high, log_scale, step)");
    const auto samples = RealVector(prhs[1], "samples");
    const auto weights = RealVector(prhs[2], "weights");
    const auto mu = RealVector(prhs[3], "mu");
    const auto sigma = RealVector(prhs[4], "sigma");
    const double internal_low = Scalar(prhs[5], "internal_low");
    const double internal_high = Scalar(prhs[6], "internal_high");
    const bool log_scale = Boolean(prhs[7], "log_scale");
    const auto step_values = RealVector(prhs[8], "step");
    if (step_values.size() != 1)
        BadArgument("step must be a real scalar or NaN");
    const double step = step_values.front();
    const bool discrete = std::isfinite(step);
    if (!(internal_low < internal_high) || weights.empty() ||
        weights.size() != mu.size() || weights.size() != sigma.size())
        BadArgument("weights, mu, and sigma must be nonempty equal-length vectors");
    if (discrete && step <= 0.0)
        BadArgument("finite step must be positive");
    if (!discrete && !std::isnan(step))
        BadArgument("step must be positive and finite, or NaN for continuous density");
    double weight_sum = 0.0;
    for (std::size_t kernel = 0; kernel < weights.size(); ++kernel) {
        if (!std::isfinite(weights[kernel]) || weights[kernel] < 0.0 ||
            !std::isfinite(mu[kernel]) || !std::isfinite(sigma[kernel]) ||
            sigma[kernel] <= 0.0)
            BadArgument("Parzen weights/mu/sigma contain invalid entries");
        weight_sum += weights[kernel];
    }
    if (!(weight_sum > 0.0))
        BadArgument("Parzen weights must contain positive mass");

    for (double sample : samples)
        if (!std::isfinite(sample) || (log_scale && sample <= 0.0) ||
            (discrete && log_scale && sample - 0.5 * step <= 0.0))
            BadArgument("samples contain values outside the numerical density domain");
    plhs[0] = RealColumn(optuna_cuda::log_pdf(samples, {{
        &weights, &mu, &sigma, nullptr, 0, internal_low, internal_high,
        step, false, log_scale}}));
    return;

    std::vector<double> denominator(weights.size());
    for (std::size_t kernel = 0; kernel < weights.size(); ++kernel)
        denominator[kernel] = OptunaLogNormalMass(
            (internal_low - mu[kernel]) / sigma[kernel],
            (internal_high - mu[kernel]) / sigma[kernel]);

    std::vector<double> result(samples.size());
    std::vector<double> components(weights.size());
    constexpr double half_log_two_pi = 0.91893853320467274178;
    for (std::size_t sample_index = 0; sample_index < samples.size(); ++sample_index) {
        const double sample = samples[sample_index];
        if (!std::isfinite(sample) || (log_scale && sample <= 0.0))
            BadArgument("samples contain values outside the numerical density domain");
        for (std::size_t kernel = 0; kernel < weights.size(); ++kernel) {
            double component = 0.0;
            if (discrete) {
                double sample_low = sample - 0.5 * step;
                double sample_high = sample + 0.5 * step;
                if (log_scale) {
                    if (!(sample_low > 0.0 && sample_high > 0.0))
                        BadArgument("log-discrete sample interval must be positive");
                    sample_low = std::log(sample_low);
                    sample_high = std::log(sample_high);
                }
                component = OptunaLogNormalMass(
                    (sample_low - mu[kernel]) / sigma[kernel],
                    (sample_high - mu[kernel]) / sigma[kernel]) -
                    denominator[kernel];
            } else {
                const double internal_sample = log_scale ? std::log(sample) : sample;
                const double z = (internal_sample - mu[kernel]) / sigma[kernel];
                component = -0.5 * z * z - half_log_two_pi -
                    std::log(sigma[kernel]) - denominator[kernel];
            }
            components[kernel] = component + std::log(weights[kernel]);
        }
        result[sample_index] = OptunaLogSumExp(components);
    }
    plhs[0] = RealColumn(result);
}

void OptunaParzenLogPdfCategorical(int nlhs, mxArray* plhs[], int nrhs,
                                   const mxArray* prhs[]) {
    CheckArity(nrhs, 4, nlhs, 1,
        "value = radia_mex('optuna.parzen.log_pdf_categorical', samples, weights, probabilities)");
    const auto samples = RealVector(prhs[1], "samples");
    const auto weights = RealVector(prhs[2], "weights");
    std::size_t kernel_count = 0, choice_count = 0;
    const auto probabilities = RealMatrix(
        prhs[3], kernel_count, choice_count, "probabilities");
    if (kernel_count == 0 || choice_count == 0 || weights.size() != kernel_count)
        BadArgument("probabilities must have one row per mixture weight");
    double weight_sum = 0.0;
    for (double weight : weights) {
        if (!std::isfinite(weight) || weight < 0.0)
            BadArgument("categorical weights must be finite and nonnegative");
        weight_sum += weight;
    }
    if (!(weight_sum > 0.0))
        BadArgument("categorical weights must contain positive mass");
    for (double probability : probabilities)
        if (!std::isfinite(probability) || probability < 0.0)
            BadArgument("categorical probabilities must be finite and nonnegative");

    for (double sample : samples)
        if (!std::isfinite(sample) || sample != std::floor(sample) ||
            sample < 1.0 || sample > static_cast<double>(choice_count))
            BadArgument("categorical samples must be one-based integer choice indices");
    plhs[0] = RealColumn(optuna_cuda::log_pdf(samples, {{
        &weights, nullptr, nullptr, &probabilities, choice_count, 0.0, 0.0,
        std::numeric_limits<double>::quiet_NaN(), true, false}}));
    return;

    std::vector<double> result(samples.size());
    std::vector<double> components(kernel_count);
    for (std::size_t sample_index = 0; sample_index < samples.size(); ++sample_index) {
        const double raw_choice = samples[sample_index];
        if (!std::isfinite(raw_choice) || raw_choice != std::floor(raw_choice) ||
            raw_choice < 1.0 || raw_choice > static_cast<double>(choice_count))
            BadArgument("categorical samples must be one-based integer choice indices");
        const std::size_t choice = static_cast<std::size_t>(raw_choice) - 1;
        for (std::size_t kernel = 0; kernel < kernel_count; ++kernel) {
            const double probability = std::max(
                probabilities[kernel * choice_count + choice],
                std::numeric_limits<double>::min());
            components[kernel] = std::log(probability) + std::log(weights[kernel]);
        }
        result[sample_index] = OptunaLogSumExp(components);
    }
    plhs[0] = RealColumn(result);
}

void OptunaRandomStateSeed(OptunaRandomStateHandle& state,
                           std::uint32_t seed) {
    state.mt[0] = seed;
    for (std::size_t index = 1; index < state.mt.size(); ++index) {
        const std::uint32_t previous = state.mt[index - 1];
        state.mt[index] = 1812433253U *
            (previous ^ (previous >> 30U)) + static_cast<std::uint32_t>(index);
    }
    state.index = state.mt.size();
    state.has_gauss = false;
    state.gauss = 0.0;
}

void OptunaRandomStateTwist(OptunaRandomStateHandle& state) {
    constexpr std::uint32_t upper = 0x80000000U;
    constexpr std::uint32_t lower = 0x7fffffffU;
    constexpr std::uint32_t matrix_a = 0x9908b0dfU;
    constexpr std::size_t offset = 397;
    for (std::size_t index = 0; index < state.mt.size() - offset; ++index) {
        const std::uint32_t joined = (state.mt[index] & upper) |
                                     (state.mt[index + 1] & lower);
        state.mt[index] = state.mt[index + offset] ^ (joined >> 1U) ^
                          ((joined & 1U) ? matrix_a : 0U);
    }
    for (std::size_t index = state.mt.size() - offset;
         index < state.mt.size() - 1; ++index) {
        const std::uint32_t joined = (state.mt[index] & upper) |
                                     (state.mt[index + 1] & lower);
        state.mt[index] = state.mt[index + offset - state.mt.size()] ^
                          (joined >> 1U) ^ ((joined & 1U) ? matrix_a : 0U);
    }
    const std::uint32_t joined = (state.mt.back() & upper) |
                                 (state.mt.front() & lower);
    state.mt.back() = state.mt[offset - 1] ^ (joined >> 1U) ^
                      ((joined & 1U) ? matrix_a : 0U);
    state.index = 0;
}

std::uint32_t OptunaRandomStateUInt32(OptunaRandomStateHandle& state) {
    if (state.index >= state.mt.size())
        OptunaRandomStateTwist(state);
    std::uint32_t value = state.mt[state.index++];
    value ^= value >> 11U;
    value ^= (value << 7U) & 0x9d2c5680U;
    value ^= (value << 15U) & 0xefc60000U;
    value ^= value >> 18U;
    return value;
}

double OptunaRandomStateUniform(OptunaRandomStateHandle& state) {
    const std::uint32_t left = OptunaRandomStateUInt32(state) >> 5U;
    const std::uint32_t right = OptunaRandomStateUInt32(state) >> 6U;
    return (static_cast<double>(left) * 67108864.0 +
            static_cast<double>(right)) / 9007199254740992.0;
}

// Wichura's AS241 inverse standard-normal CDF.  The coefficient ordering
// below follows the original Horner form and retains full double precision.
double OptunaNormalInverse(double probability) {
    const double q = probability - 0.5;
    if (std::abs(q) <= 0.425) {
        const double r = 0.180625 - q * q;
        const double numerator = (((((((
            2509.0809287301226727 * r + 33430.575583588128105) * r +
            67265.770927008700853) * r + 45921.953931549871457) * r +
            13731.693765509461125) * r + 1971.5909503065514427) * r +
            133.14166789178437745) * r + 3.387132872796366608);
        const double denominator = (((((((
            5226.495278852854561 * r + 28729.085735721942674) * r +
            39307.89580009271061) * r + 21213.794301586595867) * r +
            5394.1960214247511077) * r + 687.18700749205790830) * r +
            42.313330701600911252) * r + 1.0);
        return q * numerator / denominator;
    }
    const double tail = q < 0.0 ? probability : 1.0 - probability;
    double r = std::sqrt(-std::log(tail));
    double value = 0.0;
    if (r <= 5.0) {
        r -= 1.6;
        const double numerator = (((((((
            7.7454501427834140764e-4 * r + 0.0227238449892691845833) * r +
            0.24178072517745061177) * r + 1.27045825245236838258) * r +
            3.64784832476320460504) * r + 5.76949722146069140550) * r +
            4.63033784615654529590) * r + 1.42343711074968357734);
        const double denominator = (((((((
            1.05075007164441684324e-9 * r +
            5.47593808499534494600e-4) * r +
            0.0151986665636164571966) * r + 0.148103976427480074590) * r +
            0.689767334985100004550) * r + 1.67638483018380384940) * r +
            2.05319162663775882187) * r + 1.0);
        value = numerator / denominator;
    } else {
        r -= 5.0;
        const double numerator = (((((((
            2.01033439929228813265e-7 * r +
            2.71155556874348757815e-5) * r +
            0.00124266094738807843860) * r + 0.0265321895265761230930) * r +
            0.29656057182850489123) * r + 1.78482653991729133580) * r +
            5.46378491116411436990) * r + 6.65790464350110377720);
        const double denominator = (((((((
            2.04426310338993978564e-15 * r +
            1.42151175831644588870e-7) * r +
            1.84631831751005468180e-5) * r +
            7.86869131145613259100e-4) * r + 0.0148753612908506148525) * r +
            0.136929880922735805310) * r + 0.599832206555887937690) * r + 1.0);
        value = numerator / denominator;
    }
    return q < 0.0 ? -value : value;
}

double OptunaParzenLogPdfAt(
    double sample, const std::vector<double>& weights,
    const std::vector<double>& mu, const std::vector<double>& sigma,
    const std::vector<double>& denominator, bool log_scale, double step) {
    constexpr double half_log_two_pi = 0.91893853320467274178;
    const bool discrete = std::isfinite(step);
    std::vector<double> components(weights.size());
    for (std::size_t kernel = 0; kernel < weights.size(); ++kernel) {
        double component = 0.0;
        if (discrete) {
            double sample_low = sample - 0.5 * step;
            double sample_high = sample + 0.5 * step;
            if (log_scale) {
                sample_low = std::log(sample_low);
                sample_high = std::log(sample_high);
            }
            component = OptunaLogNormalMass(
                (sample_low - mu[kernel]) / sigma[kernel],
                (sample_high - mu[kernel]) / sigma[kernel]) -
                denominator[kernel];
        } else {
            const double internal_sample = log_scale ? std::log(sample) : sample;
            const double z = (internal_sample - mu[kernel]) / sigma[kernel];
            component = -0.5 * z * z - half_log_two_pi -
                std::log(sigma[kernel]) - denominator[kernel];
        }
        components[kernel] = component + std::log(weights[kernel]);
    }
    return OptunaLogSumExp(components);
}

void OptunaTPEBestNumerical(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    CheckArity(nrhs, 15, nlhs, 1,
        "value = radia_mex('optuna.tpe.best_numerical', handle, count, below_weights, below_mu, below_sigma, above_weights, above_mu, above_sigma, internal_low, internal_high, low, high, log_scale, step)");
    OptunaRandomStateHandle& state = OptunaRandomState(Handle(prhs[1]));
    const int count = PositiveInteger(prhs[2], "count");
    const auto below_weights = RealVector(prhs[3], "below_weights");
    const auto below_mu = RealVector(prhs[4], "below_mu");
    const auto below_sigma = RealVector(prhs[5], "below_sigma");
    const auto above_weights = RealVector(prhs[6], "above_weights");
    const auto above_mu = RealVector(prhs[7], "above_mu");
    const auto above_sigma = RealVector(prhs[8], "above_sigma");
    const double internal_low = Scalar(prhs[9], "internal_low");
    const double internal_high = Scalar(prhs[10], "internal_high");
    const double low = Scalar(prhs[11], "low");
    const double high = Scalar(prhs[12], "high");
    const bool log_scale = Boolean(prhs[13], "log_scale");
    const auto step_values = RealVector(prhs[14], "step");
    if (step_values.size() != 1)
        BadArgument("step must be a real scalar or NaN");
    const double step = step_values.front();
    const bool discrete = std::isfinite(step);
    if (!(internal_low < internal_high) || !(low < high) ||
        below_weights.empty() || above_weights.empty() ||
        below_weights.size() != below_mu.size() ||
        below_weights.size() != below_sigma.size() ||
        above_weights.size() != above_mu.size() ||
        above_weights.size() != above_sigma.size() ||
        (discrete && step <= 0.0) || (!discrete && !std::isnan(step)))
        BadArgument("invalid fused numerical TPE estimator");

    auto validate = [](const std::vector<double>& weights,
                       const std::vector<double>& mu,
                       const std::vector<double>& sigma) {
        double total = 0.0;
        for (std::size_t index = 0; index < weights.size(); ++index) {
            if (!std::isfinite(weights[index]) || weights[index] < 0.0 ||
                !std::isfinite(mu[index]) || !std::isfinite(sigma[index]) ||
                sigma[index] <= 0.0)
                BadArgument("fused numerical TPE estimator contains invalid entries");
            total += weights[index];
        }
        if (!(total > 0.0))
            BadArgument("fused numerical TPE weights must contain positive mass");
    };
    validate(below_weights, below_mu, below_sigma);
    validate(above_weights, above_mu, above_sigma);

    auto denominator = [&](const std::vector<double>& mu,
                           const std::vector<double>& sigma) {
        std::vector<double> values(mu.size());
        for (std::size_t index = 0; index < mu.size(); ++index)
            values[index] = OptunaLogNormalMass(
                (internal_low - mu[index]) / sigma[index],
                (internal_high - mu[index]) / sigma[index]);
        return values;
    };
    const auto below_denominator = denominator(below_mu, below_sigma);
    const auto above_denominator = denominator(above_mu, above_sigma);

    RandomStateCheckpoint checkpoint(state);
    std::vector<double> cumulative(below_weights.size());
    std::partial_sum(below_weights.begin(), below_weights.end(), cumulative.begin());
    cumulative.back() = 1.0;
    std::vector<std::size_t> components(static_cast<std::size_t>(count));
    for (std::size_t& component : components) {
        const double uniform = OptunaRandomStateUniform(state);
        component = static_cast<std::size_t>(std::lower_bound(
            cumulative.begin(), cumulative.end(), uniform) - cumulative.begin());
    }

    constexpr double inverse_sqrt_two = 0.70710678118654752440;
    const double probability_epsilon = std::numeric_limits<double>::epsilon();
    std::vector<double> gpu_candidates(static_cast<std::size_t>(count));
    double best_value = std::numeric_limits<double>::quiet_NaN();
    double best_acquisition = -std::numeric_limits<double>::infinity();
    for (int index = 0; index < count; ++index) {
        const std::size_t component = components[static_cast<std::size_t>(index)];
        const double mean = below_mu[component];
        const double deviation = below_sigma[component];
        const double a = (internal_low - mean) / deviation;
        const double b = (internal_high - mean) / deviation;
        const double cdf_low = 0.5 * std::erfc(-a * inverse_sqrt_two);
        const double cdf_high = 0.5 * std::erfc(-b * inverse_sqrt_two);
        double probability = cdf_low + OptunaRandomStateUniform(state) *
            (cdf_high - cdf_low);
        probability = std::min(std::max(probability, probability_epsilon),
                               1.0 - probability_epsilon);
        const double z = OptunaNormalInverse(probability);
        double candidate = mean + deviation * z;
        if (log_scale)
            candidate = std::exp(candidate);
        if (discrete)
            candidate = low + std::round((candidate - low) / step) * step;
        candidate = std::min(std::max(candidate, low), high);
        gpu_candidates[static_cast<std::size_t>(index)] = candidate;
    }
    const auto best = optuna_cuda::best(gpu_candidates,
        {{&below_weights, &below_mu, &below_sigma, nullptr, 0,
          internal_low, internal_high, step, false, log_scale}},
        {{&above_weights, &above_mu, &above_sigma, nullptr, 0,
          internal_low, internal_high, step, false, log_scale}});
    best_value = gpu_candidates[best];
    checkpoint.accepted = true;
    plhs[0] = mxCreateDoubleScalar(best_value);
}

struct OptunaJointEstimator {
    std::shared_ptr<optuna_cuda::PreparedModel> device_model;
    bool categorical = false;
    std::vector<double> weights;
    std::vector<double> mu;
    std::vector<double> sigma;
    std::vector<double> denominator;
    std::vector<double> probabilities;
    std::size_t choice_count = 0;
    double internal_low = 0.0;
    double internal_high = 0.0;
    double low = 0.0;
    double high = 0.0;
    double step = std::numeric_limits<double>::quiet_NaN();
    bool log_scale = false;
};

const mxArray* OptunaRequiredField(const mxArray* value, const char* field) {
    if (!mxIsStruct(value) || mxGetNumberOfElements(value) != 1)
        BadArgument("joint TPE estimators must be scalar structs");
    const mxArray* result = mxGetField(value, 0, field);
    if (result == nullptr)
        BadArgument(std::string("joint TPE estimator is missing field '") +
                    field + "'");
    return result;
}

OptunaJointEstimator OptunaParseJointEstimator(const mxArray* value,
                                                bool categorical) {
    OptunaJointEstimator result;
    result.categorical = categorical;
    result.weights = RealVector(OptunaRequiredField(value, "weights"), "weights");
    if (result.weights.empty())
        BadArgument("joint TPE estimator weights must not be empty");
    if (categorical) {
        std::size_t rows = 0;
        result.probabilities = RealMatrix(
            OptunaRequiredField(value, "probabilities"), rows,
            result.choice_count, "probabilities");
        if (rows != result.weights.size() || result.choice_count == 0)
            BadArgument("joint categorical probabilities have invalid shape");
        for (std::size_t row = 0; row < rows; ++row) {
            double total = 0.0;
            for (std::size_t choice = 0; choice < result.choice_count; ++choice) {
                const double probability = result.probabilities[row * result.choice_count + choice];
                if (!std::isfinite(probability) || probability < 0.0)
                    BadArgument("joint categorical probabilities must be finite and nonnegative");
                total += probability;
            }
            if (std::abs(total - 1.0) > 1e-8)
                BadArgument("joint categorical probability rows must sum to one");
        }
    } else {
        result.mu = RealVector(OptunaRequiredField(value, "mu"), "mu");
        result.sigma = RealVector(OptunaRequiredField(value, "sigma"), "sigma");
        result.internal_low = Scalar(
            OptunaRequiredField(value, "internal_low"), "internal_low");
        result.internal_high = Scalar(
            OptunaRequiredField(value, "internal_high"), "internal_high");
        result.low = Scalar(OptunaRequiredField(value, "low"), "low");
        result.high = Scalar(OptunaRequiredField(value, "high"), "high");
        result.log_scale = Boolean(
            OptunaRequiredField(value, "log"), "log_scale");
        const auto step_values = RealVector(
            OptunaRequiredField(value, "step"), "step");
        if (step_values.size() != 1)
            BadArgument("joint numerical step must be a scalar or NaN");
        result.step = step_values.front();
        if ((!std::isnan(result.step) && (!std::isfinite(result.step) || result.step <= 0.0)) ||
            (result.log_scale && (result.low <= 0.0 ||
                (std::isfinite(result.step) && result.low - 0.5 * result.step <= 0.0))))
            BadArgument("joint numerical estimator has invalid step or logarithmic domain");
        if (result.mu.size() != result.weights.size() ||
            result.sigma.size() != result.weights.size() ||
            !(result.internal_low < result.internal_high) ||
            !(result.low < result.high))
            BadArgument("joint numerical estimator has invalid shape or bounds");
        result.denominator.resize(result.weights.size());
        for (std::size_t kernel = 0; kernel < result.weights.size(); ++kernel) {
            if (!std::isfinite(result.mu[kernel]) ||
                !std::isfinite(result.sigma[kernel]) || !(result.sigma[kernel] > 0.0))
                BadArgument("joint numerical mu/sigma must be finite with positive sigma");
            result.denominator[kernel] = OptunaLogNormalMass(
                (result.internal_low - result.mu[kernel]) / result.sigma[kernel],
                (result.internal_high - result.mu[kernel]) / result.sigma[kernel]);
        }
    }
    double total = 0.0;
    for (double weight : result.weights) {
        if (!std::isfinite(weight) || weight < 0.0)
            BadArgument("joint TPE weights must be finite and nonnegative");
        total += weight;
    }
    if (!(total > 0.0))
        BadArgument("joint TPE weights must contain positive mass");
    return result;
}

double OptunaJointComponentLogPdf(const OptunaJointEstimator& estimator,
                                  double sample, std::size_t kernel) {
    if (estimator.categorical) {
        const std::size_t choice = static_cast<std::size_t>(sample) - 1;
        return std::log(std::max(
            estimator.probabilities[kernel * estimator.choice_count + choice],
            std::numeric_limits<double>::min()));
    }
    constexpr double half_log_two_pi = 0.91893853320467274178;
    if (std::isfinite(estimator.step)) {
        double sample_low = sample - 0.5 * estimator.step;
        double sample_high = sample + 0.5 * estimator.step;
        if (estimator.log_scale) {
            sample_low = std::log(sample_low);
            sample_high = std::log(sample_high);
        }
        return OptunaLogNormalMass(
            (sample_low - estimator.mu[kernel]) / estimator.sigma[kernel],
            (sample_high - estimator.mu[kernel]) / estimator.sigma[kernel]) -
            estimator.denominator[kernel];
    }
    const double internal_sample = estimator.log_scale ? std::log(sample) : sample;
    const double z = (internal_sample - estimator.mu[kernel]) /
                     estimator.sigma[kernel];
    return -0.5 * z * z - half_log_two_pi -
           std::log(estimator.sigma[kernel]) - estimator.denominator[kernel];
}

std::vector<double> OptunaBestJoint(
    OptunaRandomStateHandle& state, int count,
    const std::vector<OptunaJointEstimator>& below,
    const std::vector<OptunaJointEstimator>& above) {
    const std::size_t dimension = below.size();
    if (dimension == 0 || above.size() != dimension)
        BadArgument("joint TPE estimator dimensions must agree");
    if (dimension > std::numeric_limits<std::size_t>::max() /
                        sizeof(double) / static_cast<std::size_t>(count))
        BadArgument("joint TPE candidate allocation size overflow");
    for (std::size_t index = 0; index < dimension; ++index) {
        if (below[index].categorical != above[index].categorical ||
            (below[index].categorical && below[index].choice_count != above[index].choice_count))
            BadArgument("joint TPE below/above distribution domains must agree");
    }
    const std::size_t below_kernels = below.front().weights.size();
    const std::size_t above_kernels = above.front().weights.size();
    for (std::size_t index = 1; index < dimension; ++index)
        if (below[index].weights.size() != below_kernels ||
            above[index].weights.size() != above_kernels)
            BadArgument("joint TPE dimensions must share mixture component counts");

    RandomStateCheckpoint checkpoint(state);
    std::vector<double> cumulative(below.front().weights.size());
    std::partial_sum(below.front().weights.begin(), below.front().weights.end(),
                     cumulative.begin());
    cumulative.back() = 1.0;
    std::vector<std::size_t> components(static_cast<std::size_t>(count));
    for (std::size_t& component : components) {
        const double uniform = OptunaRandomStateUniform(state);
        component = static_cast<std::size_t>(std::lower_bound(
            cumulative.begin(), cumulative.end(), uniform) - cumulative.begin());
    }
    std::vector<double> candidates(static_cast<std::size_t>(count) * dimension);
    for (std::size_t dim = 0; dim < dimension; ++dim) {
        if (!below[dim].categorical)
            continue;
        for (int row = 0; row < count; ++row) {
            const std::size_t component = components[static_cast<std::size_t>(row)];
            const double uniform = OptunaRandomStateUniform(state);
            double running = 0.0;
            std::size_t choice = 0;
            for (; choice + 1 < below[dim].choice_count; ++choice) {
                running += below[dim].probabilities[
                    component * below[dim].choice_count + choice];
                if (running >= uniform)
                    break;
            }
            candidates[static_cast<std::size_t>(row) * dimension + dim] =
                static_cast<double>(choice + 1);
        }
    }
    constexpr double inverse_sqrt_two = 0.70710678118654752440;
    const double probability_epsilon = std::numeric_limits<double>::epsilon();
    for (std::size_t dim = 0; dim < dimension; ++dim) {
        if (below[dim].categorical)
            continue;
        for (int row = 0; row < count; ++row) {
            const std::size_t component = components[static_cast<std::size_t>(row)];
            const double mean = below[dim].mu[component];
            const double deviation = below[dim].sigma[component];
            const double a = (below[dim].internal_low - mean) / deviation;
            const double b = (below[dim].internal_high - mean) / deviation;
            const double cdf_low = 0.5 * std::erfc(-a * inverse_sqrt_two);
            const double cdf_high = 0.5 * std::erfc(-b * inverse_sqrt_two);
            double probability = cdf_low + OptunaRandomStateUniform(state) *
                (cdf_high - cdf_low);
            probability = std::min(std::max(probability, probability_epsilon),
                                   1.0 - probability_epsilon);
            double candidate = mean + deviation * OptunaNormalInverse(probability);
            if (below[dim].log_scale)
                candidate = std::exp(candidate);
            if (std::isfinite(below[dim].step))
                candidate = below[dim].low + std::round(
                    (candidate - below[dim].low) / below[dim].step) *
                    below[dim].step;
            candidates[static_cast<std::size_t>(row) * dimension + dim] =
                std::min(std::max(candidate, below[dim].low), below[dim].high);
        }
    }

    auto views = [](const std::vector<OptunaJointEstimator>& source) {
        std::vector<optuna_cuda::Estimator> result;
        for (const auto& e : source)
            result.push_back({&e.weights, &e.mu, &e.sigma, &e.probabilities,
                e.choice_count, e.internal_low, e.internal_high, e.step,
                e.categorical, e.log_scale, e.device_model.get(), &e.denominator});
        return result;
    };
    const auto best = optuna_cuda::best(candidates, views(below), views(above));
    checkpoint.accepted = true;
    std::vector<double> result(dimension);
    std::copy(candidates.begin() + best * dimension,
              candidates.begin() + (best + 1) * dimension, result.begin());
    return result;
}

void OptunaTPEBestJoint(int nlhs, mxArray* plhs[], int nrhs,
                        const mxArray* prhs[]) {
    CheckArity(nrhs, 6, nlhs, 1,
        "values = radia_mex('optuna.tpe.best_joint', handle, count, categorical, below_estimators, above_estimators)");
    OptunaRandomStateHandle& state = OptunaRandomState(Handle(prhs[1]));
    const int count = PositiveInteger(prhs[2], "count");
    const mxArray* categorical = prhs[3];
    const mxArray* below_cells = prhs[4];
    const mxArray* above_cells = prhs[5];
    if (!mxIsLogical(categorical) || !mxIsCell(below_cells) ||
        !mxIsCell(above_cells))
        BadArgument("joint TPE categorical flags must be logical and estimators cells");
    const std::size_t dimension = mxGetNumberOfElements(categorical);
    if (dimension == 0 || mxGetNumberOfElements(below_cells) != dimension ||
        mxGetNumberOfElements(above_cells) != dimension)
        BadArgument("joint TPE estimator dimensions must agree");
    const auto* categorical_data = mxGetLogicals(categorical);
    std::vector<OptunaJointEstimator> below(dimension), above(dimension);
    for (std::size_t index = 0; index < dimension; ++index) {
        below[index] = OptunaParseJointEstimator(
            mxGetCell(below_cells, index), categorical_data[index]);
        above[index] = OptunaParseJointEstimator(
            mxGetCell(above_cells, index), categorical_data[index]);
    }
    plhs[0] = RealRow(OptunaBestJoint(state, count, below, above));
}

std::vector<double> OptunaMixtureWeights(
    std::size_t observation_count, double prior_weight,
    const std::vector<double>& supplied) {
    if (!std::isfinite(prior_weight) || prior_weight < 0.0)
        BadArgument("joint TPE prior weight must be finite and nonnegative");
    if (observation_count == 0)
        return {1.0};
    std::vector<double> result(observation_count + 1, 1.0);
    if (!supplied.empty()) {
        if (supplied.size() != observation_count)
            BadArgument("joint TPE observation weights have invalid length");
        std::copy(supplied.begin(), supplied.end(), result.begin());
    } else if (observation_count >= 25) {
        const std::size_t ramp_count = observation_count - 25;
        if (ramp_count == 1) {
            result[0] = 1.0 / static_cast<double>(observation_count);
        } else if (ramp_count > 1) {
            const double start = 1.0 / static_cast<double>(observation_count);
            for (std::size_t index = 0; index < ramp_count; ++index)
                result[index] = start + (1.0 - start) *
                    static_cast<double>(index) /
                    static_cast<double>(ramp_count - 1);
        }
    }
    result.back() = prior_weight;
    double total = 0.0;
    for (double weight : result) {
        if (!std::isfinite(weight) || weight < 0.0)
            BadArgument("joint TPE observation weights are invalid");
        total += weight;
    }
    if (!(total > 0.0))
        BadArgument("joint TPE weights must contain positive mass");
    for (double& weight : result)
        weight /= total;
    return result;
}

namespace {

// NumPy 2.5 uses the x86-simd-sort four-lane bitonic network for the
// small arrays passed to np.argsort by Optuna's Parzen estimator.  Equal
// discrete observations can receive different neighbour bandwidths, so a
// stable_sort is observably different from the upstream 5.0 algorithm.
// Reproduce the scalar equivalent of that network for the <=256-observation
// hot path.  Larger arrays have the same numerical ordering unless values are
// exactly tied; their fallback remains deterministic.
void OptunaNumpySortStage(std::array<double, 4>& keys,
                          std::array<std::size_t, 4>& indices,
                          const std::array<int, 4>& permutation,
                          unsigned maximum_mask) {
    const auto input_keys = keys;
    const auto input_indices = indices;
    for (std::size_t lane = 0; lane < 4; ++lane) {
        const std::size_t other = static_cast<std::size_t>(permutation[lane]);
        const bool choose_maximum =
            (maximum_mask & (1U << static_cast<unsigned>(lane))) != 0;
        const bool choose_first = choose_maximum
            ? input_keys[lane] >= input_keys[other]
            : input_keys[lane] <= input_keys[other];
        keys[lane] = choose_first ? input_keys[lane] : input_keys[other];
        indices[lane] = choose_first
            ? input_indices[lane] : input_indices[other];
    }
}

void OptunaNumpySortFour(std::array<double, 4>& keys,
                         std::array<std::size_t, 4>& indices) {
    OptunaNumpySortStage(keys, indices, {1, 0, 3, 2}, 0xAU);
    OptunaNumpySortStage(keys, indices, {3, 2, 1, 0}, 0xCU);
    OptunaNumpySortStage(keys, indices, {1, 0, 3, 2}, 0xAU);
}

void OptunaNumpyMergeFour(std::array<double, 4>& keys,
                          std::array<std::size_t, 4>& indices) {
    OptunaNumpySortStage(keys, indices, {2, 3, 0, 1}, 0xCU);
    OptunaNumpySortStage(keys, indices, {1, 0, 3, 2}, 0xAU);
}

void OptunaNumpyCompareExchange(std::array<double, 4>& lower_keys,
                                std::array<std::size_t, 4>& lower_indices,
                                std::array<double, 4>& upper_keys,
                                std::array<std::size_t, 4>& upper_indices) {
    for (std::size_t lane = 0; lane < 4; ++lane) {
        if (lower_keys[lane] <= upper_keys[lane])
            continue;
        std::swap(lower_keys[lane], upper_keys[lane]);
        std::swap(lower_indices[lane], upper_indices[lane]);
    }
}

void OptunaNumpyMergeBlocks(
    std::vector<std::array<double, 4>>& keys,
    std::vector<std::array<std::size_t, 4>>& indices,
    std::size_t start, std::size_t count) {
    if (count == 2) {
        std::reverse(keys[start + 1].begin(), keys[start + 1].end());
        std::reverse(indices[start + 1].begin(), indices[start + 1].end());
        OptunaNumpyCompareExchange(keys[start], indices[start],
                                   keys[start + 1], indices[start + 1]);
        std::reverse(keys[start + 1].begin(), keys[start + 1].end());
        std::reverse(indices[start + 1].begin(), indices[start + 1].end());
    } else {
        for (std::size_t offset = 0; offset < count / 2; ++offset) {
            const std::size_t lower = start + offset;
            const std::size_t upper = start + count - offset - 1;
            std::reverse(keys[upper].begin(), keys[upper].end());
            std::reverse(indices[upper].begin(), indices[upper].end());
            OptunaNumpyCompareExchange(keys[lower], indices[lower],
                                       keys[upper], indices[upper]);
            std::reverse(keys[upper].begin(), keys[upper].end());
            std::reverse(indices[upper].begin(), indices[upper].end());
        }
    }
    for (std::size_t width = count / 2; width >= 2; width /= 2) {
        for (std::size_t base = start; base < start + count; base += width)
            for (std::size_t offset = 0; offset < width / 2; ++offset)
                OptunaNumpyCompareExchange(
                    keys[base + offset], indices[base + offset],
                    keys[base + offset + width / 2],
                    indices[base + offset + width / 2]);
    }
    for (std::size_t block = start; block < start + count; ++block)
        OptunaNumpyMergeFour(keys[block], indices[block]);
}

std::vector<std::size_t> OptunaNumpyArgsort(
    const std::vector<double>& values) {
    const std::size_t count = values.size();
    std::vector<std::size_t> order(count);
    std::iota(order.begin(), order.end(), 0);
    if (count <= 1)
        return order;
    if (count > 256) {
        std::stable_sort(order.begin(), order.end(),
            [&](std::size_t left, std::size_t right) {
                return values[left] < values[right];
            });
        return order;
    }

    std::size_t block_count = 1;
    while (block_count * 4 < count)
        block_count *= 2;
    std::vector<std::array<double, 4>> keys(block_count);
    std::vector<std::array<std::size_t, 4>> indices(block_count);
    const std::size_t padding_index =
        std::numeric_limits<std::size_t>::max();
    for (std::size_t position = 0; position < block_count * 4; ++position) {
        const std::size_t block = position / 4;
        const std::size_t lane = position % 4;
        keys[block][lane] = position < count
            ? values[position] : std::numeric_limits<double>::infinity();
        indices[block][lane] = position < count ? position : padding_index;
    }
    for (std::size_t block = 0; block < block_count; ++block)
        OptunaNumpySortFour(keys[block], indices[block]);
    for (std::size_t width = 2; width <= block_count; width *= 2)
        for (std::size_t start = 0; start < block_count; start += width)
            OptunaNumpyMergeBlocks(keys, indices, start, width);
    for (std::size_t position = 0; position < count; ++position)
        order[position] = indices[position / 4][position % 4];
    return order;
}

}  // namespace

std::vector<OptunaJointEstimator> OptunaBuildJointEstimators(
    const std::vector<double>& observations, std::size_t observation_count,
    std::size_t dimension, const mxLogical* categorical,
    const std::vector<double>& lows, const std::vector<double>& highs,
    const mxLogical* log_scale, const std::vector<double>& steps,
    const std::vector<int>& choice_counts, double prior_weight,
    bool magic_clip, bool endpoints,
    const std::vector<double>& observation_weights, bool host_sampling = true) {
    std::vector<OptunaJointEstimator> result(dimension);
    const auto weights = OptunaMixtureWeights(
        observation_count, prior_weight, observation_weights);
    if (optuna_cuda::use_gpu()) {
        std::vector<optuna_cuda::ModelDimension> settings(dimension);
        std::memset(settings.data(), 0, dimension * sizeof(optuna_cuda::ModelDimension));
        for (std::size_t dim = 0; dim < dimension; ++dim) {
            auto& estimator = result[dim];
            estimator.categorical = categorical[dim]; estimator.weights = weights;
            estimator.low = lows[dim]; estimator.high = highs[dim];
            estimator.step = steps[dim]; estimator.log_scale = log_scale[dim];
            estimator.choice_count = choice_counts[dim] > 0 ? static_cast<std::size_t>(choice_counts[dim]) : 0;
            if (estimator.categorical && !estimator.choice_count)
                BadArgument("joint TPE categorical choice count must be positive");
            const bool discrete = std::isfinite(estimator.step);
            const double low = discrete ? estimator.low - 0.5 * estimator.step : estimator.low;
            const double high = discrete ? estimator.high + 0.5 * estimator.step : estimator.high;
            if (!estimator.categorical && (!(low < high) || (estimator.log_scale && low <= 0.0)))
                BadArgument("joint TPE numerical support is invalid");
            estimator.internal_low = estimator.log_scale ? std::log(low) : low;
            estimator.internal_high = estimator.log_scale ? std::log(high) : high;
            const double span = estimator.internal_high - estimator.internal_low;
            const double minimum = magic_clip ? span / std::min<double>(100.0,
                static_cast<double>(observation_count) + 2.0) :
                std::numeric_limits<double>::epsilon();
            auto& setting = settings[dim];
            setting.low = estimator.internal_low; setting.high = estimator.internal_high;
            setting.step = estimator.step; setting.minimum = minimum;
            setting.choices = estimator.choice_count; setting.categorical = estimator.categorical;
            setting.log_scale = estimator.log_scale;
            for (std::size_t row = 0; row < observation_count; ++row) {
                const double raw = observations[row * dimension + dim];
                if (!std::isfinite(raw) || (estimator.categorical ?
                    (raw != std::floor(raw) || raw < 1.0 || raw > static_cast<double>(estimator.choice_count)) :
                    (estimator.log_scale && raw <= 0.0)))
                    BadArgument("TPE observation is invalid");
            }
        }
        auto model = optuna_cuda::build_model(observations, weights, settings,
            prior_weight, endpoints, host_sampling);
        for (std::size_t dim = 0; dim < dimension; ++dim) {
            result[dim].device_model = model.device;
            result[dim].mu = std::move(model.mu[dim]);
            result[dim].sigma = std::move(model.sigma[dim]);
            result[dim].probabilities = std::move(model.probabilities[dim]);


        }
        return result;
    }
    for (std::size_t dim = 0; dim < dimension; ++dim) {
        auto& estimator = result[dim];
        estimator.categorical = categorical[dim];
        estimator.weights = weights;
        if (estimator.categorical) {
            if (choice_counts[dim] < 1)
                BadArgument("joint TPE categorical choice count must be positive");
            estimator.choice_count = static_cast<std::size_t>(choice_counts[dim]);
            estimator.probabilities.assign(
                (observation_count + 1) * estimator.choice_count,
                observation_count == 0
                    ? 1.0 / static_cast<double>(estimator.choice_count)
                    : prior_weight / static_cast<double>(observation_count + 1));
            if (observation_count > 0) {
                for (std::size_t row = 0; row < observation_count; ++row) {
                    const double raw = observations[row * dimension + dim];
                    if (!std::isfinite(raw) || raw != std::floor(raw) ||
                        raw < 1.0 || raw > static_cast<double>(estimator.choice_count))
                        BadArgument("joint TPE categorical observation is invalid");
                    estimator.probabilities[
                        row * estimator.choice_count +
                        static_cast<std::size_t>(raw - 1.0)] += 1.0;
                }
                for (std::size_t row = 0; row <= observation_count; ++row) {
                    double total = 0.0;
                    for (std::size_t choice = 0; choice < estimator.choice_count;
                         ++choice)
                        total += estimator.probabilities[
                            row * estimator.choice_count + choice];
                    if (total > 0.0) {
                        for (std::size_t choice = 0;
                             choice < estimator.choice_count; ++choice)
                            estimator.probabilities[
                                row * estimator.choice_count + choice] /= total;
                    } else {
                        for (std::size_t choice = 0;
                             choice < estimator.choice_count; ++choice)
                            estimator.probabilities[
                                row * estimator.choice_count + choice] =
                                1.0 / static_cast<double>(estimator.choice_count);
                    }
                }
            }
            continue;
        }
        estimator.low = lows[dim];
        estimator.high = highs[dim];
        estimator.log_scale = log_scale[dim];
        estimator.step = steps[dim];
        const bool discrete = std::isfinite(estimator.step);
        const double support_low = discrete
            ? estimator.low - 0.5 * estimator.step : estimator.low;
        const double support_high = discrete
            ? estimator.high + 0.5 * estimator.step : estimator.high;
        if (!(support_low < support_high) ||
            (estimator.log_scale && support_low <= 0.0))
            BadArgument("joint TPE numerical support is invalid");
        estimator.internal_low = estimator.log_scale
            ? std::log(support_low) : support_low;
        estimator.internal_high = estimator.log_scale
            ? std::log(support_high) : support_high;
        const double span = estimator.internal_high - estimator.internal_low;
        double minimum = std::numeric_limits<double>::epsilon();
        if (magic_clip)
            minimum = span / std::min<double>(
                100.0, 1.0 + static_cast<double>(observation_count + 1));
        estimator.mu.resize(observation_count);
        for (std::size_t row = 0; row < observation_count; ++row) {
            const double raw = observations[row * dimension + dim];
            estimator.mu[row] = estimator.log_scale ? std::log(raw) : raw;
        }
        estimator.sigma.resize(observation_count);
        if (observation_count > 0) {
            const auto order = OptunaNumpyArgsort(estimator.mu);
            std::vector<double> sorted(observation_count);
            std::vector<double> sorted_sigma(observation_count);
            for (std::size_t index = 0; index < observation_count; ++index)
                sorted[index] = estimator.mu[order[index]];
            for (std::size_t index = 0; index < observation_count; ++index) {
                const double left = index == 0 ? estimator.internal_low
                                               : sorted[index - 1];
                const double right = index + 1 == observation_count
                    ? estimator.internal_high : sorted[index + 1];
                sorted_sigma[index] = std::max(
                    sorted[index] - left, right - sorted[index]);
            }
            if (!endpoints && observation_count >= 2) {
                sorted_sigma.front() = sorted[1] - sorted[0];
                sorted_sigma.back() =
                    sorted[observation_count - 1] -
                    sorted[observation_count - 2];
            }
            for (std::size_t index = 0; index < observation_count; ++index)
                estimator.sigma[order[index]] = std::min(
                    std::max(sorted_sigma[index], minimum), span);
        }
        estimator.mu.push_back(0.5 *
            (estimator.internal_low + estimator.internal_high));
        estimator.sigma.push_back(span);
        estimator.denominator.resize(observation_count + 1);
        for (std::size_t kernel = 0; kernel <= observation_count; ++kernel)
            estimator.denominator[kernel] = OptunaLogNormalMass(
                (estimator.internal_low - estimator.mu[kernel]) /
                    estimator.sigma[kernel],
                (estimator.internal_high - estimator.mu[kernel]) /
                    estimator.sigma[kernel]);
    }
    return result;
}

void OptunaTPEBestJointObservations(int nlhs, mxArray* plhs[], int nrhs,
                                    const mxArray* prhs[]) {
    CheckArity(nrhs, 16, nlhs, 1,
        "values = radia_mex('optuna.tpe.best_joint_observations', handle, count, categorical, lows, highs, log_scale, steps, choice_counts, good, bad, prior_weight, magic_clip, endpoints, below_weights, above_weights)");
    OptunaRandomStateHandle& state = OptunaRandomState(Handle(prhs[1]));
    const int count = PositiveInteger(prhs[2], "count");
    const mxArray* categorical_array = prhs[3];
    const mxArray* log_array = prhs[6];
    if (!mxIsLogical(categorical_array) || !mxIsLogical(log_array))
        BadArgument("joint TPE categorical and log flags must be logical");
    const std::size_t dimension = mxGetNumberOfElements(categorical_array);
    if (dimension == 0 || mxGetNumberOfElements(log_array) != dimension)
        BadArgument("joint TPE distribution dimensions must agree");
    const auto lows = FixedRealVector(prhs[4], dimension, "lows");
    const auto highs = FixedRealVector(prhs[5], dimension, "highs");
    const auto steps = FixedRealVector(prhs[7], dimension, "steps");
    const auto choice_counts = IntegerVector(prhs[8], "choice_counts");
    if (choice_counts.size() != dimension)
        BadArgument("joint TPE choice counts have invalid length");
    std::size_t good_count = 0, good_dimension = 0;
    std::size_t bad_count = 0, bad_dimension = 0;
    const auto good = RealMatrix(prhs[9], good_count, good_dimension, "good");
    const auto bad = RealMatrix(prhs[10], bad_count, bad_dimension, "bad");
    if (good_dimension != dimension || bad_dimension != dimension)
        BadArgument("joint TPE observation matrices have invalid width");
    const double prior_weight = Scalar(prhs[11], "prior_weight");
    const bool magic_clip = Boolean(prhs[12], "magic_clip");
    const bool endpoints = Boolean(prhs[13], "endpoints");
    const auto below_weights = RealVector(prhs[14], "below_weights");
    const auto above_weights = RealVector(prhs[15], "above_weights");
    const auto* categorical = mxGetLogicals(categorical_array);
    const auto* log_scale = mxGetLogicals(log_array);
    const auto below = OptunaBuildJointEstimators(
        good, good_count, dimension, categorical, lows, highs, log_scale,
        steps, choice_counts, prior_weight, magic_clip, endpoints,
        below_weights);
    const auto above = OptunaBuildJointEstimators(
        bad, bad_count, dimension, categorical, lows, highs, log_scale,
        steps, choice_counts, prior_weight, magic_clip, endpoints,
        above_weights, false);
    plhs[0] = RealRow(OptunaBestJoint(state, count, below, above));
}

OptunaJointEstimator OptunaBuildUnivariateNumericalEstimator(
    const std::vector<double>& observations, double low, double high,
    bool log_scale, double step, double prior_weight, bool magic_clip,
    bool endpoints, const std::vector<double>& observation_weights, bool host_sampling = true) {
    if (optuna_cuda::use_gpu()) {
        const mxLogical categorical = false, logarithmic = log_scale;
        return OptunaBuildJointEstimators(observations, observations.size(), 1,
            &categorical, {low}, {high}, &logarithmic, {step}, {0}, prior_weight,
            magic_clip, endpoints, observation_weights, host_sampling).front();
    }
    OptunaJointEstimator estimator;
    estimator.weights = OptunaMixtureWeights(
        observations.size(), prior_weight, observation_weights);
    estimator.low = low;
    estimator.high = high;
    estimator.log_scale = log_scale;
    estimator.step = step;
    const bool discrete = std::isfinite(step);
    const double support_low = discrete ? low - 0.5 * step : low;
    const double support_high = discrete ? high + 0.5 * step : high;
    if (!(support_low < support_high) || (log_scale && support_low <= 0.0))
        BadArgument("univariate TPE numerical support is invalid");
    estimator.internal_low = log_scale ? std::log(support_low) : support_low;
    estimator.internal_high = log_scale ? std::log(support_high) : support_high;
    const std::size_t n = observations.size();
    std::vector<double> mus(n);
    for (std::size_t index = 0; index < n; ++index) {
        if (!std::isfinite(observations[index]) ||
            (log_scale && observations[index] <= 0.0))
            BadArgument("univariate TPE observation is invalid");
        mus[index] = log_scale ? std::log(observations[index]) : observations[index];
    }
    std::vector<double> sigmas(n);
    if (n > 0) {
        const auto order = OptunaNumpyArgsort(mus);
        std::vector<double> sorted(n), sorted_sigma(n);
        for (std::size_t index = 0; index < n; ++index)
            sorted[index] = mus[order[index]];
        for (std::size_t index = 0; index < n; ++index) {
            const double left = index == 0 ? estimator.internal_low
                                           : sorted[index - 1];
            const double right = index + 1 == n ? estimator.internal_high
                                                 : sorted[index + 1];
            sorted_sigma[index] = std::max(sorted[index] - left,
                                           right - sorted[index]);
        }
        if (!endpoints && n >= 2) {
            sorted_sigma.front() = sorted[1] - sorted[0];
            sorted_sigma.back() = sorted[n - 1] - sorted[n - 2];
        }
        const double span = estimator.internal_high - estimator.internal_low;
        double minimum = std::numeric_limits<double>::epsilon();
        if (magic_clip)
            minimum = span / std::min<double>(
                100.0, 1.0 + static_cast<double>(n + 1));
        for (std::size_t index = 0; index < n; ++index)
            sigmas[order[index]] = std::min(
                std::max(sorted_sigma[index], minimum), span);
    }
    estimator.mu = mus;
    estimator.sigma = sigmas;
    estimator.mu.push_back(0.5 *
        (estimator.internal_low + estimator.internal_high));
    estimator.sigma.push_back(estimator.internal_high - estimator.internal_low);
    estimator.denominator.resize(n + 1);
    for (std::size_t kernel = 0; kernel <= n; ++kernel)
        estimator.denominator[kernel] = OptunaLogNormalMass(
            (estimator.internal_low - estimator.mu[kernel]) /
                estimator.sigma[kernel],
            (estimator.internal_high - estimator.mu[kernel]) /
                estimator.sigma[kernel]);
    return estimator;
}

void OptunaTPEBestNumericalObservations(int nlhs, mxArray* plhs[], int nrhs,
                                        const mxArray* prhs[]) {
    CheckArity(nrhs, 14, nlhs, 1,
        "value = radia_mex('optuna.tpe.best_numerical_observations', handle, count, good, bad, low, high, log_scale, step, prior_weight, magic_clip, endpoints, below_weights, above_weights)");
    OptunaRandomStateHandle& state = OptunaRandomState(Handle(prhs[1]));
    const int count = PositiveInteger(prhs[2], "count");
    const auto good = RealVector(prhs[3], "good");
    const auto bad = RealVector(prhs[4], "bad");
    const double low = Scalar(prhs[5], "low");
    const double high = Scalar(prhs[6], "high");
    const bool log_scale = Boolean(prhs[7], "log_scale");
    const auto step_values = RealVector(prhs[8], "step");
    if (step_values.size() != 1)
        BadArgument("univariate TPE step must be a scalar or NaN");
    const double prior_weight = Scalar(prhs[9], "prior_weight");
    const bool magic_clip = Boolean(prhs[10], "magic_clip");
    const bool endpoints = Boolean(prhs[11], "endpoints");
    const auto below_weights = RealVector(prhs[12], "below_weights");
    const auto above_weights = RealVector(prhs[13], "above_weights");
    std::vector<OptunaJointEstimator> below = {
        OptunaBuildUnivariateNumericalEstimator(
            good, low, high, log_scale, step_values.front(), prior_weight,
            magic_clip, endpoints, below_weights)};
    std::vector<OptunaJointEstimator> above = {
        OptunaBuildUnivariateNumericalEstimator(
            bad, low, high, log_scale, step_values.front(), prior_weight,
            magic_clip, endpoints, above_weights, false)};
    plhs[0] = mxCreateDoubleScalar(
        OptunaBestJoint(state, count, below, above).front());
}

void OptunaTPEHistoryReset(int nlhs, mxArray*[], int nrhs,
                           const mxArray* prhs[]) {
    CheckArity(nrhs, 2, nlhs, 0,
        "radia_mex('optuna.tpe.history.reset', handle)");
    auto& state = OptunaRandomState(Handle(prhs[1]));
    state.tpe_history.clear();
    state.tpe_history_rows.clear();
    state.tpe_chronological.clear();
}

void OptunaTPEHistoryAppendComplete(int nlhs, mxArray*[], int nrhs,
                                    const mxArray* prhs[]) {
    CheckArity(nrhs, 6, nlhs, 0,
        "radia_mex('optuna.tpe.history.append_complete', handle, trial_number, objective, distribution_ids, values)");
    auto& state = OptunaRandomState(Handle(prhs[1]));
    const int trial_number = NonnegativeInteger(prhs[2], "trial_number");
    const double objective = Scalar(prhs[3], "objective");
    if (!std::isfinite(objective))
        BadArgument("completed TPE history objective must be finite");
    const auto distribution_ids = IntegerVector(
        prhs[4], "distribution_ids");
    const auto values = RealVector(prhs[5], "values");
    if (distribution_ids.size() != values.size())
        BadArgument("TPE history distribution IDs and values must agree");
    OptunaTPEHistoryTrial trial;
    trial.trial_number = trial_number;
    trial.objective = objective;
    for (std::size_t index = 0; index < distribution_ids.size(); ++index) {
        const int identifier = distribution_ids[index];
        if (identifier < 1 || !std::isfinite(values[index]) ||
            trial.parameters.count(identifier) != 0)
            BadArgument(
                "TPE history distribution IDs must be unique positive integers with finite values");
        trial.parameters.emplace(identifier, values[index]);
    }
    const auto found = state.tpe_history_rows.find(trial_number);
    if (found == state.tpe_history_rows.end()) {
        const std::size_t row = state.tpe_history.size();
        state.tpe_history.push_back(std::move(trial));
        state.tpe_history_rows.emplace(trial_number, row);
        if (state.tpe_chronological.empty() ||
            state.tpe_history[state.tpe_chronological.back()].trial_number < trial_number) {
            state.tpe_chronological.push_back(row);
        } else {
            const auto position = std::lower_bound(
                state.tpe_chronological.begin(), state.tpe_chronological.end(), trial_number,
                [&](std::size_t index, int number) {
                    return state.tpe_history[index].trial_number < number;
                });
            state.tpe_chronological.insert(position, row);
        }
    } else {
        state.tpe_history[found->second] = std::move(trial);
    }
}

void OptunaTPEBestGroupedHistory(int nlhs, mxArray* plhs[], int nrhs,
                                 const mxArray* prhs[]) {
    CheckArity(nrhs, 17, nlhs, 1,
        "values = radia_mex('optuna.tpe.best_grouped_history', handle, count, group_offsets, distribution_ids, categorical, lows, highs, log_scale, steps, choice_counts, minimize, gamma, max_good, prior_weight, magic_clip, endpoints)");
    auto& state = OptunaRandomState(Handle(prhs[1]));
    const int count = PositiveInteger(prhs[2], "count");
    const auto group_offsets = IntegerVector(prhs[3], "group_offsets");
    const auto distribution_ids = IntegerVector(
        prhs[4], "distribution_ids");
    const std::size_t dimension = distribution_ids.size();
    if (dimension == 0 || group_offsets.size() < 2 ||
        group_offsets.front() != 0 ||
        group_offsets.back() != static_cast<int>(dimension))
        BadArgument("grouped TPE offsets do not span the distributions");
    for (std::size_t index = 1; index < group_offsets.size(); ++index)
        if (group_offsets[index] <= group_offsets[index - 1])
            BadArgument("grouped TPE groups must be nonempty and ordered");
    for (int identifier : distribution_ids)
        if (identifier < 1)
            BadArgument("grouped TPE distribution IDs must be positive");
    const mxArray* categorical_array = prhs[5];
    const mxArray* log_array = prhs[8];
    if (!mxIsLogical(categorical_array) || !mxIsLogical(log_array) ||
        mxGetNumberOfElements(categorical_array) != dimension ||
        mxGetNumberOfElements(log_array) != dimension)
        BadArgument("grouped TPE categorical and log flags are invalid");
    const auto* categorical = mxGetLogicals(categorical_array);
    const auto* log_scale = mxGetLogicals(log_array);
    const auto lows = FixedRealVector(prhs[6], dimension, "lows");
    const auto highs = FixedRealVector(prhs[7], dimension, "highs");
    const auto steps = FixedRealVector(prhs[9], dimension, "steps");
    const auto choice_counts = IntegerVector(prhs[10], "choice_counts");
    if (choice_counts.size() != dimension)
        BadArgument("grouped TPE choice counts have invalid length");
    const bool minimize = Boolean(prhs[11], "minimize");
    const double gamma = Scalar(prhs[12], "gamma");
    if (!(gamma > 0.0 && gamma <= 1.0))
        BadArgument("grouped TPE gamma must be in (0, 1]");
    const int max_good = PositiveInteger(prhs[13], "max_good");
    const double prior_weight = Scalar(prhs[14], "prior_weight");
    const bool magic_clip = Boolean(prhs[15], "magic_clip");
    const bool endpoints = Boolean(prhs[16], "endpoints");
    if (state.tpe_history.empty())
        BadArgument("grouped TPE history has no completed trials");

    const auto& chronological = state.tpe_chronological;
    auto ranked = chronological;
    const int finished = static_cast<int>(ranked.size());
    const int good_count = std::max(
        1, std::min(finished, std::min(
            max_good, static_cast<int>(std::ceil(gamma * finished)))));
    // Rustuna's partial-selection design avoids sorting the whole history.
    // Unlike an arbitrary unstable tie, use trial number as a secondary key
    // to preserve Optuna's stable chronological split and RNG observations.
    // Reference: optuna/rustuna ebb5e6a, rustuna_sampler/src/tpe/sampler.rs.
    if (good_count < finished) {
        std::nth_element(ranked.begin(), ranked.begin() + good_count, ranked.end(),
            [&](std::size_t left, std::size_t right) {
                const auto& lhs = state.tpe_history[left];
                const auto& rhs = state.tpe_history[right];
                if (lhs.objective == rhs.objective)
                    return lhs.trial_number < rhs.trial_number;
                return minimize ? lhs.objective < rhs.objective : lhs.objective > rhs.objective;
            });
    }
    std::vector<unsigned char> globally_good(state.tpe_history.size(), 0);
    for (int index = 0; index < good_count; ++index)
        globally_good[ranked[static_cast<std::size_t>(index)]] = 1;

    RandomStateCheckpoint checkpoint(state);
    std::vector<double> result(dimension);
    for (std::size_t group = 0; group + 1 < group_offsets.size(); ++group) {
        const std::size_t begin = static_cast<std::size_t>(group_offsets[group]);
        const std::size_t end = static_cast<std::size_t>(group_offsets[group + 1]);
        const std::size_t group_dimension = end - begin;
        std::vector<double> good;
        std::vector<double> bad;
        good.reserve(static_cast<std::size_t>(good_count) * group_dimension);
        bad.reserve((state.tpe_history.size() - good_count) * group_dimension);
        std::vector<double> row_values(group_dimension);
        for (std::size_t history_index : chronological) {
            const auto& trial = state.tpe_history[history_index];
            bool present = true;
            for (std::size_t dim = begin; dim < end; ++dim) {
                const auto parameter = trial.parameters.find(distribution_ids[dim]);
                if (parameter == trial.parameters.end()) {
                    present = false;
                    break;
                }
                row_values[dim - begin] = parameter->second;
            }
            if (!present)
                continue;
            auto& target = globally_good[history_index] ? good : bad;
            target.insert(target.end(), row_values.begin(), row_values.end());
        }
        const std::size_t group_observations =
            (good.size() + bad.size()) / group_dimension;
        if (group_observations == 0)
            BadArgument("grouped TPE group has no compatible observations");
        std::vector<double> group_lows(lows.begin() + begin, lows.begin() + end);
        std::vector<double> group_highs(
            highs.begin() + begin, highs.begin() + end);
        std::vector<double> group_steps(
            steps.begin() + begin, steps.begin() + end);
        std::vector<int> group_choices(
            choice_counts.begin() + begin, choice_counts.begin() + end);
        const auto below = OptunaBuildJointEstimators(
            good, good.size() / group_dimension, group_dimension,
            categorical + begin, group_lows, group_highs, log_scale + begin,
            group_steps, group_choices, prior_weight, magic_clip, endpoints, {});
        const auto above = OptunaBuildJointEstimators(
            bad, bad.size() / group_dimension, group_dimension,
            categorical + begin, group_lows, group_highs, log_scale + begin,
            group_steps, group_choices, prior_weight, magic_clip, endpoints, {}, false);
        const auto candidate = OptunaBestJoint(state, count, below, above);
        std::copy(candidate.begin(), candidate.end(), result.begin() + begin);
    }
    plhs[0] = RealRow(result);
    checkpoint.accepted = true;
}

std::uint32_t OptunaRandomStateInterval(OptunaRandomStateHandle& state,
                                        std::uint32_t maximum) {
    std::uint32_t mask = maximum;
    mask |= mask >> 1U;
    mask |= mask >> 2U;
    mask |= mask >> 4U;
    mask |= mask >> 8U;
    mask |= mask >> 16U;
    while (true) {
        const std::uint32_t value = OptunaRandomStateUInt32(state) & mask;
        if (value <= maximum)
            return value;
    }
}

mxArray* OptunaRandomStateSnapshot(const OptunaRandomStateHandle& state) {
    const char* fields[] = {"schema", "mt", "index", "has_gauss", "gauss"};
    mxArray* result = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(result, 0, "schema",
               mxCreateString("numpy.randomstate.mt19937.v1"));
    mxArray* mt = mxCreateNumericMatrix(1, state.mt.size(), mxUINT32_CLASS,
                                        mxREAL);
    auto* mt_data = static_cast<std::uint32_t*>(mxGetData(mt));
    std::copy(state.mt.begin(), state.mt.end(), mt_data);
    mxSetField(result, 0, "mt", mt);
    mxSetField(result, 0, "index",
               mxCreateDoubleScalar(static_cast<double>(state.index)));
    mxSetField(result, 0, "has_gauss", mxCreateLogicalScalar(state.has_gauss));
    mxSetField(result, 0, "gauss", mxCreateDoubleScalar(state.gauss));
    return result;
}

void OptunaRandomStateRestoreValue(OptunaRandomStateHandle& state,
                                   const mxArray* value) {
    if (!mxIsStruct(value) || mxGetNumberOfElements(value) != 1)
        BadArgument("random-state snapshot must be a scalar struct");
    const mxArray* schema = mxGetField(value, 0, "schema");
    const mxArray* mt = mxGetField(value, 0, "mt");
    const mxArray* index = mxGetField(value, 0, "index");
    const mxArray* has_gauss = mxGetField(value, 0, "has_gauss");
    const mxArray* gauss = mxGetField(value, 0, "gauss");
    if (schema == nullptr || mt == nullptr || index == nullptr ||
        has_gauss == nullptr || gauss == nullptr ||
        Text(schema, "snapshot.schema") != "numpy.randomstate.mt19937.v1" ||
        !mxIsUint32(mt) || mxIsComplex(mt) ||
        mxGetNumberOfElements(mt) != state.mt.size())
        BadArgument("invalid NumPy RandomState snapshot");
    const double raw_index = Scalar(index, "snapshot.index");
    if (raw_index < 0.0 || raw_index > static_cast<double>(state.mt.size()) ||
        raw_index != std::floor(raw_index))
        BadArgument("snapshot.index must be an integer in [0, 624]");
    const auto* mt_data = static_cast<const std::uint32_t*>(mxGetData(mt));
    std::copy(mt_data, mt_data + state.mt.size(), state.mt.begin());
    state.index = static_cast<std::size_t>(raw_index);
    state.has_gauss = Boolean(has_gauss, "snapshot.has_gauss");
    state.gauss = Scalar(gauss, "snapshot.gauss");
}

void OptunaRandomStateCommand(const std::string& command, int nlhs,
                              mxArray* plhs[], int nrhs,
                              const mxArray* prhs[]) {
    if (command == "optuna.random_state.create") {
        CheckArity(nrhs, 2, nlhs, 1,
            "handle = radia_mex('optuna.random_state.create', seed)");
        const double raw_seed = Scalar(prhs[1], "seed");
        if (raw_seed < 0.0 || raw_seed != std::floor(raw_seed) ||
            raw_seed > static_cast<double>(std::numeric_limits<std::uint32_t>::max()))
            BadArgument("seed must be an integer in the uint32 range");
        auto state = std::make_unique<OptunaRandomStateHandle>();
        OptunaRandomStateSeed(*state, static_cast<std::uint32_t>(raw_seed));
        plhs[0] = Uint64Output(RegisterOptunaRandomState(std::move(state)));
        return;
    }

    const bool count_command = command == "optuna.random_state.rand" ||
        command == "optuna.random_state.randn" ||
        command == "optuna.random_state.randperm";
    CheckArity(nrhs, command == "optuna.random_state.randi" ? 4 :
                       ((command == "optuna.random_state.restore" || count_command) ? 3 : 2),
               nlhs,
               command == "optuna.random_state.restore" ||
               command == "optuna.random_state.destroy" ? 0 : 1,
               "invalid optuna.random_state command arguments");
    const std::uint64_t handle = Handle(prhs[1]);
    if (command == "optuna.random_state.destroy") {
        DestroyOptunaRandomState(handle);
        return;
    }
    OptunaRandomStateHandle& state = OptunaRandomState(handle);
    if (command == "optuna.random_state.snapshot") {
        plhs[0] = OptunaRandomStateSnapshot(state);
        return;
    }
    if (command == "optuna.random_state.restore") {
        OptunaRandomStateRestoreValue(state, prhs[2]);
        return;
    }
    if (command == "optuna.random_state.rand" ||
        command == "optuna.random_state.randn" ||
        command == "optuna.random_state.randperm") {
        const int count = NonnegativeInteger(prhs[2], "count");
        std::vector<double> result(static_cast<std::size_t>(count));
        if (command == "optuna.random_state.rand") {
            for (double& value : result)
                value = OptunaRandomStateUniform(state);
            plhs[0] = RealColumn(result);
            return;
        }
        if (command == "optuna.random_state.randn") {
            for (double& value : result) {
                if (state.has_gauss) {
                    value = state.gauss;
                    state.has_gauss = false;
                    state.gauss = 0.0;
                    continue;
                }
                double left = 0.0, right = 0.0, radius_squared = 2.0;
                while (radius_squared >= 1.0 || radius_squared == 0.0) {
                    left = 2.0 * OptunaRandomStateUniform(state) - 1.0;
                    right = 2.0 * OptunaRandomStateUniform(state) - 1.0;
                    radius_squared = left * left + right * right;
                }
                const double factor = std::sqrt(
                    -2.0 * std::log(radius_squared) / radius_squared);
                state.gauss = factor * left;
                state.has_gauss = true;
                value = factor * right;
            }
            plhs[0] = RealColumn(result);
            return;
        }
        std::iota(result.begin(), result.end(), 1.0);
        for (std::size_t right = result.size(); right > 1; --right) {
            const std::size_t left = static_cast<std::size_t>(
                OptunaRandomStateInterval(state,
                    static_cast<std::uint32_t>(right - 1)));
            std::swap(result[right - 1], result[left]);
        }
        plhs[0] = RealRow(result);
        return;
    }
    if (command == "optuna.random_state.randi") {
        const double raw_maximum = Scalar(prhs[2], "maximum");
        if (raw_maximum < 1.0 || raw_maximum != std::floor(raw_maximum) ||
            raw_maximum > static_cast<double>(std::numeric_limits<std::uint32_t>::max()))
            BadArgument("maximum must be a positive integer in the uint32 range");
        const int count = NonnegativeInteger(prhs[3], "count");
        std::vector<double> result(static_cast<std::size_t>(count));
        const std::uint32_t maximum = static_cast<std::uint32_t>(raw_maximum - 1.0);
        for (double& value : result)
            value = static_cast<double>(OptunaRandomStateInterval(state, maximum)) + 1.0;
        plhs[0] = RealColumn(result);
        return;
    }
    BadArgument("unknown optuna.random_state command");
}

void ApiInfo(int nlhs, mxArray* plhs[], int nrhs) {
    CheckArity(nrhs, 1, nlhs, 1, "info = optuna_mex('api.info')");
    const char* fields[] = {
        "api_version", "backend", "command_count", "handle_count", "build_route", "taskmanager", "gpu_operations", "host_operations"
        , "execution", "host_parallel_runtime", "allocator", "density_launches", "selection_launches", "host_to_device_bytes",
        "device_to_host_bytes", "peak_temporary_bytes", "cuda_runtime", "cuda_driver", "cuda_device", "execution_mode", "cpu_density_calls", "cpu_selection_calls", "cpu_parallel_density_calls", "cpu_max_density_threads", "device_allocation_calls", "workspace_reuses", "avoided_upload_bytes", "cached_device_bytes", "model_build_launches", "model_reuses", "preparation_launches", "model_download_bytes", "cpu_preparation_calls", "cpu_cached_normalizer_values", "model_radix_sorts", "peak_sort_temporary_bytes"
    };
    plhs[0] = mxCreateStructMatrix(1, 1, sizeof(fields) / sizeof(fields[0]), fields);
    mxSetField(plhs[0], 0, "build_route", mxCreateString("nvcc-shared-library"));
    mxSetField(plhs[0], 0, "taskmanager", mxCreateLogicalScalar(false));
    mxSetField(plhs[0], 0, "gpu_operations", mxCreateString("TPE bandwidth/category model construction, density preparation and candidate evaluation/selection"));
    mxSetField(plhs[0], 0, "host_operations", mxCreateString("Management, ordered RNG, candidate generation, Pareto and Sobol"));
    mxSetField(plhs[0], 0, "api_version", mxCreateDoubleScalar(1.0));
    mxSetField(plhs[0], 0, "backend", mxCreateString("optuna_mex"));
    mxSetField(plhs[0], 0, "command_count", mxCreateDoubleScalar(22.0));
    const auto stats = optuna_cuda::statistics();
    mxSetField(plhs[0], 0, "execution_mode", mxCreateString(optuna_cuda::use_gpu() ? "gpu" : "cpu"));
    mxSetField(plhs[0], 0, "cpu_density_calls", mxCreateDoubleScalar(static_cast<double>(stats.cpu_density_calls)));
    mxSetField(plhs[0], 0, "cpu_selection_calls", mxCreateDoubleScalar(static_cast<double>(stats.cpu_selection_calls)));
    mxSetField(plhs[0], 0, "cpu_parallel_density_calls", mxCreateDoubleScalar(static_cast<double>(stats.cpu_parallel_density_calls)));
    mxSetField(plhs[0], 0, "cpu_max_density_threads", mxCreateDoubleScalar(stats.cpu_max_density_threads));
    mxSetField(plhs[0], 0, "device_allocation_calls", mxCreateDoubleScalar(static_cast<double>(stats.device_allocation_calls)));
    mxSetField(plhs[0], 0, "workspace_reuses", mxCreateDoubleScalar(static_cast<double>(stats.workspace_reuses)));
    mxSetField(plhs[0], 0, "avoided_upload_bytes", mxCreateDoubleScalar(static_cast<double>(stats.avoided_upload_bytes)));
    mxSetField(plhs[0], 0, "cached_device_bytes", mxCreateDoubleScalar(static_cast<double>(stats.cached_device_bytes)));
    mxSetField(plhs[0], 0, "model_build_launches", mxCreateDoubleScalar(static_cast<double>(stats.model_build_launches)));
    mxSetField(plhs[0], 0, "model_reuses", mxCreateDoubleScalar(static_cast<double>(stats.model_reuses)));
    mxSetField(plhs[0], 0, "model_radix_sorts", mxCreateDoubleScalar(static_cast<double>(stats.model_radix_sorts)));
    mxSetField(plhs[0], 0, "peak_sort_temporary_bytes", mxCreateDoubleScalar(static_cast<double>(stats.peak_sort_temporary_bytes)));
    mxSetField(plhs[0], 0, "preparation_launches", mxCreateDoubleScalar(static_cast<double>(stats.preparation_launches)));
    mxSetField(plhs[0], 0, "model_download_bytes", mxCreateDoubleScalar(static_cast<double>(stats.model_download_bytes)));
    mxSetField(plhs[0], 0, "cpu_preparation_calls", mxCreateDoubleScalar(static_cast<double>(stats.cpu_preparation_calls)));
    mxSetField(plhs[0], 0, "cpu_cached_normalizer_values", mxCreateDoubleScalar(static_cast<double>(stats.cpu_cached_normalizer_values)));
    mxSetField(plhs[0], 0, "execution", mxCreateString("matlab-mex-shared-cuda"));
    mxSetField(plhs[0], 0, "host_parallel_runtime", mxCreateString("openmp"));
    mxSetField(plhs[0], 0, "allocator", mxCreateString(
        "runtime-adapter-owned"));
    mxSetField(plhs[0], 0, "density_launches", mxCreateDoubleScalar(static_cast<double>(stats.density_launches)));
    mxSetField(plhs[0], 0, "selection_launches", mxCreateDoubleScalar(static_cast<double>(stats.selection_launches)));
    mxSetField(plhs[0], 0, "host_to_device_bytes", mxCreateDoubleScalar(static_cast<double>(stats.host_to_device_bytes)));
    mxSetField(plhs[0], 0, "device_to_host_bytes", mxCreateDoubleScalar(static_cast<double>(stats.device_to_host_bytes)));
    mxSetField(plhs[0], 0, "peak_temporary_bytes", mxCreateDoubleScalar(static_cast<double>(stats.peak_temporary_bytes)));
    mxSetField(plhs[0], 0, "cuda_runtime", mxCreateDoubleScalar(stats.runtime_version));
    mxSetField(plhs[0], 0, "cuda_driver", mxCreateDoubleScalar(stats.driver_version));
    mxSetField(plhs[0], 0, "cuda_device", mxCreateDoubleScalar(stats.device));
    std::lock_guard<std::mutex> guard(registry_mutex);
    mxSetField(plhs[0], 0, "handle_count", mxCreateDoubleScalar(
        static_cast<double>(optuna_random_state_registry.size())));
}

void OptunaSobolPoints(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    CheckArity(nrhs, 3, nlhs, 1,
               "points = optuna_mex('optuna.sobol.points', directions, sample_ids)");
    const mxArray* direction_value = prhs[1];
    const mxArray* sample_value = prhs[2];
    if (!mxIsUint32(direction_value) || mxIsComplex(direction_value) ||
        mxGetNumberOfDimensions(direction_value) != 2 ||
        mxGetN(direction_value) != 32)
        BadArgument("directions must be a dimension-by-32 uint32 matrix");
    if (!mxIsUint32(sample_value) || mxIsComplex(sample_value))
        BadArgument("sample_ids must be a uint32 array");
    const std::size_t dimension = mxGetM(direction_value);
    const std::size_t count = mxGetNumberOfElements(sample_value);
    if (dimension == 0 || dimension > 21201)
        BadArgument("Sobol dimension must be between 1 and 21201");
    const auto* directions = static_cast<const std::uint32_t*>(
        mxGetData(direction_value));
    const auto* sample_ids = static_cast<const std::uint32_t*>(
        mxGetData(sample_value));
    plhs[0] = mxCreateDoubleMatrix(count, dimension, mxREAL);
    double* points = mxGetDoubles(plhs[0]);
    constexpr double scale = 1.0 / 4294967296.0;
    for (std::size_t dim = 0; dim < dimension; ++dim) {
        std::uint32_t previous_sample = 0;
        std::uint32_t previous_integer = 0;
        for (std::size_t row = 0; row < count; ++row) {
            const std::uint32_t sample = sample_ids[row];
            std::uint32_t integer;
            if (row != 0 &&
                previous_sample != std::numeric_limits<std::uint32_t>::max() &&
                sample == previous_sample + 1U) {
                // Consecutive Sobol points differ by one Gray-code bit.
                // QMC batches are contiguous, so update in O(1) instead of
                // recomputing every active bit for every point.
                unsigned bit = 0;
                std::uint32_t value = sample;
                while ((value & 1U) == 0U) {
                    value >>= 1U;
                    ++bit;
                }
                integer = previous_integer ^ directions[dim +
                    static_cast<std::size_t>(bit) * dimension];
            } else {
                std::uint32_t gray = sample ^ (sample >> 1U);
                integer = 0;
                unsigned bit = 0;
                while (gray != 0) {
                    if ((gray & 1U) != 0)
                        integer ^= directions[dim +
                            static_cast<std::size_t>(bit) * dimension];
                    gray >>= 1U;
                    ++bit;
                }
            }
            points[row + dim * count] = static_cast<double>(integer) * scale;
            previous_sample = sample;
            previous_integer = integer;
        }
    }
}

void Dispatch(const std::string& command, int nlhs, mxArray* plhs[], int nrhs,
              const mxArray* prhs[]) {
    if (command == "api.info") {
        ApiInfo(nlhs, plhs, nrhs);
        return;
    }
    if (command == "api.execution") {
        if ((nrhs != 1 && nrhs != 2) || nlhs != 1)
            BadArgument("usage: mode = optuna_mex('api.execution' [, 'cpu'|'gpu'])");
        if (nrhs == 2) {
            const auto mode = Text(prhs[1], "execution mode");
            if (mode != "cpu" && mode != "gpu") BadArgument("execution mode must be cpu or gpu");
            optuna_cuda::set_gpu(mode == "gpu");
        }
        plhs[0] = mxCreateString(optuna_cuda::use_gpu() ? "gpu" : "cpu");
        return;
    }
    if (command == "api.commands") {
        CheckArity(nrhs, 1, nlhs, 1,
                   "commands = optuna_mex('api.commands')");
        plhs[0] = Commands();
        return;
    }
    if (command == "optuna.pareto.rank_crowding") {
        OptunaParetoRankCrowding(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "optuna.parzen.log_pdf_numerical") {
        OptunaParzenLogPdfNumerical(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "optuna.parzen.log_pdf_categorical") {
        OptunaParzenLogPdfCategorical(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "optuna.tpe.best_numerical") {
        OptunaTPEBestNumerical(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "optuna.tpe.best_joint") {
        OptunaTPEBestJoint(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "optuna.tpe.best_joint_observations") {
        OptunaTPEBestJointObservations(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "optuna.tpe.best_numerical_observations") {
        OptunaTPEBestNumericalObservations(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "optuna.tpe.history.reset") {
        OptunaTPEHistoryReset(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "optuna.tpe.history.append_complete") {
        OptunaTPEHistoryAppendComplete(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "optuna.tpe.best_grouped_history") {
        OptunaTPEBestGroupedHistory(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command == "optuna.sobol.points") {
        OptunaSobolPoints(nlhs, plhs, nrhs, prhs);
        return;
    }
    if (command.rfind("optuna.random_state.", 0) == 0) {
        OptunaRandomStateCommand(command, nlhs, plhs, nrhs, prhs);
        return;
    }
    BadArgument("unknown optuna_mex command: " + command);
}

} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    try {
        if (nrhs < 1)
            BadArgument("first input must be a command character vector");
        EnsureExitHandler();
        Dispatch(Text(prhs[0], "command"), nlhs, plhs, nrhs, prhs);
    } catch (const std::exception& exception) {
        mexErrMsgIdAndTxt("radia:optuna:mex:Exception", "%s",
                         exception.what());
    } catch (...) {
        mexErrMsgIdAndTxt("radia:optuna:mex:UnknownException",
                         "unknown C++ exception");
    }
}
