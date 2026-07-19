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
#include "rad_hdiv_vim.h"
#include "rad_hybrid_vim_schur.h"

#include <core/taskmanager.hpp>
#include <flags.hpp>
#include <hcurlhofespace.hpp>
#include <hdivhofespace.hpp>
#include <meshaccess.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
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
    std::unique_ptr<HACApKChargeGram> manager;
};

std::mutex registry_mutex;
std::unordered_map<std::uint64_t, std::unique_ptr<EnergyStopMaterial>> energy_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<BEMHandle>> bem_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<PEECHandle>> peec_registry;
std::unordered_map<std::uint64_t, std::unique_ptr<ChargeGramHandle>> charge_gram_registry;
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
    while (lock_count > 0) {
        mexUnlock();
        --lock_count;
    }
}

std::uint64_t RegisterBEM(std::unique_ptr<BEMHandle> bem) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    while (next_handle == 0 || energy_registry.count(next_handle) != 0 ||
           bem_registry.count(next_handle) != 0 || peec_registry.count(next_handle) != 0 ||
           charge_gram_registry.count(next_handle) != 0)
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
    while (next_handle == 0 || energy_registry.count(next_handle) != 0 ||
           bem_registry.count(next_handle) != 0 || peec_registry.count(next_handle) != 0 ||
           charge_gram_registry.count(next_handle) != 0)
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
    while (next_handle == 0 || energy_registry.count(next_handle) != 0 ||
           bem_registry.count(next_handle) != 0 || peec_registry.count(next_handle) != 0 ||
           charge_gram_registry.count(next_handle) != 0)
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
    while (next_handle == 0 || energy_registry.count(next_handle) != 0 ||
           bem_registry.count(next_handle) != 0 || peec_registry.count(next_handle) != 0 ||
           charge_gram_registry.count(next_handle) != 0)
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

mxArray* Commands() {
    static const char* names[] = {
        "api.info", "api.commands", "taskmanager.probe", "ngsolve.space_info",
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
        "biot_savart.h_segments_complex", "biot_savart.a_segments_complex",
        "biot_savart.a_triangles_complex", "biot_savart.b_triangles_complex",
        "bem.assemble_sldl", "bem.assemble_sldl_p2",
        "hacapk.bem.create", "hacapk.bem.destroy", "hacapk.bem.build",
        "hacapk.bem.matvec", "hacapk.bem.info", "hacapk.peec.create",
        "hacapk.peec.destroy", "hacapk.peec.build", "hacapk.peec.matvec",
        "hacapk.peec.info", "hacapk.charge_gram.create_monopole",
        "hacapk.charge_gram.destroy", "hacapk.charge_gram.build",
        "hacapk.charge_gram.matvec", "hacapk.charge_gram.matvec_transpose",
        "hacapk.charge_gram.matvec_sym", "hacapk.charge_gram.entry",
        "hacapk.charge_gram.info",
        "radia.ObjHexahedron", "radia.ObjTetrahedron", "radia.ObjWedge",
        "radia.ObjPyramid", "radia.ObjCnt", "radia.MatLin",
        "radia.MatSatIsoTab", "radia.MatApl", "radia.Solve", "radia.Fld",
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
                                                           charge_gram_registry.size())));
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
    holder->manager = std::make_unique<HACApKChargeGram>(
        std::move(centroids), std::move(measures), std::move(self_energy));
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
    const int n = holder.manager->GetNDOF();
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
    const int n = holder.manager->GetNDOF();
    const int i = PositiveInteger(prhs[2], "i");
    const int j = PositiveInteger(prhs[3], "j");
    if (i > n || j > n)
        BadArgument("charge entry indices are out of range");
    plhs[0] = mxCreateDoubleScalar(
        holder.manager->GetInteractionMatrixElement(i - 1, j - 1));
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

void Dispatch(const std::string& command, int nlhs, mxArray* plhs[], int nrhs,
              const mxArray* prhs[]) {
    if (command == "api.info")
        ApiInfo(nlhs, plhs, nrhs);
    else if (command == "api.commands") {
        CheckArity(nrhs, 1, nlhs, 1, "commands = radia_mex('api.commands')");
        plhs[0] = Commands();
    } else if (command == "taskmanager.probe")
        TaskProbe(nlhs, plhs, nrhs, prhs);
    else if (command == "ngsolve.space_info")
        SpaceInfo(nlhs, plhs, nrhs, prhs);
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
    else if (command == "hacapk.charge_gram.create_monopole")
        ChargeGramCreateMonopole(nlhs, plhs, nrhs, prhs);
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
    else if (command == "hacapk.charge_gram.info")
        ChargeGramInfo(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.ObjHexahedron" ||
             command == "radia.ObjTetrahedron" ||
             command == "radia.ObjWedge" || command == "radia.ObjPyramid")
        RadiaPolyhedron(command, nlhs, plhs, nrhs, prhs);
    else if (command == "radia.ObjCnt")
        RadiaObjCnt(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.MatLin")
        RadiaMatLin(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.MatSatIsoTab")
        RadiaMatSatIsoTab(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.MatApl")
        RadiaMatApl(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.Solve")
        RadiaSolve(nlhs, plhs, nrhs, prhs);
    else if (command == "radia.Fld")
        RadiaFld(nlhs, plhs, nrhs, prhs);
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
