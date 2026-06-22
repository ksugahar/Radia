/**
 * @file radia_pybind.cpp
 * @brief Radia Python bindings using pybind11 (NGSolve-style)
 *
 * This file provides pybind11 bindings for all Radia C API functions.
 * It is the sole Python binding implementation (legacy radpy_pyapi.cpp removed 2026-01).
 *
 * The module is named _radia_pybind and is imported by src/radia/__init__.py.
 *
 * @author Claude Opus 4.5
 * @date 2026-01-25
 */

// NGSolve headers MUST be included FIRST, before Radia headers.
// Radia's radentry.h defines EXP as __declspec(dllexport) which conflicts
// with the EXP enum value in NGSolve's evalfunc.hpp.
#include <fem.hpp>
#include <comp.hpp>
#include <python_ngstd.hpp>

// Temporarily undefine EXP if NGSolve headers left it undefined
// (radentry.h will redefine it)
#ifdef EXP
#undef EXP
#endif

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>

#include <vector>
#include <array>
#include <complex>
#include <cmath>
#include <string>
#include <stdexcept>
#include <memory>
#include <unordered_map>
#include <optional>

// Radia core headers (after NGSolve to avoid EXP macro conflict)
#include "radentry.h"
#include "rad_constants.h"
#include "rad_highorder_nodes.h"
#include "rad_hacapk_peec.h"  // HACApK PEEC adapter (manager + sanity check)

// HACApK PCA cluster strategy setter/getter (C symbols).
// Forward-declared at file scope so the lambdas in PYBIND11_MODULE can call them.
extern "C" {
    void cHACApK_set_cluster_strategy(int strategy);
    int  cHACApK_get_cluster_strategy(void);
    double cHACApK_harith_self_test(int depth, int n_per_block);
    double cHACApK_harith_self_test_rk(int n_per_block, int rk_rank);
    double cHACApK_harith_self_test_addmul_rkrk(int m, int n, int inner,
                                                  int kA, int kB, int kC);
    double cHACApK_harith_self_test_rk_deep(int n_per_block, int rk_rank);
    double cHACApK_harith_self_test_mixed_sibling(int nb_small);
    void cHACApK_hlu_get_timings(double *out_t_decomp, double *out_t_solve,
                                  long *out_n_dense_lu, long *out_n_dense_gemm);
    void cHACApK_hlu_set_trunc_tol(double tol);
    double cHACApK_hlu_get_trunc_tol(void);
    void cHACApK_hlu_get_materialize_stats(long *out_n_calls, long *out_n_elems);
    void cHACApK_hlu_get_materialize_split(long *out_internal, long *out_leaf);
    double cHACApK_hlu_run_on_hacapk(void *leafmtxp_void, void *control_void,
                                     const double *x_orig, const double *y_orig, int nffc);
    void* cHACApK_hlu_factor_leafmtxp(void* leafmtxp_void, void* control_void, int nffc);
    int   cHACApK_hlu_apply(void* root_void, void* control_void, const double* r, double* z, int nd);
    void  cHACApK_hlu_free_factors(void* root_void);
    void cHACApK_hlu_get_mixed_breakdown(long *out_addmul9, long *out_lln9, long *out_run9);
    void cHACApK_hlu_set_parallel(int on);
    int  cHACApK_hlu_get_parallel(void);
    void cHACApK_hlu_set_par_cutoff(long c);
    int  chacapk_max_threads(void);
    void cHACApK_hlu_set_accum_cap(int c);
    int  cHACApK_hlu_get_accum_cap(void);
    double cHACApK_harith_self_test_mixed_sibling_nonuniform(int n1, int n2, int m1, int m3);
    double cHACApK_harith_self_test_mixed_sibling_via_conversion(int nb_small);
    double cHACApK_harith_self_test_depth3_asymmetric(int nb_tiny);
    double cHACApK_harith_self_test_radia_exact(void);
    double cHACApK_harith_self_test_radia_exact_diag(double diag_boost);
    double cHACApK_harith_self_test_radia_exact_with_matrix(
        const double *A_full, const double *b);
}
#include "rad_hacapk_bem.h"   // HACApK scalar BEM adapter (Laplace SL/DL Galerkin)
#include "rad_bem_galerkin.h" // Fast Galerkin SL/DL assembler
#include "rad_biot_savart_filaments.h" // Fast Biot-Savart H/A from finite-segment filaments
#include "rad_biot_savart_surface.h"   // Fast Biot-Savart B/A from triangulated surface
#include "rad_equivalence_source.h"    // Stratton-Chu equivalence-theorem reconstruction
#include "rad_average_field.h" // Closed-form cuboid average B (Wakao Part 6 §7)
#include "rad_stream_function.h" // (ACA+)+TSVD stream-function coil solver
#include "rad_peec_matrices.h"  // PEECMatrixBuilder for filament input
#include "rad_hdiv_vim.h"        // Symmetric HDiv-type VIM demag operator (N = B^T G B)
#include "rad_hacapk_hdiv.h"     // HACApK H-matrix for the HDiv-type VIM demag operator
#include <core/taskmanager.hpp>  // ngcore::ParallelFor / TaskManager (HDiv-VIM batched field, obs-parallel)

namespace py = pybind11;
using namespace pybind11::literals;

// ============================================================================
// Helper Functions
// ============================================================================

namespace {

/**
 * @brief Convert Python list/array to std::vector<double>
 */
std::vector<double> to_vector(const py::object& obj) {
    if (py::isinstance<py::array>(obj)) {
        auto arr = obj.cast<py::array_t<double>>();
        auto buf = arr.request();
        double* ptr = static_cast<double*>(buf.ptr);
        return std::vector<double>(ptr, ptr + buf.size);
    } else if (py::isinstance<py::list>(obj)) {
        return obj.cast<std::vector<double>>();
    } else if (py::isinstance<py::tuple>(obj)) {
        py::tuple tup = obj.cast<py::tuple>();
        std::vector<double> result;
        result.reserve(tup.size());
        for (size_t i = 0; i < tup.size(); i++) {
            result.push_back(tup[i].cast<double>());
        }
        return result;
    }
    throw std::runtime_error("Expected list, tuple, or numpy array");
}

/**
 * @brief Convert 2D Python list to flat array with vertex data
 * Returns vertex data and count
 */
std::pair<std::vector<double>, int> to_vertex_array(const py::list& vertices) {
    int nv = static_cast<int>(vertices.size());
    std::vector<double> flat;
    flat.reserve(nv * 3);

    for (const auto& v : vertices) {
        auto coord = to_vector(v.cast<py::object>());
        if (coord.size() != 3) {
            throw std::runtime_error("Each vertex must have 3 coordinates");
        }
        flat.insert(flat.end(), coord.begin(), coord.end());
    }

    return {flat, nv};
}

/**
 * @brief Convert C++ error code to Python exception if needed
 */
void check_error(int errCode) {
    if (errCode != 0) {
        char errText[2048] = {0};
        RadErrGetText(errText, errCode);
        throw std::runtime_error(std::string("Radia error: ") + errText);
    }
}

} // anonymous namespace


// ============================================================================
// Object Creation Functions
// ============================================================================

namespace radia_objects {

/**
 * @brief Create rectangular permanent magnet
 *
 * @param center Center point [x, y, z]
 * @param dimensions Block dimensions [Lx, Ly, Lz]
 * @param magnetization Magnetization vector [Mx, My, Mz] in A/m
 * @return Object handle
 */
int ObjRecMag(py::array_t<double> center,
              py::array_t<double> dimensions,
              py::array_t<double> magnetization) {

    auto c = center.unchecked<1>();
    auto d = dimensions.unchecked<1>();
    auto m = magnetization.unchecked<1>();

    if (c.size() != 3 || d.size() != 3 || m.size() != 3) {
        throw std::runtime_error("center, dimensions, and magnetization must have 3 elements");
    }

    double P[3] = {c(0), c(1), c(2)};
    double L[3] = {d(0), d(1), d(2)};
    double M[3] = {m(0), m(1), m(2)};

    int handle = 0;
    int err = RadObjRecMag(&handle, P, L, M);
    check_error(err);

    return handle;
}

/**
 * @brief Compute cross product of two 3D vectors
 */
static void cross3(const double* a, const double* b, double* result) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

/**
 * @brief Compute dot product of two 3D vectors
 */
static double dot3(const double* a, const double* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * @brief Normalize a 3D vector in place, return magnitude
 */
static double normalize3(double* v) {
    double mag = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (mag > 1e-15) {
        v[0] /= mag;
        v[1] /= mag;
        v[2] /= mag;
    }
    return mag;
}

/**
 * @brief Create hexahedral element from 8 vertices
 *
 * Face ordering follows ELF kkh array convention: Face 0-5
 *
 * ELF assigns faces based on TOPOLOGY (vertex connectivity from kkh array).
 * Each face is defined by the 4 vertices specified in the kkh array.
 *
 * @param vertices 8 vertices in CHEXA ordering
 * @param magnetization Magnetization vector [Mx, My, Mz] in A/m
 * @return Object handle
 */
int ObjHexahedron(py::list vertices, py::array_t<double> magnetization) {
    if (py::len(vertices) != 8) {
        throw std::runtime_error("Hexahedron requires exactly 8 vertices");
    }

    auto [flat_verts, nv] = to_vertex_array(vertices);
    auto m = magnetization.unchecked<1>();

    if (m.size() != 3) {
        throw std::runtime_error("magnetization must have 3 elements");
    }

    double M[3] = {m(0), m(1), m(2)};

    // Extract vertex coordinates (0-indexed)
    double v[8][3];
    for (int i = 0; i < 8; i++) {
        v[i][0] = flat_verts[i * 3 + 0];
        v[i][1] = flat_verts[i * 3 + 1];
        v[i][2] = flat_verts[i * 3 + 2];
    }

    // Netgen hexahedron vertex convention (1-indexed):
    //      5----8      Vertices 1-4: z- (bottom)
    //     /|   /|      Vertices 5-8: z+ (top)
    //    6----7 |
    //    | 1--|-4      Face vertex order: CCW when viewed from outside
    //    |/   |/       -> outward normal by right-hand rule
    //    2----3
    //
    // Face definitions (Netgen/ELF convention, 1-indexed):
    //   Face 1 (index 0): z- = 1-2-3-4  (bottom)
    //   Face 2 (index 1): x+ = 2-6-7-3  (right)
    //   Face 3 (index 2): y- = 1-5-6-2  (front)
    //   Face 4 (index 3): x- = 1-4-8-5  (left)
    //   Face 5 (index 4): y+ = 3-7-8-4  (back)
    //   Face 6 (index 5): z+ = 5-8-7-6  (top)
    // NOTE: 1-indexed vertex references (Radia C++ core convention).
    // v[] array above is 0-indexed for coordinate access.
    // These indices are passed directly to PolyhedronDLL which expects 1-indexed.
    static const int NETGEN_FACES[6][4] = {
        {1, 2, 3, 4},  // Face 0: z- (bottom)
        {2, 6, 7, 3},  // Face 1: x+ (right)
        {1, 5, 6, 2},  // Face 2: y- (front)
        {1, 4, 8, 5},  // Face 3: x- (left)
        {3, 7, 8, 4},  // Face 4: y+ (back)
        {5, 8, 7, 6}   // Face 5: z+ (top)
    };

    // Build flat face array in Netgen order
    std::vector<int> flatFaces;
    int faceLengths[6] = {4, 4, 4, 4, 4, 4};

    for (int face_idx = 0; face_idx < 6; face_idx++) {
        for (int v_idx = 0; v_idx < 4; v_idx++) {
            flatFaces.push_back(NETGEN_FACES[face_idx][v_idx]);
        }
    }

    int handle = 0;
    double M_LinCoef[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double J[3] = {0, 0, 0};
    double J_LinCoef[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    int err = RadObjPolyhdr(&handle, flat_verts.data(), nv,
                           flatFaces.data(), faceLengths, 6,
                           M, M_LinCoef, J, J_LinCoef);
    check_error(err);

    return handle;
}

/**
 * @brief Create tetrahedral element from 4 vertices
 *
 * @param vertices 4 vertices
 * @param magnetization Magnetization vector [Mx, My, Mz] in A/m
 * @return Object handle
 */
int ObjTetrahedron(py::list vertices, py::array_t<double> magnetization) {
    if (py::len(vertices) != 4) {
        throw std::runtime_error("Tetrahedron requires exactly 4 vertices");
    }

    auto [flat_verts, nv] = to_vertex_array(vertices);
    auto m = magnetization.unchecked<1>();

    if (m.size() != 3) {
        throw std::runtime_error("magnetization must have 3 elements");
    }

    double M[3] = {m(0), m(1), m(2)};

    // Tetrahedron face definitions (1-indexed vertices for each triangular face)
    // Following right-hand rule for outward normals
    static const int TET_FACES[4][3] = {
        {1, 3, 2},  // Base (opposite vertex 4)
        {1, 2, 4},  // Face opposite vertex 3
        {2, 3, 4},  // Face opposite vertex 1
        {3, 1, 4}   // Face opposite vertex 2
    };

    // Build flat face array and face length array
    std::vector<int> flatFaces;
    int faceLengths[4] = {3, 3, 3, 3};

    for (int f = 0; f < 4; f++) {
        for (int v = 0; v < 3; v++) {
            flatFaces.push_back(TET_FACES[f][v]);
        }
    }

    int handle = 0;
    double M_LinCoef[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double J[3] = {0, 0, 0};
    double J_LinCoef[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    int err = RadObjPolyhdr(&handle, flat_verts.data(), nv,
                           flatFaces.data(), faceLengths, 4,
                           M, M_LinCoef, J, J_LinCoef);
    check_error(err);

    return handle;
}

/**
 * @brief Create wedge/prism element from 6 vertices
 *
 * Wedge element with triangular top and bottom faces.
 * 5 faces total: 2 triangular (top/bottom) + 3 quadrilateral (sides)
 * 5 DOF for MSC method.
 *
 * Vertex convention (ELF MMB6T compatible):
 *      3-----5
 *     /|    /|
 *    / |   / |
 *   4-----+  |      Vertices 0-2: z- (bottom triangle)
 *   |  0--|--2      Vertices 3-5: z+ (top triangle)
 *   | /   | /       Bottom: v0, v1, v2 (CCW when viewed from below)
 *   |/    |/        Top: v3, v4, v5 (CCW when viewed from above)
 *   1-----+
 *
 * @param vertices 6 vertices in wedge order
 * @param magnetization Magnetization vector [Mx, My, Mz] in A/m
 * @return Object handle
 */
int ObjWedge(py::list vertices, py::array_t<double> magnetization) {
    if (py::len(vertices) != 6) {
        throw std::runtime_error("Wedge requires exactly 6 vertices");
    }

    auto [flat_verts, nv] = to_vertex_array(vertices);
    auto m = magnetization.unchecked<1>();

    if (m.size() != 3) {
        throw std::runtime_error("magnetization must have 3 elements");
    }

    double M[3] = {m(0), m(1), m(2)};

    // Wedge face definitions (1-indexed vertices)
    // Face 0: bottom triangle (v0, v2, v1) - outward normal points -z
    // Face 1: top triangle (v3, v4, v5) - outward normal points +z
    // Face 2: quad side (v0, v1, v4, v3)
    // Face 3: quad side (v1, v2, v5, v4)
    // Face 4: quad side (v2, v0, v3, v5)

    // Build flat face array (1-indexed for RadObjPolyhdr)
    std::vector<int> flatFaces;
    int faceLengths[5] = {3, 3, 4, 4, 4};

    // Face 0: bottom triangle (CCW from below = CW from above)
    flatFaces.push_back(1);  // v0
    flatFaces.push_back(3);  // v2
    flatFaces.push_back(2);  // v1

    // Face 1: top triangle (CCW from above)
    flatFaces.push_back(4);  // v3
    flatFaces.push_back(5);  // v4
    flatFaces.push_back(6);  // v5

    // Face 2: front quad (v0, v1, v4, v3)
    flatFaces.push_back(1);  // v0
    flatFaces.push_back(2);  // v1
    flatFaces.push_back(5);  // v4
    flatFaces.push_back(4);  // v3

    // Face 3: right quad (v1, v2, v5, v4)
    flatFaces.push_back(2);  // v1
    flatFaces.push_back(3);  // v2
    flatFaces.push_back(6);  // v5
    flatFaces.push_back(5);  // v4

    // Face 4: left quad (v2, v0, v3, v5)
    flatFaces.push_back(3);  // v2
    flatFaces.push_back(1);  // v0
    flatFaces.push_back(4);  // v3
    flatFaces.push_back(6);  // v5

    int handle = 0;
    double M_LinCoef[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double J[3] = {0, 0, 0};
    double J_LinCoef[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    int err = RadObjPolyhdr(&handle, flat_verts.data(), nv,
                           flatFaces.data(), faceLengths, 5,
                           M, M_LinCoef, J, J_LinCoef);
    check_error(err);

    return handle;
}

/**
 * @brief Create container for objects
 *
 * @param objects List of object handles
 * @return Container handle
 */
int ObjCnt(py::list objects) {
    std::vector<int> handles;
    for (const auto& obj : objects) {
        handles.push_back(obj.cast<int>());
    }

    int handle = 0;
    int err = RadObjCnt(&handle, handles.data(), static_cast<int>(handles.size()));
    check_error(err);

    return handle;
}

/**
 * @brief Create background field source with callback
 *
 * @param callback Python function that takes [x,y,z] and returns [Bx,By,Bz] in Tesla
 * @return Object handle
 */
int ObjBckg(py::function callback) {
    // Store callback as PyObject for C API
    PyObject* py_callback = callback.ptr();
    Py_INCREF(py_callback);

    int handle = 0;
    int err = RadObjBckgCF(&handle, py_callback);
    check_error(err);

    return handle;
}

} // namespace radia_objects


// ============================================================================
// Field Computation Functions
// ============================================================================

namespace radia_field {

/**
 * @brief Single-point field computation (internal helper)
 */
static py::object FldSingle(int obj, const std::string& field_type, const double* coords) {
    double result[6] = {0};
    int nResult = 0;

    std::string ft = field_type;
    if (ft == "phi") ft = "p";
    char* id = const_cast<char*>(ft.c_str());
    int err = RadFld(result, &nResult, obj, id, const_cast<double*>(coords), 1);
    check_error(err);

    if (nResult == 1) {
        return py::float_(result[0]);
    } else {
        py::array_t<double> arr(nResult);
        auto buf = arr.mutable_unchecked<1>();
        for (int i = 0; i < nResult; i++) {
            buf(i) = result[i];
        }
        return arr;
    }
}

/**
 * @brief Batch field computation for B or H (internal helper)
 *
 * Returns only the requested component (B or H), not both.
 */
static py::array_t<double> FldBatchBorH(int obj, const std::string& field_type,
                                         int n_points, double* pts) {
    std::vector<double> B_out(n_points * 3);
    std::vector<double> H_out(n_points * 3);

    {
        py::gil_scoped_release release;
        int err = RadFldBatch(B_out.data(), H_out.data(), n_points, pts, obj);
        check_error(err);
    }

    // Select requested field
    double* src = (field_type == "h") ? H_out.data() : B_out.data();

    py::array_t<double> result({n_points, 3});
    auto buf = result.mutable_unchecked<2>();
    for (int i = 0; i < n_points; i++) {
        buf(i, 0) = src[i * 3 + 0];
        buf(i, 1) = src[i * 3 + 1];
        buf(i, 2) = src[i * 3 + 2];
    }
    return result;
}

/**
 * @brief Unified field computation at single point or multiple points
 *
 * Auto-detects single point (shape (3,)) vs batch (shape (N,3)).
 * For batch mode, returns only the requested field type as an array.
 *
 * @param obj Object handle
 * @param field_type Field type: "b", "h", "a", "phi", "m", "bx", "by", "bz", etc.
 * @param points Evaluation point(s): [x,y,z] or array of shape (N,3)
 * @return Single point: scalar or [Fx,Fy,Fz]. Batch: array of shape (N,3) or (N,)
 */
py::object Fld(int obj, const std::string& field_type, py::array_t<double> points) {
    py::buffer_info buf = points.request();

    // --- Single point: shape (3,) ---
    if (buf.ndim == 1) {
        if (buf.shape[0] != 3) {
            throw std::runtime_error("Single point must have 3 coordinates");
        }
        double* p = static_cast<double*>(buf.ptr);
        return FldSingle(obj, field_type, p);
    }

    // --- Batch: shape (N, 3) ---
    if (buf.ndim != 2 || buf.shape[1] != 3) {
        throw std::runtime_error("points must be shape (3,) for single point or (N,3) for batch");
    }

    int n_points = static_cast<int>(buf.shape[0]);
    double* pts = static_cast<double*>(buf.ptr);

    // Dispatch by field type
    if (field_type == "b" || field_type == "h") {
        return FldBatchBorH(obj, field_type, n_points, pts);
    }

    if (field_type == "a") {
        std::vector<double> A_out(n_points * 3);
        {
            py::gil_scoped_release release;
            int err = RadFldA(A_out.data(), n_points, pts, obj);
            check_error(err);
        }
        py::array_t<double> result({n_points, 3});
        auto r = result.mutable_unchecked<2>();
        for (int i = 0; i < n_points; i++) {
            r(i, 0) = A_out[i * 3 + 0];
            r(i, 1) = A_out[i * 3 + 1];
            r(i, 2) = A_out[i * 3 + 2];
        }
        return result;
    }

    if (field_type == "phi") {
        std::vector<double> phi_out(n_points);
        {
            py::gil_scoped_release release;
            int err = RadFldPhi(phi_out.data(), n_points, pts, obj);
            check_error(err);
        }
        py::array_t<double> result(n_points);
        auto r = result.mutable_unchecked<1>();
        for (int i = 0; i < n_points; i++) {
            r(i) = phi_out[i];
        }
        return result;
    }

    // Fallback: per-point evaluation for "m", "bx", "by", "bz", etc.
    // Determine result dimension from first point
    double first_result[6] = {0};
    int nResult = 0;
    {
        std::string ft = field_type;
        if (ft == "phi") ft = "p";
        char* id = const_cast<char*>(ft.c_str());
        int err = RadFld(first_result, &nResult, obj, id, pts, 1);
        check_error(err);
    }

    if (nResult == 1) {
        // Scalar field (bx, by, bz, etc.)
        py::array_t<double> result(n_points);
        auto r = result.mutable_unchecked<1>();
        r(0) = first_result[0];
        for (int i = 1; i < n_points; i++) {
            double res[6] = {0};
            int nr = 0;
            std::string ft = field_type;
            char* id = const_cast<char*>(ft.c_str());
            int err = RadFld(res, &nr, obj, id, pts + i * 3, 1);
            check_error(err);
            r(i) = res[0];
        }
        return result;
    } else {
        // Vector field (m, etc.)
        py::array_t<double> result({n_points, 3});
        auto r = result.mutable_unchecked<2>();
        for (int j = 0; j < nResult && j < 3; j++) r(0, j) = first_result[j];
        for (int i = 1; i < n_points; i++) {
            double res[6] = {0};
            int nr = 0;
            std::string ft = field_type;
            char* id = const_cast<char*>(ft.c_str());
            int err = RadFld(res, &nr, obj, id, pts + i * 3, 1);
            check_error(err);
            for (int j = 0; j < nr && j < 3; j++) r(i, j) = res[j];
        }
        return result;
    }
}

/**
 * @brief Export field to VTS file
 *
 * @param obj Object handle
 * @param filename Output filename
 * @param x_range [xmin, xmax]
 * @param y_range [ymin, ymax]
 * @param z_range [zmin, zmax]
 * @param nx Number of points in x
 * @param ny Number of points in y
 * @param nz Number of points in z
 * @param include_B Include B field
 * @param include_H Include H field
 * @param unit_scale Coordinate scale factor
 * @return Filename
 */

} // namespace radia_field


// ============================================================================
// Material Functions
// ============================================================================

namespace radia_material {

/**
 * @brief Create linear magnetic material
 *
 * @param mu_r Relative permeability (isotropic) or [mu_par, mu_perp]
 * @param easy_axis Optional easy axis for anisotropic materials
 * @return Material handle
 */
int MatLin(py::object mu_r, py::object easy_axis = py::none()) {
    int handle = 0;

    if (py::isinstance<py::float_>(mu_r) || py::isinstance<py::int_>(mu_r)) {
        // Isotropic material: convert mu_r to chi (susceptibility)
        double mu_r_val = mu_r.cast<double>();
        if (mu_r_val < 1.0) {
            throw std::runtime_error("mu_r must be >= 1.0 (relative permeability)");
        }
        double ksi = mu_r_val - 1.0;  // chi = mu_r - 1
        int err = RadMatLinIso(&handle, ksi);
        check_error(err);
    } else {
        // Anisotropic material: convert [mu_par, mu_perp] to [chi_par, chi_perp]
        auto mu = to_vector(mu_r);
        if (mu.size() != 2) {
            throw std::runtime_error("For anisotropic material, mu_r must be [mu_par, mu_perp]");
        }
        if (mu[0] < 1.0 || mu[1] < 1.0) {
            throw std::runtime_error("mu_r values must be >= 1.0 (relative permeability)");
        }

        if (easy_axis.is_none()) {
            throw std::runtime_error("Anisotropic material requires easy_axis");
        }

        auto axis = to_vector(easy_axis);
        if (axis.size() != 3) {
            throw std::runtime_error("easy_axis must have 3 components");
        }

        double Ksi[2] = {mu[0] - 1.0, mu[1] - 1.0};  // chi = mu_r - 1
        int err = RadMatLinAniso(&handle, Ksi, axis.data());
        check_error(err);
    }

    return handle;
}

/**
 * @brief Create nonlinear isotropic material from B-H curve
 *
 * @param bh_data List of [H, B] pairs
 * @return Material handle
 */
int MatSatIsoTab(py::list bh_data) {
    int np = static_cast<int>(py::len(bh_data));
    std::vector<double> flat_data;
    flat_data.reserve(np * 2);

    for (const auto& pair : bh_data) {
        auto p = to_vector(pair.cast<py::object>());
        if (p.size() != 2) {
            throw std::runtime_error("Each B-H pair must have 2 elements [H, B]");
        }
        flat_data.push_back(p[0]);  // H
        flat_data.push_back(p[1]);  // B
    }

    int handle = 0;
    int err = RadMatSatIsoTab(&handle, flat_data.data(), np);
    check_error(err);

    return handle;
}

/**
 * @brief Apply material to object
 *
 * @param obj Object handle
 * @param mat Material handle
 * @return Result object handle
 */
int MatApl(int obj, int mat) {
    int result = 0;
    int err = RadMatApl(&result, obj, mat);
    check_error(err);
    return result;
}

// NOTE: MatSIBC removed (2026-02-13). Use Python SIBC (scipy Bessel) instead.

} // namespace radia_material


// ============================================================================
// Solver Functions
// ============================================================================

namespace radia_solver {

/**
 * @brief Solve magnetostatic problem
 *
 * @param obj Object or container handle
 * @param prec Convergence precision
 * @param max_iter Maximum iterations
 * @param method Solver method: 0=LU, 1=BiCGSTAB, 2=HACApK
 * @param image Image symmetry string (e.g., "+x", "-z", "+x-z") or empty
 * @return Result tuple [residual, ?, ?, iterations]
 */
py::tuple Solve(int obj, double prec, int max_iter, int method, const std::string& image = "") {
    double D[4] = {0};
    int n = 0;

    // Note: GIL release disabled because OpenMP threads in the solver
    // don't have proper Python thread state, causing crashes.
    const char* img = image.empty() ? nullptr : image.c_str();
    int err = RadSolve(D, &n, obj, prec, max_iter, method, img);
    check_error(err);

    return py::make_tuple(D[0], D[1], D[2], D[3]);
}

/**
 * @brief Build interaction matrix without solving
 *
 * @param obj Object or container handle
 * @param image Image symmetry string (e.g., "+x", "-z", "+x-z") or empty
 * @return Interaction matrix handle (for GetInteractMatrix)
 */
int BuildMatrix(int obj, const std::string& image = "") {
    int handle = 0;
    const char* img = image.empty() ? nullptr : image.c_str();
    int err = RadBuildMatrix(&handle, obj, img);
    check_error(err);
    return handle;
}

/**
 * @brief Get solve statistics
 * @return Dictionary with timing and iteration counts
 */
py::dict GetSolveStats() {
    double stats[20] = {0};
    int n = 0;

    int err = RadGetSolveStats(stats, &n);
    if (err != 0 || n == 0) {
        return py::dict();
    }

    py::dict result;
    if (n >= 4) {
        result["t_matrix_build"] = stats[0];
        result["t_linear_solve"] = stats[1];
        result["linear_iterations"] = static_cast<int>(stats[2]);
        result["nonl_iterations"] = static_cast<int>(stats[3]);
    }
    if (n >= 6) {
        result["taskmanager_enabled"] = (stats[4] > 0.5);
        result["num_threads"] = static_cast<int>(stats[5]);
    }
    if (n >= 7) {
        result["t_lu_decomp"] = stats[6];
    }
    if (n >= 8) {
        result["t_hmatrix_build"] = stats[7];
    }
    if (n >= 9) {
        result["deflation_cycles"] = static_cast<int>(stats[8]);
    }
    if (n >= 10) {
        result["deflation_alpha"] = stats[9];
    }

    return result;
}

/**
 * @brief Set HACApK parameters
 *
 * @param eps ACA+ tolerance
 * @param leaf_size Minimum cluster size
 * @param eta Admissibility parameter
 */
void SetHACApKParams(double eps, int leaf_size, double eta);  // Forward declaration

/**
 * @brief Set BiCGSTAB tolerance
 * @param tol Convergence tolerance
 */
void SetBiCGSTABTol(double tol);  // Forward declaration

// PreRelax REMOVED (2026-01-31) - Use BuildMatrix() instead
// The new API is: rad.BuildMatrix(obj, image='+x-z')
// or simply: rad.Solve(obj, prec, maxiter, method, image='+x-z')

/**
 * @brief Get interaction matrix as numpy array
 *
 * @param intrc_handle Interaction handle from BuildMatrix
 * @return Tuple (matrix as 2D numpy array, dof)
 */
py::tuple GetInteractMatrix(int intrc_handle) {
    // First call to get DOF only
    int dof = 0;
    int err = RadGetInteractMatrix(nullptr, &dof, intrc_handle);
    check_error(err);

    if (dof <= 0) {
        throw std::runtime_error("No interaction matrix built");
    }

    // Allocate and get matrix
    std::vector<double> matrix_data(static_cast<size_t>(dof) * dof);
    err = RadGetInteractMatrix(matrix_data.data(), &dof, intrc_handle);
    check_error(err);

    // Create numpy array from row-major FlatInteract data
    py::array_t<double> result({dof, dof});
    auto r = result.mutable_unchecked<2>();

    // FlatInteract is stored ROW-MAJOR: A[i][j] at index [i * dof + j]
    // This matches C/NumPy convention, so direct copy
    for (int i = 0; i < dof; i++) {
        for (int j = 0; j < dof; j++) {
            r(i, j) = matrix_data[i * dof + j];  // Row-major: A[target][source]
        }
    }

    return py::make_tuple(result, dof);
}

/**
 * @brief Get the yano-MSC cell-graph cycle (loop) basis as a numpy array
 * @param intrc_handle Interaction handle from BuildMatrix
 * @return Tuple (L as 2D numpy array (dof x nLoop), nLoop)
 */
py::tuple GetLoopBasis(int intrc_handle) {
    int nLoop = 0, dof = 0;
    int err = RadGetLoopBasis(nullptr, &nLoop, &dof, intrc_handle);
    check_error(err);
    if (dof <= 0) {
        throw std::runtime_error("No interaction matrix built");
    }
    int ncol = nLoop > 0 ? nLoop : 0;
    py::array_t<double> result({dof, ncol});
    if (nLoop > 0) {
        std::vector<double> data((size_t)dof * (size_t)nLoop);
        err = RadGetLoopBasis(data.data(), &nLoop, &dof, intrc_handle);
        check_error(err);
        auto r = result.mutable_unchecked<2>();
        // L is ROW-MAJOR: L[d * nLoop + c]
        for (int i = 0; i < dof; i++) {
            for (int c = 0; c < nLoop; c++) {
                r(i, c) = data[(size_t)i * (size_t)nLoop + (size_t)c];
            }
        }
    }
    return py::make_tuple(result, nLoop);
}

/**
 * @brief Get per-DOF hex face geometry as a numpy array
 * @param intrc_handle Interaction handle from BuildMatrix
 * @return numpy array (dof x 11): [elem_local, area, cx,cy,cz, nx,ny,nz(outward), ecx,ecy,ecz] per DOF
 */
py::array_t<double> GetFaceGeom(int intrc_handle) {
    int dof = 0;
    int err = RadGetFaceGeom(nullptr, &dof, intrc_handle);
    check_error(err);
    if (dof <= 0) {
        throw std::runtime_error("No interaction matrix built");
    }
    py::array_t<double> result({dof, 11});
    std::vector<double> data((size_t)dof * 11);
    err = RadGetFaceGeom(data.data(), &dof, intrc_handle);
    check_error(err);
    auto r = result.mutable_unchecked<2>();
    // G is ROW-MAJOR: G[d * 11 + j]
    for (int i = 0; i < dof; i++)
        for (int j = 0; j < 11; j++)
            r(i, j) = data[(size_t)i * 11 + (size_t)j];
    return result;
}

/**
 * @brief Per-hex centroid demag field + gradient functionals (the moment-formulation kernel)
 * @param intrc_handle Interaction handle from BuildMatrix
 * @return numpy array (nHex, 9, dof): comp k (Hx,Hy,Hz, gxx,gyy,gzz,gxy,gxz,gyz), source DOF g
 */
py::array_t<double> GetCentroidFieldGrad(int intrc_handle) {
    int nHex = 0, dof = 0;
    int err = RadGetCentroidFieldGrad(nullptr, &nHex, &dof, intrc_handle);
    check_error(err);
    if (nHex <= 0 || dof <= 0) {
        throw std::runtime_error("No hex DOFs in interaction matrix");
    }
    py::array_t<double> result({nHex, 9, dof});
    std::vector<double> data((size_t)nHex * 9 * dof);
    err = RadGetCentroidFieldGrad(data.data(), &nHex, &dof, intrc_handle);
    check_error(err);
    auto r = result.mutable_unchecked<3>();
    // C is ROW-MAJOR: C[(h*9 + k)*dof + g]
    for (int h = 0; h < nHex; h++)
        for (int k = 0; k < 9; k++)
            for (int g = 0; g < dof; g++)
                r(h, k, g) = data[((size_t)h * 9 + (size_t)k) * dof + (size_t)g];
    return result;
}

// Moment-yano system matrix A (dof x dof, row-major) + rhs (dof) for uniform linear chi + uniform applied
// field (hx,hy,hz).  Step-1 verification of the EIEM2 -> moment-yano upgrade (vs examples/vim prototype).
py::tuple BuildMomentSystem(int intrc_handle, double chi, double hx, double hy, double hz) {
    double Happ[3] = {hx, hy, hz};
    int dof = 0;
    int err = RadBuildMomentSystem(chi, Happ, nullptr, nullptr, &dof, intrc_handle);
    check_error(err);
    if (dof <= 0) throw std::runtime_error("No DOFs in interaction matrix");
    py::array_t<double> Aout({dof, dof});
    py::array_t<double> rout(dof);
    err = RadBuildMomentSystem(chi, Happ, (double*)Aout.request().ptr, (double*)rout.request().ptr, &dof, intrc_handle);
    check_error(err);
    return py::make_tuple(Aout, rout, dof);
}

/**
 * @brief Densify the actual HACApK (ACA+) operator as a numpy array
 * @param intrc_handle Interaction handle from BuildMatrix
 * @return Tuple (matrix as 2D numpy array, dof)
 */
py::tuple HMatrixDensify(int intrc_handle) {
    int dof = 0;
    int err = RadHMatrixDensify(nullptr, &dof, intrc_handle);
    check_error(err);

    if (dof <= 0) {
        throw std::runtime_error("HMatrixDensify: no operator (HACApK unavailable or build failed)");
    }

    std::vector<double> matrix_data(static_cast<size_t>(dof) * dof);
    err = RadHMatrixDensify(matrix_data.data(), &dof, intrc_handle);
    check_error(err);

    py::array_t<double> result({dof, dof});
    auto r = result.mutable_unchecked<2>();
    for (int i = 0; i < dof; i++) {
        for (int j = 0; j < dof; j++) {
            r(i, j) = matrix_data[i * dof + j];  // row-major A[target][source]
        }
    }

    return py::make_tuple(result, dof);
}

/**
 * @brief Phase 4: H-LU smoke test on a real Radia HACApK tree.
 *
 * Builds a fresh HACApK manager from the interaction handle, computes
 * y = A * x via MatVec (in HACApK format), converts the leafmtxp to our
 * internal column-major layout, factors via cHACApK_hlu_decomp, solves
 * x_solved = A^-1 y, and returns max relative error |x_solved - x| / max|x|.
 *
 * Returns:
 *   < 0    : driver failed (see negative sentinels in cHACApK_hlu_run_on_hacapk)
 *   ~ 1e-10..1e-15 : success at near-machine precision
 *   >  1e-3  : the H-matrix had rk leaves whose recompression lost accuracy
 *              (Phase 3.5 tolerance is 1e-14 default; ACA-truncated leaves
 *              add to that floor proportional to aca_eps)
 */
double HLUTestOnHACApK(int intrc_handle) {
    return RadHLUTestOnHACApK(intrc_handle);
}

/**
 * @brief Phase 4 debug: materialize the post-convert tree as a dense matrix.
 * Returns (A_perm, lod) where A_perm is in permuted ordering.
 */
py::tuple HLUDebugMaterialize(int intrc_handle) {
    // Pre-query: build manager to get dof count
    // Use a hack: allocate generous buffer (say up to 10000 dof for diagnostics)
    int max_dof = 10000;
    std::vector<double> A_buf((size_t)max_dof * max_dof);
    std::vector<int> lod_buf(max_dof);
    int nd = 0;
    int ok = RadHLUDebugMaterialize(intrc_handle, A_buf.data(), lod_buf.data(), &nd);
    if (!ok || nd <= 0) {
        throw std::runtime_error("HLUDebugMaterialize failed");
    }
    py::array_t<double> A_result({nd, nd});
    auto Ar = A_result.mutable_unchecked<2>();
    // A_buf is column-major nd x nd; convert to row-major for numpy
    for (int i = 0; i < nd; i++)
        for (int j = 0; j < nd; j++)
            Ar(i, j) = A_buf[(size_t)i + (size_t)j * (size_t)nd];
    py::array_t<int> lod_result(nd);
    auto lr = lod_result.mutable_unchecked<1>();
    for (int i = 0; i < nd; i++) lr(i) = lod_buf[i];
    return py::make_tuple(A_result, lod_result);
}


} // namespace radia_solver


// ============================================================================
// Utility Functions
// ============================================================================

namespace radia_utility {

/**
 * @brief Delete object
 * @param obj Object handle
 */
void UtiDel(int obj) {
    int n = 0;
    int err = RadUtiDel(&n, obj);
    check_error(err);
}

/**
 * @brief Delete all objects
 */
void UtiDelAll() {
    int n = 0;
    int err = RadUtiDelAll(&n);
    check_error(err);
}

/**
 * @brief Get object geometry volume
 * @param obj Object handle
 * @return Volume
 */
double ObjGeoVol(int obj) {
    double v = 0;
    int err = RadObjGeoVol(&v, obj);
    check_error(err);
    return v;
}

/**
 * @brief Get number of degrees of freedom
 * @param obj Object handle
 * @return Number of DOF
 */
int ObjDegFre(int obj) {
    int num = 0;
    int err = RadObjDegFre(&num, obj);
    check_error(err);
    return num;
}

} // namespace radia_utility


// ============================================================================
// Additional Object Creation Functions
// ============================================================================

namespace radia_objects_ext {

/**
 * @brief Create extruded polygon block
 */
int ObjThckPgn(double xc, double lx, py::list polygon, const std::string& axis, py::array_t<double> magnetization) {
    std::vector<double> flat_verts;
    for (const auto& pt : polygon) {
        auto coord = to_vector(pt.cast<py::object>());
        if (coord.size() != 2) {
            throw std::runtime_error("Each polygon point must have 2 coordinates");
        }
        flat_verts.insert(flat_verts.end(), coord.begin(), coord.end());
    }
    int nv = static_cast<int>(polygon.size());

    auto m = magnetization.unchecked<1>();
    double M[3] = {m(0), m(1), m(2)};

    char a = axis[0];

    int handle = 0;
    int err = RadObjThckPgn(&handle, xc, lx, flat_verts.data(), nv, a, M);
    check_error(err);
    return handle;
}

/**
 * @brief Create cylindrical magnet
 */
int ObjCylMag(py::array_t<double> center, double r, double h, int nseg,
              const std::string& axis, py::array_t<double> magnetization) {
    auto c = center.unchecked<1>();
    auto m = magnetization.unchecked<1>();
    double P[3] = {c(0), c(1), c(2)};
    double M[3] = {m(0), m(1), m(2)};
    char a = axis[0];

    int handle = 0;
    int err = RadObjCylMag(&handle, P, r, h, nseg, a, M);
    check_error(err);
    return handle;
}

/**
 * @brief Create rectangular current block
 */
int ObjRecCur(py::array_t<double> center, py::array_t<double> dimensions, py::array_t<double> current_density) {
    auto c = center.unchecked<1>();
    auto d = dimensions.unchecked<1>();
    auto j = current_density.unchecked<1>();
    double P[3] = {c(0), c(1), c(2)};
    double L[3] = {d(0), d(1), d(2)};
    double J[3] = {j(0), j(1), j(2)};

    int handle = 0;
    int err = RadObjRecCur(&handle, P, L, J);
    check_error(err);
    return handle;
}

/**
 * @brief Create arc current coil
 */
int ObjArcCur(py::array_t<double> center, py::array_t<double> radii, py::array_t<double> phi,
              double h, int nseg, const std::string& man_auto, const std::string& axis, double j) {
    auto c = center.unchecked<1>();
    auto r = radii.unchecked<1>();
    auto p = phi.unchecked<1>();
    double P[3] = {c(0), c(1), c(2)};
    double R[2] = {r(0), r(1)};
    double Phi[2] = {p(0), p(1)};
    char ma = man_auto[0];
    char a = axis[0];

    int handle = 0;
    int err = RadObjArcCur(&handle, P, R, Phi, h, nseg, ma, a, j);
    check_error(err);
    return handle;
}

/**
 * @brief Create racetrack coil
 */
int ObjRaceTrk(py::array_t<double> center, py::array_t<double> radii, py::array_t<double> lengths,
               double h, int nseg, const std::string& man_auto, const std::string& axis, double j) {
    auto c = center.unchecked<1>();
    auto r = radii.unchecked<1>();
    auto l = lengths.unchecked<1>();
    double P[3] = {c(0), c(1), c(2)};
    double R[2] = {r(0), r(1)};
    double L[2] = {l(0), l(1)};
    char ma = man_auto[0];
    char a = axis[0];

    int handle = 0;
    int err = RadObjRaceTrk(&handle, P, R, L, h, nseg, ma, a, j);
    check_error(err);
    return handle;
}

/**
 * @brief Create filament current
 */
int ObjFlmCur(py::list points, double current) {
    auto [flat_pts, np] = to_vertex_array(points);

    int handle = 0;
    int err = RadObjFlmCur(&handle, flat_pts.data(), np, current);
    check_error(err);
    return handle;
}

/**
 * @brief Add objects to container
 */
void ObjAddToCnt(int cnt, py::list objects) {
    std::vector<int> handles;
    for (const auto& obj : objects) {
        handles.push_back(obj.cast<int>());
    }
    int err = RadObjAddToCnt(cnt, handles.data(), static_cast<int>(handles.size()));
    check_error(err);
}

/**
 * @brief Get container size
 */
int ObjCntSize(int cnt) {
    int n = 0;
    int err = RadObjCntSize(&n, cnt);
    check_error(err);
    return n;
}

/**
 * @brief Get container contents
 */
py::list ObjCntStuf(int cnt) {
    int n = 0;
    RadObjCntSize(&n, cnt);

    std::vector<int> objs(n);
    int err = RadObjCntStuf(objs.data(), cnt);
    check_error(err);

    py::list result;
    for (int i = 0; i < n; i++) {
        result.append(objs[i]);
    }
    return result;
}

/**
 * @brief Duplicate object
 */
int ObjDpl(int obj, const std::string& opt = "") {
    int handle = 0;
    char opt_str[256];
    strcpy(opt_str, opt.c_str());
    int err = RadObjDpl(&handle, obj, opt_str);
    check_error(err);
    return handle;
}

/**
 * @brief Get object magnetization
 * For a single object: returns dict with 'center' and 'magnetization' tuples.
 * For a container: returns list of (center_tuple, magnetization_tuple) tuples.
 */
py::object ObjM(int obj) {
    // Get DOF count to determine buffer size
    // Each element outputs 6 doubles (3 center + 3 magnetization)
    // DOF per element is 3 or 6, so numElements <= numDOF/3
    // Buffer needed: numElements * 6 <= numDOF * 2
    int numDOF = 0;
    RadObjDegFre(&numDOF, obj);
    int bufSize = (numDOF > 0) ? (numDOF * 2 + 6) : 6;

    std::vector<double> M(bufSize, 0.0);
    int arMesh[20] = {0};
    int err = RadObjM(M.data(), arMesh, obj);
    check_error(err);

    // arMesh[0] = NumDims, arMesh[1..NumDims] = Dims
    // Single element: Dims = {3, 2, 1}, NumDims = 3
    // N elements: Dims = {3, 2, N}, NumDims = 3
    int numDims = arMesh[0];
    int numPoints = 1;
    if(numDims >= 3) numPoints = arMesh[3];

    if(numPoints <= 1) {
        // Single object
        py::dict result;
        result["center"] = py::make_tuple(M[0], M[1], M[2]);
        result["magnetization"] = py::make_tuple(M[3], M[4], M[5]);
        return result;
    } else {
        // Container: return list of (center, magnetization) tuples
        py::list result;
        for(int i = 0; i < numPoints; i++) {
            py::tuple center = py::make_tuple(M[i*6], M[i*6+1], M[i*6+2]);
            py::tuple magn = py::make_tuple(M[i*6+3], M[i*6+4], M[i*6+5]);
            result.append(py::make_tuple(center, magn));
        }
        return result;
    }
}

/**
 * @brief Set object magnetization
 */
void ObjSetM(int obj, py::array_t<double> magnetization) {
    auto m = magnetization.unchecked<1>();
    double M[3] = {m(0), m(1), m(2)};
    int err = RadObjSetM(obj, M);
    check_error(err);
}

/**
 * @brief Scale current
 */
void ObjScaleCur(int obj, double scale) {
    int err = RadObjScaleCur(obj, scale);
    check_error(err);
}

} // namespace radia_objects_ext


// ============================================================================
// Transformation Functions
// ============================================================================

namespace radia_transform {

int TrfTrsl(py::array_t<double> vector) {
    auto v = vector.unchecked<1>();
    double V[3] = {v(0), v(1), v(2)};

    int handle = 0;
    int err = RadTrfTrsl(&handle, V);
    check_error(err);
    return handle;
}

int TrfRot(py::array_t<double> point, py::array_t<double> vector, double phi) {
    auto p = point.unchecked<1>();
    auto v = vector.unchecked<1>();
    double P[3] = {p(0), p(1), p(2)};
    double V[3] = {v(0), v(1), v(2)};

    int handle = 0;
    int err = RadTrfRot(&handle, P, V, phi);
    check_error(err);
    return handle;
}

// TrfPlSym REMOVED (2026-01-31) - Use Image symmetry instead
// Use: Solve(..., image="+x") or BuildMatrix(obj, image="+x")

int TrfInv() {
    int handle = 0;
    int err = RadTrfInv(&handle);
    check_error(err);
    return handle;
}

int TrfCmbL(int orig_trf, int trf) {
    int handle = 0;
    int err = RadTrfCmbL(&handle, orig_trf, trf);
    check_error(err);
    return handle;
}

int TrfCmbR(int orig_trf, int trf) {
    int handle = 0;
    int err = RadTrfCmbR(&handle, orig_trf, trf);
    check_error(err);
    return handle;
}

// TrfMlt REMOVED (2026-01-31) - Use Image symmetry instead
// TrfMlt shared DOFs between original and virtual elements, which is incorrect for MSC 6DOF hexahedra
// Use: Solve(..., image="+x-z") or BuildMatrix(obj, image="+x-z")

int TrfOrnt(int obj, int trf) {
    int handle = 0;
    int err = RadTrfOrnt(&handle, obj, trf);
    check_error(err);
    return handle;
}

// TrfZerPara REMOVED (2026-01-31) - Use Image symmetry instead
// TrfZerPerp REMOVED (2026-01-31) - Use Image symmetry instead
// Use: Solve(..., image="+x-z") or BuildMatrix(obj, image="+x-z")

} // namespace radia_transform


// ============================================================================
// Additional Material Functions
// ============================================================================

namespace radia_material_ext {

int MatPM(double Br, double Hc, py::array_t<double> easy_axis) {
    auto ea = easy_axis.unchecked<1>();
    double EA[3] = {ea(0), ea(1), ea(2)};

    int handle = 0;
    int err = RadMatPM(&handle, Br, Hc, EA);
    check_error(err);
    return handle;
}

int MatMagFixed(py::array_t<double> magnetization) {
    auto m = magnetization.unchecked<1>();
    double M[3] = {m(0), m(1), m(2)};

    int handle = 0;
    int err = RadMatMagFixed(&handle, M);
    check_error(err);
    return handle;
}

int MatMagLinear(double Br, double Hc, py::array_t<double> easy_axis) {
    auto ea = easy_axis.unchecked<1>();
    double EA[3] = {ea(0), ea(1), ea(2)};

    int handle = 0;
    int err = RadMatMagLinear(&handle, Br, Hc, EA);
    check_error(err);
    return handle;
}

int MatMagCurve(py::list bh_data, py::array_t<double> easy_axis) {
    std::vector<double> flat_data;
    for (const auto& pair : bh_data) {
        auto p = pair.cast<py::list>();
        flat_data.push_back(p[0].cast<double>());
        flat_data.push_back(p[1].cast<double>());
    }
    int np = static_cast<int>(bh_data.size());

    auto ea = easy_axis.unchecked<1>();
    double EA[3] = {ea(0), ea(1), ea(2)};

    int handle = 0;
    int err = RadMatMagCurve(&handle, flat_data.data(), np, EA);
    check_error(err);
    return handle;
}

int MatSatIsoFrm(py::list params) {
    double KsiMs1[2] = {0, 0};
    double KsiMs2[2] = {0, 0};
    double KsiMs3[2] = {0, 0};

    int n = static_cast<int>(params.size());
    if (n >= 1) {
        auto p1 = params[0].cast<py::list>();
        KsiMs1[0] = p1[0].cast<double>();
        KsiMs1[1] = p1[1].cast<double>();
    }
    if (n >= 2) {
        auto p2 = params[1].cast<py::list>();
        KsiMs2[0] = p2[0].cast<double>();
        KsiMs2[1] = p2[1].cast<double>();
    }
    if (n >= 3) {
        auto p3 = params[2].cast<py::list>();
        KsiMs3[0] = p3[0].cast<double>();
        KsiMs3[1] = p3[1].cast<double>();
    }

    int handle = 0;
    int err = RadMatSatIsoFrm(&handle, KsiMs1, KsiMs2, KsiMs3);
    check_error(err);
    return handle;
}

int MatSatAniso(py::list data_par, py::list data_per) {
    std::vector<double> flat_par, flat_per;
    for (const auto& pair : data_par) {
        auto p = pair.cast<py::list>();
        flat_par.push_back(p[0].cast<double>());
        flat_par.push_back(p[1].cast<double>());
    }
    for (const auto& pair : data_per) {
        auto p = pair.cast<py::list>();
        flat_per.push_back(p[0].cast<double>());
        flat_per.push_back(p[1].cast<double>());
    }

    int handle = 0;
    int err = RadMatSatAniso(&handle, flat_par.data(), static_cast<int>(data_par.size()),
                             flat_per.data(), static_cast<int>(data_per.size()));
    check_error(err);
    return handle;
}

int MatSatLamTab(py::list mh_data, double packing_factor, py::array_t<double> normal) {
    std::vector<double> flat_data;
    for (const auto& pair : mh_data) {
        auto p = pair.cast<py::list>();
        flat_data.push_back(p[0].cast<double>());
        flat_data.push_back(p[1].cast<double>());
    }
    int np = static_cast<int>(mh_data.size());

    auto n = normal.unchecked<1>();
    double N[3] = {n(0), n(1), n(2)};

    int handle = 0;
    int err = RadMatSatLamTab(&handle, flat_data.data(), np, packing_factor, N);
    check_error(err);
    return handle;
}

py::array_t<double> MatMvsH(int obj, const std::string& component, py::array_t<double> h_field) {
    auto h = h_field.unchecked<1>();
    double H[3] = {h(0), h(1), h(2)};

    double M[3] = {0};
    int nM = 3;
    char comp[16];
    strcpy(comp, component.c_str());

    int err = RadMatMvsH(M, &nM, obj, comp, H);
    check_error(err);

    py::array_t<double> result(nM);
    auto r = result.mutable_unchecked<1>();
    for (int i = 0; i < nM; i++) {
        r(i) = M[i];
    }
    return result;
}

/**
 * @brief Create energy-based vector hysteresis material with table-based shape functions
 *
 * @param K Number of partial polarizations (play operators)
 * @param chi Pinning strengths chi_k [A/m], array of K values
 * @param f_k_tables List of K tuples (r_array, f_array), shape function tables
 * @param eps Regularization parameter for smoothed norm (default 1e-8)
 * @return Material handle
 */
int MatEnergyHysteresis(int K, py::array_t<double> chi,
                         py::list f_k_tables, double eps) {
    auto a_chi = chi.unchecked<1>();

    if (a_chi.shape(0) != K) {
        throw std::runtime_error("chi array must have length K");
    }
    if ((int)f_k_tables.size() != K) {
        throw std::runtime_error("f_k_tables must have length K");
    }

    // Copy chi
    std::vector<double> v_chi(K);
    for (int k = 0; k < K; k++) {
        v_chi[k] = a_chi(k);
    }

    // Extract table data and flatten for C interface
    std::vector<int> table_sizes(K);
    std::vector<double> r_flat, f_flat;

    for (int k = 0; k < K; k++) {
        py::tuple tab = f_k_tables[k].cast<py::tuple>();
        if (tab.size() != 2) {
            throw std::runtime_error("Each f_k_tables entry must be a (r, f) tuple");
        }
        py::array_t<double> r_arr = tab[0].cast<py::array_t<double>>();
        py::array_t<double> f_arr = tab[1].cast<py::array_t<double>>();
        auto r = r_arr.unchecked<1>();
        auto f = f_arr.unchecked<1>();

        if (r.shape(0) != f.shape(0)) {
            throw std::runtime_error("r and f arrays must have same length");
        }

        int n = (int)r.shape(0);
        table_sizes[k] = n;
        for (int i = 0; i < n; i++) {
            r_flat.push_back(r(i));
            f_flat.push_back(f(i));
        }
    }

    int handle = 0;
    int err = RadMatEnergyHysteresis(&handle, K, v_chi.data(),
                                      r_flat.data(), f_flat.data(),
                                      table_sizes.data(), eps);
    check_error(err);
    return handle;
}

/**
 * @brief Create a direct B-input play hysteresis material.
 *
 * Unlike energy-based hysteresis, the play model evaluates B->H directly
 * in O(K) without Newton iteration. Shape functions f_k can be negative.
 * Inverse (H->B) uses Newton with analytical Jacobian.
 *
 * @param K Number of play operators
 * @param eta Play thresholds [Tesla], array of K values
 * @param f_k_tables List of K tuples (r_array, f_array), shape function tables
 * @return Material handle
 */
int MatPlayHysteresis(int K, py::array_t<double> eta,
                       py::list f_k_tables) {
    auto a_eta = eta.unchecked<1>();

    if (a_eta.shape(0) != K) {
        throw std::runtime_error("eta array must have length K");
    }
    if ((int)f_k_tables.size() != K) {
        throw std::runtime_error("f_k_tables must have length K");
    }

    // Copy eta
    std::vector<double> v_eta(K);
    for (int k = 0; k < K; k++) {
        v_eta[k] = a_eta(k);
    }

    // Extract table data and flatten for C interface
    std::vector<int> table_sizes(K);
    std::vector<double> r_flat, f_flat;

    for (int k = 0; k < K; k++) {
        py::tuple tab = f_k_tables[k].cast<py::tuple>();
        if (tab.size() != 2) {
            throw std::runtime_error("Each f_k_tables entry must be a (r, f) tuple");
        }
        py::array_t<double> r_arr = tab[0].cast<py::array_t<double>>();
        py::array_t<double> f_arr = tab[1].cast<py::array_t<double>>();
        auto r = r_arr.unchecked<1>();
        auto f = f_arr.unchecked<1>();

        if (r.shape(0) != f.shape(0)) {
            throw std::runtime_error("r and f arrays must have same length");
        }

        int n = (int)r.shape(0);
        table_sizes[k] = n;
        for (int i = 0; i < n; i++) {
            r_flat.push_back(r(i));
            f_flat.push_back(f(i));
        }
    }

    int handle = 0;
    int err = RadMatPlayHysteresis(&handle, K, v_eta.data(),
                                    r_flat.data(), f_flat.data(),
                                    table_sizes.data());
    check_error(err);
    return handle;
}

/**
 * @brief Save the internal state of an energy hysteresis material.
 *
 * Returns the state as a flat numpy array. Use with MatHysRestoreState()
 * to save/restore state during Picard iteration for FEM hysteresis solves.
 */
py::array_t<double> MatHysSaveState(int mat) {
    // First query the size
    int len = 0;
    int err = RadMatHysSaveState(mat, nullptr, &len);
    if (err < 0) throw std::runtime_error("MatHysSaveState: invalid material handle");

    // Allocate and fill
    py::array_t<double> state(len);
    auto s = state.mutable_unchecked<1>();
    err = RadMatHysSaveState(mat, s.mutable_data(0), &len);
    if (err < 0) throw std::runtime_error("MatHysSaveState: save failed");
    return state;
}

/**
 * @brief Restore the internal state of an energy hysteresis material.
 */
void MatHysRestoreState(int mat, py::array_t<double> state) {
    auto s = state.unchecked<1>();
    int err = RadMatHysRestoreState(mat, s.data(0), (int)s.shape(0));
    if (err < 0) throw std::runtime_error("MatHysRestoreState: restore failed");
}

/**
 * @brief Commit the current state for the next time step.
 *
 * After Picard iteration converges, call this to commit the converged
 * state as the reference for the next quasi-static step.
 */
void MatHysCommitState(int mat) {
    int err = RadMatHysCommitState(mat);
    if (err < 0) throw std::runtime_error("MatHysCommitState: commit failed");
}

/**
 * Get the reversible reluctivity nu_rev for energy-based decomposition.
 * H = nu_rev * B + H_irr(B)
 */
double MatHysGetNuRev(int mat) {
    double nu_rev = 0;
    int err = RadMatHysGetNuRev(mat, &nu_rev);
    if (err < 0) throw std::runtime_error("MatHysGetNuRev: not a Play hysteresis material");
    return nu_rev;
}

/**
 * Compute irreversible field H_irr = H(B) - nu_rev * B.
 * For Hantila solver: constant matrix uses nu_rev, residual uses H_irr.
 */
py::array_t<double> MatHysIrreversible(int mat, py::array_t<double> B) {
    auto b = B.unchecked<1>();
    if (b.shape(0) != 3) throw std::runtime_error("B must have 3 components");
    double pB[3] = {b(0), b(1), b(2)};
    double pHirr[3] = {0, 0, 0};
    int err = RadMatHysIrreversible(mat, pB, pHirr);
    if (err < 0) throw std::runtime_error("MatHysIrreversible: not a Play hysteresis material");
    auto result = py::array_t<double>(3);
    auto r = result.mutable_unchecked<1>();
    r(0) = pHirr[0]; r(1) = pHirr[1]; r(2) = pHirr[2];
    return result;
}

} // namespace radia_material_ext


// ============================================================================
// Additional Solver Functions
// ============================================================================

// Opt-in "improved yano-type" surface-charge (MSC) kernel flag (pyramid-cloud + pyramid-centroid eval);
// defined in rad_polyhedron.cpp (global scope), read by the yano-MSC matrix-build paths in
// rad_interaction.cpp.  Declared here at GLOBAL scope so the in-namespace SolverConfig/GetSolverConfig
// usages resolve to the global definition (not a namespace-scoped symbol).
extern bool g_yano_pyramid_cloud;
extern bool g_yano_no_center_charge;
extern double g_yano_eval_alpha;
extern bool g_yano_moment;

namespace radia_solver_ext {

py::tuple SolveNonl(int obj, double prec, int max_iter, int method, int nonl_method, const std::string& image = "") {
    double D[4] = {0};
    int n = 4;

    // Note: GIL release disabled because OpenMP threads in the solver
    // don't have proper Python thread state, causing crashes.
    const char* img = image.empty() ? nullptr : image.c_str();
    int err = RadSolveNonl(D, &n, obj, prec, max_iter, method, nonl_method, img);
    check_error(err);

    return py::make_tuple(D[0], D[1], D[2], D[3]);
}

void SetHACApKParams(double eps, int leaf_size, double eta) {
    int n = 0;
    int err = RadSetHACApKParams(&n, eps, leaf_size, eta);
    check_error(err);
}

void SetHMatrixEpsilon(double eps) {
    int n = 0;
    int err = RadSetHMatrixEpsilon(&n, eps);
    check_error(err);
}

py::dict GetHACApKStats() {
    double dOut[20] = {0};
    int nOut = 0;

    int err = RadGetHACApKStats(dOut, &nOut);
    if (err != 0 || nOut == 0) {
        return py::none();
    }

    py::dict result;
    result["n_lowrank"] = static_cast<int>(dOut[0]);
    result["n_dense"] = static_cast<int>(dOut[1]);
    result["max_rank"] = static_cast<int>(dOut[2]);
    result["n_leaves"] = static_cast<int>(dOut[3]);
    result["n_dof"] = static_cast<int>(dOut[4]);
    result["compression"] = dOut[5];
    result["build_time"] = dOut[6];
    result["hmatrix_build_time"] = dOut[7];
    result["linear_iterations"] = static_cast<int>(dOut[9]);
    result["memory_mb"] = dOut[10];
    result["dense_memory_mb"] = dOut[11];
    return result;
}

void SetBiCGSTABTol(double tol) {
    int n = 0;
    int err = RadSetBiCGSTABTol(&n, tol);
    check_error(err);
}

double GetBiCGSTABTol() {
    double tol = 0;
    int err = RadGetBiCGSTABTol(&tol);
    check_error(err);
    return tol;
}

void SetRelaxParam(double relax) {
    int n = 0;
    int err = RadSetRelaxParam(&n, relax);
    check_error(err);
}

double GetRelaxParam() {
    double relax = 0;
    int err = RadGetRelaxParam(&relax);
    check_error(err);
    return relax;
}

void SetKeepMagnetization(bool keep) {
    int n = 0;
    int err = RadSetKeepMagnetization(&n, keep ? 1 : 0);
    check_error(err);
}

bool GetKeepMagnetization() {
    int keep = 0;
    int err = RadGetKeepMagnetization(&keep);
    check_error(err);
    return keep != 0;
}

void SetNewtonMethod(bool use_newton) {
    int n = 0;
    int err = RadSetNewtonMethod(&n, use_newton ? 1 : 0);
    check_error(err);
}

bool GetNewtonMethod() {
    int use_newton = 0;
    int err = RadGetNewtonMethod(&use_newton);
    check_error(err);
    return use_newton != 0;
}

void SetNewtonDamping(bool enabled = true, int max_iter = 5, double min_omega = 0.01) {
    int n = 0;
    int err = RadSetNewtonDamping(&n, enabled ? 1 : 0, max_iter, min_omega);
    check_error(err);
}

py::dict GetNewtonDampingStats() {
    int enabled = 0;
    int max_iter = 0;
    double min_omega = 0.0;
    int err = RadGetNewtonDampingStats(&enabled, &max_iter, &min_omega);
    check_error(err);

    py::dict stats;
    stats["enabled"] = (enabled != 0);
    stats["max_iter"] = max_iter;
    stats["min_omega"] = min_omega;
    return stats;
}

void SetBInputNewton(bool enabled) {
    int n = 0;
    int err = RadSetBInputNewton(&n, enabled ? 1 : 0);
    check_error(err);
}

bool GetBInputNewton() {
    int enabled = 0;
    int err = RadGetBInputNewton(&enabled);
    check_error(err);
    return enabled != 0;
}

void SetBInputHantila(bool enabled) {
    int n = 0;
    int err = RadSetBInputHantila(&n, enabled ? 1 : 0);
    check_error(err);
}

bool GetBInputHantila() {
    int enabled = 0;
    int err = RadGetBInputHantila(&enabled);
    check_error(err);
    return enabled != 0;
}

void SetHantilaAlpha(double alpha) {
    int n = 0;
    int err = RadSetHantilaAlpha(&n, alpha);
    check_error(err);
}

double GetHantilaAlpha() {
    double alpha = 0.0;
    int err = RadGetHantilaAlpha(&alpha);
    check_error(err);
    return alpha;
}

void SetHantilaRelax(double relax) {
    int n = 0;
    int err = RadSetHantilaRelax(&n, relax);
    check_error(err);
}

double GetHantilaRelax() {
    double relax = 0.0;
    int err = RadGetHantilaRelax(&relax);
    check_error(err);
    return relax;
}

// SetIMASymmetry, BuildIMAMatrix REMOVED (2026-01-31)
// Use BuildMatrix(obj, image="+x-z") or Solve(obj, ..., image="+x-z") instead

// ---- Unified SolverConfig / GetSolverConfig ----
// (g_yano_pyramid_cloud is declared at global scope above the namespace; usages below resolve to it.)

void SolverConfig(py::kwargs kwargs) {
    // HACApK parameters
    if (kwargs.contains("hacapk_eps") || kwargs.contains("hacapk_leaf") || kwargs.contains("hacapk_eta")) {
        double eps = kwargs.contains("hacapk_eps") ? kwargs["hacapk_eps"].cast<double>() : -1;
        int leaf = kwargs.contains("hacapk_leaf") ? kwargs["hacapk_leaf"].cast<int>() : -1;
        double eta = kwargs.contains("hacapk_eta") ? kwargs["hacapk_eta"].cast<double>() : -1;
        // Get current values for unset params
        if (eps < 0 || leaf < 0 || eta < 0) {
            double cur_eps = 1e-4; int cur_leaf = 10; double cur_eta = 2.0;
            // Read current from stats if available
            double dOut[20] = {0}; int nOut = 0;
            RadGetHACApKStats(dOut, &nOut);
            // Defaults used if not previously set
            if (eps < 0) eps = cur_eps;
            if (leaf < 0) leaf = cur_leaf;
            if (eta < 0) eta = cur_eta;
        }
        SetHACApKParams(eps, leaf, eta);
    }

    if (kwargs.contains("hmatrix_eps")) {
        SetHMatrixEpsilon(kwargs["hmatrix_eps"].cast<double>());
    }

    if (kwargs.contains("bicgstab_tol")) {
        SetBiCGSTABTol(kwargs["bicgstab_tol"].cast<double>());
    }

    if (kwargs.contains("relax_param")) {
        SetRelaxParam(kwargs["relax_param"].cast<double>());
    }

    if (kwargs.contains("newton_method")) {
        SetNewtonMethod(kwargs["newton_method"].cast<bool>());
    }

    if (kwargs.contains("newton_damping") || kwargs.contains("newton_damping_max_iter") || kwargs.contains("newton_damping_min_omega")) {
        bool enabled = kwargs.contains("newton_damping") ? kwargs["newton_damping"].cast<bool>() : true;
        int max_iter = kwargs.contains("newton_damping_max_iter") ? kwargs["newton_damping_max_iter"].cast<int>() : 5;
        double min_omega = kwargs.contains("newton_damping_min_omega") ? kwargs["newton_damping_min_omega"].cast<double>() : 0.01;
        SetNewtonDamping(enabled, max_iter, min_omega);
    }

    if (kwargs.contains("b_input_newton")) {
        SetBInputNewton(kwargs["b_input_newton"].cast<bool>());
    }

    if (kwargs.contains("b_input_hantila")) {
        SetBInputHantila(kwargs["b_input_hantila"].cast<bool>());
    }

    if (kwargs.contains("hantila_alpha")) {
        SetHantilaAlpha(kwargs["hantila_alpha"].cast<double>());
    }

    if (kwargs.contains("hantila_relax")) {
        SetHantilaRelax(kwargs["hantila_relax"].cast<double>());
    }

    if (kwargs.contains("keep_magnetization")) {
        SetKeepMagnetization(kwargs["keep_magnetization"].cast<bool>());
    }

    // Improved yano-type MSC kernel (opt-in research flag): pyramid-cloud + pyramid-centroid eval.
    // Default false == the historical EIEM2 single-point kernel (goldens unchanged).  Stage 1: supported
    // for the dense LU/BiCGSTAB yano-MSC path (no HACApK method 2, no IMA image symmetry).
    if (kwargs.contains("yano_pyramid_cloud")) {
        g_yano_pyramid_cloud = kwargs["yano_pyramid_cloud"].cast<bool>();
    }
    // Research flag: drop the element-center cancellation charge (raw collocation).  For studying the
    // "div(B)=0 via Lagrange instead of a center charge" question.  Pure-hex dense/BiCGSTAB matrix path.
    if (kwargs.contains("yano_no_center_charge")) {
        g_yano_no_center_charge = kwargs["yano_no_center_charge"].cast<bool>();
    }
    // Research: override the EIEM2 collocation-point alpha (eval = a*FaceCenter + (1-a)*center; -1 = 0.5).
    if (kwargs.contains("yano_eval_alpha")) {
        g_yano_eval_alpha = kwargs["yano_eval_alpha"].cast<double>();
    }
    if (kwargs.contains("yano_moment")) {
        g_yano_moment = kwargs["yano_moment"].cast<bool>();
    }
}

py::dict GetSolverConfig() {
    py::dict config;

    // BiCGSTAB tolerance
    { double tol = 1e-4;
      RadGetBiCGSTABTol(&tol);
      config["bicgstab_tol"] = tol; }

    // Relaxation parameter
    { double relax = 0.0;
      RadGetRelaxParam(&relax);
      config["relax_param"] = relax; }

    // Keep magnetization
    config["keep_magnetization"] = GetKeepMagnetization();

    // Newton method
    { int use_newton = 0;
      RadGetNewtonMethod(&use_newton);
      config["newton_method"] = (use_newton != 0); }

    // Newton damping
    { int enabled = 0; int max_iter = 5; double min_omega = 0.01;
      RadGetNewtonDampingStats(&enabled, &max_iter, &min_omega);
      config["newton_damping"] = (enabled != 0);
      config["newton_damping_max_iter"] = max_iter;
      config["newton_damping_min_omega"] = min_omega; }

    // B-input Newton for hysteresis
    { int b_input = 0;
      RadGetBInputNewton(&b_input);
      config["b_input_newton"] = (b_input != 0); }

    // B-input Hantila for hysteresis
    { int b_hantila = 0;
      RadGetBInputHantila(&b_hantila);
      config["b_input_hantila"] = (b_hantila != 0); }

    { double alpha = 0.0;
      RadGetHantilaAlpha(&alpha);
      config["hantila_alpha"] = alpha; }

    { double relax = 0.0;
      RadGetHantilaRelax(&relax);
      config["hantila_relax"] = relax; }

    // HACApK stats (if available)
    { double dOut[20] = {0}; int nOut = 0;
      RadGetHACApKStats(dOut, &nOut);
      if (nOut > 0) {
          py::dict stats;
          stats["n_lowrank"] = static_cast<int>(dOut[0]);
          stats["n_dense"] = static_cast<int>(dOut[1]);
          stats["max_rank"] = static_cast<int>(dOut[2]);
          stats["n_leaves"] = static_cast<int>(dOut[3]);
          stats["n_dof"] = static_cast<int>(dOut[4]);
          stats["compression"] = dOut[5];
          stats["build_time"] = dOut[6];
          stats["hmatrix_build_time"] = dOut[7];
          stats["linear_iterations"] = static_cast<int>(dOut[9]);
          stats["memory_mb"] = dOut[10];
          stats["dense_memory_mb"] = dOut[11];
          config["hacapk_stats"] = stats;
      } }

    config["yano_pyramid_cloud"] = g_yano_pyramid_cloud;
    config["yano_no_center_charge"] = g_yano_no_center_charge;
    config["yano_eval_alpha"] = g_yano_eval_alpha;
    config["yano_moment"] = g_yano_moment;

    return config;
}

} // namespace radia_solver_ext


// ============================================================================
// Additional Field Functions
// ============================================================================

namespace radia_field_ext {

// FldEnr (energy-based) REMOVED (Phase C, 2026-04-16)

py::array_t<double> FldFrc(int obj, int shape) {
    double f[6] = {0};
    int nf = 6;

    int err = RadFldFrc(f, &nf, obj, shape);
    check_error(err);

    py::array_t<double> result(nf);
    auto r = result.mutable_unchecked<1>();
    for (int i = 0; i < nf; i++) {
        r(i) = f[i];
    }
    return result;
}

int FldFrcShpRtg(py::array_t<double> center, py::array_t<double> dimensions) {
    auto c = center.unchecked<1>();
    auto d = dimensions.unchecked<1>();
    double P[3] = {c(0), c(1), c(2)};
    double W[2] = {d(0), d(1)};

    int handle = 0;
    int err = RadFldFrcShpRtg(&handle, P, W);
    check_error(err);
    return handle;
}


void FldCmpCrt(double prcB, double prcA, double prcBInt, double prcFrc, double prcTrjCrd, double prcTrjAng) {
    int n = 0;
    int err = RadFldCmpCrt(&n, prcB, prcA, prcBInt, prcFrc, prcTrjCrd, prcTrjAng);
    check_error(err);
}

void FldLenRndSw(const std::string& on_off) {
    int n = 0;
    char opt[16];
    strcpy(opt, on_off.c_str());
    int err = RadFldLenRndSw(&n, opt);
    check_error(err);
}

} // namespace radia_field_ext


// Replaced by Python-based PEEC topology solver (peec_topology.py, peec_coupled.py)


// ============================================================================
// Additional Utility Functions
// ============================================================================

namespace radia_utility_ext {

double UtiVer() {
    double ver = 0;
    int err = RadUtiVer(&ver);
    check_error(err);
    return ver;
}

} // namespace radia_utility_ext


// ============================================================================
// NGSolve CoefficientFunction: RadiaFieldCF
// ============================================================================

namespace ngfem
{

class RadiaFieldCF : public CoefficientFunction
{
public:
	int radia_obj;
	std::string field_type;

	// Coordinate transformation
	double origin[3];
	double u_axis[3];
	double v_axis[3];
	double w_axis[3];
	bool use_transform;

	// Computation settings
	std::optional<double> precision;

	// Point cache for batch evaluation
	mutable std::unordered_map<uint64_t, std::array<double,3>> point_cache_;
	mutable bool use_cache_;
	double cache_tolerance_;
	mutable size_t cache_hits_;
	mutable size_t cache_misses_;

	// Cached Radia module
	mutable py::module_ rad_module_;

	RadiaFieldCF(int obj, const std::string& ftype = "b",
	             std::optional<std::vector<double>> opt_origin = std::nullopt,
	             std::optional<std::vector<double>> opt_u = std::nullopt,
	             std::optional<std::vector<double>> opt_v = std::nullopt,
	             std::optional<std::vector<double>> opt_w = std::nullopt,
	             std::optional<double> opt_precision = std::nullopt,
	             const std::string& units = "m")
	    : CoefficientFunction(ftype == "phi" ? 1 : 3),
	      radia_obj(obj), field_type(ftype), use_transform(false),
	      precision(opt_precision),
	      use_cache_(false), cache_tolerance_(1e-10), cache_hits_(0), cache_misses_(0)
	{
		if (field_type != "b" && field_type != "h" &&
		    field_type != "a" && field_type != "m" && field_type != "phi") {
			throw std::invalid_argument(
				"Invalid field_type. Must be 'b', 'h', 'a', 'm', or 'phi'");
		}
		if (units != "m") {
			throw std::invalid_argument(
				"RadiaField requires units='m'. Radia always uses meters.");
		}

		origin[0] = 0; origin[1] = 0; origin[2] = 0;
		u_axis[0] = 1; u_axis[1] = 0; u_axis[2] = 0;
		v_axis[0] = 0; v_axis[1] = 1; v_axis[2] = 0;
		w_axis[0] = 0; w_axis[1] = 0; w_axis[2] = 1;

		auto apply_vec = [this](const std::optional<std::vector<double>>& opt,
		                        double dst[3], bool do_normalize) {
			if (!opt.has_value()) return;
			const auto& v = opt.value();
			if (v.size() != 3)
				throw std::invalid_argument("Vector must have 3 components");
			dst[0] = v[0]; dst[1] = v[1]; dst[2] = v[2];
			if (do_normalize) normalize(dst);
			use_transform = true;
		};

		apply_vec(opt_origin, origin, false);
		apply_vec(opt_u, u_axis, true);
		apply_vec(opt_v, v_axis, true);
		apply_vec(opt_w, w_axis, true);

		py::gil_scoped_acquire acquire;
		rad_module_ = py::module_::import("radia");

		if (precision.has_value()) {
			double prec = precision.value();
			std::string prec_str = "PrcB->" + std::to_string(prec) +
			                       ",PrcA->" + std::to_string(prec) +
			                       ",PrcH->" + std::to_string(prec) +
			                       ",PrcM->" + std::to_string(prec);
			rad_module_.attr("FldCmpPrc")(prec_str);
		}
	}

private:
	void normalize(double vec[3]) {
		double norm = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		if (norm < 1e-12)
			throw std::invalid_argument("Cannot normalize zero vector");
		vec[0] /= norm; vec[1] /= norm; vec[2] /= norm;
	}

	double dot(const double a[3], const double b[3]) const {
		return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	}

	uint64_t hash_point(double x, double y, double z) const {
		int64_t ix = static_cast<int64_t>(x / cache_tolerance_);
		int64_t iy = static_cast<int64_t>(y / cache_tolerance_);
		int64_t iz = static_cast<int64_t>(z / cache_tolerance_);
		uint64_t hash = 14695981039346656037ULL;
		hash ^= static_cast<uint64_t>(ix); hash *= 1099511628211ULL;
		hash ^= static_cast<uint64_t>(iy); hash *= 1099511628211ULL;
		hash ^= static_cast<uint64_t>(iz); hash *= 1099511628211ULL;
		return hash;
	}

	void transform_to_local(const double p_global[3], double p_local[3]) const {
		if (use_transform) {
			double p_t[3] = {p_global[0]-origin[0], p_global[1]-origin[1], p_global[2]-origin[2]};
			p_local[0] = dot(u_axis, p_t);
			p_local[1] = dot(v_axis, p_t);
			p_local[2] = dot(w_axis, p_t);
		} else {
			p_local[0] = p_global[0]; p_local[1] = p_global[1]; p_local[2] = p_global[2];
		}
	}

	void transform_to_global(const double f_local[3], double f_global[3]) const {
		if (use_transform) {
			f_global[0] = u_axis[0]*f_local[0] + v_axis[0]*f_local[1] + w_axis[0]*f_local[2];
			f_global[1] = u_axis[1]*f_local[0] + v_axis[1]*f_local[1] + w_axis[1]*f_local[2];
			f_global[2] = u_axis[2]*f_local[0] + v_axis[2]*f_local[1] + w_axis[2]*f_local[2];
		} else {
			f_global[0] = f_local[0]; f_global[1] = f_local[1]; f_global[2] = f_local[2];
		}
	}

public:
	void PrepareCache(py::list points_list) {
		py::gil_scoped_acquire acquire;
		size_t npts = points_list.size();
		point_cache_.clear();
		cache_hits_ = 0;
		cache_misses_ = 0;

		if (npts == 0) { use_cache_ = false; return; }

		py::array_t<double> pts_arr({(py::ssize_t)npts, (py::ssize_t)3});
		auto pts_buf = pts_arr.mutable_unchecked<2>();
		std::vector<double> globals(npts * 3);

		for (size_t i = 0; i < npts; i++) {
			py::list pt = points_list[i].cast<py::list>();
			double x = pt[0].cast<double>();
			double y = pt[1].cast<double>();
			double z = pt[2].cast<double>();
			globals[i*3] = x; globals[i*3+1] = y; globals[i*3+2] = z;
			double p_global[3] = {x, y, z};
			double p_local[3];
			transform_to_local(p_global, p_local);
			pts_buf(i, 0) = p_local[0];
			pts_buf(i, 1) = p_local[1];
			pts_buf(i, 2) = p_local[2];
		}

		// Use unified Fld(obj, field_type, points_array)
		py::object fld_result = rad_module_.attr("Fld")(radia_obj, field_type, pts_arr);

		if (field_type == "phi") {
			py::array_t<double> phi_arr = fld_result.cast<py::array_t<double>>();
			auto phi = phi_arr.unchecked<1>();
			for (size_t i = 0; i < npts; i++) {
				uint64_t hash = hash_point(globals[i*3], globals[i*3+1], globals[i*3+2]);
				point_cache_[hash] = {phi(i), 0.0, 0.0};
			}
		} else {
			py::array_t<double> fld_arr = fld_result.cast<py::array_t<double>>();
			auto fld = fld_arr.unchecked<2>();
			for (size_t i = 0; i < npts; i++) {
				double f_local[3] = {fld(i, 0), fld(i, 1), fld(i, 2)};
				double f_global[3];
				transform_to_global(f_local, f_global);
				uint64_t hash = hash_point(globals[i*3], globals[i*3+1], globals[i*3+2]);
				point_cache_[hash] = {f_global[0], f_global[1], f_global[2]};
			}
		}
		use_cache_ = true;
	}

	void ClearCache() {
		point_cache_.clear();
		use_cache_ = false;
		cache_hits_ = 0;
		cache_misses_ = 0;
	}

	py::dict GetCacheStats() const {
		py::dict stats;
		stats["enabled"] = use_cache_;
		stats["size"] = point_cache_.size();
		stats["hits"] = cache_hits_;
		stats["misses"] = cache_misses_;
		double total = cache_hits_ + cache_misses_;
		stats["hit_rate"] = (total > 0) ? (cache_hits_ / total) : 0.0;
		return stats;
	}

	// VoxelCoefficient generation for trajectory calculations
	py::object AsVoxelCF(py::object mesh, int resolution) const {
		py::gil_scoped_acquire acquire;

		py::module_ ngsolve = py::module_::import("ngsolve");
		py::object VoxelCoefficient = ngsolve.attr("VoxelCoefficient");
		py::object CF = ngsolve.attr("CF");
		py::module_ np = py::module_::import("numpy");

		// Get bounding box from mesh
		py::object ngmesh = mesh.attr("ngmesh");
		py::tuple bbox = ngmesh.attr("bounding_box").cast<py::tuple>();
		py::object pmin_obj = bbox[0];
		py::object pmax_obj = bbox[1];
		double pmin[3], pmax[3];
		for (int i = 0; i < 3; i++) {
			pmin[i] = pmin_obj[py::int_(i)].cast<double>();
			pmax[i] = pmax_obj[py::int_(i)].cast<double>();
		}
		double max_dim = 0;
		for (int i = 0; i < 3; i++)
			max_dim = std::max(max_dim, pmax[i] - pmin[i]);
		double margin = 0.01 * max_dim;
		for (int i = 0; i < 3; i++) {
			pmin[i] -= margin;
			pmax[i] += margin;
		}

		int nx = resolution, ny = resolution, nz = resolution;
		size_t total = (size_t)nx * ny * nz;

		// Generate grid points
		py::array_t<double> pts_arr({(py::ssize_t)total, (py::ssize_t)3});
		auto pts = pts_arr.mutable_unchecked<2>();
		size_t idx = 0;
		for (int ix = 0; ix < nx; ix++) {
			double x = pmin[0] + (pmax[0] - pmin[0]) * ix / (nx - 1);
			for (int iy = 0; iy < ny; iy++) {
				double y = pmin[1] + (pmax[1] - pmin[1]) * iy / (ny - 1);
				for (int iz = 0; iz < nz; iz++) {
					double z = pmin[2] + (pmax[2] - pmin[2]) * iz / (nz - 1);
					double p_global[3] = {x, y, z};
					double p_local[3];
					transform_to_local(p_global, p_local);
					pts(idx, 0) = p_local[0];
					pts(idx, 1) = p_local[1];
					pts(idx, 2) = p_local[2];
					idx++;
				}
			}
		}

		py::tuple start = py::make_tuple(pmin[0], pmin[1], pmin[2]);
		py::tuple end = py::make_tuple(pmax[0], pmax[1], pmax[2]);

		// Use unified Fld(obj, field_type, points_array)
		py::object fld_result = rad_module_.attr("Fld")(radia_obj, field_type, pts_arr);

		if (field_type == "phi") {
			py::array_t<double> phi_arr = fld_result.cast<py::array_t<double>>();
			py::object data = phi_arr.attr("reshape")(nx, ny, nz);
			data = np.attr("ascontiguousarray")(data.attr("transpose")(2, 1, 0));
			return VoxelCoefficient(start, end, data, "linear"_a = true);
		}

		// Vector fields
		py::array_t<double> field_arr = fld_result.cast<py::array_t<double>>();

		{
			// Vector: 3 components, each (nz, ny, nx)
			auto fld = field_arr.unchecked<2>();
			py::list cfs;
			for (int comp = 0; comp < 3; comp++) {
				py::array_t<double> comp_data({(py::ssize_t)total});
				auto cd = comp_data.mutable_unchecked<1>();
				for (size_t i = 0; i < total; i++) {
					double f_local[3] = {fld(i, 0), fld(i, 1), fld(i, 2)};
					double f_global[3];
					transform_to_global(f_local, f_global);
					cd(i) = f_global[comp];
				}
				py::object data = comp_data.attr("reshape")(nx, ny, nz);
				data = np.attr("ascontiguousarray")(data.attr("transpose")(2, 1, 0));
				cfs.append(VoxelCoefficient(start, end, data, "linear"_a = true));
			}
			return CF(py::tuple(cfs));
		}
	}

	virtual ~RadiaFieldCF() {}

	// Scalar evaluation for 'phi'
	virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
	{
		if (field_type != "phi") return 0.0;

		auto pnt = mip.GetPoint();
		int dim = pnt.Size();
		double p_global[3] = {pnt[0], (dim>=2)?pnt[1]:0.0, (dim>=3)?pnt[2]:0.0};
		double p_local[3];
		transform_to_local(p_global, p_local);

		double phi_value = 0.0;
		{
			py::gil_scoped_acquire acquire;
			try {
				py::list coords;
				coords.append(p_local[0]); coords.append(p_local[1]); coords.append(p_local[2]);
				// Use 'p' (not 'phi') to get scalar-only return (faster)
				py::object result = rad_module_.attr("Fld")(radia_obj, "p", coords);
				phi_value = result.cast<double>();
			} catch (std::exception&) {
				return 0.0;
			}
		}
		return phi_value;
	}

	// Single-point vector evaluation
	virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
	                      FlatVector<> result) const override
	{
		auto pnt = mip.GetPoint();
		int dim = pnt.Size();
		double p_global[3] = {pnt[0], (dim>=2)?pnt[1]:0.0, (dim>=3)?pnt[2]:0.0};

		if (field_type == "phi") {
			result(0) = Evaluate(mip);
			return;
		}

		// Check cache
		if (use_cache_) {
			uint64_t hash = hash_point(p_global[0], p_global[1], p_global[2]);
			auto it = point_cache_.find(hash);
			if (it != point_cache_.end()) {
				cache_hits_++;
				result(0) = it->second[0];
				result(1) = it->second[1];
				result(2) = it->second[2];
				return;
			}
			cache_misses_++;
		}

		double p_local[3];
		transform_to_local(p_global, p_local);

		double f_local[3];
		{
			py::gil_scoped_acquire acquire;
			try {
				py::list coords;
				coords.append(p_local[0]); coords.append(p_local[1]); coords.append(p_local[2]);
				py::object field_result = rad_module_.attr("Fld")(radia_obj, field_type, coords);
				f_local[0] = field_result[py::int_(0)].cast<double>();
				f_local[1] = field_result[py::int_(1)].cast<double>();
				f_local[2] = field_result[py::int_(2)].cast<double>();
			} catch (std::exception&) {
				result(0) = 0.0; result(1) = 0.0; result(2) = 0.0;
				return;
			}
		}

		double f_global[3];
		transform_to_global(f_local, f_global);
		result(0) = f_global[0]; result(1) = f_global[1]; result(2) = f_global[2];
	}

	// Batch evaluation (called by NGSolve for integration rules)
	virtual void Evaluate(const BaseMappedIntegrationRule& mir,
	                      BareSliceMatrix<> result) const override
	{
		size_t npts = mir.Size();

		py::gil_scoped_acquire acquire;

		try {
			// Try cache first
			if (use_cache_) {
				bool all_cached = true;
				for (size_t i = 0; i < npts; i++) {
					auto pnt = mir[i].GetPoint();
					int dim = pnt.Size();
					double p[3] = {pnt[0], (dim>=2)?pnt[1]:0.0, (dim>=3)?pnt[2]:0.0};
					uint64_t hash = hash_point(p[0], p[1], p[2]);
					auto it = point_cache_.find(hash);
					if (it != point_cache_.end()) {
						cache_hits_++;
						result(i,0) = it->second[0];
						result(i,1) = it->second[1];
						result(i,2) = it->second[2];
					} else {
						cache_misses_++;
						all_cached = false;
						break;
					}
				}
				if (all_cached) return;
			}

			// Build numpy array of local coordinates
			py::array_t<double> pts_arr({(py::ssize_t)npts, (py::ssize_t)3});
			auto pts_buf = pts_arr.mutable_unchecked<2>();

			for (size_t i = 0; i < npts; i++) {
				auto pnt = mir[i].GetPoint();
				int dim = pnt.Size();
				double p_global[3] = {pnt[0], (dim>=2)?pnt[1]:0.0, (dim>=3)?pnt[2]:0.0};
				double p_local[3];
				transform_to_local(p_global, p_local);
				pts_buf(i, 0) = p_local[0];
				pts_buf(i, 1) = p_local[1];
				pts_buf(i, 2) = p_local[2];
			}

			// Batch compute via unified Fld(obj, field_type, points_array)
			py::object fld_result = rad_module_.attr("Fld")(radia_obj, field_type, pts_arr);

			if (field_type == "phi") {
				// Scalar result: shape (N,)
				py::array_t<double> phi_arr = fld_result.cast<py::array_t<double>>();
				auto phi = phi_arr.unchecked<1>();
				for (size_t i = 0; i < npts; i++) {
					result(i, 0) = phi(i);
				}
			} else {
				// Vector result: shape (N, 3)
				py::array_t<double> field_arr = fld_result.cast<py::array_t<double>>();
				auto fld = field_arr.unchecked<2>();
				for (size_t i = 0; i < npts; i++) {
					double f_local[3] = {fld(i, 0), fld(i, 1), fld(i, 2)};
					double f_global[3];
					transform_to_global(f_local, f_global);
					result(i, 0) = f_global[0];
					result(i, 1) = f_global[1];
					result(i, 2) = f_global[2];
				}
			}
		} catch (std::exception&) {
			size_t dim = (field_type == "phi") ? 1 : 3;
			for (size_t i = 0; i < npts; i++) {
				for (size_t j = 0; j < dim; j++)
					result(i, j) = 0.0;
			}
		}
	}
};

} // namespace ngfem


// ============================================================================
// HDiv-type VIM demag operator (N = B^T G B) -- thin binding for golden testing.
// Returns the assembled N and charge map B as flat row-major lists; the Python golden
// test (tests/feec/) reshapes and checks symmetry + loops-field-null (loops = ker B).
// ============================================================================
namespace radia_hdivvim {
py::dict HDivVimAssemble(int nx, int ny, int nz, int nsub, double distort) {
    rad_hdiv::Mesh m = rad_hdiv::BuildStructuredRT0(nx, ny, nz, 1.0, distort);
    std::vector<double> B, N, M_mass;
    int n_charge = 0, n_bnd = 0;
    rad_hdiv::AssembleChargeMap(m, B, n_charge, n_bnd);
    rad_hdiv::AssembleN(m, N, nsub);   // nsub=0 centroid-monopole (fast); >=1 accurate sub-point
    rad_hdiv::AssembleMass(m, M_mass);
    py::dict d;
    d["nf"]       = m.n_face();
    d["n_cell"]   = m.n_cell;
    d["n_charge"] = n_charge;
    d["n_bnd"]    = n_bnd;
    d["N"]        = N;        // row-major (nf x nf)
    d["B"]        = B;        // row-major (n_charge x nf)
    d["M_mass"]   = M_mass;   // row-major (nf x nf), RT0 HDiv mass
    return d;
}


// Build the HDiv-type VIM demag operator N as a HACApK H-matrix and PROBE it against the dense
// reference (rad_hdiv::AssembleN) -- the production verification that the symmetric operator now
// scales as an H-matrix.  Self-contained (the C++ holds both the H-matrix and the dense N, and
// returns the comparison): for a few deterministic probe vectors x, compares H-matvec(x) vs
// N_dense*x (matvec_relerr) and checks H-matvec symmetry (x^T N y == y^T N x).  Returns the
// H-matrix stats (compression / n_lowrank / build_time) too.  Golden sizes only (dense N formed).
py::dict HDivVimHMatrixProbe(int nx, int ny, int nz, int nsub, double distort,
                             double eps, int leaf, double eta) {
    rad_hdiv::Mesh m = rad_hdiv::BuildStructuredRT0(nx, ny, nz, 1.0, distort);
    const int nf = m.n_face();
    std::vector<double> Nd;
    rad_hdiv::AssembleN(m, Nd, nsub);                 // dense reference (same nsub)

    RadHACApKHDivManager mgr(nx, ny, nz, 1.0, distort, nsub);
    RadHACApKParams prm;
    prm.aca_eps = eps; prm.leaf_size = leaf; prm.eta = eta; prm.print_level = 0;
    bool ok = mgr.BuildHMatrix(prm);

    py::dict d;
    d["ok"] = ok;
    d["nf"] = nf;
    if (!ok) return d;

    const RadHACApKStats& st = mgr.GetStats();
    d["n_dof"]           = st.n_dof;
    d["n_leaves"]        = st.n_leaves;
    d["n_lowrank"]       = st.n_lowrank;
    d["n_dense"]         = st.n_dense;
    d["max_rank"]        = st.max_rank;
    d["compression"]     = st.compression;
    d["build_time"]      = st.build_time;
    d["memory_mb"]       = st.memory_mb;
    d["dense_memory_mb"] = st.dense_memory_mb;

    auto densemv = [&](const std::vector<double>& x, std::vector<double>& y) {
        y.assign(nf, 0.0);
        for (int i = 0; i < nf; ++i) {
            const double* Nr = &Nd[(size_t)i * nf];
            double s = 0.0;
            for (int j = 0; j < nf; ++j) s += Nr[j] * x[j];
            y[i] = s;
        }
    };

    // matvec relerr vs dense N over a few deterministic probe vectors
    double max_rel = 0.0;
    for (int k = 0; k < 4; ++k) {
        std::vector<double> x(nf);
        for (int f = 0; f < nf; ++f)
            x[f] = std::sin(0.7 * (f + 1) + 1.3 * k) + 0.3 * std::cos(0.21 * (f + 1) * (k + 1));
        std::vector<double> yd, yh(nf, 0.0);
        densemv(x, yd);
        mgr.MatVec(x, yh);
        double num = 0.0, den = 0.0;
        for (int f = 0; f < nf; ++f) {
            num = std::max(num, std::fabs(yh[f] - yd[f]));
            den = std::max(den, std::fabs(yd[f]));
        }
        if (den > 0.0) max_rel = std::max(max_rel, num / den);
    }
    d["matvec_relerr"] = max_rel;

    // symmetry of the H-matvec: |x^T (N y) - y^T (N x)| / |x^T (N y)|
    std::vector<double> xv(nf), yv(nf), Nx(nf, 0.0), Ny(nf, 0.0);
    for (int f = 0; f < nf; ++f) { xv[f] = std::sin(0.3 * (f + 1)); yv[f] = std::cos(0.17 * (f + 1)); }
    mgr.MatVec(xv, Nx);
    mgr.MatVec(yv, Ny);
    double xtNy = 0.0, ytNx = 0.0;
    for (int f = 0; f < nf; ++f) { xtNy += xv[f] * Ny[f]; ytNx += yv[f] * Nx[f]; }
    d["symmetry_relerr"] = std::fabs(xtNy - ytNx) / (std::fabs(xtNy) + 1e-300);
    return d;
}

// SYSTEM-A H-LU on the HDiv-VIM operator: build A = M_mass + chi*N as a HACApK H-matrix on the FACE
// DOFs (SetSystemMode), then factor via cHACApK_hlu_decomp and round-trip A^-1(A x_orig)=x_orig.
// This is the HDiv analogue of HLUTestOnHACApK (which runs on the COLLOCATION matrix): it proves the
// materialize-free H-LU is a scalable DIRECT solve / strong preconditioner for the soft-iron material
// system, ON the actual HDiv-VIM operator.  Returns round-trip rel err + the materialize split
// (n_internal must be 0 = no cubic-driving densification) + decomp time + stats.
#ifdef RADIA_USE_HACAPK
py::dict HDivVimHLUProbe(int nx, int ny, int nz, int nsub, double distort,
                         double chi, double eps, int leaf, double eta, double trunc_tol) {
    py::dict d;
    RadHACApKHDivManager mgr(nx, ny, nz, 1.0, distort, nsub);
    mgr.SetSystemMode(chi);                       // store A = M_mass + chi*N (BEFORE build)
    RadHACApKParams prm;
    prm.aca_eps = eps; prm.leaf_size = leaf; prm.eta = eta; prm.print_level = 0;
    if (!mgr.BuildHMatrix(prm)) { d["ok"] = false; return d; }
    const int nf = mgr.GetNDOF();
    d["ok"] = true; d["ndof"] = nf;

    // deterministic x_orig; y = A x_orig via the H-matvec (the H-matrix IS A in system mode)
    std::vector<double> x_orig((size_t)nf), y_orig((size_t)nf, 0.0);
    unsigned long long seed = 13579ULL;
    for (int i = 0; i < nf; ++i) {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        x_orig[i] = (double)((seed >> 33) & 0x7fffffff) / 2147483647.0 - 0.5;
    }
    mgr.MatVec(x_orig, y_orig);

    cHACApK_hlu_set_trunc_tol(trunc_tol);
    double err = cHACApK_hlu_run_on_hacapk(mgr.GetLeafmtxp(), mgr.GetLcontrol(),
                                           x_orig.data(), y_orig.data(), 1);  // nffc=1 (one DOF/face)
    d["roundtrip_relerr"] = err;
    double t_decomp = 0, t_solve = 0; long n_lu = 0, n_gemm = 0;
    cHACApK_hlu_get_timings(&t_decomp, &t_solve, &n_lu, &n_gemm);
    d["t_decomp_sec"] = t_decomp; d["t_solve_sec"] = t_solve;
    long n_internal = 0, n_leaf = 0;
    cHACApK_hlu_get_materialize_split(&n_internal, &n_leaf);
    d["n_internal"] = n_internal; d["n_leaf"] = n_leaf;
    const RadHACApKStats& st = mgr.GetStats();
    d["compression"] = st.compression; d["max_rank"] = st.max_rank; d["build_time"] = st.build_time;
    return d;
}

// Phase-2 UNSTRUCTURED (real tet mesh) system-A H-LU: build A = M_mass + chi*N on the RT0 face DOFs of
// an arbitrary mesh (geometry from radia.vim build_demag), H-LU, round-trip A^-1(A x)=x.  This is the
// production-path operator (the C-type tet mesh); proves the materialize-free H-LU is a scalable
// direct-solve / Mprec for the soft-iron material system on unstructured geometry.
py::dict HDivVimTetHLUProbe(std::vector<double> face_centroids, double chi,
                            std::vector<int> face_charge, std::vector<double> face_coef,
                            std::vector<int> mI, std::vector<int> mJ, std::vector<double> mV,
                            std::vector<double> cell_verts, std::vector<double> face_verts, int n_el,
                            double eps, int leaf, double eta, double trunc_tol, double gram_near_factor) {
    py::dict d;
    RadHACApKHDivSystemTet mgr(std::move(face_centroids), chi, std::move(face_charge), std::move(face_coef),
                               std::move(mI), std::move(mJ), std::move(mV),
                               std::move(cell_verts), std::move(face_verts), n_el, gram_near_factor);
    RadHACApKParams prm; prm.aca_eps = eps; prm.leaf_size = leaf; prm.eta = eta; prm.print_level = 0;
    if (!mgr.BuildHMatrix(prm)) { d["ok"] = false; return d; }
    const int nf = mgr.GetNDOF();
    d["ok"] = true; d["ndof"] = nf;
    std::vector<double> x_orig((size_t)nf), y_orig((size_t)nf, 0.0);
    unsigned long long seed = 13579ULL;
    for (int i = 0; i < nf; ++i) {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        x_orig[i] = (double)((seed >> 33) & 0x7fffffff) / 2147483647.0 - 0.5;
    }
    mgr.MatVec(x_orig, y_orig);
    cHACApK_hlu_set_trunc_tol(trunc_tol);
    double err = cHACApK_hlu_run_on_hacapk(mgr.GetLeafmtxp(), mgr.GetLcontrol(),
                                           x_orig.data(), y_orig.data(), 1);
    d["roundtrip_relerr"] = err;
    double t_decomp = 0, t_solve = 0; long n_lu = 0, n_gemm = 0;
    cHACApK_hlu_get_timings(&t_decomp, &t_solve, &n_lu, &n_gemm);
    d["t_decomp_sec"] = t_decomp; d["t_solve_sec"] = t_solve;
    long n_internal = 0, n_leaf = 0;
    cHACApK_hlu_get_materialize_split(&n_internal, &n_leaf);
    d["n_internal"] = n_internal; d["n_leaf"] = n_leaf;
    const RadHACApKStats& st = mgr.GetStats();
    d["compression"] = st.compression; d["max_rank"] = st.max_rank; d["build_time"] = st.build_time;
    return d;
}

// Persistent factored H-LU(A) for the unstructured HDiv-VIM material system: build A = M_mass + chi*N on
// the RT0 face DOFs, H-LU-factor it ONCE, then apply A^-1 repeatedly (per GMRES iteration).  This is the
// production Mprec for the soft-iron demag solve -- A^-1 A ~ I -> N-INDEPENDENT outer iteration count
// (vs the M_mass^-1 preconditioner whose count grows with N).  factor in the ctor, apply via solve(),
// free in the dtor.
struct PyHDivVimTetSolver {
    std::unique_ptr<RadHACApKHDivSystemTet> mgr;
    void* root = nullptr;
    int   nd   = 0;
    PyHDivVimTetSolver(std::vector<double> face_centroids, double chi,
                       std::vector<int> face_charge, std::vector<double> face_coef,
                       std::vector<int> mI, std::vector<int> mJ, std::vector<double> mV,
                       std::vector<double> cell_verts, std::vector<double> face_verts, int n_el,
                       double eps, int leaf, double eta, double trunc_tol, double gram_near_factor) {
        mgr.reset(new RadHACApKHDivSystemTet(std::move(face_centroids), chi,
                  std::move(face_charge), std::move(face_coef), std::move(mI), std::move(mJ), std::move(mV),
                  std::move(cell_verts), std::move(face_verts), n_el, gram_near_factor));
        RadHACApKParams prm; prm.aca_eps = eps; prm.leaf_size = leaf; prm.eta = eta; prm.print_level = 0;
        if (!mgr->BuildHMatrix(prm)) throw std::runtime_error("HDivVimTetSolver: H-matrix build failed");
        cHACApK_hlu_set_trunc_tol(trunc_tol);
        nd = mgr->GetNDOF();
        root = cHACApK_hlu_factor_leafmtxp(mgr->GetLeafmtxp(), mgr->GetLcontrol(), 1);
        if (!root) throw std::runtime_error("HDivVimTetSolver: H-LU factor failed (singular block?)");
    }
    ~PyHDivVimTetSolver() { if (root) cHACApK_hlu_free_factors(root); }
    std::vector<double> solve(const std::vector<double>& rhs) const {
        if ((int)rhs.size() != nd) throw std::runtime_error("HDivVimTetSolver.solve: rhs size mismatch");
        std::vector<double> z((size_t)nd, 0.0);
        int rc = cHACApK_hlu_apply(root, mgr->GetLcontrol(), rhs.data(), z.data(), nd);
        if (rc != 0) throw std::runtime_error("HDivVimTetSolver.solve: H-LU apply failed");
        return z;
    }
    int ndof() const { return nd; }
};
#endif

// M2 verification probes: the C++ analytic charge-Gram potentials vs the dense Python reference.
double TriPotentialProbe(const std::vector<double>& V, const std::vector<double>& r) {
    double Vv[3][3], rr[3];
    for (int i = 0; i < 3; ++i) { rr[i] = r[i]; for (int k = 0; k < 3; ++k) Vv[i][k] = V[3*i+k]; }
    return rad_hdiv::TriPotential(Vv, rr);
}
double PhiTetProbe(const std::vector<double>& V, const std::vector<double>& P) {
    double Vv[4][3], PP[3];
    for (int i = 0; i < 3; ++i) PP[i] = P[i];
    for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) Vv[i][k] = V[3*i+k];
    return rad_hdiv::PhiTet(Vv, PP);
}
std::vector<double> TriFieldProbe(const std::vector<double>& V, const std::vector<double>& r) {
    double Vv[3][3], rr[3], out[3];
    for (int i = 0; i < 3; ++i) { rr[i] = r[i]; for (int k = 0; k < 3; ++k) Vv[i][k] = V[3*i+k]; }
    rad_hdiv::TriField(Vv, rr, out);
    return {out[0], out[1], out[2]};
}
std::vector<double> TetFieldProbe(const std::vector<double>& V, const std::vector<double>& P) {
    double Vv[4][3], PP[3], out[3];
    for (int i = 0; i < 3; ++i) PP[i] = P[i];
    for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) Vv[i][k] = V[3*i+k];
    rad_hdiv::TetField(Vv, PP, out);
    return {out[0], out[1], out[2]};
}
// ---- degree-1/2 polynomial-charge field kernel probes (entry-by-entry validation vs Python) ----
std::vector<double> TriMoment1Probe(const std::vector<double>& V, const std::vector<double>& r) {
    double Vv[3][3], rr[3], out[3];
    for (int i = 0; i < 3; ++i) { rr[i] = r[i]; for (int k = 0; k < 3; ++k) Vv[i][k] = V[3*i+k]; }
    rad_hdiv::TriMoment1(Vv, rr, out);
    return {out[0], out[1], out[2]};
}
std::vector<double> TriMoment2Probe(const std::vector<double>& V, const std::vector<double>& r) {
    double Vv[3][3], rr[3], out[3][3];
    for (int i = 0; i < 3; ++i) { rr[i] = r[i]; for (int k = 0; k < 3; ++k) Vv[i][k] = V[3*i+k]; }
    rad_hdiv::TriMoment2(Vv, rr, out);
    return {out[0][0],out[0][1],out[0][2], out[1][0],out[1][1],out[1][2], out[2][0],out[2][1],out[2][2]};
}
std::vector<double> TetMoment1Probe(const std::vector<double>& V, const std::vector<double>& r) {
    double Vv[4][3], rr[3], out[3];
    for (int i = 0; i < 3; ++i) rr[i] = r[i];
    for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) Vv[i][k] = V[3*i+k];
    rad_hdiv::TetMoment1(Vv, rr, out);
    return {out[0], out[1], out[2]};
}
std::vector<double> TetVolFieldLinearProbe(const std::vector<double>& V, const std::vector<double>& r,
                                           double rho0, const std::vector<double>& g) {
    double Vv[4][3], rr[3], gg[3], out[3];
    for (int i = 0; i < 3; ++i) { rr[i] = r[i]; gg[i] = g[i]; }
    for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) Vv[i][k] = V[3*i+k];
    rad_hdiv::TetVolFieldLinear(Vv, rr, rho0, gg, out);
    return {out[0], out[1], out[2]};
}
std::vector<double> TetVolFieldQuadraticProbe(const std::vector<double>& V, const std::vector<double>& r,
                                              double rho0, const std::vector<double>& g,
                                              const std::vector<double>& Q) {
    double Vv[4][3], rr[3], gg[3], QQ[3][3], out[3];
    for (int i = 0; i < 3; ++i) { rr[i] = r[i]; gg[i] = g[i]; for (int k=0;k<3;++k) QQ[i][k]=Q[3*i+k]; }
    for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) Vv[i][k] = V[3*i+k];
    rad_hdiv::TetVolFieldQuadratic(Vv, rr, rho0, gg, QQ, out);
    return {out[0], out[1], out[2]};
}
std::vector<double> LinTriFieldProbe(const std::vector<double>& V, const std::vector<double>& r,
                                     double sigma0, const std::vector<double>& s) {
    double Vv[3][3], rr[3], ss[3], out[3];
    for (int i = 0; i < 3; ++i) { rr[i] = r[i]; ss[i] = s[i]; for (int k = 0; k < 3; ++k) Vv[i][k] = V[3*i+k]; }
    rad_hdiv::LinTriField(Vv, rr, sigma0, ss, out);
    return {out[0], out[1], out[2]};
}
std::vector<double> QuadTriFieldProbe(const std::vector<double>& V, const std::vector<double>& r,
                                      double sigma0, const std::vector<double>& s,
                                      const std::vector<double>& S) {
    double Vv[3][3], rr[3], ss[3], SS[3][3], out[3];
    for (int i = 0; i < 3; ++i) { rr[i] = r[i]; ss[i] = s[i]; for (int k=0;k<3;++k) SS[i][k]=S[3*i+k]; }
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) Vv[i][k] = V[3*i+k];
    rad_hdiv::QuadTriField(Vv, rr, sigma0, ss, SS, out);
    return {out[0], out[1], out[2]};
}

// Batched HDiv-VIM demag field: the analytic charge field of ALL volume (linear-rho tet) + surface
// (quadratic-sigma tri) sources summed at ALL obs points, ngcore::ParallelFor OVER OBS (each obs is an
// independent sum -> no race).  This is the order_p / assemble_demag_field hot loop moved out of Python:
// one C++ call instead of O(n_obs x n_src) per-pair pybind crossings, and it honours the CALLER's
// `with TaskManager()` (bare ParallelFor falls back to serial when there is no active TaskManager -- per
// the caller-wraps policy, this function does NOT create its own RegionTaskManager).
//   vol  : n_vol  blocks of 16 doubles  [v0x..v3z (12), rho0 (1), gx,gy,gz (3)]
//   surf : n_surf blocks of 22 doubles  [v0x..v2z (9),  sigma0 (1), sx,sy,sz (3), S00..S22 (9 row-major)]
//   obs  : n_obs * 3.  Returns n_obs * 3, NO 1/(4pi) (applied in Python, matching the per-pair kernels).
std::vector<double> HDivDemagFieldBatch(const std::vector<double>& vol, const std::vector<double>& surf,
                                        const std::vector<double>& obs) {
    const long n_vol  = (long)(vol.size()  / 16);
    const long n_surf = (long)(surf.size() / 22);
    const long n_obs  = (long)(obs.size()  / 3);
    std::vector<double> out((size_t)n_obs * 3, 0.0);
    ngcore::ParallelFor(ngcore::IntRange((size_t)n_obs), [&](size_t jj) {
        const long j = (long)jj;
        double r[3] = { obs[3*j], obs[3*j+1], obs[3*j+2] };
        double acc[3] = {0.0, 0.0, 0.0};
        for (long e = 0; e < n_vol; ++e) {
            const double* b = &vol[(size_t)e * 16];
            double V[4][3]; for (int i=0;i<4;++i) for(int k=0;k<3;++k) V[i][k] = b[3*i+k];
            double g[3] = { b[13], b[14], b[15] };
            double o[3]; rad_hdiv::TetVolFieldLinear(V, r, b[12], g, o);
            acc[0]+=o[0]; acc[1]+=o[1]; acc[2]+=o[2];
        }
        for (long e = 0; e < n_surf; ++e) {
            const double* b = &surf[(size_t)e * 22];
            double V[3][3]; for (int i=0;i<3;++i) for(int k=0;k<3;++k) V[i][k] = b[3*i+k];
            double s[3] = { b[10], b[11], b[12] };
            double S[3][3]; for (int i=0;i<3;++i) for(int k=0;k<3;++k) S[i][k] = b[13 + 3*i + k];
            double o[3]; rad_hdiv::QuadTriField(V, r, b[9], s, S, o);
            acc[0]+=o[0]; acc[1]+=o[1]; acc[2]+=o[2];
        }
        out[3*j+0] = acc[0]; out[3*j+1] = acc[1]; out[3*j+2] = acc[2];
    });
    return out;
}
} // namespace radia_hdivvim

// ============================================================================
// Module Definition
// ============================================================================

PYBIND11_MODULE(_radia_pybind, m) {
    m.doc() = R"pbdoc(
        Radia - 3D Magnetostatics Library (pybind11 bindings)

        This module provides Python bindings for Radia using pybind11.
        It follows NGSolve design patterns for clean, efficient bindings.

        Example:
            import radia as rad

            # Create rectangular magnet (all coordinates in meters)
            magnet = rad.ObjRecMag([0,0,0], [0.04, 0.04, 0.02], [0, 0, 954930])

            # Compute field
            B = rad.Fld(magnet, 'b', [0.05, 0, 0])
    )pbdoc";

    // Version info
    m.attr("__version__") = "1.4.0";

    // ========================================================================
    // HDiv-type VIM (symmetric demag operator) -- golden-test entry
    // ========================================================================
    m.def("_hdiv_vim_assemble", &radia_hdivvim::HDivVimAssemble,
          py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("nsub") = 0, py::arg("distort") = 0.0,
          R"pbdoc(
              Assemble the symmetric HDiv-type VIM demag operator N = B^T G B on a structured
              nx*ny*nz hex grid (RT0 faces).  nsub=0 (default) uses the fast centroid-monopole Gram;
              nsub>=1 uses accurate sub-point quadrature (positive semi-definite N, true demag
              factors -- O(n_charge^2 nsub^6), golden sizes only).  Returns a dict with nf, n_cell,
              n_charge, n_bnd, and the row-major flat lists N (nf x nf), B (n_charge x nf), M_mass
              (nf x nf).  For testing: symmetry, loops field-null (ker B), demag factors (N, M_mass).
          )pbdoc");

    m.def("_hdiv_tri_potential", &radia_hdivvim::TriPotentialProbe, py::arg("V"), py::arg("r"),
          "M2 verify: Wilton triangle 1/r potential INT_T 1/|r-r'| dA' (V = 9 flat doubles = 3 verts, "
          "r = 3).  Pure 1/r integral (no 1/4pi).  Should match radia.vim._core.tri_potential.");
    m.def("_hdiv_phi_tet", &radia_hdivvim::PhiTetProbe, py::arg("V"), py::arg("P"),
          "M2 verify: tet Newtonian potential INT_tet 1/|P-r'| dV' (V = 12 flat doubles = 4 verts, "
          "P = 3) via the divergence theorem.  Should match radia.vim._core.phi_tet.");
    m.def("_hdiv_tri_field", &radia_hdivvim::TriFieldProbe, py::arg("V"), py::arg("r"),
          "Wilton triangle FIELD INT_T (r-r')/|r-r'|^3 dA' (V = 9 flat doubles, r = 3) -> 3-vector, no "
          "1/4pi.  = -grad TriPotential.  Should match radia.vim.flat_triangle_charge_field.");
    m.def("_hdiv_tet_field", &radia_hdivvim::TetFieldProbe, py::arg("V"), py::arg("P"),
          "Tet volume-charge FIELD INT_tet (P-r')/|P-r'|^3 dV' (V = 12 flat doubles, P = 3) -> 3-vector, "
          "no 1/4pi.  = -grad PhiTet.  Should match radia.vim.tet_self_volume_field * 4pi.");
    m.def("_hdiv_tri_moment1", &radia_hdivvim::TriMoment1Probe, py::arg("V"), py::arg("r"),
          "Surface first moment INT_T r'/R dS' (V=9, r=3) -> 3-vector.  == triangle_potential_moment.");
    m.def("_hdiv_tri_moment2", &radia_hdivvim::TriMoment2Probe, py::arg("V"), py::arg("r"),
          "Surface second moment INT_T r'(x)r'/R dS' (V=9, r=3) -> 9 (row-major 3x3).  == triangle_potential_moment2.");
    m.def("_hdiv_tet_moment1", &radia_hdivvim::TetMoment1Probe, py::arg("V"), py::arg("r"),
          "Volume first moment INT_V r'/R dV' (V=12, r=3) -> 3-vector.  == tet_newtonian_moment.");
    m.def("_hdiv_tet_volfield_linear", &radia_hdivvim::TetVolFieldLinearProbe,
          py::arg("V"), py::arg("r"), py::arg("rho0"), py::arg("g"),
          "Linear volume-charge field (V=12, r=3, rho0 scalar, g=3) -> 3-vector.  == tet_volume_field_linear.");
    m.def("_hdiv_tet_volfield_quadratic", &radia_hdivvim::TetVolFieldQuadraticProbe,
          py::arg("V"), py::arg("r"), py::arg("rho0"), py::arg("g"), py::arg("Q"),
          "Quadratic volume-charge field (V=12, r=3, rho0, g=3, Q=9 row-major) -> 3-vector.  == tet_volume_field_quadratic.");
    m.def("_hdiv_lin_tri_field", &radia_hdivvim::LinTriFieldProbe,
          py::arg("V"), py::arg("r"), py::arg("sigma0"), py::arg("s"),
          "Linear surface-charge field (V=9, r=3, sigma0, s=3) -> 3-vector.  == linear_triangle_charge_field.");
    m.def("_hdiv_quad_tri_field", &radia_hdivvim::QuadTriFieldProbe,
          py::arg("V"), py::arg("r"), py::arg("sigma0"), py::arg("s"), py::arg("S"),
          "Quadratic surface-charge field (V=9, r=3, sigma0, s=3, S=9 row-major) -> 3-vector.  == quadratic_triangle_charge_field.");
    m.def("_hdiv_demag_field_batch", &radia_hdivvim::HDivDemagFieldBatch,
          py::arg("vol"), py::arg("surf"), py::arg("obs"),
          "Batched HDiv-VIM demag field, ngcore::ParallelFor OVER OBS (honours the caller's TaskManager; "
          "serial fallback if none).  vol = n_vol x 16 [verts12, rho0, g3], surf = n_surf x 22 "
          "[verts9, sigma0, s3, S9 row-major], obs = n_obs x 3 -> n_obs x 3 (NO 1/4pi).  == the per-pair "
          "sum of _hdiv_tet_volfield_linear + _hdiv_quad_tri_field, the assemble_demag_field hot loop in C++.");


    m.def("_hdiv_vim_hmatrix_probe", &radia_hdivvim::HDivVimHMatrixProbe,
          py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("nsub") = 0, py::arg("distort") = 0.0,
          py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0,
          R"pbdoc(
              Build the HDiv-type VIM demag operator N = B^T G B as a HACApK H-matrix (DOFs = RT0
              faces; ACA entry function = on-demand charge-cluster Coulomb sum) and probe it against
              the dense rad_hdiv::AssembleN(nsub) reference.  Returns a dict: ok, nf, n_dof, n_leaves,
              n_lowrank, n_dense, max_rank, compression (H-bytes/dense-bytes), build_time, memory_mb,
              dense_memory_mb, matvec_relerr (max over deterministic probe vectors of ||H x - N x||/||N x||),
              and symmetry_relerr (|x^T N y - y^T N x| from the H-matvec).  Golden sizes only (forms dense N).
          )pbdoc");

#ifdef RADIA_USE_HACAPK
    m.def("_hdiv_vim_hlu_probe", &radia_hdivvim::HDivVimHLUProbe,
          py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("nsub") = 0, py::arg("distort") = 0.0,
          py::arg("chi") = 999.0, py::arg("eps") = 1e-5, py::arg("leaf") = 32, py::arg("eta") = 2.0,
          py::arg("trunc_tol") = 1e-4,
          R"pbdoc(
              SYSTEM-A H-LU on the HDiv-VIM operator: build A = M_mass + chi*N as a HACApK H-matrix on
              the RT0 face DOFs, factor via the materialize-free H-LU, and round-trip A^-1(A x)=x.
              Returns: ok, ndof, roundtrip_relerr, t_decomp_sec, t_solve_sec, n_internal (cubic-driver
              materialize, must be 0), n_leaf, compression, max_rank, build_time.  Proves the H-LU is a
              scalable direct-solve / strong preconditioner for the soft-iron material system on the
              actual HDiv operator (not just collocation).
          )pbdoc");

    m.def("_hdiv_vim_tet_hlu_probe", &radia_hdivvim::HDivVimTetHLUProbe,
          py::arg("face_centroids"), py::arg("chi"), py::arg("face_charge"), py::arg("face_coef"),
          py::arg("mI"), py::arg("mJ"), py::arg("mV"),
          py::arg("cell_verts"), py::arg("face_verts"), py::arg("n_el"),
          py::arg("eps") = 1e-5, py::arg("leaf") = 32, py::arg("eta") = 2.0,
          py::arg("trunc_tol") = 1e-8, py::arg("gram_near_factor") = 2.0,
          R"pbdoc(
              UNSTRUCTURED (real tet mesh) system-A H-LU: build A = M_mass + chi*N on the RT0 face DOFs
              from radia.vim build_demag geometry (face centroids; per-face <=2 charge map
              face_charge/face_coef; mass COO mI/mJ/mV; cell_verts/face_verts/n_el for the embedded
              analytic charge Gram), H-LU it, round-trip A^-1(A x)=x.  Returns ok, ndof, roundtrip_relerr,
              t_decomp_sec, n_internal (must be 0), n_leaf, compression, max_rank.  The production-path
              operator (the C-type); proves the materialize-free H-LU scales on unstructured geometry.
          )pbdoc");
#endif

    // Build-once HACApK H-matrix operator for the HDiv-type VIM demag operator -- the production
    // API: build N as an H-matrix ONCE, then drive a scalable symmetric solve (MINRES) over the
    // O(N log N) matvec.  matvec(x)=N x; apply_system(x,inv_chi)=inv_chi M_mass x - N x (the
    // material system operator); diag_system(inv_chi) is the Jacobi preconditioner diagonal.
    py::class_<RadHACApKHDivManager>(m, "_HDivVimHMatrix")
        .def(py::init([](int nx, int ny, int nz, int nsub, double distort,
                         double eps, int leaf, double eta) {
                 auto mgr = std::unique_ptr<RadHACApKHDivManager>(
                     new RadHACApKHDivManager(nx, ny, nz, 1.0, distort, nsub));
                 RadHACApKParams p;
                 p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                 if (!mgr->BuildHMatrix(p)) throw std::runtime_error("HDiv VIM H-matrix build failed");
                 return mgr;
             }),
             py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("nsub") = 0, py::arg("distort") = 0.0,
             py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0,
             "Build the HDiv-type VIM demag operator N as a HACApK H-matrix on a structured "
             "nx*ny*nz RT0 hex grid (nsub Gram quadrature, optional distortion; ACA eps/leaf/eta).")
        .def("ndof", [](RadHACApKHDivManager& s) { return s.GetNDOF(); })
        .def("matvec", [](RadHACApKHDivManager& s, const std::vector<double>& x) {
                 std::vector<double> y((size_t)s.GetNDOF(), 0.0);
                 s.MatVec(x, y);
                 return y;
             }, py::arg("x"), "N x (the O(N log N) H-matvec).")
        .def("apply_system", [](RadHACApKHDivManager& s, const std::vector<double>& x, double inv_chi) {
                 std::vector<double> y;
                 s.ApplySystem(x, inv_chi, y);
                 return y;
             }, py::arg("x"), py::arg("inv_chi"),
             "A x = inv_chi * (M_mass x) - (N x), the symmetric-indefinite material system operator.")
        .def("diag_system", [](RadHACApKHDivManager& s, double inv_chi) { return s.DiagSystem(inv_chi); },
             py::arg("inv_chi"), "diag(A) = inv_chi M_mass_ff - N_ff (Jacobi preconditioner).")
        .def("stats", [](RadHACApKHDivManager& s) {
                 const RadHACApKStats& st = s.GetStats();
                 py::dict d;
                 d["n_dof"] = st.n_dof; d["n_leaves"] = st.n_leaves; d["n_lowrank"] = st.n_lowrank;
                 d["n_dense"] = st.n_dense; d["max_rank"] = st.max_rank;
                 d["compression"] = st.compression; d["build_time"] = st.build_time;
                 d["memory_mb"] = st.memory_mb; d["dense_memory_mb"] = st.dense_memory_mb;
                 return d;
             }, "H-matrix stats dict (n_dof, n_leaves, n_lowrank, compression, build_time, ...).");

#ifdef RADIA_USE_HACAPK
    // Persistent factored H-LU(A) preconditioner for the unstructured HDiv-VIM material system.
    // Factor A = M_mass + chi*N ONCE (ctor), apply A^-1 per GMRES iteration (solve) -> N-independent
    // outer iters (the production Mprec replacing M_mass^-1).
    py::class_<radia_hdivvim::PyHDivVimTetSolver>(m, "_HDivVimTetSolver")
        .def(py::init([](std::vector<double> face_centroids, double chi,
                         std::vector<int> face_charge, std::vector<double> face_coef,
                         std::vector<int> mI, std::vector<int> mJ, std::vector<double> mV,
                         std::vector<double> cell_verts, std::vector<double> face_verts, int n_el,
                         double eps, int leaf, double eta, double trunc_tol, double gram_near_factor) {
                 return std::unique_ptr<radia_hdivvim::PyHDivVimTetSolver>(
                     new radia_hdivvim::PyHDivVimTetSolver(std::move(face_centroids), chi,
                         std::move(face_charge), std::move(face_coef), std::move(mI), std::move(mJ),
                         std::move(mV), std::move(cell_verts), std::move(face_verts), n_el,
                         eps, leaf, eta, trunc_tol, gram_near_factor));
             }),
             py::arg("face_centroids"), py::arg("chi"), py::arg("face_charge"), py::arg("face_coef"),
             py::arg("mI"), py::arg("mJ"), py::arg("mV"),
             py::arg("cell_verts"), py::arg("face_verts"), py::arg("n_el"),
             py::arg("eps") = 1e-5, py::arg("leaf") = 32, py::arg("eta") = 2.0,
             py::arg("trunc_tol") = 1e-6, py::arg("gram_near_factor") = 2.0,
             "Build + H-LU-factor A = M_mass + chi*N on the RT0 face DOFs (factor once).")
        .def("solve", &radia_hdivvim::PyHDivVimTetSolver::solve, py::arg("rhs"),
             "Apply A^-1 to rhs (per-GMRES-iteration preconditioner application).")
        .def("ndof", &radia_hdivvim::PyHDivVimTetSolver::ndof);
#endif

    // Charge-charge Coulomb Gram G as a HACApK H-matrix -- the UNSTRUCTURED / general-mesh path.
    // Charges (cell rho + boundary-face sigma) extracted from ANY RT0 mesh (e.g. NGSolve tet
    // HDiv(0)); pass charge centroids/measures + the caller-computed diagonal self-energies.  The
    // demag operator N = B^T G B is applied as B^T (matvec(B m)) with B the sparse charge map.
    py::class_<RadHACApKChargeGram>(m, "_ChargeGramHMatrix")
        .def(py::init([](std::vector<double> centroids, std::vector<double> measures,
                         std::vector<double> self_energy, double eps, int leaf, double eta) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(
                     new RadHACApKChargeGram(std::move(centroids), std::move(measures),
                                             std::move(self_energy)));
                 RadHACApKParams p;
                 p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                 if (!mgr->BuildHMatrix(p)) throw std::runtime_error("charge Gram H-matrix build failed");
                 return mgr;
             }),
             py::arg("centroids"), py::arg("measures"), py::arg("self_energy"),
             py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0,
             "Build the n_charge x n_charge Coulomb Gram G as a HACApK H-matrix over the charge "
             "centroids (G[a!=b] = meas_a meas_b/(4pi r), G[a][a] = self_energy[a]).")
        .def(py::init([](std::vector<double> cell_verts, std::vector<double> face_verts,
                         int n_el, double eps, int leaf, double eta, double near_factor) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(
                     new RadHACApKChargeGram(std::move(cell_verts), std::move(face_verts), n_el, near_factor));
                 RadHACApKParams p;
                 p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                 if (!mgr->BuildHMatrix(p)) throw std::runtime_error("analytic charge Gram H-matrix build failed");
                 return mgr;
             }),
             py::arg("cell_verts"), py::arg("face_verts"), py::arg("n_el"),
             py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0, py::arg("near_factor") = 1e30,
             "ANALYTIC mode (M2b): build the EXACT charge Gram as a HACApK H-matrix from per-charge "
             "geometry (cell_verts [n_el*12] tets, face_verts [n_bf*9] triangles). Entry = analytic "
             "PhiTet/TriPotential inner x outer quadrature (matches dense build_demag(analytic_gram=True)). "
             "near_factor (default 1e30 = all-analytic) gives the NEAR/FAR build speedup: pass ~2 to use "
             "the cheap centroid-monopole for far pairs, analytic only for near (non-uniform-M) pairs.")
        .def(py::init([](std::vector<double> cell_tris, std::vector<int> cell_troff,
                         std::vector<double> cell_cent, std::vector<double> cell_meas,
                         std::vector<double> face_tris, std::vector<int> face_troff,
                         std::vector<double> face_cent, std::vector<double> face_meas,
                         int n_el, double eps, int leaf, double eta, double near_factor) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(
                     new RadHACApKChargeGram(std::move(cell_tris), std::move(cell_troff),
                                             std::move(cell_cent), std::move(cell_meas),
                                             std::move(face_tris), std::move(face_troff),
                                             std::move(face_cent), std::move(face_meas), n_el, near_factor));
                 RadHACApKParams p;
                 p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                 if (!mgr->BuildHMatrix(p)) throw std::runtime_error("polytope charge Gram H-matrix build failed");
                 return mgr;
             }),
             py::arg("cell_tris"), py::arg("cell_troff"), py::arg("cell_cent"), py::arg("cell_meas"),
             py::arg("face_tris"), py::arg("face_troff"), py::arg("face_cent"), py::arg("face_meas"),
             py::arg("n_el"), py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0,
             py::arg("near_factor") = 1e30,
             "POLYTOPE mode (hex/wedge cells + quad faces): the EXACT analytic charge Gram for any "
             "flat-faced convex cell, the triangulation supplied from Python.  cell_tris is a flat "
             "triangle soup (9 doubles/tri) of all cells' convex-hull triangles, cell_troff [n_el+1] the "
             "CSR offsets (in triangles), cell_cent [n_el*3] the vertex-mean centroid (fan apex + outward "
             "normal ref), cell_meas [n_el] the cell volume; face_tris/face_troff/face_cent/face_meas the "
             "boundary faces' sub-triangles (quad->2).  Entry = divergence-theorem polytope potential x "
             "centroid-fan / Dunavant outer quadrature (matches dense analytic_charge_gram polytope path). "
             "near_factor (default 1e30 = all-analytic); pass ~2 for the NEAR/FAR build speedup.")
        .def(py::init([](std::vector<double> cell_verts, std::vector<double> face_verts, int n_el,
                         std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
                         std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
                         std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
                         std::vector<double> ref_tet_pts_lo, std::vector<double> ref_tet_w_lo,
                         std::vector<double> ref_tri_pts_lo, std::vector<double> ref_tri_w_lo,
                         double ho_far_factor,
                         std::vector<double> ref_tet_pts_in, std::vector<double> ref_tet_w_in,
                         std::vector<double> ref_tri_pts_in, std::vector<double> ref_tri_w_in,
                         double eps, int leaf, double eta) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(
                     new RadHACApKChargeGram(std::move(cell_verts), std::move(face_verts), n_el,
                                             std::move(charge_host), std::move(charge_kind),
                                             std::move(charge_expo), std::move(ref_tet_pts),
                                             std::move(ref_tet_w), std::move(ref_tri_pts), std::move(ref_tri_w),
                                             std::move(ref_tet_pts_lo), std::move(ref_tet_w_lo),
                                             std::move(ref_tri_pts_lo), std::move(ref_tri_w_lo), ho_far_factor,
                                             std::move(ref_tet_pts_in), std::move(ref_tet_w_in),
                                             std::move(ref_tri_pts_in), std::move(ref_tri_w_in)));
                 RadHACApKParams p;
                 p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                 if (!mgr->BuildHMatrix(p)) throw std::runtime_error("high-order charge Gram H-matrix build failed");
                 return mgr;
             }),
             py::arg("cell_verts"), py::arg("face_verts"), py::arg("n_el"),
             py::arg("charge_host"), py::arg("charge_kind"), py::arg("charge_expo"),
             py::arg("ref_tet_pts"), py::arg("ref_tet_w"), py::arg("ref_tri_pts"), py::arg("ref_tri_w"),
             py::arg("ref_tet_pts_lo") = std::vector<double>{}, py::arg("ref_tet_w_lo") = std::vector<double>{},
             py::arg("ref_tri_pts_lo") = std::vector<double>{}, py::arg("ref_tri_w_lo") = std::vector<double>{},
             py::arg("ho_far_factor") = 1e30,
             py::arg("ref_tet_pts_in") = std::vector<double>{}, py::arg("ref_tet_w_in") = std::vector<double>{},
             py::arg("ref_tri_pts_in") = std::vector<double>{}, py::arg("ref_tri_w_in") = std::vector<double>{},
             py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0,
             "HIGH-ORDER (order-p) mode: POLYNOMIAL charges (monomial basis per host). charge_host[c]/"
             "charge_kind[c] (0=cell,1=face)/charge_expo[3c] define each charge; ref_tet_pts[nqt*3]/ref_tet_w "
             "(sum 1/6) + ref_tri_pts[nqr*2]/ref_tri_w (sum 1/2) are the reference Gauss-Duffy rules. Entry = "
             "monomial-weighted outer quad x the subtraction inner potential (matches dense build_demag_highorder). "
             "ref_*_lo + ho_far_factor (<inf) enable the accuracy-preserving NEAR/FAR adaptive quadrature: far "
             "pairs (|c_a-c_b| > ho_far_factor*(size_a+size_b)) use the cheap LOW-quad plain double-Gauss.")
        .def("ndof", [](RadHACApKChargeGram& s) { return s.GetNDOF(); })
        .def("matvec", [](RadHACApKChargeGram& s, const std::vector<double>& x) {
                 std::vector<double> y((size_t)s.GetNDOF(), 0.0);
                 s.MatVec(x, y);
                 return y;
             }, py::arg("x"), "G q (the O(N log N) Gram H-matvec).")
        .def("solve_linear_material",
             [](RadHACApKChargeGram& s,
                std::vector<int> B_indptr, std::vector<int> B_indices, std::vector<double> B_data,
                int n_face, std::vector<int> mI, std::vector<int> mJ, std::vector<double> mV,
                double inv_chi, std::vector<double> prec, std::vector<double> rhs,
                double tol, int maxit) {
                 int iters = 0;
                 std::vector<double> m = s.SolveLinearMaterial(B_indptr, B_indices, B_data, n_face,
                                                               mI, mJ, mV, inv_chi, prec, rhs,
                                                               tol, maxit, iters);
                 py::dict d; d["m"] = m; d["iters"] = iters; return d;
             },
             py::arg("B_indptr"), py::arg("B_indices"), py::arg("B_data"), py::arg("n_face"),
             py::arg("mI"), py::arg("mJ"), py::arg("mV"), py::arg("inv_chi"), py::arg("prec"),
             py::arg("rhs"), py::arg("tol") = 1e-9, py::arg("maxit") = 5000,
             "M3: solve the SPD HDiv-VIM linear material system ((1/chi)M_mass + B^T G B) m = rhs by "
             "Jacobi-preconditioned CG (G applied as the charge-Gram H-matvec). Returns {m, iters}.")
        .def("solve_material_minres",
             [](RadHACApKChargeGram& s,
                std::vector<int> B_indptr, std::vector<int> B_indices, std::vector<double> B_data,
                int n_face, std::vector<int> mI, std::vector<int> mJ, std::vector<double> mV,
                double inv_chi, std::vector<double> prec, std::vector<double> rhs,
                double tol, int maxit) {
                 int iters = 0;
                 std::vector<double> m = s.SolveMaterialMINRES(B_indptr, B_indices, B_data, n_face,
                                                               mI, mJ, mV, inv_chi, prec, rhs,
                                                               tol, maxit, iters);
                 py::dict d; d["m"] = m; d["iters"] = iters; return d;
             },
             py::arg("B_indptr"), py::arg("B_indices"), py::arg("B_data"), py::arg("n_face"),
             py::arg("mI"), py::arg("mJ"), py::arg("mV"), py::arg("inv_chi"), py::arg("prec"),
             py::arg("rhs"), py::arg("tol") = 1e-8, py::arg("maxit") = 5000,
             "mu_r-INDEPENDENT material solve: Jacobi-preconditioned MINRES for the SYMMETRIC INDEFINITE "
             "A m = rhs, A = inv_chi*M_mass - B^T G B (G via the analytic charge-Gram H-matvec, "
             "HACApK-parallel under ngcore::RegionTaskManager). Returns {m, iters}.")
        .def("solve_nonlinear_picard",
             [](RadHACApKChargeGram& s,
                std::vector<int> B_indptr, std::vector<int> B_indices, std::vector<double> B_data,
                int n_face, std::vector<int> mI, std::vector<int> mJ, std::vector<double> mV,
                std::vector<double> Mmass_diag, std::vector<double> N_diag, std::vector<double> mu,
                double denom, double chi0, double Msat, double H0,
                int picard_iters, double cg_tol, int cg_maxit) {
                 auto r = s.SolveNonlinearPicard(B_indptr, B_indices, B_data, n_face, mI, mJ, mV,
                                                 Mmass_diag, N_diag, mu, denom, chi0, Msat, H0,
                                                 picard_iters, cg_tol, cg_maxit);
                 py::dict d;
                 d["m"] = r.m; d["Mavg"] = r.Mavg; d["chi"] = r.chi; d["Dscal"] = r.Dscal;
                 d["iters"] = r.iters;
                 return d;
             },
             py::arg("B_indptr"), py::arg("B_indices"), py::arg("B_data"), py::arg("n_face"),
             py::arg("mI"), py::arg("mJ"), py::arg("mV"), py::arg("Mmass_diag"), py::arg("N_diag"),
             py::arg("mu"), py::arg("denom"), py::arg("chi0"), py::arg("Msat"), py::arg("H0"),
             py::arg("picard_iters") = 100, py::arg("cg_tol") = 1e-10, py::arg("cg_maxit") = 5000,
             "M3 (nonlinear): scalar-chi Picard solve of the isotropic nonlinear demag M=Mof(H0-Dscal M) "
             "entirely in C++ (each step a SolveLinearMaterial + closed-form chi update; G via the analytic "
             "H-matvec). Returns {m, Mavg, chi, Dscal, iters}.")
        .def("stats", [](RadHACApKChargeGram& s) {
                 const RadHACApKStats& st = s.GetStats();
                 py::dict d;
                 d["n_dof"] = st.n_dof; d["n_leaves"] = st.n_leaves; d["n_lowrank"] = st.n_lowrank;
                 d["n_dense"] = st.n_dense; d["max_rank"] = st.max_rank;
                 d["compression"] = st.compression; d["build_time"] = st.build_time;
                 d["memory_mb"] = st.memory_mb; d["dense_memory_mb"] = st.dense_memory_mb;
                 return d;
             }, "H-matrix stats dict.");

    // ========================================================================
    // Object Creation
    // ========================================================================

    m.def("ObjRecMag", &radia_objects::ObjRecMag,
          py::arg("center"), py::arg("dimensions"), py::arg("magnetization"),
          R"pbdoc(
              Create rectangular permanent magnet.

              Args:
                  center: Center point [x, y, z]
                  dimensions: Block dimensions [Lx, Ly, Lz]
                  magnetization: Magnetization vector [Mx, My, Mz] in A/m

              Returns:
                  Object handle

              Example:
                  magnet = rad.ObjRecMag([0,0,0], [0.04, 0.04, 0.02], [0, 0, 954930])
          )pbdoc");

    m.def("ObjHexahedron", &radia_objects::ObjHexahedron,
          py::arg("vertices"), py::arg("magnetization"),
          R"pbdoc(
              Create hexahedral element from 8 vertices.

              Vertices should be ordered as: bottom face (vertices 0-3) counterclockwise,
              then top face (vertices 4-7) counterclockwise.

              Args:
                  vertices: List of 8 vertex coordinates [[x,y,z], ...]
                  magnetization: Magnetization vector [Mx, My, Mz] in A/m

              Returns:
                  Object handle
          )pbdoc");

    m.def("ObjTetrahedron", &radia_objects::ObjTetrahedron,
          py::arg("vertices"), py::arg("magnetization"),
          R"pbdoc(
              Create tetrahedral element from 4 vertices.

              Args:
                  vertices: List of 4 vertex coordinates [[x,y,z], ...]
                  magnetization: Magnetization vector [Mx, My, Mz] in A/m

              Returns:
                  Object handle
          )pbdoc");

    m.def("ObjWedge", &radia_objects::ObjWedge,
          py::arg("vertices"), py::arg("magnetization"),
          R"pbdoc(
              Create wedge/prism element from 6 vertices.

              Wedge element with triangular top and bottom faces.
              5 faces total: 2 triangular (top/bottom) + 3 quadrilateral (sides).
              5 DOF for MSC method.

              Vertex ordering (ELF MMB6T compatible):
                  Bottom triangle: vertices 0, 1, 2 (CCW when viewed from below)
                  Top triangle: vertices 3, 4, 5 (CCW when viewed from above)
                  v3-v5 are directly above v0-v2 respectively.

              Args:
                  vertices: List of 6 vertex coordinates [[x,y,z], ...]
                  magnetization: Magnetization vector [Mx, My, Mz] in A/m

              Returns:
                  Object handle
          )pbdoc");

    m.def("ObjCnt", &radia_objects::ObjCnt,
          py::arg("objects"),
          R"pbdoc(
              Create container for objects.

              Args:
                  objects: List of object handles

              Returns:
                  Container handle
          )pbdoc");

    m.def("ObjBckg", &radia_objects::ObjBckg,
          py::arg("callback"),
          R"pbdoc(
              Create background field source with callback.

              The callback receives [x,y,z] in current units and must return
              [Bx, By, Bz] in Tesla.

              Args:
                  callback: Python function(point) -> [Bx, By, Bz]

              Returns:
                  Object handle

              Example:
                  # Uniform field
                  ext = rad.ObjBckg(lambda p: [0, 0, 0.1])  # 0.1 T in z

                  # Quadrupole field
                  def quad(p):
                      return [10*p[1], 10*p[0], 0]
                  ext = rad.ObjBckg(quad)
          )pbdoc");

    // ========================================================================
    // Field Computation
    // ========================================================================

    m.def("Fld", &radia_field::Fld,
          py::arg("obj"), py::arg("field_type"), py::arg("points"),
          R"pbdoc(
              Compute field at single point or multiple points.

              Auto-detects single vs batch from input shape:
              - Shape (3,): single point evaluation
              - Shape (N, 3): batch evaluation (TaskManager-parallelized)

              Args:
                  obj: Object handle
                  field_type: "b", "h", "a", "phi", "m", "bx", "by", "bz", etc.
                  points: [x,y,z] for single point, or array of shape (N,3) for batch

              Returns:
                  Single point: scalar or array [Fx, Fy, Fz]
                  Batch: array of shape (N, 3) for vector fields, (N,) for scalar fields

              Examples:
                  B = rad.Fld(obj, "b", [0, 0, 0.1])        # Single point -> [Bx, By, Bz]
                  B = rad.Fld(obj, "b", points_Nx3)          # Batch -> (N, 3) array
                  phi = rad.Fld(obj, "phi", points_Nx3)      # Batch scalar -> (N,) array
          )pbdoc");


    // ========================================================================
    // Materials
    // ========================================================================

    m.def("MatLin", &radia_material::MatLin,
          py::arg("mu_r"), py::arg("easy_axis") = py::none(),
          R"pbdoc(
              Create linear magnetic material.

              Args:
                  mu_r: Relative permeability (scalar for isotropic,
                        [mu_par, mu_perp] for anisotropic)
                  easy_axis: [ex, ey, ez] for anisotropic materials

              Returns:
                  Material handle

              Example:
                  mat = rad.MatLin(1000)  # Isotropic mu_r=1000
                  mat = rad.MatLin([5000, 100], [0, 0, 1])  # Anisotropic
          )pbdoc");

    m.def("MatSatIsoTab", &radia_material::MatSatIsoTab,
          py::arg("bh_data"),
          R"pbdoc(
              Create nonlinear isotropic material from B-H curve.

              Args:
                  bh_data: List of [H, B] pairs where H is in A/m and B in Tesla

              Returns:
                  Material handle
          )pbdoc");

    m.def("MatApl", &radia_material::MatApl,
          py::arg("obj"), py::arg("mat"),
          R"pbdoc(
              Apply material to object.

              Args:
                  obj: Object handle
                  mat: Material handle

              Returns:
                  Result object handle
          )pbdoc");

    // NOTE: MatSIBC binding removed (2026-02-13). Use Python SIBC instead.

    // ========================================================================
    // Solver
    // ========================================================================

    m.def("Solve", &radia_solver::Solve,
          py::arg("obj"), py::arg("prec"), py::arg("max_iter"), py::arg("method") = 1,
          py::arg("image") = "",
          R"pbdoc(
              Solve magnetostatic problem.

              Args:
                  obj: Object or container handle
                  prec: Convergence precision (e.g., 0.001 = 0.1%)
                  max_iter: Maximum iterations
                  method: Solver method (0=LU, 1=BiCGSTAB, 2=HACApK)
                  image: Image symmetry string (e.g., "+x", "-z", "+x-z")
                         +x: X-axis symmetric mirror (default)
                         -z: Z-axis antisymmetric mirror
                         +x-z: X symmetric + Z antisymmetric (ELF quarter model)

              Returns:
                  Tuple (residual, ?, ?, iterations)

              Solver Selection:
                  - Small (<500 elements): method=0 (LU)
                  - Medium (500-5000): method=1 (BiCGSTAB)
                  - Large (>5000): method=2 (HACApK)
          )pbdoc");

    m.def("BuildMatrix", &radia_solver::BuildMatrix,
          py::arg("obj"), py::arg("image") = "",
          R"pbdoc(
              Build interaction matrix without solving.

              This allows inspection of the matrix before solving. The matrix is
              cached for subsequent Solve() calls with the same object and image.

              Args:
                  obj: Object or container handle
                  image: Image symmetry string (e.g., "+x", "-z", "+x-z")

              Returns:
                  Interaction matrix handle (for GetInteractMatrix)

              Example:
                  >>> handle = rad.BuildMatrix(model, image='+x-z')
                  >>> matrix, dof = rad.GetInteractMatrix(handle)
                  >>> rad.Solve(model, 0.0001, 100, 0)  # Uses cached matrix
          )pbdoc");

    m.def("GetSolveStats", &radia_solver::GetSolveStats,
          R"pbdoc(
              Get solve statistics from last Solve() call.

              Returns:
                  Dictionary with:
                    - t_matrix_build: Matrix construction time [s]
                    - t_linear_solve: Linear solver time [s]
                    - linear_iterations: BiCGSTAB iterations
                    - nonl_iterations: Nonlinear iterations
          )pbdoc");

    // PreRelax REMOVED (2026-01-31) - Use BuildMatrix() instead

    m.def("GetInteractMatrix", &radia_solver::GetInteractMatrix,
          py::arg("intrc_handle"),
          R"pbdoc(
              Get interaction matrix as numpy array.

              Args:
                  intrc_handle: Interaction handle from BuildMatrix()

              Returns:
                  Tuple (matrix, dof) where matrix is (dof x dof) numpy array
          )pbdoc");

    m.def("GetLoopBasis", &radia_solver::GetLoopBasis,
          py::arg("intrc_handle"),
          R"pbdoc(
              Get the yano-MSC cell-graph cycle (loop) basis L of an interaction matrix.

              The loops are the field-null subspace of the surface-charge collocation
              operator N (== HDiv ker(B), the cell/internal-face cycle space).  They sit
              at eigenvalue 1/chi in A = (1/chi) I - N, so cond(A) ~ mu_r; deflating L
              makes the high-mu_r solve bounded and mu_r-independent.  Built geometry-only
              (no SVD): face-center matching -> cell adjacency -> fundamental cycles.

              Args:
                  intrc_handle: Interaction handle from BuildMatrix()

              Returns:
                  Tuple (L, nLoop) where L is a (dof x nLoop) numpy array whose columns
                  span ker(N) (verify ||N @ L|| ~ 0).
          )pbdoc");

    m.def("GetFaceGeom", &radia_solver::GetFaceGeom,
          py::arg("intrc_handle"),
          R"pbdoc(
              Get per-DOF hex face geometry of an interaction matrix, in the matrix DOF order.

              Rows align 1:1 with GetInteractMatrix.  Each row (stride 11) is
              [elem_local, area, cx,cy,cz, nx,ny,nz(outward), ecx,ecy,ecz]: the face area,
              face centroid, outward unit normal, and owning-element center.  Non-hex DOFs
              (tet/wedge) get elem_local = -1 and zeros.  Lets Python form the div(B)=0
              constraint (Sum_f area_f*sigma_f = 0 per element), the uniform-field RHS
              (n . H), and the dipole moment (Sum_f sigma_f*area_f*(c_f - c_e)).

              Args:
                  intrc_handle: Interaction handle from BuildMatrix()

              Returns:
                  numpy array (dof x 11).
          )pbdoc");

    m.def("GetCentroidFieldGrad", &radia_solver::GetCentroidFieldGrad,
          py::arg("intrc_handle"),
          R"pbdoc(
              Per-hex centroid demag field + gradient functionals -- the kernel of the
              parameter-free yano-MSC moment formulation (replaces the eval-point alpha and
              the finite-difference conditioning noise; see examples/vim/
              yano_moment_analytic_selfterm.py).

              For each hex element, the demag field H and gradient gradH at the element
              CENTROID as linear functionals of every source DOF charge:
                SELF face  -> bare charged-face field (interior centroid, finite, no center
                              charge -- patch-test exact: single cube -> demag N=1/3);
                MUTUAL face -> yano dipole layer = bare face - area*(point @ source center)
                              (finite distance -> singularity-free).
              Field convention H = (1/4pi) int sigma (r-r')/|r-r'|^3 dA'.

              Args:
                  intrc_handle: Interaction handle from BuildMatrix()

              Returns:
                  numpy array (nHex, 9, dof): component k in
                  (Hx,Hy,Hz, gxx,gyy,gzz,gxy,gxz,gyz), source DOF g.  Use C[e,0:3,:] as the
                  centroid field functional F0 and C[e,3:9,:] as the symmetric gradient Ginv.
          )pbdoc");

    m.def("BuildMomentSystem", &radia_solver::BuildMomentSystem,
          py::arg("intrc_handle"), py::arg("chi"), py::arg("hx"), py::arg("hy"), py::arg("hz"),
          R"pbdoc(
              Parameter-free MOMENT-yano system matrix A (dof x dof, row-major) + rhs (dof) for a
              UNIFORM linear material (chi) in a UNIFORM applied field (hx,hy,hz).  Per hex (6 face-
              charge DOF): 3 dipole + 1 monopole + 2 diagonal-quadrupole rows = moment of sigma matched
              to chi*{H,gradH}(centroid) via GetCentroidFieldGrad.  Step-1 verification of the
              EIEM2 -> moment-yano upgrade (matches examples/vim/yano_moment_iter_scaling.py::build).

              Returns: (A, rhs, dof) -- A is (dof, dof), rhs is (dof,).
          )pbdoc");

    m.def("HMatrixDensify", &radia_solver::HMatrixDensify,
          py::arg("intrc_handle"),
          R"pbdoc(
              Densify the actual HACApK (ACA+) system operator.

              Builds the MSC H-matrix for the interaction handle and applies it
              to unit vectors, returning the dense A = -N + diag(1/chi) in the
              original DOF ordering. Use to validate the H-matrix against the
              exact dense matrix (eigenvalues, deflation).

              Args:
                  intrc_handle: Interaction handle from BuildMatrix()
              Returns:
                  Tuple (matrix, dof) where matrix is (dof x dof) numpy array
          )pbdoc");

    m.def("HLUDebugMaterialize", &radia_solver::HLUDebugMaterialize,
          py::arg("intrc_handle"),
          R"pbdoc(
              Phase 4 debug: materialize the post-convert tree as dense.
              Returns (A_perm, lod) where A_perm is in HACApK's permuted
              ordering and lod[i] = 0-based original index of permuted
              position i. Compare A_perm against P A_orig P^T to verify
              the conversion + tree-building is correct.
          )pbdoc");

    m.def("HLUSetTruncTol", [](double tol) { cHACApK_hlu_set_trunc_tol(tol); },
          py::arg("tol"),
          R"pbdoc(
              Set the rk-leaf SVD recompression tolerance for H-LU. KEY
              speed/accuracy knob: 1e-14 = machine precision (ranks grow,
              slow); 1e-4 = ACA accuracy (ranks low, intended O(N log^2 N)).
          )pbdoc");
    m.def("HLUGetTruncTol", []() { return cHACApK_hlu_get_trunc_tol(); });

    m.def("HLULastTimings", []() -> py::dict {
        double t_decomp = 0, t_solve = 0; long n_lu = 0, n_gemm = 0;
        cHACApK_hlu_get_timings(&t_decomp, &t_solve, &n_lu, &n_gemm);
        py::dict d;
        d["t_decomp_sec"] = t_decomp;
        d["t_solve_sec"] = t_solve;
        d["n_dense_lu"] = n_lu;
        d["n_dense_gemm"] = n_gemm;
        return d;
    },
    R"pbdoc(
        Return timing + op counts from the most recent H-LU decomp/solve
        (cHACApK_hlu_decomp + cHACApK_hlu_solve_vec). Call right after
        rad.HLUTestOnHACApK() to get the factorization and solve wall times.
    )pbdoc");

    m.def("HLUMaterializeStats", []() -> py::dict {
        long n_calls = 0, n_elems = 0, n_internal = 0, n_leaf = 0;
        cHACApK_hlu_get_materialize_stats(&n_calls, &n_elems);
        cHACApK_hlu_get_materialize_split(&n_internal, &n_leaf);
        py::dict d;
        d["n_calls"] = n_calls;
        d["n_elems"] = n_elems;
        d["n_internal"] = n_internal;  /* internal-subtree densifications (the cubic driver; must be 0) */
        d["n_leaf"] = n_leaf;          /* benign dense-leaf data copies for leaf-level dgemms */
        return d;
    },
    R"pbdoc(
        Return the Phase 3.6 materialize-fallback profile from the most
        recent H-LU decomp: n_calls (how many mixed leaf+internal nodes
        were densified) and n_elems (total densified matrix entries, a
        proxy for peak scratch memory). Both should be 0 for balanced
        (power-of-2 element-count) trees; large values indicate the
        materialize-and-redo path dominating on unbalanced trees.
    )pbdoc");

    m.def("HLUSetParallel", [](bool on) { cHACApK_hlu_set_parallel(on ? 1 : 0); },
          py::arg("on"),
          "Enable/disable block-level parallelism (ngcore TaskManager) in the "
          "H-LU h_addmul recursion. Default on.");
    m.def("HLUGetParallel", []() { return cHACApK_hlu_get_parallel() != 0; });
    m.def("HLUSetParCutoff", [](long c) { cHACApK_hlu_set_par_cutoff(c); },
          py::arg("cutoff"),
          "Minimum C-block area (rows*cols) above which h_addmul parallelizes "
          "its output blocks. Below this it runs serial (avoids tiny tasks).");
    m.def("HLUMaxThreads", []() { return chacapk_max_threads(); },
          "ngcore TaskManager max threads available to the H-LU.");
    m.def("HLUSetAccumCap", [](int c) { cHACApK_hlu_set_accum_cap(c); },
          py::arg("cap"),
          "Accumulator (lazy recompression) rank cap: an rk leaf accumulates "
          "low-rank updates by column-append and recompresses only when its "
          "rank exceeds 'cap' (plus a final flush). cap=0 disables (recompress "
          "every update = previous behavior). Default 64.");
    m.def("HLUGetAccumCap", []() { return cHACApK_hlu_get_accum_cap(); });


    m.def("HLUMixedBreakdown", []() -> py::dict {
        long a[9] = {0}, l[9] = {0}, r[9] = {0};
        cHACApK_hlu_get_mixed_breakdown(a, l, r);
        auto pack = [](long *v) {
            py::dict d;  const char *kn[3] = {"internal","rk","dense"};
            for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
                d[(std::string(kn[i])+"x"+kn[j]).c_str()] = v[i*3+j];
            return d;
        };
        py::dict out;
        out["addmul"] = pack(a);  /* keyed kindA x kindB */
        out["lln"] = pack(l);     /* kindL x kindX */
        out["run"] = pack(r);     /* kindU x kindX */
        return out;
    },
    R"pbdoc(
        Breakdown of which operand-kind combinations triggered the H-LU
        materialize fallback in the most recent decomp, for h_addmul / the
        two htrsm directions. Each value is a count keyed "<kindA>x<kindB>"
        (kind in {internal, rk, dense}). Tells whether the rk-factored
        mixed-multiply optimization (rk operand) or sub-view split (dense
        operand) is the dominant case to optimize.
    )pbdoc");

    m.def("HLUTestOnHACApK", &radia_solver::HLUTestOnHACApK,
          py::arg("intrc_handle"),
          R"pbdoc(
              Phase 4 smoke test: run H-LU on a real Radia HACApK tree.

              Builds a fresh HACApK MSC manager for the given interaction
              handle, computes y = A * x via MatVec on a deterministic test
              vector, converts the leafmtxp leaves to internal column-major
              layout, factors via cHACApK_hlu_decomp, solves x_solved = A^-1 y,
              and returns the max relative error
              max|x_solved - x| / max|x| in the original DOF ordering.

              Returns:
                  float  -- < 0 if internal failure (negative sentinel from
                            cHACApK_hlu_run_on_hacapk), otherwise the
                            max-elementwise relative error vs the test vector.
                            Expected ~ aca_eps + machine-precision rounding.
          )pbdoc");

    m.def("GetClusterStrategy", []() -> int {
        return cHACApK_get_cluster_strategy();
    },
    "Return the current HACApK cluster strategy (0=BBOX, 1=PCA).");

    // -- H-LU self-test (Phase 1/2 dense-leaf) --
    m.def("HLUSelfTest", [](int depth, int n_per_block) -> double {
        return cHACApK_harith_self_test(depth, n_per_block);
    },
    py::arg("depth") = 1, py::arg("n_per_block") = 100,
    R"pbdoc(
        Self-test the H-LU dense-leaf implementation on a synthetic
        (2^depth) x (2^depth) block-tree of dense leaves. Returns the max
        relative error between the H-LU solve and a reference LAPACK dgesv
        solve on the same matrix. Expected < 1e-10 for diag-dominant random.
        Total matrix size N = (2^depth) * n_per_block; depth=2 means 16 leaves,
        depth=3 means 64 leaves -- the latter exercises deep recursion.
    )pbdoc");

    // -- H-LU rk-aware self-test (Phase 3 partial: rk off-diag at depth=1) --
    m.def("HLUSelfTestRk", [](int n_per_block, int rk_rank) -> double {
        return cHACApK_harith_self_test_rk(n_per_block, rk_rank);
    },
    py::arg("n_per_block") = 100, py::arg("rk_rank") = 5,
    R"pbdoc(
        Phase 3 partial validation: build a 2x2 block-tree with DENSE diagonal
        leaves (random diag-dominant) and explicit-rank RK off-diagonal leaves
        (A_ij = U_ij V_ij^T from random U_ij, V_ij of rank rk_rank), then
        run cHACApK_hlu_decomp + solve and compare against LAPACK dgesv.

        Exercises the Phase 3 rk paths:
          - htrsm_lln(L=dense, X=rk),  htrsm_run(U=dense, X=rk)
          - h_addmul(rk*rk -> dense) on the trailing update
          - hmatvec_subtract on rk leaves in forward/backward sweeps

        Returns max relative error (should be ~ machine precision).
    )pbdoc");

    // -- Phase 3.5 unit test: h_addmul rk*rk -> rk + recompression --
    m.def("HLUSelfTestAddmulRkRk", [](int m_, int n_, int inner_,
                                        int kA, int kB, int kC) -> double {
        return cHACApK_harith_self_test_addmul_rkrk(m_, n_, inner_, kA, kB, kC);
    },
    py::arg("m") = 64, py::arg("n") = 64, py::arg("inner") = 64,
    py::arg("kA") = 5, py::arg("kB") = 5, py::arg("kC") = 5,
    R"pbdoc(
        Phase 3.5 unit test: rk(A) * rk(B) -> rk(C) h_addmul with ACA recompression.
        Builds three random rk leaves of given ranks, computes the dense ground
        truth (U_c V_c^T + alpha A B), runs h_addmul, and verifies the result by
        comparing C's new dense reconstruction to the truth. Expected ~1e-13.
    )pbdoc");

    m.def("HLUSelfTestRadiaExactWithMatrix",
          [](py::array_t<double, py::array::c_style | py::array::forcecast> A,
             py::array_t<double, py::array::c_style | py::array::forcecast> b) -> double {
        if (A.ndim() != 2 || A.shape(0) != 162 || A.shape(1) != 162)
            throw std::runtime_error("A must be 162x162");
        if (b.ndim() != 1 || b.shape(0) != 162)
            throw std::runtime_error("b must be length 162");
        /* Convert row-major numpy to column-major for the test. */
        std::vector<double> A_colmajor(162 * 162);
        auto Ar = A.unchecked<2>();
        for (int j = 0; j < 162; j++)
            for (int i = 0; i < 162; i++)
                A_colmajor[(size_t)i + (size_t)j * 162] = Ar(i, j);
        return cHACApK_harith_self_test_radia_exact_with_matrix(
            A_colmajor.data(), b.data());
    },
    py::arg("A"), py::arg("b"),
    R"pbdoc(
        Phase 4 debug: same Radia tree shape but with caller-provided
        162x162 matrix A and 162 RHS b. Use to test whether H-LU works
        on the EXACT real Radia matrix data, isolating the bug from
        the HACApK build/permutation path.
    )pbdoc");

    m.def("HLUSelfTestRadiaExactDiag", [](double diag_boost) -> double {
        return cHACApK_harith_self_test_radia_exact_diag(diag_boost);
    },
    py::arg("diag_boost") = 2.0,
    R"pbdoc(
        Phase 4 debug: Radia-shape with adjustable diag_boost. Sweep low
        values (1.5, 1.2, 1.05) to find no-pivot LU stability threshold.
        At boost ~ row_sum (mildly dominant), no-pivot can lose stability.
    )pbdoc");

    // -- Phase 4 debug: EXACT mimic of Radia nx=3 leaf=10 sizes --
    m.def("HLUSelfTestRadiaExact", []() -> double {
        return cHACApK_harith_self_test_radia_exact();
    },
    R"pbdoc(
        Phase 4 debug: hardcoded mimic of Radia nx=3 leaf=10 tree shape
        + EXACT cluster sizes (108/54 at root, 72/36 at TL, 48/24 at TL.TL).
        Synthetic random matrix, no permutation, no layout conversion.
        If this fails, the bug is in the recursive H-LU with rectangular
        non-uniform leaves. If it passes, the bug is in the permutation
        or layout conversion path.
    )pbdoc");

    // -- Phase 4 debug: depth-3 asymmetric (mimics Radia exact shape) --
    m.def("HLUSelfTestDepth3Asymmetric", [](int nb_tiny) -> double {
        return cHACApK_harith_self_test_depth3_asymmetric(nb_tiny);
    },
    py::arg("nb_tiny") = 3,
    R"pbdoc(
        Phase 4 debug: synthetic test mimicking Radia's nx=3 leaf=10 tree
        shape exactly. 10 leaves, 3 internal nodes, depth 3. Only the
        TL.TL sub-block goes deep (4 leaves at depth 3); everything else
        is leaves at depth 1 or 2.

        If mixed_sibling (depth 2) passes but this fails, the bug is in
        the deeper recursion with mixed leaf+internal across multiple
        levels.
    )pbdoc");

    // -- Phase 4 debug: mixed_sibling via HACApK row-major -> convert path --
    m.def("HLUSelfTestMixedSiblingViaConversion", [](int nb_small) -> double {
        return cHACApK_harith_self_test_mixed_sibling_via_conversion(nb_small);
    },
    py::arg("nb_small") = 5,
    R"pbdoc(
        Phase 4 debug: same shape as HLUSelfTestMixedSibling but leaves are
        built in HACApK row-major format then transposed before H-LU.
        Isolates the conversion path from the H-LU correctness.
    )pbdoc");

    // -- Phase 4 debug: non-uniform mixed-sibling test --
    m.def("HLUSelfTestMixedSiblingNonUniform", [](int n1, int n2, int m1, int m3) -> double {
        return cHACApK_harith_self_test_mixed_sibling_nonuniform(n1, n2, m1, m3);
    },
    py::arg("n1") = 5, py::arg("n2") = 7, py::arg("m1") = 2, py::arg("m3") = 3,
    R"pbdoc(
        Phase 4 debug: non-uniform-size mixed-sibling test. Reproduces
        HACApK's asymmetric element-count splits (13 -> 6+7 etc.) in a
        synthetic setting. Root splits (n1, n2), TL splits (m1, n1-m1),
        BR splits (m3, n2-m3). All sub-leaves have different shapes.

        If uniform passes and this fails, the bug is in non-uniform leaf
        sub-views (sub-view offset calc in Phase 3.6 mixed cases or
        materialize_node_as_dense).
    )pbdoc");

    // -- Phase 4 debug: mixed-sibling test (mimics Radia tree shape) --
    m.def("HLUSelfTestMixedSibling", [](int nb_small) -> double {
        return cHACApK_harith_self_test_mixed_sibling(nb_small);
    },
    py::arg("nb_small") = 5,
    R"pbdoc(
        Phase 4 debug: synthetic test that mimics Radia's nx=3 leaf=10 tree
        shape (10 dense leaves, 3 internal nodes, depth 3) using uniform
        leaf sizes (nb_small per side). Root has 2x2 children: TL & BR are
        internal (each containing 2x2 sub-leaves), TR & BL are leaves.

        If this passes while real Radia trees fail, the bug is in
        non-uniform leaf sizes. If it fails, the bug is in mixed-sibling
        recursion itself.
    )pbdoc");

    // -- Phase 3.5 integration test: depth=2 H-LU with rk off-diagonals --
    m.def("HLUSelfTestRkDeep", [](int n_per_block, int rk_rank) -> double {
        return cHACApK_harith_self_test_rk_deep(n_per_block, rk_rank);
    },
    py::arg("n_per_block") = 30, py::arg("rk_rank") = 5,
    R"pbdoc(
        Phase 3.5 integration test: depth=2 H-LU (4x4 leaf grid) with
        DENSE diagonal leaves and RK off-diagonal leaves of rank rk_rank.
        Exercises the full Phase 1 + 2 + 3 partial + 3.5 pipeline including
        rk(A)*rk(B) -> rk(C) trailing updates inside off-diagonal sub-blocks
        plus ACA recompression. Returns max relative error vs LAPACK dgesv.
    )pbdoc");


    // ========================================================================
    // Utility
    // ========================================================================

    m.def("UtiDel", &radia_utility::UtiDel,
          py::arg("obj"),
          "Delete object.");

    m.def("UtiDelAll", &radia_utility::UtiDelAll,
          "Delete all objects.");

    m.def("ObjGeoVol", &radia_utility::ObjGeoVol,
          py::arg("obj"),
          "Get object geometry volume.");

    m.def("ObjDegFre", &radia_utility::ObjDegFre,
          py::arg("obj"),
          "Get number of degrees of freedom.");

    // ========================================================================
    // Extended Object Creation
    // ========================================================================

    m.def("ObjThckPgn", &radia_objects_ext::ObjThckPgn,
          py::arg("xc"), py::arg("lx"), py::arg("polygon"), py::arg("axis"), py::arg("magnetization"),
          R"pbdoc(
              Create extruded polygon block.

              Args:
                  xc: Center position along extrusion axis
                  lx: Length along extrusion axis
                  polygon: List of 2D polygon vertices [[y1,z1], [y2,z2], ...]
                  axis: Extrusion axis ("x", "y", or "z")
                  magnetization: Magnetization vector [Mx, My, Mz] in A/m

              Returns:
                  Object handle
          )pbdoc");

    m.def("ObjCylMag", &radia_objects_ext::ObjCylMag,
          py::arg("center"), py::arg("r"), py::arg("h"), py::arg("nseg"),
          py::arg("axis"), py::arg("magnetization"),
          R"pbdoc(
              Create cylindrical magnet.

              Args:
                  center: Center point [x, y, z]
                  r: Radius
                  h: Height
                  nseg: Number of segments (polygonal approximation)
                  axis: Cylinder axis ("x", "y", or "z")
                  magnetization: Magnetization vector [Mx, My, Mz] in A/m

              Returns:
                  Object handle
          )pbdoc");

    m.def("ObjRecCur", &radia_objects_ext::ObjRecCur,
          py::arg("center"), py::arg("dimensions"), py::arg("current_density"),
          R"pbdoc(
              Create rectangular current block.

              Args:
                  center: Center point [x, y, z]
                  dimensions: Block dimensions [Lx, Ly, Lz]
                  current_density: Current density vector [Jx, Jy, Jz] in A/m^2

              Returns:
                  Object handle
          )pbdoc");

    m.def("ObjArcCur", &radia_objects_ext::ObjArcCur,
          py::arg("center"), py::arg("radii"), py::arg("phi"), py::arg("h"),
          py::arg("nseg"), py::arg("man_auto"), py::arg("axis"), py::arg("j"),
          R"pbdoc(
              Create arc current coil.

              Args:
                  center: Center point [x, y, z]
                  radii: Inner and outer radii [r_in, r_out]
                  phi: Start and end angles [phi1, phi2] in radians
                  h: Height
                  nseg: Number of segments
                  man_auto: "man" (manual) or "auto" segment distribution
                  axis: Coil axis ("x", "y", or "z")
                  j: Current density in A/m^2

              Returns:
                  Object handle
          )pbdoc");

    m.def("ObjRaceTrk", &radia_objects_ext::ObjRaceTrk,
          py::arg("center"), py::arg("radii"), py::arg("lengths"), py::arg("h"),
          py::arg("nseg"), py::arg("man_auto"), py::arg("axis"), py::arg("j"),
          R"pbdoc(
              Create racetrack coil.

              Args:
                  center: Center point [x, y, z]
                  radii: Inner and outer radii [r_in, r_out]
                  lengths: Straight section lengths [Lx, Ly]
                  h: Height
                  nseg: Number of segments in arcs
                  man_auto: "man" (manual) or "auto" segment distribution
                  axis: Coil axis ("x", "y", or "z")
                  j: Current density in A/m^2

              Returns:
                  Object handle
          )pbdoc");

    m.def("ObjFlmCur", &radia_objects_ext::ObjFlmCur,
          py::arg("points"), py::arg("current"),
          R"pbdoc(
              Create filament current (wire).

              Args:
                  points: List of 3D points defining the wire path [[x,y,z], ...]
                  current: Current in Amperes

              Returns:
                  Object handle
          )pbdoc");

    m.def("ObjAddToCnt", &radia_objects_ext::ObjAddToCnt,
          py::arg("cnt"), py::arg("objects"),
          "Add objects to existing container.");

    m.def("ObjCntSize", &radia_objects_ext::ObjCntSize,
          py::arg("cnt"),
          "Get number of objects in container.");

    m.def("ObjCntStuf", &radia_objects_ext::ObjCntStuf,
          py::arg("cnt"),
          "Get list of object handles in container.");

    m.def("ObjDpl", &radia_objects_ext::ObjDpl,
          py::arg("obj"), py::arg("opt") = "",
          "Duplicate object.");

    m.def("ObjM", &radia_objects_ext::ObjM,
          py::arg("obj"),
          R"pbdoc(
              Get object magnetization.

              Returns:
                  For single object: dict with 'center' and 'magnetization' tuples.
                  For container: list of (center, magnetization) tuples.
          )pbdoc");

    m.def("ObjSetM", &radia_objects_ext::ObjSetM,
          py::arg("obj"), py::arg("magnetization"),
          "Set object magnetization.");

    m.def("ObjScaleCur", &radia_objects_ext::ObjScaleCur,
          py::arg("obj"), py::arg("scale"),
          "Scale current in object by factor.");

    // ========================================================================
    // Transformations
    // ========================================================================

    m.def("TrfTrsl", &radia_transform::TrfTrsl,
          py::arg("vector"),
          R"pbdoc(
              Create translation transformation.

              Args:
                  vector: Translation vector [vx, vy, vz]

              Returns:
                  Transformation handle
          )pbdoc");

    m.def("TrfRot", &radia_transform::TrfRot,
          py::arg("point"), py::arg("vector"), py::arg("phi"),
          R"pbdoc(
              Create rotation transformation.

              Args:
                  point: Point on rotation axis [x, y, z]
                  vector: Rotation axis vector [vx, vy, vz]
                  phi: Rotation angle in radians

              Returns:
                  Transformation handle
          )pbdoc");

    // TrfPlSym REMOVED (2026-01-31) - Use IMA symmetry instead

    m.def("TrfInv", &radia_transform::TrfInv,
          "Create inversion transformation (point reflection through origin).");

    m.def("TrfCmbL", &radia_transform::TrfCmbL,
          py::arg("orig_trf"), py::arg("trf"),
          "Combine transformations: orig_trf = trf * orig_trf (left multiply).");

    m.def("TrfCmbR", &radia_transform::TrfCmbR,
          py::arg("orig_trf"), py::arg("trf"),
          "Combine transformations: orig_trf = orig_trf * trf (right multiply).");

    // TrfMlt REMOVED (2026-01-31) - Use IMA symmetry instead
    // The shared-DOF approach was fundamentally incompatible with MSC 6DOF hexahedra

    m.def("TrfOrnt", &radia_transform::TrfOrnt,
          py::arg("obj"), py::arg("trf"),
          "Apply transformation to orient object (modifies in place).");

    // TrfZerPara REMOVED (2026-01-31) - Use IMA symmetry instead
    // TrfZerPerp REMOVED (2026-01-31) - Use IMA symmetry instead

    // ========================================================================
    // Extended Materials
    // ========================================================================

    m.def("MatPM", &radia_material_ext::MatPM,
          py::arg("Br"), py::arg("Hc"), py::arg("easy_axis"),
          R"pbdoc(
              Create permanent magnet material.

              Args:
                  Br: Remanent field [T]
                  Hc: Coercive force [A/m]
                  easy_axis: Easy axis direction [ex, ey, ez]

              Returns:
                  Material handle
          )pbdoc");

    m.def("MatMagFixed", &radia_material_ext::MatMagFixed,
          py::arg("magnetization"),
          R"pbdoc(
              Create fixed magnetization material.

              Args:
                  magnetization: Fixed magnetization [Mx, My, Mz] in A/m

              Returns:
                  Material handle
          )pbdoc");

    m.def("MatMagLinear", &radia_material_ext::MatMagLinear,
          py::arg("Br"), py::arg("Hc"), py::arg("easy_axis"),
          R"pbdoc(
              Create linear demagnetization permanent magnet material.

              Args:
                  Br: Remanent field [T]
                  Hc: Coercive force [A/m]
                  easy_axis: Easy axis direction [ex, ey, ez]

              Returns:
                  Material handle
          )pbdoc");

    m.def("MatMagCurve", &radia_material_ext::MatMagCurve,
          py::arg("bh_data"), py::arg("easy_axis"),
          R"pbdoc(
              Create permanent magnet material with B-H demagnetization curve.

              Args:
                  bh_data: List of [H, B] pairs (second quadrant)
                  easy_axis: Easy axis direction [ex, ey, ez]

              Returns:
                  Material handle
          )pbdoc");

    m.def("MatSatIsoFrm", &radia_material_ext::MatSatIsoFrm,
          py::arg("params"),
          R"pbdoc(
              Create nonlinear isotropic material from formula.

              Args:
                  params: List of parameter pairs [[ksi1, Ms1], [ksi2, Ms2], [ksi3, Ms3]]

              Returns:
                  Material handle
          )pbdoc");

    m.def("MatSatAniso", &radia_material_ext::MatSatAniso,
          py::arg("data_par"), py::arg("data_per"),
          R"pbdoc(
              Create nonlinear anisotropic material.

              Args:
                  data_par: M-H curve parallel to easy axis [[H, M], ...]
                  data_per: M-H curve perpendicular to easy axis [[H, M], ...]

              Returns:
                  Material handle
          )pbdoc");

    m.def("MatSatLamTab", &radia_material_ext::MatSatLamTab,
          py::arg("mh_data"), py::arg("packing_factor"), py::arg("normal"),
          R"pbdoc(
              Create laminated nonlinear material.

              Args:
                  mh_data: M-H curve [[H, M], ...]
                  packing_factor: Lamination packing factor (0-1)
                  normal: Lamination normal vector [nx, ny, nz]

              Returns:
                  Material handle
          )pbdoc");

    m.def("MatEnergyHysteresis", &radia_material_ext::MatEnergyHysteresis,
          py::arg("K"), py::arg("chi"), py::arg("f_k_tables"),
          py::arg("eps") = 1e-8,
          R"pbdoc(
              Create energy-based vector hysteresis material (Type 5).

              Refactor 2026-05-XX: Type 5 is now a thin subclass of Type 6
              (radTPlayHysteresisMaterial). The original Schur-complement
              Newton implementation (Henrotte-Egger formulation with
              per-hysteron J_k as variables) is structurally incompatible
              with Potter B-input shape functions (sign-indefinite per-k
              Hessian). The current implementation uses Type 6's 3D Newton
              with analytical Jacobian on F(B)=H(B)-H_target=0 in O(K) per
              element, while preserving the energy formulation
              W(B) = (1/2) nu_rev |B|^2 + sum_k G_k(|p_k|) as a property.

              Args:
                  K: Number of partial polarizations (play operators)
                  chi: Play thresholds chi_k = eta_k [Tesla], array of K
                      values. The "chi" naming is retained for API
                      compatibility; numerically these are the Play
                      thresholds eta_k from Hane-Sugahara identification.
                  f_k_tables: List of K tuples (r_array, f_array), shape
                      function tables.
                      r_array: |p| grid points [0, r_max], monotonically increasing.
                      f_array: f_k(|p|) shape function values (can be
                          negative for k >= 1 by Potter loop closure).
                  eps: Retained for ABI compatibility, unused. The
                      Bergqvist regularization |x|_eps had no role in the
                      rev/irrev separated formulation.

              Returns:
                  Material handle
          )pbdoc");

    m.def("MatPlayHysteresis", &radia_material_ext::MatPlayHysteresis,
          py::arg("K"), py::arg("eta"), py::arg("f_k_tables"),
          R"pbdoc(
              Create direct B-input play hysteresis material.

              Unlike energy-based hysteresis, this model evaluates B->H directly
              in O(K) without Newton iteration. Shape functions f_k can be negative
              (sign-unconstrained). Inverse (H->B) uses Newton with analytical Jacobian.

              Args:
                  K: Number of play operators
                  eta: Play thresholds [Tesla], array of K values
                  f_k_tables: List of K tuples (r_array, f_array), shape function tables.
                      r_array: |p| grid points [0, r_max], monotonically increasing.
                      f_array: f_k(|p|) shape function values (can be negative).

              Returns:
                  Material handle
          )pbdoc");

    m.def("MatMvsH", &radia_material_ext::MatMvsH,
          py::arg("obj"), py::arg("component"), py::arg("h_field"),
          R"pbdoc(
              Get magnetization from material at given H field.

              Args:
                  obj: Material handle
                  component: "mx", "my", "mz", or "abs" for magnitude
                  h_field: H field [Hx, Hy, Hz] in A/m

              Returns:
                  Magnetization component(s)
          )pbdoc");

    m.def("MatHysSaveState", &radia_material_ext::MatHysSaveState,
          py::arg("mat"),
          R"pbdoc(
              Save internal state of a hysteresis material (Energy or Play model).

              Returns the full state (prev, pinning, current vectors for all K
              operators) as a flat numpy array. Use with MatHysRestoreState() to
              save/restore state during Picard iteration.

              Args:
                  mat: Material handle from MatEnergyHysteresis() or MatPlayHysteresis()

              Returns:
                  numpy array of state values (length K*9)
          )pbdoc");

    m.def("MatHysRestoreState", &radia_material_ext::MatHysRestoreState,
          py::arg("mat"), py::arg("state"),
          R"pbdoc(
              Restore internal state of a hysteresis material (Energy or Play model).

              Args:
                  mat: Material handle from MatEnergyHysteresis() or MatPlayHysteresis()
                  state: State array from MatHysSaveState()
          )pbdoc");

    m.def("MatHysCommitState", &radia_material_ext::MatHysCommitState,
          py::arg("mat"),
          R"pbdoc(
              Commit current state for the next time step.

              After Picard iteration converges, call this to commit the converged
              state as the reference for the next quasi-static step.

              Args:
                  mat: Material handle from MatEnergyHysteresis() or MatPlayHysteresis()
          )pbdoc");

    m.def("MatHysGetNuRev", &radia_material_ext::MatHysGetNuRev,
          py::arg("mat"),
          R"pbdoc(
              Get reversible reluctivity nu_rev for energy-based decomposition.

              H(B) = nu_rev * B + H_irr(B, history)

              nu_rev is computed automatically at material construction as the
              maximum dH/dB on the virgin curve. Used by Hantila solver for
              the constant stiffness matrix (LU factored once).

              Args:
                  mat: Material handle from MatPlayHysteresis()

              Returns:
                  float: nu_rev in A/m/T
          )pbdoc");

    m.def("MatHysIrreversible", &radia_material_ext::MatHysIrreversible,
          py::arg("mat"), py::arg("B"),
          R"pbdoc(
              Compute irreversible field H_irr(B) = H(B) - nu_rev * B.

              For Hantila polarization method:
              - Constant matrix: nu_rev (LU factored once)
              - Nonlinear residual: H_irr (updated each iteration)

              Args:
                  mat: Material handle from MatPlayHysteresis()
                  B: numpy array [Bx, By, Bz] in Tesla

              Returns:
                  numpy array [Hx_irr, Hy_irr, Hz_irr] in A/m
          )pbdoc");

    // ========================================================================
    // Extended Solver Functions
    // ========================================================================

    m.def("SolveNonl", &radia_solver_ext::SolveNonl,
          py::arg("obj"), py::arg("prec"), py::arg("max_iter"),
          py::arg("method"), py::arg("nonl_method"), py::arg("image") = "",
          R"pbdoc(
              Solve with specific nonlinear iteration method.

              Args:
                  obj: Object handle
                  prec: Convergence precision
                  max_iter: Maximum iterations
                  method: Linear solver (0=LU, 1=BiCGSTAB, 2=HACApK)
                  nonl_method: Nonlinear method
                  image: Image symmetry string (e.g., "+x", "-z", "+x-z")

              Returns:
                  Tuple (residual, max_M, avg_M, iterations)
          )pbdoc");

    m.def("SolverConfig", &radia_solver_ext::SolverConfig,
          R"pbdoc(
              Configure solver parameters (unified API).

              All parameters are optional keyword arguments. Only specified
              parameters are changed; others retain their current values.

              Keyword Args:
                  hacapk_eps (float): H-matrix ACA tolerance (default: 1e-4)
                  hacapk_leaf (int): H-matrix minimum cluster size (default: 10)
                  hacapk_eta (float): H-matrix admissibility parameter (default: 2.0)
                  hmatrix_eps (float): H-matrix field evaluation epsilon
                  bicgstab_tol (float): BiCGSTAB convergence tolerance (default: 1e-4)
                  relax_param (float): Under-relaxation (0=full step, <1=damped)
                  newton_method (bool): True=Newton-Raphson, False=Picard (default)
                  newton_damping (bool): Enable Newton line search damping
                  newton_damping_max_iter (int): Max line search iterations (default: 5)
                  newton_damping_min_omega (float): Minimum omega (default: 0.01)
                  b_input_newton (bool): Enable B-input Newton for hysteresis (default: False)
                  b_input_hantila (bool): Enable B-input Hantila for hysteresis (default: False)
                  hantila_alpha (float): Hantila polarization parameter (0=auto, default: 0)
                  hantila_relax (float): Hantila under-relaxation (0=full step, default: 0)

              Example:
                  rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
                  rad.SolverConfig(bicgstab_tol=1e-6, relax_param=0.3)
                  rad.SolverConfig(newton_method=True, newton_damping=True)
                  rad.SolverConfig(b_input_newton=True)  # For energy-based hysteresis
                  rad.SolverConfig(b_input_hantila=True)  # Hantila polarization (faster)
          )pbdoc");

    m.def("GetSolverConfig", &radia_solver_ext::GetSolverConfig,
          R"pbdoc(
              Get current solver configuration.

              Returns:
                  Dictionary with all solver parameters:
                  - bicgstab_tol, relax_param, newton_method, b_input_newton, b_input_hantila
                  - newton_damping, newton_damping_max_iter, newton_damping_min_omega
                  - hacapk_stats (if H-matrix solve has been performed)
          )pbdoc");

    // Image symmetry functions REMOVED (2026-01-31)
    // SetIMASymmetry, BuildIMAMatrix, etc. are replaced by the unified API:
    //   rad.Solve(obj, prec, maxiter, method, image='+x-z')
    //   rad.BuildMatrix(obj, image='+x-z')
    // The 'image' parameter specifies mirror symmetry: "+x", "-z", "+x-z", etc.

    // ========================================================================
    // Extended Field Functions
    // ========================================================================

    // FldEnr binding REMOVED (Phase C, 2026-04-16)

    m.def("FldFrc", &radia_field_ext::FldFrc,
          py::arg("obj"), py::arg("shape"),
          R"pbdoc(
              Compute force on object using Maxwell stress tensor.

              Args:
                  obj: Object handle
                  shape: Integration shape handle

              Returns:
                  Array [Fx, Fy, Fz, Tx, Ty, Tz] (force and torque)
          )pbdoc");

    m.def("FldFrcShpRtg", &radia_field_ext::FldFrcShpRtg,
          py::arg("center"), py::arg("dimensions"),
          R"pbdoc(
              Create rectangular integration shape for force calculation.

              Args:
                  center: Center point [x, y, z]
                  dimensions: Rectangle dimensions [wx, wy]

              Returns:
                  Shape handle
          )pbdoc");


    m.def("FldCmpCrt", &radia_field_ext::FldCmpCrt,
          py::arg("prcB"), py::arg("prcA"), py::arg("prcBInt"),
          py::arg("prcFrc"), py::arg("prcTrjCrd"), py::arg("prcTrjAng"),
          "Set field computation precision criteria.");

    m.def("FldLenRndSw", &radia_field_ext::FldLenRndSw,
          py::arg("on_off"),
          "Switch field lens/rendering mode ('on' or 'off').");

    // ========================================================================
    // Extended Field Functions
    // ========================================================================

    // FldLst - Field along line segment
    m.def("FldLst", [](int obj, const std::string& field_id,
                       const py::list& p1, const py::list& p2,
                       int np, const std::string& arg_opt = "noarg",
                       double start = 0.0) -> py::object {
        auto v1 = to_vector(p1.cast<py::object>());
        auto v2 = to_vector(p2.cast<py::object>());
        if (v1.size() != 3 || v2.size() != 3) {
            throw std::runtime_error("Points must have 3 coordinates");
        }

        // Determine output size based on arg_opt
        bool with_arg = (arg_opt == "arg");
        int n_out_per_point = with_arg ? 4 : 3;  // x,y,z or s,x,y,z or scalar

        // Check if field_id is a single component or vector
        bool is_vector = field_id.empty() || field_id == "b" || field_id == "h" ||
                         field_id == "a" || field_id == "m";
        int values_per_point = is_vector ? 3 : 1;
        int total_out = np * (with_arg ? (1 + values_per_point) : values_per_point);

        std::vector<double> result(total_out);
        int nB = 0;

        char id_cstr[16];
        strncpy(id_cstr, field_id.c_str(), 15);
        id_cstr[15] = '\0';

        char arg_cstr[16];
        strncpy(arg_cstr, arg_opt.c_str(), 15);
        arg_cstr[15] = '\0';

        int err = RadFldLst(result.data(), &nB, obj, id_cstr, v1.data(), v2.data(),
                            np, arg_cstr, start);
        check_error(err);

        // Return as list
        py::list out;
        if (with_arg) {
            int stride = 1 + values_per_point;
            for (int i = 0; i < np; ++i) {
                py::list row;
                for (int j = 0; j < stride; ++j) {
                    row.append(result[i * stride + j]);
                }
                out.append(row);
            }
        } else {
            if (is_vector) {
                for (int i = 0; i < np; ++i) {
                    py::list row;
                    for (int j = 0; j < 3; ++j) {
                        row.append(result[i * 3 + j]);
                    }
                    out.append(row);
                }
            } else {
                for (int i = 0; i < np; ++i) {
                    out.append(result[i]);
                }
            }
        }
        return out;
    },
    py::arg("obj"), py::arg("field_id"), py::arg("p1"), py::arg("p2"),
    py::arg("np"), py::arg("arg_opt") = "noarg", py::arg("start") = 0.0,
    R"pbdoc(
        Compute field along line segment.

        Args:
            obj: Object handle
            field_id: Field type ("bx", "by", "bz", "b", "hx", etc.)
            p1: Start point [x, y, z]
            p2: End point [x, y, z]
            np: Number of points
            arg_opt: "arg" to include position, "noarg" for values only
            start: Starting position for arg

        Returns:
            List of field values or [position, values] pairs
    )pbdoc");

    // FldInt - Field integral
    m.def("FldInt", [](int obj, const std::string& inf_fin,
                       const std::string& field_id,
                       const py::list& p1, const py::list& p2) -> py::object {
        auto v1 = to_vector(p1.cast<py::object>());
        auto v2 = to_vector(p2.cast<py::object>());
        if (v1.size() != 3 || v2.size() != 3) {
            throw std::runtime_error("Points must have 3 coordinates");
        }

        double result[3] = {0, 0, 0};
        int nf = 0;

        char inf_cstr[16], id_cstr[16];
        strncpy(inf_cstr, inf_fin.c_str(), 15);
        strncpy(id_cstr, field_id.c_str(), 15);
        inf_cstr[15] = id_cstr[15] = '\0';

        int err = RadFldInt(result, &nf, obj, inf_cstr, id_cstr, v1.data(), v2.data());
        check_error(err);

        if (nf == 1) {
            return py::cast(result[0]);
        } else {
            py::list out;
            for (int i = 0; i < nf; ++i) {
                out.append(result[i]);
            }
            return out;
        }
    },
    py::arg("obj"), py::arg("inf_fin"), py::arg("field_id"),
    py::arg("p1"), py::arg("p2"),
    R"pbdoc(
        Compute field integral along line.

        Args:
            obj: Object handle
            inf_fin: "inf" for infinite line, "fin" for finite segment
            field_id: Integral component ("ibx", "iby", "ibz", or "")
            p1: Start point [x, y, z]
            p2: End point [x, y, z]

        Returns:
            Field integral value(s) in Tesla*mm
    )pbdoc");

    // ObjCenFld - Object center and field
    m.def("ObjCenFld", [](int obj, const std::string& field_type) -> py::dict {
        double B[6] = {0};  // center + field
        int arMesh[3] = {1, 1, 1};

        char type_char = field_type.empty() ? 'B' : toupper(field_type[0]);

        int err = RadObjCenFld(B, arMesh, obj, type_char);
        check_error(err);

        py::dict result;
        result["center"] = py::make_tuple(B[0], B[1], B[2]);
        result["field"] = py::make_tuple(B[3], B[4], B[5]);
        return result;
    },
    py::arg("obj"), py::arg("field_type") = "B",
    R"pbdoc(
        Get object center point and field at that point.

        Args:
            obj: Object handle
            field_type: "A", "B", "H", "J", or "M"

        Returns:
            Dictionary with 'center' and 'field' tuples
    )pbdoc");

    // Use CERN Xsuite/Xtrack for GPU-accelerated beam tracking.

    // FldEnrFrc / FldEnrTrq (energy-based force/torque) REMOVED (Phase C, 2026-04-16)

    // FldCmpPrc - Set computation precision (string options version)
    m.def("FldCmpPrc", [](const std::string& opt) {
        char opt_cstr[256];
        strncpy(opt_cstr, opt.c_str(), 255);
        opt_cstr[255] = '\0';

        int n = 0;
        int err = RadFldCmpPrc(&n, opt_cstr);
        check_error(err);
    },
    py::arg("opt"),
    R"pbdoc(
        Set computation precision using options string.

        Args:
            opt: Options string like "PrcB->1e-6,PrcA->1e-5,..."
    )pbdoc");

    // FldLenTol - Length tolerance
    m.def("FldLenTol", [](double abs_val, double rel_val, double zero_val = 0.0) {
        int n = 0;
        int err = RadFldLenTol(&n, abs_val, rel_val, zero_val);
        check_error(err);
    },
    py::arg("abs_val"), py::arg("rel_val"), py::arg("zero_val") = 0.0,
    R"pbdoc(
        Set length tolerance for randomization.

        Args:
            abs_val: Absolute tolerance
            rel_val: Relative tolerance
            zero_val: Zero threshold (default 0)
    )pbdoc");

    // ========================================================================
    // Extended Object Creation Functions
    // ========================================================================

    // ObjMltExtPgn - Multiple extruded polygons
    m.def("ObjMltExtPgn", [](const py::list& slices,
                              const py::list& magnetization = py::list()) -> int {
        // Parse slices: [[[[x11,y11],[x12,y12],...],z1], ...]
        int ns = static_cast<int>(py::len(slices));
        if (ns < 2) {
            throw std::runtime_error("At least 2 slices required");
        }

        std::vector<double> flatVert;
        std::vector<int> slicesLen;
        std::vector<double> attitudes;

        for (const auto& slice : slices) {
            py::list s = slice.cast<py::list>();
            if (py::len(s) != 2) {
                throw std::runtime_error("Each slice must be [polygon, altitude]");
            }

            py::list polygon = s[0].cast<py::list>();
            double z = s[1].cast<double>();

            attitudes.push_back(z);
            slicesLen.push_back(static_cast<int>(py::len(polygon)));

            for (const auto& pt : polygon) {
                auto coord = to_vector(pt.cast<py::object>());
                if (coord.size() != 2) {
                    throw std::runtime_error("Each 2D point must have 2 coordinates");
                }
                flatVert.push_back(coord[0]);
                flatVert.push_back(coord[1]);
            }
        }

        double M[3] = {0, 0, 0};
        if (py::len(magnetization) >= 3) {
            auto m = to_vector(magnetization.cast<py::object>());
            M[0] = m[0]; M[1] = m[1]; M[2] = m[2];
        }

        int n = 0;
        int err = RadObjMltExtPgn(&n, flatVert.data(), slicesLen.data(),
                                   attitudes.data(), ns, M);
        check_error(err);
        return n;
    },
    py::arg("slices"), py::arg("magnetization") = py::list(),
    R"pbdoc(
        Create polyhedron from extruded polygon slices.

        Args:
            slices: List of [polygon_2d, altitude] pairs
            magnetization: [Mx, My, Mz] in A/m

        Returns:
            Object handle
    )pbdoc");

    // ObjMltExtRtg - Multiple extruded rectangles
    m.def("ObjMltExtRtg", [](const py::list& slices,
                              const py::list& magnetization = py::list()) -> int {
        int ns = static_cast<int>(py::len(slices));
        if (ns < 2) {
            throw std::runtime_error("At least 2 slices required");
        }

        std::vector<double> flatCenPts;
        std::vector<double> flatRtgSizes;

        for (const auto& slice : slices) {
            py::list s = slice.cast<py::list>();
            if (py::len(s) != 2) {
                throw std::runtime_error("Each slice must be [[x,y,z], [wx,wy]]");
            }

            auto center = to_vector(s[0].cast<py::object>());
            auto size = to_vector(s[1].cast<py::object>());

            if (center.size() != 3 || size.size() != 2) {
                throw std::runtime_error("Center must be [x,y,z], size must be [wx,wy]");
            }

            flatCenPts.insert(flatCenPts.end(), center.begin(), center.end());
            flatRtgSizes.insert(flatRtgSizes.end(), size.begin(), size.end());
        }

        double M[3] = {0, 0, 0};
        if (py::len(magnetization) >= 3) {
            auto m = to_vector(magnetization.cast<py::object>());
            M[0] = m[0]; M[1] = m[1]; M[2] = m[2];
        }

        int n = 0;
        int err = RadObjMltExtRtg(&n, flatCenPts.data(), flatRtgSizes.data(), ns, M);
        check_error(err);
        return n;
    },
    py::arg("slices"), py::arg("magnetization") = py::list(),
    R"pbdoc(
        Create polyhedron from rectangular slices.

        Args:
            slices: List of [[x,y,z], [wx,wy]] pairs
            magnetization: [Mx, My, Mz] in A/m

        Returns:
            Object handle
    )pbdoc");

    // ObjMltExtTri - disabled legacy Triangle-based API
    m.def("ObjMltExtTri", [](double xc, double lx,
                              const py::list& vertices,
                              const py::list& subdiv,
                              const std::string& axis = "x",
                              const py::list& magnetization = py::list(),
                              const std::string& opt = "") -> int {
        (void)xc;
        (void)lx;
        (void)vertices;
        (void)subdiv;
        (void)axis;
        (void)magnetization;
        (void)opt;
        throw std::runtime_error(
            "ObjMltExtTri is disabled because Radia no longer bundles Triangle. "
            "Use the Netgen/Cubit mesh workflow instead.");
    },
    py::arg("xc"), py::arg("lx"), py::arg("vertices"), py::arg("subdiv"),
    py::arg("axis") = "x", py::arg("magnetization") = py::list(),
    py::arg("opt") = "",
    R"pbdoc(
        Legacy triangulated extruded polygon API.

        This function is disabled because Radia no longer bundles Triangle.
        Use the Netgen/Cubit mesh workflow instead.

        Args:
            xc: Center position in extrusion direction
            lx: Thickness in extrusion direction
            vertices: 2D polygon vertices [[y1,z1], [y2,z2], ...]
            subdiv: Subdivision params [[k1,q1], [k2,q2], ...]
            axis: Extrusion axis ("x", "y", or "z")
            magnetization: [Mx, My, Mz] in A/m
            opt: Options string

        Returns:
            Object handle
    )pbdoc");

    // ObjArcPgnMag - Arc polygon magnet
    m.def("ObjArcPgnMag", [](const py::list& center,
                              const std::string& axis,
                              const py::list& vertices,
                              const py::list& phi_range,
                              int nseg,
                              const std::string& sym_nosym = "nosym",
                              const py::list& magnetization = py::list()) -> int {
        auto p = to_vector(center.cast<py::object>());
        if (p.size() != 2) {
            throw std::runtime_error("Center must have 2 coordinates [r, z]");
        }

        auto phi = to_vector(phi_range.cast<py::object>());
        if (phi.size() != 2) {
            throw std::runtime_error("Phi range must be [phi_min, phi_max]");
        }

        int nv = static_cast<int>(py::len(vertices));
        std::vector<double> flatVert;
        for (const auto& pt : vertices) {
            auto coord = to_vector(pt.cast<py::object>());
            if (coord.size() != 2) {
                throw std::runtime_error("Each vertex must have 2 coordinates");
            }
            flatVert.push_back(coord[0]);
            flatVert.push_back(coord[1]);
        }

        double M[3] = {0, 0, 0};
        if (py::len(magnetization) >= 3) {
            auto m = to_vector(magnetization.cast<py::object>());
            M[0] = m[0]; M[1] = m[1]; M[2] = m[2];
        }

        char a = axis.empty() ? 'z' : axis[0];
        char sym_no = (sym_nosym == "sym") ? 's' : 'n';

        int n = 0;
        int err = RadObjArcPgnMag(&n, p.data(), a, flatVert.data(), nv,
                                   phi.data(), nseg, sym_no, M);
        check_error(err);
        return n;
    },
    py::arg("center"), py::arg("axis"), py::arg("vertices"),
    py::arg("phi_range"), py::arg("nseg"),
    py::arg("sym_nosym") = "nosym", py::arg("magnetization") = py::list(),
    R"pbdoc(
        Create arc polygon magnet.

        Args:
            center: Axis position [x, y] or [r, z]
            axis: Rotation axis ("x", "y", or "z")
            vertices: 2D cross-section [[r1,z1], [r2,z2], ...]
            phi_range: [phi_min, phi_max] in radians
            nseg: Number of azimuthal segments
            sym_nosym: "sym" or "nosym"
            magnetization: [Mx, My, Mz] in A/m

        Returns:
            Object handle
    )pbdoc");

    // ========================================================================
    // Extended Material Functions
    // ========================================================================

    // MatSatLamFrm - Laminated saturable material (formula)
    m.def("MatSatLamFrm", [](const py::list& ksi_ms1,
                              const py::list& ksi_ms2,
                              const py::list& ksi_ms3,
                              double packing,
                              const py::list& normal) -> int {
        auto k1 = to_vector(ksi_ms1.cast<py::object>());
        auto k2 = to_vector(ksi_ms2.cast<py::object>());
        auto k3 = to_vector(ksi_ms3.cast<py::object>());
        auto n = to_vector(normal.cast<py::object>());

        if (k1.size() < 2) k1.resize(2, 0.0);
        if (k2.size() < 2) k2.resize(2, 0.0);
        if (k3.size() < 2) k3.resize(2, 0.0);
        if (n.size() != 3) {
            throw std::runtime_error("Normal must have 3 components");
        }

        int mat = 0;
        int err = RadMatSatLamFrm(&mat, k1.data(), k2.data(), k3.data(), packing, n.data());
        check_error(err);
        return mat;
    },
    py::arg("ksi_ms1"), py::arg("ksi_ms2"), py::arg("ksi_ms3"),
    py::arg("packing"), py::arg("normal"),
    R"pbdoc(
        Create laminated nonlinear material using formula.

        Args:
            ksi_ms1: [ksi1, ms1] first term
            ksi_ms2: [ksi2, ms2] second term
            ksi_ms3: [ksi3, ms3] third term
            packing: Packing factor (0 < p <= 1)
            normal: Lamination normal [nx, ny, nz]

        Returns:
            Material handle
    )pbdoc");

    // NOTE: Extended Conductor (CndHexahedron, CndWire, CndSpiral, CndDefinePort,
    // CndImpedanceSweep) and Coupled Solver (CplMagSetConductor, CplMagSweep)

    // ========================================================================
    // Extended Utility Functions
    // ========================================================================

    m.def("UtiVer", &radia_utility_ext::UtiVer,
          "Get Radia library version number.");

    // UtiDmp / UtiDmpPrs REMOVED (Phase B1, 2026-04-15) -
    // .rad save/load is no longer supported. Reconstruct objects from
    // user scripts instead.

    // ========================================================================
    // NGSolve CoefficientFunction: RadiaField
    // ========================================================================

    py::class_<ngfem::RadiaFieldCF,
               std::shared_ptr<ngfem::RadiaFieldCF>,
               ngfem::CoefficientFunction>(m, "RadiaField")
        .def(py::init<int, const std::string&,
                      std::optional<std::vector<double>>,
                      std::optional<std::vector<double>>,
                      std::optional<std::vector<double>>,
                      std::optional<std::vector<double>>,
                      std::optional<double>,
                      const std::string&>(),
             py::arg("radia_obj"),
             py::arg("field_type") = "b",
             py::arg("origin") = py::none(),
             py::arg("u_axis") = py::none(),
             py::arg("v_axis") = py::none(),
             py::arg("w_axis") = py::none(),
             py::arg("precision") = py::none(),
             py::arg("units") = "m",
             R"pbdoc(
                 NGSolve CoefficientFunction for Radia field evaluation.

                 Creates a CoefficientFunction that can be used directly with
                 GridFunction.Set() for field projection onto FE spaces.

                 Args:
                     radia_obj: Radia object handle
                     field_type: 'b', 'h', 'a', 'm', or 'phi'
                     origin: Translation [x,y,z] in meters
                     u_axis, v_axis, w_axis: Local coordinate axes (auto-normalized)
                     precision: Computation precision in Tesla
                     units: Must be 'm' (meters)

                 Example:
                     B_cf = rad.RadiaField(magnet, 'b')
                     gf = GridFunction(HDiv(mesh, order=2))
                     gf.Set(B_cf)

                     # VoxelCoefficient for trajectory:
                     B_voxel = B_cf.as_voxel_cf(mesh, resolution=61)
             )pbdoc")
        .def_readonly("radia_obj", &ngfem::RadiaFieldCF::radia_obj)
        .def_readonly("field_type", &ngfem::RadiaFieldCF::field_type)
        .def_readonly("use_transform", &ngfem::RadiaFieldCF::use_transform)
        .def_readonly("precision", &ngfem::RadiaFieldCF::precision)
        .def("PrepareCache", &ngfem::RadiaFieldCF::PrepareCache,
             py::arg("points"),
             "Pre-cache field values at given points for fast gf.Set()")
        .def("ClearCache", &ngfem::RadiaFieldCF::ClearCache,
             "Clear cached field values")
        .def("GetCacheStats", &ngfem::RadiaFieldCF::GetCacheStats,
             "Get cache statistics: enabled, size, hits, misses, hit_rate")
        .def("as_voxel_cf", &ngfem::RadiaFieldCF::AsVoxelCF,
             py::arg("mesh"), py::arg("resolution") = 41,
             R"pbdoc(
                 Create VoxelCoefficient for fast repeated evaluation.

                 Pre-computes field on a regular grid with trilinear interpolation.
                 Ideal for trajectory calculations where evaluation speed matters.

                 Args:
                     mesh: NGSolve mesh (for bounding box)
                     resolution: Grid points per dimension (default: 41)

                 Returns:
                     NGSolve CoefficientFunction (VoxelCoefficient-based)
             )pbdoc");

    // ================================================================
    // High-order BND node computation for GMSH export (C++ accelerated)
    // ================================================================
    m.def("_compute_ho_bnd_nodes", [](std::shared_ptr<ngcomp::MeshAccess> ma, int curve_order) {
        auto result = radia::ComputeHighOrderBndNodes(ma, curve_order);

        // Convert nodes to numpy array (N, 3)
        size_t nn = result.nodes.size();
        py::array_t<double> nodes_arr({(py::ssize_t)nn, (py::ssize_t)3});
        auto na = nodes_arr.mutable_unchecked<2>();
        for (size_t i = 0; i < nn; i++) {
            na(i, 0) = result.nodes[i][0];
            na(i, 1) = result.nodes[i][1];
            na(i, 2) = result.nodes[i][2];
        }

        py::dict out;
        out["nodes"] = nodes_arr;
        out["elem_conn"] = py::cast(result.elem_conn);
        out["elem_materials"] = py::cast(result.elem_materials);
        out["elem_gmsh_types"] = py::cast(result.elem_gmsh_types);
        out["elem_orig_idx"] = py::cast(result.elem_orig_idx);
        out["n_vertices"] = result.n_vertices;
        return out;
    }, py::arg("mesh"), py::arg("curve_order"),
    R"pbdoc(
        Compute high-order BND node positions for GMSH export (C++ accelerated).

        Evaluates mesh.GetTrafo() at GMSH Lagrange reference points for all
        BND elements. Edge nodes are cached across shared edges.

        Replaces the Python loop in gmsh_post_export.py (~1000x faster).

        Args:
            mesh: NGSolve Mesh (must have Curve(p) applied)
            curve_order: polynomial order (1..5)

        Returns:
            dict with keys: nodes (N,3 array), elem_conn (list of lists),
            elem_materials (list of str), elem_gmsh_types (list of int),
            elem_orig_idx (list of int), n_vertices (int)
    )pbdoc");

    // ========================================================================
    // HACApK PEEC adapter sanity check (Step 3 of HACApK-PEEC integration)
    // ========================================================================
    m.def("_TestPEECHACApKSanity", &RadHACApKPEECSanityCheck, py::arg("n_filaments"),
          R"pbdoc(
              Internal sanity check for RadHACApKPEECManager.

              Builds N parallel filaments along +z spaced 1 cm apart, forms
              the dense Ruehli L matrix and the HACApK L matrix, runs MatVec
              with a deterministic test vector in both, and returns the max
              relative error (HACApK vs dense). Used during development of
              the HACApK-PEEC adapter; should return a value below aca_eps
              (internally set to 1e-8).

              Args:
                  n_filaments: number of filaments (>= 2)

              Returns:
                  float: max relative error of HACApK MatVec vs dense L @ x.
                  Negative values indicate failure (-1 bad arg, -2 build
                  failed, -3 ndof mismatch).
          )pbdoc");

    // ========================================================================
    // HACApK PEEC manager (Step 4 stage 1: pybind class exposure)
    // ========================================================================
    //
    // Self-contained wrapper that owns its own PEECMatrixBuilder and
    // RadHACApKPEECManager. Inputs are numpy arrays of filament geometry
    // (no dependency on PyPEECBuilder which lives in a separate pybind
    // module). The Python-side complex BiCGSTAB combines two real MatVec
    // calls per complex iteration to handle (R + jωL + Zs) systems —
    // see hacapk_peec_prima_plan.md "Option A".
    //
    py::class_<radia::PEECMatrixBuilder>(m, "_PEECBuilderInternal");
    py::class_<RadHACApKPEECManager>(m, "_HACApKPEECManagerInternal");

    class PyHACApKPEECManager {
    public:
        PyHACApKPEECManager(py::array_t<double, py::array::c_style | py::array::forcecast> centers,
                            py::array_t<double, py::array::c_style | py::array::forcecast> directions,
                            py::array_t<double, py::array::c_style | py::array::forcecast> lengths,
                            py::array_t<double, py::array::c_style | py::array::forcecast> widths,
                            py::array_t<double, py::array::c_style | py::array::forcecast> heights,
                            py::array_t<double, py::array::c_style | py::array::forcecast> sigmas)
        {
            auto c = centers.unchecked<2>();
            auto d = directions.unchecked<2>();
            auto l = lengths.unchecked<1>();
            auto w = widths.unchecked<1>();
            auto h = heights.unchecked<1>();
            auto s = sigmas.unchecked<1>();
            const int64_t n = static_cast<int64_t>(c.shape(0));
            if (c.shape(1) != 3) throw std::invalid_argument("centers must have shape (N, 3)");
            if (d.shape(0) != n || d.shape(1) != 3) throw std::invalid_argument("directions must have shape (N, 3)");
            if (l.shape(0) != n) throw std::invalid_argument("lengths must have shape (N,)");
            if (w.shape(0) != n) throw std::invalid_argument("widths must have shape (N,)");
            if (h.shape(0) != n) throw std::invalid_argument("heights must have shape (N,)");
            if (s.shape(0) != n) throw std::invalid_argument("sigmas must have shape (N,)");

            for (int64_t i = 0; i < n; ++i) {
                radia::PEECSegment seg;
                seg.center = TVector3d(c(i, 0), c(i, 1), c(i, 2));
                seg.direction = TVector3d(d(i, 0), d(i, 1), d(i, 2));
                seg.length = l(i);
                seg.width = w(i);
                seg.height = h(i);
                seg.sigma = s(i);
                builder_.AddSegment(seg);
            }
            manager_ = std::make_unique<RadHACApKPEECManager>(builder_);
        }

        bool BuildHMatrix(double aca_eps, int leaf_size, double eta, int max_rank, int print_level) {
            RadHACApKParams p = RadHACApKPEECDefaultParams();
            if (aca_eps > 0)   p.aca_eps   = aca_eps;
            if (leaf_size > 0) p.leaf_size = leaf_size;
            if (eta > 0)       p.eta       = eta;
            if (max_rank > 0)  p.max_rank  = max_rank;
            p.print_level = print_level;
            return manager_->BuildHMatrix(p);
        }

        py::array_t<double> MatVec(py::array_t<double, py::array::c_style | py::array::forcecast> x) {
            const int n = manager_->GetNDOF();
            auto xb = x.unchecked<1>();
            if ((int)xb.shape(0) != n) {
                throw std::invalid_argument("x size must equal NDOF (= number of filaments)");
            }
            std::vector<double> xv(n), yv(n);
            for (int i = 0; i < n; ++i) xv[i] = xb(i);
            manager_->MatVec(xv, yv);
            py::array_t<double> y(n);
            auto yb = y.mutable_unchecked<1>();
            for (int i = 0; i < n; ++i) yb(i) = yv[i];
            return y;
        }

        int GetNDOF() const { return manager_->GetNDOF(); }
        bool IsValid() const { return manager_->IsValid(); }

        py::dict GetStats() const {
            const auto& s = manager_->GetStats();
            py::dict d;
            d["n_dof"] = s.n_dof;
            d["n_leaves"] = s.n_leaves;
            d["n_lowrank"] = s.n_lowrank;
            d["n_dense"] = s.n_dense;
            d["max_rank"] = s.max_rank;
            d["compression"] = s.compression;
            d["build_time"] = s.build_time;
            d["memory_mb"] = s.memory_mb;
            d["dense_memory_mb"] = s.dense_memory_mb;
            return d;
        }

    private:
        radia::PEECMatrixBuilder builder_;
        std::unique_ptr<RadHACApKPEECManager> manager_;
    };

    py::class_<PyHACApKPEECManager>(m, "HACApKPEECManager",
        R"pbdoc(
            HACApK adapter for PEEC filament inductance.

            Builds the H-matrix for the real symmetric L (Ruehli mutual
            inductance) of a set of straight filaments. The frequency-
            dependent system (R + jωL + Zs) is assembled in Python by
            combining real MatVec(L) calls with diagonal R/Zs terms
            (Option A pattern; see docs).

            Args:
                centers: (N, 3) filament center coordinates [m]
                directions: (N, 3) filament unit direction vectors
                lengths: (N,) filament lengths [m]
                widths: (N,) filament widths [m]
                heights: (N,) filament heights [m]
                sigmas: (N,) filament conductivities [S/m]
        )pbdoc")
        .def(py::init<py::array_t<double, py::array::c_style | py::array::forcecast>,
                       py::array_t<double, py::array::c_style | py::array::forcecast>,
                       py::array_t<double, py::array::c_style | py::array::forcecast>,
                       py::array_t<double, py::array::c_style | py::array::forcecast>,
                       py::array_t<double, py::array::c_style | py::array::forcecast>,
                       py::array_t<double, py::array::c_style | py::array::forcecast>>(),
             py::arg("centers"), py::arg("directions"), py::arg("lengths"),
             py::arg("widths"), py::arg("heights"), py::arg("sigmas"))
        .def("BuildHMatrix", &PyHACApKPEECManager::BuildHMatrix,
             py::arg("aca_eps") = -1.0, py::arg("leaf_size") = -1,
             py::arg("eta") = -1.0, py::arg("max_rank") = -1,
             py::arg("print_level") = 0,
             R"pbdoc(
                 Build the H-matrix for L. Negative arguments use the
                 PEEC-tuned defaults (aca_eps=1e-4, leaf_size=128, eta=3.0,
                 max_rank=400). Returns True on success.
             )pbdoc")
        .def("MatVec", &PyHACApKPEECManager::MatVec, py::arg("x"),
             R"pbdoc(
                 Real matvec y = L * x. Both x and y are length-NDOF
                 (= n_filaments) double arrays.
             )pbdoc")
        .def("GetNDOF", &PyHACApKPEECManager::GetNDOF)
        .def("IsValid", &PyHACApKPEECManager::IsValid)
        .def("GetStats", &PyHACApKPEECManager::GetStats,
             "Return dict with n_dof, n_leaves, n_lowrank, n_dense, "
             "max_rank, compression, build_time [s], memory_mb, dense_memory_mb.");

    // ========================================================================
    // HACApK scalar BEM adapter (Phase 1.4 of in-tree SIBC HACApK pipeline)
    // ========================================================================
    //
    // Wraps RadHACApKBEMManager: takes a pre-computed dense Galerkin matrix
    // (built in Python via radia.bem.sibc_hacapk.assemble_SL_dense or
    // assemble_DL_dense) and wraps it as a HACApK H-matrix for ACA
    // compression and fast O(N log N) MatVec.
    //
    py::class_<RadHACApKBEMManager>(m, "_HACApKBEMManagerInternal");

    class PyHACApKBEMManager {
    public:
        PyHACApKBEMManager(py::array_t<double, py::array::c_style | py::array::forcecast> coords,
                            py::array_t<double, py::array::c_style | py::array::forcecast> entries)
            : coords_(coords)   // hold a Python ref so the buffer outlives manager_
            , entries_(entries)
        {
            auto c = coords_.unchecked<2>();
            auto e = entries_.unchecked<2>();
            const int64_t n = static_cast<int64_t>(c.shape(0));
            if (c.shape(1) != 3)
                throw std::invalid_argument("coords must have shape (N, 3)");
            if (e.shape(0) != n || e.shape(1) != n)
                throw std::invalid_argument("entries must have shape (N, N) matching coords");

            const double* c_ptr = static_cast<const double*>(coords_.data());
            const double* e_ptr = static_cast<const double*>(entries_.data());
            manager_ = std::make_unique<RadHACApKBEMManager>(c_ptr, e_ptr,
                                                              static_cast<int>(n));
        }

        bool BuildHMatrix(double aca_eps, int leaf_size, double eta,
                          int max_rank, int print_level) {
            RadHACApKParams p = RadHACApKBEMDefaultParams();
            if (aca_eps > 0)   p.aca_eps   = aca_eps;
            if (leaf_size > 0) p.leaf_size = leaf_size;
            if (eta > 0)       p.eta       = eta;
            if (max_rank > 0)  p.max_rank  = max_rank;
            p.print_level = print_level;
            return manager_->BuildHMatrix(p);
        }

        py::array_t<double> MatVec(py::array_t<double,
                                                py::array::c_style |
                                                py::array::forcecast> x) {
            const int n = manager_->GetNDOF();
            auto xb = x.unchecked<1>();
            if ((int)xb.shape(0) != n) {
                throw std::invalid_argument(
                    "x size must equal NDOF (= number of vertices)");
            }
            std::vector<double> xv(n), yv(n);
            for (int i = 0; i < n; ++i) xv[i] = xb(i);
            manager_->MatVec(xv, yv);
            py::array_t<double> y(n);
            auto yb = y.mutable_unchecked<1>();
            for (int i = 0; i < n; ++i) yb(i) = yv[i];
            return y;
        }

        int GetNDOF() const { return manager_->GetNDOF(); }
        bool IsValid() const { return manager_->IsValid(); }

        py::dict GetStats() const {
            const auto& s = manager_->GetStats();
            py::dict d;
            d["n_dof"] = s.n_dof;
            d["n_leaves"] = s.n_leaves;
            d["n_lowrank"] = s.n_lowrank;
            d["n_dense"] = s.n_dense;
            d["max_rank"] = s.max_rank;
            d["compression"] = s.compression;
            d["build_time"] = s.build_time;
            d["memory_mb"] = s.memory_mb;
            d["dense_memory_mb"] = s.dense_memory_mb;
            return d;
        }

    private:
        py::array_t<double, py::array::c_style | py::array::forcecast> coords_;
        py::array_t<double, py::array::c_style | py::array::forcecast> entries_;
        std::unique_ptr<RadHACApKBEMManager> manager_;
    };

    py::class_<PyHACApKBEMManager>(m, "HACApKBEMManager",
        R"pbdoc(
            HACApK adapter for a scalar Galerkin BEM matrix (Laplace SL/DL).

            Wraps a pre-computed dense (N, N) entry table as an
            ACA-compressed H-matrix.  Build the dense table once in Python
            via radia.bem.sibc_hacapk.assemble_SL_dense / assemble_DL_dense
            and hand it here for storage compression and fast MatVec.

            After BuildHMatrix() succeeds, the dense table can be discarded
            (still held internally by the manager via numpy refcounting,
            but the H-matrix is what MatVec uses).

            Args:
                coords: (N, 3) DOF coordinates [m] (= vertex coordinates
                        for P1 Galerkin)
                entries: (N, N) dense matrix, row-major.
        )pbdoc")
        .def(py::init<py::array_t<double, py::array::c_style | py::array::forcecast>,
                       py::array_t<double, py::array::c_style | py::array::forcecast>>(),
             py::arg("coords"), py::arg("entries"))
        .def("BuildHMatrix", &PyHACApKBEMManager::BuildHMatrix,
             py::arg("aca_eps") = -1.0, py::arg("leaf_size") = -1,
             py::arg("eta") = -1.0, py::arg("max_rank") = -1,
             py::arg("print_level") = 0,
             R"pbdoc(
                 Build the H-matrix.  Negative arguments use the BEM-tuned
                 defaults (aca_eps=1e-6, leaf_size=64, eta=2.0, max_rank=400).
                 Returns True on success.
             )pbdoc")
        .def("MatVec", &PyHACApKBEMManager::MatVec, py::arg("x"),
             R"pbdoc(
                 Real matvec y = M * x.  Both x and y are length-N double
                 arrays where N = number of vertices = NDOF.
             )pbdoc")
        .def("GetNDOF", &PyHACApKBEMManager::GetNDOF)
        .def("IsValid", &PyHACApKBEMManager::IsValid)
        .def("GetStats", &PyHACApKBEMManager::GetStats,
             "Return dict with n_dof, n_leaves, n_lowrank, n_dense, "
             "max_rank, compression, build_time [s], memory_mb, dense_memory_mb.");

    // ========================================================================
    // Fast C++ Galerkin SL/DL assembler -- Phase 1.9
    // ========================================================================
    m.def("_AssembleSLDL_Galerkin",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> verts,
           py::array_t<int64_t, py::array::c_style | py::array::forcecast> tris,
           py::array_t<double, py::array::c_style | py::array::forcecast> p2_nodes,
           int regular_quad_degree,
           int singular_n_q,
           int n_threads)
        {
            auto vb = verts.unchecked<2>();
            auto tb = tris.unchecked<2>();
            auto pb = p2_nodes.unchecked<3>();
            const int n_v = static_cast<int>(vb.shape(0));
            const int n_t = static_cast<int>(tb.shape(0));
            if (vb.shape(1) != 3) throw std::invalid_argument("verts must be (N, 3)");
            if (tb.shape(1) != 3) throw std::invalid_argument("tris must be (N_t, 3)");
            if (pb.shape(0) != n_t || pb.shape(1) != 6 || pb.shape(2) != 3)
                throw std::invalid_argument("p2_nodes must be (N_t, 6, 3)");

            py::array_t<double> SL({n_v, n_v});
            py::array_t<double> DL({n_v, n_v});
            // Direct pointers into numpy buffers
            const double* verts_ptr = static_cast<const double*>(verts.data());
            const int64_t* tris_ptr = static_cast<const int64_t*>(tris.data());
            const double* p2_ptr = static_cast<const double*>(p2_nodes.data());
            double* SL_ptr = static_cast<double*>(SL.mutable_data());
            double* DL_ptr = static_cast<double*>(DL.mutable_data());

            {
                py::gil_scoped_release rel;
                radia::bem::AssembleSLDL(
                    verts_ptr, n_v,
                    tris_ptr, n_t,
                    p2_ptr,
                    regular_quad_degree,
                    singular_n_q,
                    n_threads,
                    SL_ptr,
                    DL_ptr);
            }
            return py::make_tuple(SL, DL);
        },
        py::arg("verts"), py::arg("tris"), py::arg("p2_nodes"),
        py::arg("regular_quad_degree") = 11,
        py::arg("singular_n_q") = 8,
        py::arg("n_threads") = 0,
        R"pbdoc(
            Build dense SL and DL Galerkin matrices on a P1 H1 surface
            mesh with optional P2-curved geometry.  Uses Sauter-Schwab
            Duffy 4D quadrature for singular pairs and tensor-product
            Gauss for regular pairs.  Internally OpenMP-parallel.

            Args:
                verts: (n_v, 3) vertex coords [m]
                tris:  (n_t, 3) int64 triangle vertex indices into verts
                p2_nodes: (n_t, 6, 3) per-tri P2 Lagrange node coords in
                          order [v0, v1, v2, mid01, mid12, mid20].
                          For flat tris: pass corner mid-points (0.5*(v0+v1)
                          etc.).  For curved geometry: pass nodes
                          extracted via mesh.GetTrafo at the 6 P2 ref pts.
                regular_quad_degree: triangle Gauss degree for non-singular
                          pairs (default 11; uses Stroud 7-pt for <=5,
                          13-pt for <=7, Duffy n*n GL for higher).
                singular_n_q: 1D Gauss-Legendre order for SS Duffy sub-cubes.
                          n_q^4 nodes per sub-cube; 2/5/6 sub-cubes for
                          vertex/edge/identical (default 8 -> ~25k pts/pair).
                n_threads: 0 = OpenMP default (typically all cores).

            Returns:
                (SL, DL): tuple of (n_v, n_v) float64 ndarrays.
                          Convention matches NGSolve.bem LaplaceSL/LaplaceDL
                          (verified bit-exact to 1e-10).
        )pbdoc");

    // ========================================================================
    // Lagrange-P2 H1 variant of the C++ Galerkin SL/DL assembler.
    // 6 basis functions per triangle: 3 vertex + 3 edge mid-point hats.
    // ========================================================================
    m.def("_AssembleSLDL_Galerkin_P2",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> verts,
           py::array_t<int64_t, py::array::c_style | py::array::forcecast> tris,
           py::array_t<double, py::array::c_style | py::array::forcecast> p2_nodes,
           py::array_t<int64_t, py::array::c_style | py::array::forcecast> dofs_per_tri,
           int n_dof,
           int regular_quad_degree,
           int singular_n_q,
           int n_threads)
        {
            auto vb = verts.unchecked<2>();
            auto tb = tris.unchecked<2>();
            auto pb = p2_nodes.unchecked<3>();
            auto db = dofs_per_tri.unchecked<2>();
            const int n_v = static_cast<int>(vb.shape(0));
            const int n_t = static_cast<int>(tb.shape(0));
            if (vb.shape(1) != 3) throw std::invalid_argument("verts must be (N, 3)");
            if (tb.shape(1) != 3) throw std::invalid_argument("tris must be (N_t, 3)");
            if (pb.shape(0) != n_t || pb.shape(1) != 6 || pb.shape(2) != 3)
                throw std::invalid_argument("p2_nodes must be (N_t, 6, 3)");
            if (db.shape(0) != n_t || db.shape(1) != 6)
                throw std::invalid_argument("dofs_per_tri must be (N_t, 6)");
            if (n_dof <= 0)
                throw std::invalid_argument("n_dof must be positive");

            py::array_t<double> SL({n_dof, n_dof});
            py::array_t<double> DL({n_dof, n_dof});
            const double* verts_ptr = static_cast<const double*>(verts.data());
            const int64_t* tris_ptr = static_cast<const int64_t*>(tris.data());
            const double* p2_ptr = static_cast<const double*>(p2_nodes.data());
            const int64_t* dofs_ptr = static_cast<const int64_t*>(dofs_per_tri.data());
            double* SL_ptr = static_cast<double*>(SL.mutable_data());
            double* DL_ptr = static_cast<double*>(DL.mutable_data());

            {
                py::gil_scoped_release rel;
                radia::bem::AssembleSLDL_P2(
                    verts_ptr, n_v,
                    tris_ptr, n_t,
                    p2_ptr,
                    dofs_ptr,
                    n_dof,
                    regular_quad_degree,
                    singular_n_q,
                    n_threads,
                    SL_ptr,
                    DL_ptr);
            }
            return py::make_tuple(SL, DL);
        },
        py::arg("verts"), py::arg("tris"), py::arg("p2_nodes"),
        py::arg("dofs_per_tri"), py::arg("n_dof"),
        py::arg("regular_quad_degree") = 11,
        py::arg("singular_n_q") = 8,
        py::arg("n_threads") = 0,
        R"pbdoc(
            Build dense SL and DL Galerkin matrices on a P2 H1 LAGRANGE
            surface basis (6 hats per triangle: 3 vertices + 3 edge
            mid-points), with optional P2-curved geometry.  Output is
            (n_dof, n_dof) where n_dof = n_v + n_unique_edges.

            Args:
                verts: (n_v, 3) corner vertex coords [m]
                tris:  (n_t, 3) int64 corner-vertex indices
                p2_nodes: (n_t, 6, 3) per-tri Lagrange node coords in
                          order [v0, v1, v2, mid01, mid12, mid20]
                dofs_per_tri: (n_t, 6) int64 global Lagrange DOF indices
                          per local basis function, same node ordering
                          as p2_nodes.
                n_dof:    output matrix size.
                regular_quad_degree, singular_n_q: same semantics as
                          _AssembleSLDL_Galerkin (P1 entry).
                n_threads: 0 = OpenMP default.

            Returns:
                (SL, DL) of shape (n_dof, n_dof) float64.
        )pbdoc");

    // ========================================================================
    // Fast finite-segment Biot-Savart H (complex per-segment currents).
    // Replaces the pure-numpy compute_phi_inc_from_filaments inner kernel.
    // ========================================================================
    m.def("_HFromSegmentsComplex",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> segs,
           py::array_t<double, py::array::c_style | py::array::forcecast> obs,
           py::array_t<double, py::array::c_style | py::array::forcecast> I_re,
           py::array_t<double, py::array::c_style | py::array::forcecast> I_im,
           int n_threads)
        {
            auto sb = segs.unchecked<3>();
            auto ob = obs.unchecked<2>();
            auto rb = I_re.unchecked<1>();
            auto ib = I_im.unchecked<1>();
            const int n_seg = static_cast<int>(sb.shape(0));
            const int n_obs = static_cast<int>(ob.shape(0));
            if (sb.shape(1) != 2 || sb.shape(2) != 3)
                throw std::invalid_argument("segs must be (N_seg, 2, 3)");
            if (ob.shape(1) != 3)
                throw std::invalid_argument("obs must be (N_obs, 3)");
            if ((int)rb.shape(0) != n_seg || (int)ib.shape(0) != n_seg)
                throw std::invalid_argument("I_re/I_im must have length N_seg");

            py::array_t<double> H_re({n_obs, 3});
            py::array_t<double> H_im({n_obs, 3});

            const double* segs_ptr = static_cast<const double*>(segs.data());
            const double* obs_ptr  = static_cast<const double*>(obs.data());
            const double* re_ptr   = static_cast<const double*>(I_re.data());
            const double* im_ptr   = static_cast<const double*>(I_im.data());
            double* hre_ptr        = static_cast<double*>(H_re.mutable_data());
            double* him_ptr        = static_cast<double*>(H_im.mutable_data());

            {
                py::gil_scoped_release rel;
                radia::bs::HFromSegmentsComplex(
                    segs_ptr, n_seg, obs_ptr, n_obs,
                    re_ptr, im_ptr,
                    hre_ptr, him_ptr,
                    n_threads);
            }
            return py::make_tuple(H_re, H_im);
        },
        py::arg("segs"), py::arg("obs"),
        py::arg("I_re"), py::arg("I_im"),
        py::arg("n_threads") = 0,
        R"pbdoc(
            Finite-segment Biot-Savart H at obs points with complex per-segment
            currents.  Drop-in replacement for the pure-numpy
            _h_segments_complex helper used by compute_phi_inc_from_filaments.

            Args:
                segs:  (N_seg, 2, 3) endpoints (p1, p2) for each segment [m]
                obs:   (N_obs, 3) observation points [m]
                I_re:  (N_seg,) real part of complex per-segment current [A]
                I_im:  (N_seg,) imag part
                n_threads: 0 = TaskManager default; > 0 = request that many.

            Returns:
                (H_re, H_im): tuple of (N_obs, 3) float64 arrays such that
                             H = H_re + 1j*H_im is the complex H [A/m].
        )pbdoc");

    // ========================================================================
    // Fast finite-segment vector potential A (complex per-segment currents).
    // Companion to _HFromSegmentsComplex.  Used by FilamentBundleAC.fld('a')
    // and by Telegen-reciprocity DeltaL = (1/I^2) * int_S J_s . A_inc dS.
    // ========================================================================
    m.def("_AFromSegmentsComplex",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> segs,
           py::array_t<double, py::array::c_style | py::array::forcecast> obs,
           py::array_t<double, py::array::c_style | py::array::forcecast> I_re,
           py::array_t<double, py::array::c_style | py::array::forcecast> I_im,
           int n_threads)
        {
            auto sb = segs.unchecked<3>();
            auto ob = obs.unchecked<2>();
            auto rb = I_re.unchecked<1>();
            auto ib = I_im.unchecked<1>();
            const int n_seg = static_cast<int>(sb.shape(0));
            const int n_obs = static_cast<int>(ob.shape(0));
            if (sb.shape(1) != 2 || sb.shape(2) != 3)
                throw std::invalid_argument("segs must be (N_seg, 2, 3)");
            if (ob.shape(1) != 3)
                throw std::invalid_argument("obs must be (N_obs, 3)");
            if ((int)rb.shape(0) != n_seg || (int)ib.shape(0) != n_seg)
                throw std::invalid_argument("I_re/I_im must have length N_seg");

            py::array_t<double> A_re({n_obs, 3});
            py::array_t<double> A_im({n_obs, 3});

            const double* segs_ptr = static_cast<const double*>(segs.data());
            const double* obs_ptr  = static_cast<const double*>(obs.data());
            const double* re_ptr   = static_cast<const double*>(I_re.data());
            const double* im_ptr   = static_cast<const double*>(I_im.data());
            double* are_ptr        = static_cast<double*>(A_re.mutable_data());
            double* aim_ptr        = static_cast<double*>(A_im.mutable_data());

            {
                py::gil_scoped_release rel;
                radia::bs::AFromSegmentsComplex(
                    segs_ptr, n_seg, obs_ptr, n_obs,
                    re_ptr, im_ptr,
                    are_ptr, aim_ptr,
                    n_threads);
            }
            return py::make_tuple(A_re, A_im);
        },
        py::arg("segs"), py::arg("obs"),
        py::arg("I_re"), py::arg("I_im"),
        py::arg("n_threads") = 0,
        R"pbdoc(
            Finite-segment vector potential A at obs points with complex
            per-segment currents.  Companion to _HFromSegmentsComplex.

            For each segment of length L from p1 to p2 with unit direction e_l,
            and obs P, with t = (P-p1) . e_l, r1 = |P-p1|, r2 = |P-p2|:

                A = (mu_0/(4*pi)) * I * log( ((L-t) + r2) / (-t + r1) ) * e_l

            Output is in T*m (SI).  Sums contributions from all segments.

            Args:
                segs:  (N_seg, 2, 3) endpoints (p1, p2) for each segment [m]
                obs:   (N_obs, 3) observation points [m]
                I_re:  (N_seg,) real part of complex per-segment current [A]
                I_im:  (N_seg,) imag part
                n_threads: 0 = TaskManager default; > 0 = request that many.

            Returns:
                (A_re, A_im): tuple of (N_obs, 3) float64 arrays such that
                             A = A_re + 1j*A_im is the complex A [T*m].
        )pbdoc");

    // ========================================================================
    // Fast surface-current Biot-Savart on triangulated surface
    // (piecewise-constant complex J_s per triangle).  Used to evaluate the
    // workpiece's induced B and A from the SIBC stream function.
    // ========================================================================
    m.def("_AFromTrianglesComplex",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> verts,
           py::array_t<double, py::array::c_style | py::array::forcecast> J_re,
           py::array_t<double, py::array::c_style | py::array::forcecast> J_im,
           py::array_t<double, py::array::c_style | py::array::forcecast> obs,
           int n_threads)
        {
            auto vb = verts.unchecked<3>();
            auto jr = J_re.unchecked<2>();
            auto ji = J_im.unchecked<2>();
            auto ob = obs.unchecked<2>();
            const int n_t   = static_cast<int>(vb.shape(0));
            const int n_obs = static_cast<int>(ob.shape(0));
            if (vb.shape(1) != 3 || vb.shape(2) != 3)
                throw std::invalid_argument("verts must be (N_t, 3, 3)");
            if ((int)jr.shape(0) != n_t || (int)ji.shape(0) != n_t
                || jr.shape(1) != 3 || ji.shape(1) != 3)
                throw std::invalid_argument("J_re/J_im must be (N_t, 3)");
            if (ob.shape(1) != 3)
                throw std::invalid_argument("obs must be (N_obs, 3)");

            py::array_t<double> A_re({n_obs, 3});
            py::array_t<double> A_im({n_obs, 3});

            const double* v_ptr = static_cast<const double*>(verts.data());
            const double* jr_ptr = static_cast<const double*>(J_re.data());
            const double* ji_ptr = static_cast<const double*>(J_im.data());
            const double* o_ptr  = static_cast<const double*>(obs.data());
            double* are = static_cast<double*>(A_re.mutable_data());
            double* aim = static_cast<double*>(A_im.mutable_data());

            {
                py::gil_scoped_release rel;
                radia::bs::AFromTrianglesComplex(
                    v_ptr, n_t, jr_ptr, ji_ptr,
                    o_ptr, n_obs, are, aim, n_threads);
            }
            return py::make_tuple(A_re, A_im);
        },
        py::arg("verts"), py::arg("J_re"), py::arg("J_im"),
        py::arg("obs"), py::arg("n_threads") = 0,
        R"pbdoc(
            Vector potential A at obs points from a triangulated surface
            with piecewise-constant complex J_s per triangle.

            Args:
                verts: (N_t, 3, 3) triangle vertices [m]
                J_re:  (N_t, 3) real part of per-triangle J_s [A/m]
                J_im:  (N_t, 3) imag part
                obs:   (N_obs, 3) observation points [m]
                n_threads: 0 = TaskManager default; > 0 = request that many.

            Returns:
                (A_re, A_im): tuple of (N_obs, 3) float64 arrays such that
                             A = A_re + 1j*A_im is the complex A [T*m].

            3-point Gauss quadrature; non-singular evaluation only.
        )pbdoc");

    m.def("_BFromTrianglesComplex",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> verts,
           py::array_t<double, py::array::c_style | py::array::forcecast> J_re,
           py::array_t<double, py::array::c_style | py::array::forcecast> J_im,
           py::array_t<double, py::array::c_style | py::array::forcecast> obs,
           int n_threads)
        {
            auto vb = verts.unchecked<3>();
            auto jr = J_re.unchecked<2>();
            auto ji = J_im.unchecked<2>();
            auto ob = obs.unchecked<2>();
            const int n_t   = static_cast<int>(vb.shape(0));
            const int n_obs = static_cast<int>(ob.shape(0));
            if (vb.shape(1) != 3 || vb.shape(2) != 3)
                throw std::invalid_argument("verts must be (N_t, 3, 3)");
            if ((int)jr.shape(0) != n_t || (int)ji.shape(0) != n_t
                || jr.shape(1) != 3 || ji.shape(1) != 3)
                throw std::invalid_argument("J_re/J_im must be (N_t, 3)");
            if (ob.shape(1) != 3)
                throw std::invalid_argument("obs must be (N_obs, 3)");

            py::array_t<double> B_re({n_obs, 3});
            py::array_t<double> B_im({n_obs, 3});

            const double* v_ptr = static_cast<const double*>(verts.data());
            const double* jr_ptr = static_cast<const double*>(J_re.data());
            const double* ji_ptr = static_cast<const double*>(J_im.data());
            const double* o_ptr  = static_cast<const double*>(obs.data());
            double* bre = static_cast<double*>(B_re.mutable_data());
            double* bim = static_cast<double*>(B_im.mutable_data());

            {
                py::gil_scoped_release rel;
                radia::bs::BFromTrianglesComplex(
                    v_ptr, n_t, jr_ptr, ji_ptr,
                    o_ptr, n_obs, bre, bim, n_threads);
            }
            return py::make_tuple(B_re, B_im);
        },
        py::arg("verts"), py::arg("J_re"), py::arg("J_im"),
        py::arg("obs"), py::arg("n_threads") = 0,
        R"pbdoc(
            Magnetic flux density B at obs points from a triangulated
            surface with piecewise-constant complex J_s per triangle.

            Args:
                verts: (N_t, 3, 3) triangle vertices [m]
                J_re:  (N_t, 3) real part of per-triangle J_s [A/m]
                J_im:  (N_t, 3) imag part
                obs:   (N_obs, 3) observation points [m]
                n_threads: 0 = TaskManager default; > 0 = request that many.

            Returns:
                (B_re, B_im): tuple of (N_obs, 3) float64 arrays such that
                             B = B_re + 1j*B_im is the complex B [T].

            3-point Gauss quadrature; non-singular evaluation only.
        )pbdoc");

    // ========================================================================
    // Equivalence-theorem near-field source -- Phase A static H reconstruction
    // (Schelkunoff/Love).  See docs/equivalence_source/CPP_DESIGN.md.
    // Replaces the Python Stratton-Chu inner loop in
    // src/radia/equivalence_source.py:evaluate_static_H().
    // ========================================================================
    m.def("_EquivalenceSourceStaticH",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> centroids,
           py::array_t<double, py::array::c_style | py::array::forcecast> normals,
           py::array_t<double, py::array::c_style | py::array::forcecast> areas,
           py::array_t<double, py::array::c_style | py::array::forcecast> H_surf,
           py::array_t<double, py::array::c_style | py::array::forcecast> obs,
           int n_threads)
        {
            auto cb = centroids.unchecked<2>();
            auto nb = normals.unchecked<2>();
            auto ab = areas.unchecked<1>();
            auto hb = H_surf.unchecked<2>();
            auto ob = obs.unchecked<2>();
            const int n_faces = static_cast<int>(cb.shape(0));
            const int n_obs   = static_cast<int>(ob.shape(0));
            if (cb.shape(1) != 3 || nb.shape(1) != 3 || hb.shape(1) != 3)
                throw std::invalid_argument("centroids/normals/H_surf must be (N, 3)");
            if ((int)nb.shape(0) != n_faces || (int)ab.shape(0) != n_faces
                || (int)hb.shape(0) != n_faces)
                throw std::invalid_argument(
                    "centroids/normals/areas/H_surf must have matching N_faces");
            if (ob.shape(1) != 3)
                throw std::invalid_argument("obs must be (N_obs, 3)");

            py::array_t<double> H_out({n_obs, 3});

            const double* c_ptr = static_cast<const double*>(centroids.data());
            const double* n_ptr = static_cast<const double*>(normals.data());
            const double* a_ptr = static_cast<const double*>(areas.data());
            const double* h_ptr = static_cast<const double*>(H_surf.data());
            const double* o_ptr = static_cast<const double*>(obs.data());
            double* out_ptr = static_cast<double*>(H_out.mutable_data());

            {
                py::gil_scoped_release rel;
                radia::eqsrc::EvaluateStaticH(
                    c_ptr, n_ptr, a_ptr, h_ptr, n_faces,
                    o_ptr, n_obs, out_ptr, n_threads);
            }
            return H_out;
        },
        py::arg("centroids"), py::arg("normals"), py::arg("areas"),
        py::arg("H_surf"), py::arg("obs"), py::arg("n_threads") = 0,
        R"pbdoc(
            Equivalence-theorem magnetostatic H reconstruction at obs
            points from per-face surface sources on a closed
            triangulated surface.

            Stratton-Chu static reduction:
                H(r) = (1/(4 pi)) * sum_faces
                        { grad(1/R) x J_s  -  (n . H_s) grad(1/R) } * dS
            with J_s = n x H_s, R = |r - centroid|,
                 grad(1/R) = -R_vec / R^3.

            Args:
                centroids: (N_faces, 3) face centroids [m]
                normals:   (N_faces, 3) OUTWARD unit normals
                areas:     (N_faces,)   face areas [m^2]
                H_surf:    (N_faces, 3) H phasor REAL part at centroid
                                         [A/m].  Static reduction; the
                                         imaginary part is zero by
                                         definition for omega = 0.
                obs:       (N_obs, 3)   observation points [m].  MUST be
                                         outside the closed surface.
                n_threads: 0 = NGSolve TaskManager default;
                           > 0 = request that many.

            Returns:
                H_out: (N_obs, 3) float64 reconstructed H [A/m].

            Phase A (omega = 0 only); Phase B will add
            _EquivalenceSourceHarmonic for full dyadic GF.
        )pbdoc");

    // ========================================================================
    // Phase B: time-harmonic dyadic Stratton-Chu reconstruction.
    // Resolves Phase 2 KNOWN_LIMITATION (66% near-field undershoot from
    // missing (1/k^2) grad-grad psi term).
    // ========================================================================
    m.def("_EquivalenceSourceHarmonic",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> centroids,
           py::array_t<double, py::array::c_style | py::array::forcecast> normals,
           py::array_t<double, py::array::c_style | py::array::forcecast> areas,
           py::array_t<double, py::array::c_style | py::array::forcecast> E_re,
           py::array_t<double, py::array::c_style | py::array::forcecast> E_im,
           py::array_t<double, py::array::c_style | py::array::forcecast> H_re,
           py::array_t<double, py::array::c_style | py::array::forcecast> H_im,
           py::array_t<double, py::array::c_style | py::array::forcecast> obs,
           double omega,
           int n_threads)
        {
            auto cb = centroids.unchecked<2>();
            auto nb = normals.unchecked<2>();
            auto ab = areas.unchecked<1>();
            auto er = E_re.unchecked<2>();
            auto ei = E_im.unchecked<2>();
            auto hr = H_re.unchecked<2>();
            auto hi = H_im.unchecked<2>();
            auto ob = obs.unchecked<2>();
            const int n_faces = static_cast<int>(cb.shape(0));
            const int n_obs   = static_cast<int>(ob.shape(0));
            if (cb.shape(1) != 3 || nb.shape(1) != 3
                || er.shape(1) != 3 || ei.shape(1) != 3
                || hr.shape(1) != 3 || hi.shape(1) != 3)
                throw std::invalid_argument("centroids/normals/E_re/E_im/H_re/H_im must be (N, 3)");
            if ((int)nb.shape(0) != n_faces || (int)ab.shape(0) != n_faces
                || (int)er.shape(0) != n_faces || (int)ei.shape(0) != n_faces
                || (int)hr.shape(0) != n_faces || (int)hi.shape(0) != n_faces)
                throw std::invalid_argument("all face arrays must share N_faces");
            if (ob.shape(1) != 3)
                throw std::invalid_argument("obs must be (N_obs, 3)");
            if (omega <= 0.0)
                throw std::invalid_argument("omega must be > 0 for harmonic; use _EquivalenceSourceStaticH for omega = 0");

            py::array_t<double> Eo_re({n_obs, 3});
            py::array_t<double> Eo_im({n_obs, 3});
            py::array_t<double> Ho_re({n_obs, 3});
            py::array_t<double> Ho_im({n_obs, 3});

            {
                py::gil_scoped_release rel;
                radia::eqsrc::EvaluateHarmonic(
                    static_cast<const double*>(centroids.data()),
                    static_cast<const double*>(normals.data()),
                    static_cast<const double*>(areas.data()),
                    static_cast<const double*>(E_re.data()),
                    static_cast<const double*>(E_im.data()),
                    static_cast<const double*>(H_re.data()),
                    static_cast<const double*>(H_im.data()),
                    n_faces,
                    static_cast<const double*>(obs.data()),
                    n_obs, omega,
                    static_cast<double*>(Eo_re.mutable_data()),
                    static_cast<double*>(Eo_im.mutable_data()),
                    static_cast<double*>(Ho_re.mutable_data()),
                    static_cast<double*>(Ho_im.mutable_data()),
                    n_threads);
            }
            return py::make_tuple(Eo_re, Eo_im, Ho_re, Ho_im);
        },
        py::arg("centroids"), py::arg("normals"), py::arg("areas"),
        py::arg("E_re"), py::arg("E_im"),
        py::arg("H_re"), py::arg("H_im"),
        py::arg("obs"), py::arg("omega"),
        py::arg("n_threads") = 0,
        R"pbdoc(
            Equivalence-theorem time-harmonic (E, H) reconstruction at
            obs points with FULL dyadic Green's function.

            Full Stratton-Chu: includes the (1/k^2) grad-grad-psi term
            that the scalar form omits.  Resolves the deep-near-field
            undershoot (Phase 2 KNOWN_LIMITATION).

            Args:
                centroids: (N_faces, 3) face centroids [m]
                normals:   (N_faces, 3) OUTWARD unit normals
                areas:     (N_faces,)   face areas [m^2]
                E_re,E_im: (N_faces, 3) E phasor at centroid [V/m]
                H_re,H_im: (N_faces, 3) H phasor at centroid [A/m]
                obs:       (N_obs, 3)   observation points [m]
                omega:     angular frequency [rad/s] (> 0; use
                           _EquivalenceSourceStaticH for omega = 0)
                n_threads: 0 = TaskManager default

            Returns:
                (E_re, E_im, H_re, H_im) tuple of (N_obs, 3) arrays.

            Per docs/equivalence_source/CPP_DESIGN.md sec 3.2.
        )pbdoc");

    // ========================================================================
    // Closed-form cuboid average B (Wakao Part 6 §7) -- Phase beta v4.22.0
    // ========================================================================
    m.def("_average_B_in_box",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> M,
           py::array_t<double, py::array::c_style | py::array::forcecast> src_min,
           py::array_t<double, py::array::c_style | py::array::forcecast> src_max,
           py::array_t<double, py::array::c_style | py::array::forcecast> tgt_min,
           py::array_t<double, py::array::c_style | py::array::forcecast> tgt_max)
        {
            if (M.size() != 3 || src_min.size() != 3 || src_max.size() != 3
                || tgt_min.size() != 3 || tgt_max.size() != 3)
                throw std::invalid_argument("M, src_min, src_max, tgt_min, tgt_max must each have length 3");
            py::array_t<double> B(3);
            radia::average_field::AverageBInBox(
                static_cast<const double*>(M.data()),
                static_cast<const double*>(src_min.data()),
                static_cast<const double*>(src_max.data()),
                static_cast<const double*>(tgt_min.data()),
                static_cast<const double*>(tgt_max.data()),
                static_cast<double*>(B.mutable_data()));
            return B;
        },
        py::arg("M"), py::arg("src_min"), py::arg("src_max"),
        py::arg("tgt_min"), py::arg("tgt_max"),
        R"pbdoc(
            Spatial average of B over a target rectangular box from
            uniform magnetisation in a source rectangular box.

            Closed form: 64-corner sum of 3-fold antiderivatives G1, G2
            (Phase alpha derivation, sympy-verified). Both source and
            target must be axis-aligned cuboids. Overlapping cuboids are
            handled via a mu_0 * M * V_overlap / V_T correction term.

            Args:
                M:        (3,) magnetisation [A/m]
                src_min:  (3,) lower corner of source cuboid [m]
                src_max:  (3,) upper corner of source cuboid [m]
                tgt_min:  (3,) lower corner of target cuboid [m]
                tgt_max:  (3,) upper corner of target cuboid [m]

            Returns:
                B: (3,) average <B_x>, <B_y>, <B_z> over target box [T]

            Caveat: the 64-corner sum suffers ULP cancellation when V_T
            is much smaller than the source box; in that regime use the
            Gauss-Legendre numerical fallback in
            radia.analytical_formulas.cuboid_average_field.
        )pbdoc");

    m.def("_average_demag_tensor",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> src_min,
           py::array_t<double, py::array::c_style | py::array::forcecast> src_max,
           py::array_t<double, py::array::c_style | py::array::forcecast> tgt_min,
           py::array_t<double, py::array::c_style | py::array::forcecast> tgt_max)
        {
            if (src_min.size() != 3 || src_max.size() != 3
                || tgt_min.size() != 3 || tgt_max.size() != 3)
                throw std::invalid_argument("src_min, src_max, tgt_min, tgt_max must each have length 3");
            py::array_t<double> A({3, 3});
            radia::average_field::AverageDemagTensor(
                static_cast<const double*>(src_min.data()),
                static_cast<const double*>(src_max.data()),
                static_cast<const double*>(tgt_min.data()),
                static_cast<const double*>(tgt_max.data()),
                static_cast<double*>(A.mutable_data()));
            return A;
        },
        py::arg("src_min"), py::arg("src_max"),
        py::arg("tgt_min"), py::arg("tgt_max"),
        R"pbdoc(
            3x3 average demag tensor A_T such that
                <B>_T = mu_0 * (A_T @ M + M * V_overlap / V_T)
            for source/target axis-aligned cuboids.
        )pbdoc");

    // ========================================================================
    // (ACA+)+TSVD generic least-norm solver (kernel-agnostic).
    // The matrix entry A(i,j) is supplied by a Python callable entry(i,j),
    // so the same machinery serves any Radia source family (coil Biot-Savart,
    // MMM/MSC magnetic-material field, ...).  ACA+ is HACApK's cHACApK_acaplus.
    // ========================================================================
    m.def("_stream_aca_tsvd",
        [](int M, int N, std::function<double(int, int)> entry,
           int modes, int kmax, double aca_eps, int method)
        {
            if (M <= 0 || N <= 0)
                throw std::invalid_argument("M and N must be positive");
            if (!entry)
                throw std::invalid_argument("entry callback must be callable");
            const radia::stream_function::Method mth =
                (method == 2) ? radia::stream_function::Method::Method2
                              : radia::stream_function::Method::Method3;
            radia::stream_function::TSVDResult r =
                radia::stream_function::ACATSVD(
                    M, N, entry, modes, kmax, aca_eps, mth);

            py::array_t<double> U({r.M, r.modes});
            py::array_t<double> S(r.modes);
            py::array_t<double> V({r.N, r.modes});
            double* pU = static_cast<double*>(U.mutable_data());
            double* pS = static_cast<double*>(S.mutable_data());
            double* pV = static_cast<double*>(V.mutable_data());
            for (size_t i = 0; i < r.U.size(); ++i) pU[i] = r.U[i];
            for (size_t i = 0; i < r.S.size(); ++i) pS[i] = r.S[i];
            for (size_t i = 0; i < r.V.size(); ++i) pV[i] = r.V[i];
            return py::make_tuple(U, S, V, r.k_aca);
        },
        py::arg("M"), py::arg("N"), py::arg("entry"),
        py::arg("modes"), py::arg("kmax"),
        py::arg("aca_eps") = 1.0e-4, py::arg("method") = 3,
        R"pbdoc(
            (ACA+)+TSVD recompressed truncated SVD of an M x N matrix A whose
            entries are supplied on demand by the callback entry(i, j) ->
            A(i,j) (0-based row i in [0,M), col j in [0,N)).  Kernel-agnostic:
            the callback may use any Radia field computation (coil Biot-Savart,
            MMM/MSC magnetic material field, ...).  ACA+ is delegated to HACApK
            (cHACApK_acaplus).  Returns (U, S, V, k_aca):
              U:     (M, modes) row-major
              S:     (modes,)
              V:     (N, modes) row-major
              k_aca: ACA+ rank found before truncation
            method=3 (default) is the improved 2-SVD path (f90 method_aca_tsvd_2);
            method=2 is the full re-SVD of both factors (f90 method_aca_tsvd_1).
        )pbdoc");
}
