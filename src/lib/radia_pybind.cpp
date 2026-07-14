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
#include <atomic>
#include <complex>
#include <cmath>
#include <string>
#include <stdexcept>
#include <memory>
#include <unordered_map>
#include <unordered_set>
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
#include "rad_hdiv_field_evaluator.h" // Persistent RT1 field source + target tree
#include "rad_hacapk_hdiv.h"     // HACApK H-matrix for the HDiv-type VIM demag operator
#include "rad_planar_charges.h"  // Shared 2D planar exterior field + Maxwell torque
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

template <typename T>
std::vector<T> to_1d_vector(
        py::array_t<T, py::array::c_style | py::array::forcecast> arr,
        const char* name) {
    auto buf = arr.request();
    if (buf.ndim != 1)
        throw std::runtime_error(std::string(name) + " must be a 1D contiguous array");
    const T* ptr = static_cast<const T*>(buf.ptr);
    return std::vector<T>(ptr, ptr + buf.size);
}

py::dict solve_timings_dict(const RadHACApKChargeGram& s) {
    py::dict timings;
    for (const auto& kv : s.LastSolveTimings()) timings[py::str(kv.first)] = kv.second;
    return timings;
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

// NOTE: The Python-facing ObjRecMag wrapper (radia_objects::ObjRecMag) was
// removed when the Python ObjRecMag constructor was retired in favour of the
// fixed-M magnet_box path. The internal C++ surface-current rectangular block
// kernel (radTRecMag) and the C API RadObjRecMag remain in use elsewhere.

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
 * 5 surface-face coefficients for fixed-magnetization field evaluation.
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
 * @brief Create pyramid element from 5 vertices
 *
 * Square-base pyramid: 1 quadrilateral base + 4 triangular sides.
 * 5 faces total.  The object constructor is retained for fixed magnetization
 * field evaluation; mesh-backed magnetic-material solves route through HDiv-VIM.
 *
 * Vertex convention (matches netgen_mesh_import.PYRAMID_FACES):
 *   v0..v3 = base quad, v4 = apex.
 *   Base   : (v0, v3, v2, v1)  -- outward normal points away from the apex
 *   Sides  : (v0,v1,v4) (v1,v2,v4) (v2,v3,v4) (v3,v0,v4)
 *
 * @param vertices 5 vertices: base v0..v3 then apex v4
 * @param magnetization Magnetization vector [Mx, My, Mz] in A/m
 * @return Object handle
 */
int ObjPyramid(py::list vertices, py::array_t<double> magnetization) {
    if (py::len(vertices) != 5) {
        throw std::runtime_error("Pyramid requires exactly 5 vertices");
    }

    auto [flat_verts, nv] = to_vertex_array(vertices);
    auto m = magnetization.unchecked<1>();

    if (m.size() != 3) {
        throw std::runtime_error("magnetization must have 3 elements");
    }

    double M[3] = {m(0), m(1), m(2)};

    // Pyramid face definitions (1-indexed vertices; matches netgen_mesh_import.PYRAMID_FACES).
    // Face 0: base quad (v0, v3, v2, v1) - outward normal points away from the apex.
    // Faces 1-4: triangular sides (each base edge + apex v4).
    std::vector<int> flatFaces;
    int faceLengths[5] = {4, 3, 3, 3, 3};

    // base quad
    flatFaces.push_back(1); flatFaces.push_back(4); flatFaces.push_back(3); flatFaces.push_back(2);
    // side triangles
    flatFaces.push_back(1); flatFaces.push_back(2); flatFaces.push_back(5);
    flatFaces.push_back(2); flatFaces.push_back(3); flatFaces.push_back(5);
    flatFaces.push_back(3); flatFaces.push_back(4); flatFaces.push_back(5);
    flatFaces.push_back(4); flatFaces.push_back(1); flatFaces.push_back(5);

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
    // stats[8]/stats[9] were the moment loop-deflation cycles/alpha; surface-charge no longer connects to HACApK,
    // so they are always 0 and no longer surfaced.
    if (n >= 11) {
        result["t_moment_fieldgrad"] = stats[10];
    }
    if (n >= 12) {
        result["t_moment_system_build"] = stats[11];
    }

    return result;
}

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
// TrfMlt shared DOFs between original and virtual elements, which is incorrect for independent face coefficients.
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

/**
 * @brief Batched Forward evaluation for the HDiv hysteresis stepping:
 *        H[i] = nu_rev * B[i] + H_irr(B[i]; states[i]) for every row.
 *
 * Each row restores its COMMITTED state into the (single) material handle and
 * evaluates from it, so the call is PURE w.r.t. the states (Irreversible plays
 * committed -> scratch only; nothing commits).  The loop is a SERIAL C++ loop
 * with the GIL released: the handle's internal buffers are shared, so rows
 * must not run concurrently -- the win is removing the n x 2 Python<->C++
 * crossings of the per-row Python loop, which dominates large-n_el stepping.
 */
py::array_t<double> MatHysForwardBatch(
        int mat,
        py::array_t<double, py::array::c_style | py::array::forcecast> B,
        py::array_t<double, py::array::c_style | py::array::forcecast> states) {
    auto b = B.unchecked<2>();
    auto s = states.unchecked<2>();
    if (b.shape(1) != 3) throw std::runtime_error("MatHysForwardBatch: B must be (n, 3)");
    if (s.shape(0) != b.shape(0))
        throw std::runtime_error("MatHysForwardBatch: states rows must match B rows");
    const py::ssize_t n = b.shape(0);
    const int slen = (int)s.shape(1);
    double nu_rev = 0;
    if (RadMatHysGetNuRev(mat, &nu_rev) < 0)
        throw std::runtime_error("MatHysForwardBatch: not a Play hysteresis material");
    py::array_t<double> H({n, (py::ssize_t)3});
    auto h = H.mutable_unchecked<2>();
    {
        py::gil_scoped_release release;
        for (py::ssize_t i = 0; i < n; ++i) {
            if (RadMatHysRestoreState(mat, s.data(i, 0), slen) < 0)
                throw std::runtime_error("MatHysForwardBatch: state restore failed (row length mismatch?)");
            double pB[3] = {b(i, 0), b(i, 1), b(i, 2)};
            double pHirr[3] = {0, 0, 0};
            if (RadMatHysIrreversible(mat, pB, pHirr) < 0)
                throw std::runtime_error("MatHysForwardBatch: irreversible evaluation failed");
            h(i, 0) = nu_rev * pB[0] + pHirr[0];
            h(i, 1) = nu_rev * pB[1] + pHirr[1];
            h(i, 2) = nu_rev * pB[2] + pHirr[2];
        }
    }
    return H;
}

/**
 * @brief Batched state commit for the HDiv hysteresis stepping: for every row,
 *        restore states[i], play to B[i] (committed -> scratch), commit
 *        (scratch -> committed), and return the new committed state row.
 */
py::array_t<double> MatHysCommitBatch(
        int mat,
        py::array_t<double, py::array::c_style | py::array::forcecast> B,
        py::array_t<double, py::array::c_style | py::array::forcecast> states) {
    auto b = B.unchecked<2>();
    auto s = states.unchecked<2>();
    if (b.shape(1) != 3) throw std::runtime_error("MatHysCommitBatch: B must be (n, 3)");
    if (s.shape(0) != b.shape(0))
        throw std::runtime_error("MatHysCommitBatch: states rows must match B rows");
    const py::ssize_t n = b.shape(0);
    const int slen = (int)s.shape(1);
    int lenq = 0;
    if (RadMatHysSaveState(mat, nullptr, &lenq) < 0)
        throw std::runtime_error("MatHysCommitBatch: invalid material handle");
    if (lenq != slen)
        throw std::runtime_error("MatHysCommitBatch: states row length does not match the material state size");
    py::array_t<double> out({n, (py::ssize_t)slen});
    auto o = out.mutable_unchecked<2>();
    {
        py::gil_scoped_release release;
        for (py::ssize_t i = 0; i < n; ++i) {
            if (RadMatHysRestoreState(mat, s.data(i, 0), slen) < 0)
                throw std::runtime_error("MatHysCommitBatch: state restore failed");
            double pB[3] = {b(i, 0), b(i, 1), b(i, 2)};
            double pHirr[3] = {0, 0, 0};
            if (RadMatHysIrreversible(mat, pB, pHirr) < 0)
                throw std::runtime_error("MatHysCommitBatch: play evaluation failed");
            if (RadMatHysCommitState(mat) < 0)
                throw std::runtime_error("MatHysCommitBatch: commit failed");
            int len = slen;
            if (RadMatHysSaveState(mat, o.mutable_data(i, 0), &len) < 0)
                throw std::runtime_error("MatHysCommitBatch: state save failed");
        }
    }
    return out;
}

} // namespace radia_material_ext


// ============================================================================
// Additional Solver Functions
// ============================================================================

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

void SolverConfig(py::kwargs kwargs) {
    static const std::unordered_set<std::string> allowed = {
        "relax_param", "newton_method",
        "newton_damping", "newton_damping_max_iter", "newton_damping_min_omega",
        "b_input_newton", "b_input_hantila", "hantila_alpha", "hantila_relax",
        "keep_magnetization",
    };
    for (auto item : kwargs) {
        std::string key = py::cast<std::string>(item.first);
        if (allowed.find(key) == allowed.end()) {
            throw std::invalid_argument("unknown SolverConfig option: " + key);
        }
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
}

py::dict GetSolverConfig() {
    py::dict config;

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


// Replaced by Python-based PEEC topology solver (peec_topology.py)


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
	mutable std::atomic<size_t> cache_hits_;
	mutable std::atomic<size_t> cache_misses_;

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

	static void CheckRadErr(int err) {
		if (err > 0)
			throw std::runtime_error(
				"RadiaField: Radia error " + std::to_string(err) +
				" during field evaluation");
	}

	// GIL-free, job-safe field evaluation via the SERIAL direct C API.
	//
	// NGSolve assembly calls CoefficientFunction::Evaluate concurrently from
	// TaskManager worker threads, i.e. from INSIDE a running ngcore job.  Two
	// things are therefore forbidden here:
	//  1. Any Python round-trip (py::gil_scoped_acquire + rad.Fld): GIL
	//     save/restore on worker threads interleaves across the job and
	//     corrupts the interpreter state.
	//  2. The PARALLEL batch entries (RadFldBatch etc.): their internal
	//     ParallelFor issues a nested CreateJob, and ngcore job state is
	//     static -- nesting corrupts the running assembly job
	//     (0xC0000374 / 0xC0000005 heap corruption, 2026-07-10 incident).
	// The RadFld*Serial C entries evaluate the (small, per-integration-rule)
	// point loop in the calling thread; parallelism stays where it belongs,
	// in NGSolve's element loop.  Single-point RadFld is also banned: it
	// round-trips results through the non-thread-safe global ioBuffer.
	//
	// pts_local: npts*3 coordinates (already CF-local frame).
	// vals: filled with npts values for "phi", npts*3 otherwise.
	void ComputeLocalField(std::vector<double>& pts_local, size_t npts,
	                       std::vector<double>& vals) const
	{
		int n = static_cast<int>(npts);
		if (field_type == "phi") {
			vals.assign(npts, 0.0);
			CheckRadErr(RadFldPhiSerial(vals.data(), n, pts_local.data(), radia_obj));
		} else if (field_type == "a") {
			vals.assign(npts * 3, 0.0);
			CheckRadErr(RadFldASerial(vals.data(), n, pts_local.data(), radia_obj));
		} else {
			std::vector<double> B(npts * 3, 0.0), H(npts * 3, 0.0);
			CheckRadErr(RadFldBatchSerial(B.data(), H.data(), n, pts_local.data(), radia_obj));
			if (field_type == "b") vals = std::move(B);
			else if (field_type == "h") vals = std::move(H);
			else {
				// "m": M = B/mu0 - H (exact identity; B and H come from one batch call)
				const double InvMu0 = 1.0 / (4. * 3.1415926535897932 * 1.e-7);
				vals.resize(npts * 3);
				for (size_t k = 0; k < npts * 3; k++) vals[k] = B[k] * InvMu0 - H[k];
			}
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
		size_t hits = cache_hits_.load(), misses = cache_misses_.load();
		stats["enabled"] = use_cache_;
		stats["size"] = point_cache_.size();
		stats["hits"] = hits;
		stats["misses"] = misses;
		double total = static_cast<double>(hits + misses);
		stats["hit_rate"] = (total > 0) ? (hits / total) : 0.0;
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

	// Scalar evaluation for 'phi'.
	// GIL-free: direct C API only (see ComputeLocalField) -- NGSolve may call
	// this from TaskManager worker threads.
	virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
	{
		if (field_type != "phi") return 0.0;

		auto pnt = mip.GetPoint();
		int dim = pnt.Size();
		double p_global[3] = {pnt[0], (dim>=2)?pnt[1]:0.0, (dim>=3)?pnt[2]:0.0};
		std::vector<double> pts(3);
		transform_to_local(p_global, pts.data());

		std::vector<double> vals;
		ComputeLocalField(pts, 1, vals);
		return vals[0];
	}

	// Single-point vector evaluation (GIL-free, direct C API)
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

		std::vector<double> pts(3);
		transform_to_local(p_global, pts.data());

		std::vector<double> vals;
		ComputeLocalField(pts, 1, vals);

		double f_global[3];
		transform_to_global(vals.data(), f_global);
		result(0) = f_global[0]; result(1) = f_global[1]; result(2) = f_global[2];
	}

	// Batch evaluation (called by NGSolve for integration rules).
	// GIL-free: NGSolve assembly invokes this concurrently from TaskManager
	// worker threads; any Python/GIL use here corrupts the interpreter (see
	// ComputeLocalField).  Errors propagate as C++ exceptions (fail fast) --
	// no silent zero-fill.
	virtual void Evaluate(const BaseMappedIntegrationRule& mir,
	                      BareSliceMatrix<> result) const override
	{
		size_t npts = mir.Size();

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

		// Local-frame coordinates
		std::vector<double> pts(npts * 3);
		for (size_t i = 0; i < npts; i++) {
			auto pnt = mir[i].GetPoint();
			int dim = pnt.Size();
			double p_global[3] = {pnt[0], (dim>=2)?pnt[1]:0.0, (dim>=3)?pnt[2]:0.0};
			transform_to_local(p_global, pts.data() + i * 3);
		}

		std::vector<double> vals;
		ComputeLocalField(pts, npts, vals);

		if (field_type == "phi") {
			for (size_t i = 0; i < npts; i++) result(i, 0) = vals[i];
		} else {
			for (size_t i = 0; i < npts; i++) {
				double f_global[3];
				transform_to_global(vals.data() + i * 3, f_global);
				result(i, 0) = f_global[0];
				result(i, 1) = f_global[1];
				result(i, 2) = f_global[2];
			}
		}
	}
};

} // namespace ngfem


// ============================================================================
// Shared 2D planar exterior field + Maxwell torque -- rad_planar_charges
// (method-agnostic: planar solvers feed a charge cloud)
// ============================================================================
namespace radia_planar_charges {
// H at observation points from a 2D point-charge cloud.  Xq (nq,2), Q (nq,), P (nP,2) -> H (nP,2)
py::array_t<double> PlanarChargeField(
        py::array_t<double, py::array::c_style | py::array::forcecast> Xq,
        py::array_t<double, py::array::c_style | py::array::forcecast> Q,
        py::array_t<double, py::array::c_style | py::array::forcecast> P) {
    auto xb = Xq.request(); auto qb = Q.request(); auto pb = P.request();
    if (xb.ndim != 2 || xb.shape[1] != 2) throw std::runtime_error("PlanarChargeField: Xq must be (nq,2)");
    if (qb.ndim != 1 || qb.shape[0] != xb.shape[0]) throw std::runtime_error("PlanarChargeField: Q must be (nq,)");
    if (pb.ndim != 2 || pb.shape[1] != 2) throw std::runtime_error("PlanarChargeField: P must be (nP,2)");
    int nq = static_cast<int>(xb.shape[0]);
    int nP = static_cast<int>(pb.shape[0]);
    py::array_t<double> H({nP, 2});
    double* hp = static_cast<double*>(H.request().ptr);
    { py::gil_scoped_release rel;
      rad_planar_charges::Field(nq, static_cast<double*>(xb.ptr), static_cast<double*>(qb.ptr),
                                nP, static_cast<double*>(pb.ptr), hp); }
    return H;
}

// Out-of-plane vector potential A_z from a 2D point-charge cloud.  Xq (nq,2), Q (nq,), P (nP,2) -> Az (nP,)
py::array_t<double> PlanarChargeAz(
        py::array_t<double, py::array::c_style | py::array::forcecast> Xq,
        py::array_t<double, py::array::c_style | py::array::forcecast> Q,
        py::array_t<double, py::array::c_style | py::array::forcecast> P) {
    auto xb = Xq.request(); auto qb = Q.request(); auto pb = P.request();
    if (xb.ndim != 2 || xb.shape[1] != 2) throw std::runtime_error("PlanarChargeAz: Xq must be (nq,2)");
    if (qb.ndim != 1 || qb.shape[0] != xb.shape[0]) throw std::runtime_error("PlanarChargeAz: Q must be (nq,)");
    if (pb.ndim != 2 || pb.shape[1] != 2) throw std::runtime_error("PlanarChargeAz: P must be (nP,2)");
    int nq = static_cast<int>(xb.shape[0]);
    int nP = static_cast<int>(pb.shape[0]);
    py::array_t<double> Az(nP);
    double* ap = static_cast<double*>(Az.request().ptr);
    { py::gil_scoped_release rel;
      rad_planar_charges::FieldAz(nq, static_cast<double*>(xb.ptr), static_cast<double*>(qb.ptr),
                                  nP, static_cast<double*>(pb.ptr), ap); }
    return Az;
}

// Maxwell-stress torque per unit length on a circle in air (body cloud + uniform applied Hext).
double PlanarMaxwellTorqueCircle(
        py::array_t<double, py::array::c_style | py::array::forcecast> Xq,
        py::array_t<double, py::array::c_style | py::array::forcecast> Q,
        double Rc, double cx, double cy, int n, double hextx, double hexty) {
    auto xb = Xq.request(); auto qb = Q.request();
    if (xb.ndim != 2 || xb.shape[1] != 2) throw std::runtime_error("PlanarMaxwellTorqueCircle: Xq must be (nq,2)");
    if (qb.ndim != 1 || qb.shape[0] != xb.shape[0]) throw std::runtime_error("PlanarMaxwellTorqueCircle: Q must be (nq,)");
    int nq = static_cast<int>(xb.shape[0]);
    double T;
    { py::gil_scoped_release rel;
      T = rad_planar_charges::MaxwellTorqueCircle(nq, static_cast<double*>(xb.ptr),
                                                  static_cast<double*>(qb.ptr), Rc, cx, cy, n, hextx, hexty); }
    return T;
}

// Maxwell-stress FORCE per unit length on a circle in air (body cloud + uniform applied Hext) -> (Fx,Fy).
py::array_t<double> PlanarMaxwellForceCircle(
        py::array_t<double, py::array::c_style | py::array::forcecast> Xq,
        py::array_t<double, py::array::c_style | py::array::forcecast> Q,
        double Rc, double cx, double cy, int n, double hextx, double hexty) {
    auto xb = Xq.request(); auto qb = Q.request();
    if (xb.ndim != 2 || xb.shape[1] != 2) throw std::runtime_error("PlanarMaxwellForceCircle: Xq must be (nq,2)");
    if (qb.ndim != 1 || qb.shape[0] != xb.shape[0]) throw std::runtime_error("PlanarMaxwellForceCircle: Q must be (nq,)");
    int nq = static_cast<int>(xb.shape[0]);
    py::array_t<double> F(2);
    double* fp = static_cast<double*>(F.request().ptr);
    { py::gil_scoped_release rel;
      rad_planar_charges::MaxwellForceCircle(nq, static_cast<double*>(xb.ptr), static_cast<double*>(qb.ptr),
                                             Rc, cx, cy, n, hextx, hexty, fp); }
    return F;
}
} // namespace radia_planar_charges

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
            magnet = rad.magnet_box([0,0,0], [0.04, 0.04, 0.02], [0, 0, 954930])

            # Compute field
            B = rad.Fld(magnet, 'b', [0.05, 0, 0])
    )pbdoc";

    // Version info
    m.attr("__version__") = "1.4.0";

    using FieldEvaluator = rad_hdiv::HDivFieldEvaluator;
    py::class_<FieldEvaluator, std::shared_ptr<FieldEvaluator>>(m, "_HDivFieldEvaluator")
        .def_static("from_tet", [](
                py::array_t<double, py::array::c_style | py::array::forcecast> volume,
                py::array_t<double, py::array::c_style | py::array::forcecast> surface,
                py::array_t<int, py::array::c_style | py::array::forcecast> image_masks,
                py::array_t<double, py::array::c_style | py::array::forcecast> image_signs,
                int leaf_size, double theta, std::size_t tree_min_sources,
                std::size_t auto_min_work, double tree_relative_tolerance,
                int probe_count) {
            auto v = volume.request(), s = surface.request();
            if (v.ndim != 2 || v.shape[1] != 16)
                throw std::runtime_error("_HDivFieldEvaluator.from_tet: volume must have shape (n,16)");
            if (s.ndim != 2 || s.shape[1] != 22)
                throw std::runtime_error("_HDivFieldEvaluator.from_tet: surface must have shape (n,22)");
            const double* vp = static_cast<const double*>(v.ptr);
            const double* sp = static_cast<const double*>(s.ptr);
            rad_hdiv::FieldEvaluatorOptions options;
            options.leaf_size = leaf_size; options.theta = theta;
            options.tree_min_sources = tree_min_sources; options.auto_min_work = auto_min_work;
            options.tree_relative_tolerance = tree_relative_tolerance; options.probe_count = probe_count;
            return FieldEvaluator::FromTet(
                std::vector<double>(vp, vp+v.size), std::vector<double>(sp, sp+s.size),
                to_1d_vector<int>(image_masks, "image_masks"),
                to_1d_vector<double>(image_signs, "image_signs"), options);
        }, py::arg("volume"), py::arg("surface"),
           py::arg("image_masks") = py::array_t<int>(0),
           py::arg("image_signs") = py::array_t<double>(0),
           py::arg("leaf_size") = 32, py::arg("theta") = 0.05,
           py::arg("tree_min_sources") = 256, py::arg("auto_min_work") = 500000000,
           py::arg("tree_relative_tolerance") = 1.0e-5, py::arg("probe_count") = 16,
           "Materialize an immutable RT1 tet/triangle source evaluator. Source arrays are copied once.")
        .def_static("from_cloud", [](
                py::array_t<double, py::array::c_style | py::array::forcecast> xyz,
                py::array_t<double, py::array::c_style | py::array::forcecast> strength,
                py::array_t<int, py::array::c_style | py::array::forcecast> image_masks,
                py::array_t<double, py::array::c_style | py::array::forcecast> image_signs,
                int leaf_size, double theta, std::size_t tree_min_sources,
                std::size_t auto_min_work, double tree_relative_tolerance,
                int probe_count) {
            auto x = xyz.request(), q = strength.request();
            if (x.ndim != 2 || x.shape[1] != 3)
                throw std::runtime_error("_HDivFieldEvaluator.from_cloud: xyz must have shape (n,3)");
            if (q.ndim != 1 || q.shape[0] != x.shape[0])
                throw std::runtime_error("_HDivFieldEvaluator.from_cloud: strength must have shape (n,)");
            const double* xp = static_cast<const double*>(x.ptr);
            const double* qp = static_cast<const double*>(q.ptr);
            rad_hdiv::FieldEvaluatorOptions options;
            options.leaf_size = leaf_size; options.theta = theta;
            options.tree_min_sources = tree_min_sources; options.auto_min_work = auto_min_work;
            options.tree_relative_tolerance = tree_relative_tolerance; options.probe_count = probe_count;
            return FieldEvaluator::FromCloud(
                std::vector<double>(xp, xp+x.size), std::vector<double>(qp, qp+q.size),
                to_1d_vector<int>(image_masks, "image_masks"),
                to_1d_vector<double>(image_signs, "image_signs"), options);
        }, py::arg("xyz"), py::arg("strength"),
           py::arg("image_masks") = py::array_t<int>(0),
           py::arg("image_signs") = py::array_t<double>(0),
           py::arg("leaf_size") = 32, py::arg("theta") = 0.05,
           py::arg("tree_min_sources") = 256, py::arg("auto_min_work") = 500000000,
           py::arg("tree_relative_tolerance") = 1.0e-5, py::arg("probe_count") = 16,
           "Materialize an immutable RT1 quadrature-cloud evaluator. Source arrays are copied once.")
        .def("field", [](const FieldEvaluator& evaluator,
                py::array_t<double, py::array::c_style | py::array::forcecast> observations,
                const std::string& algorithm) {
            auto input = observations.request();
            if (input.ndim != 2 || input.shape[1] != 3)
                throw std::runtime_error("_HDivFieldEvaluator.field: observations must have shape (n,3)");
            const std::size_t count = static_cast<std::size_t>(input.shape[0]);
            py::array_t<double> output({static_cast<py::ssize_t>(count), py::ssize_t(3)});
            auto output_buffer = output.request();
            const double* in = static_cast<const double*>(input.ptr);
            double* out = static_cast<double*>(output_buffer.ptr);
            const auto selected = FieldEvaluator::ParseAlgorithm(algorithm);
            {
                py::gil_scoped_release release;
                evaluator.Evaluate(in, count, out, selected);
            }
            return output;
        }, py::arg("observations"), py::arg("algorithm") = "auto",
           "Evaluate all physical and IMA sources in one TaskManager-parallel call; returns NO 1/(4pi).")
        .def("candidate_algorithm_for", [](const FieldEvaluator& evaluator, std::size_t n_observations) {
            return std::string(FieldEvaluator::AlgorithmName(evaluator.AlgorithmFor(n_observations)));
        }, "Return the work-threshold candidate; auto still probes accuracy and measured cost.")
        .def("last_algorithm", [](const FieldEvaluator& evaluator) {
            return std::string(FieldEvaluator::AlgorithmName(evaluator.LastAlgorithm()));
        }, "Return the algorithm selected by the most recent field call.")
        .def("stats", [](const FieldEvaluator& evaluator) {
            double lower[3], upper[3]; evaluator.Bounds(lower, upper);
            py::dict result;
            result["source_count"] = evaluator.SourceCount();
            result["image_count"] = evaluator.ImageCount();
            result["tree_nodes"] = evaluator.TreeNodeCount();
            result["leaf_size"] = evaluator.LeafSize();
            result["theta"] = evaluator.Theta();
            result["tree_min_sources"] = evaluator.TreeMinSources();
            result["auto_min_work"] = evaluator.AutoMinWork();
            result["tree_relative_tolerance"] = evaluator.TreeRelativeTolerance();
            result["probe_count"] = evaluator.ProbeCount();
            result["bounds_min"] = py::make_tuple(lower[0], lower[1], lower[2]);
            result["bounds_max"] = py::make_tuple(upper[0], upper[1], upper[2]);
            return result;
        });


    // Charge-charge Coulomb Gram G as a HACApK H-matrix -- the UNSTRUCTURED / general-mesh path.
    // Charges (cell rho + boundary-face sigma) extracted from an HDiv mesh; pass charge
    // centroids/measures + the caller-computed diagonal self-energies.  The
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
                         int n_el, double eps, int leaf, double eta, double near_factor,
                         std::vector<int> image_masks, std::vector<double> image_signs, bool build,
                         int far_quad) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(
                     new RadHACApKChargeGram(std::move(cell_verts), std::move(face_verts), n_el, near_factor,
                                             std::move(image_masks), std::move(image_signs), far_quad));
                 if (build) {
                     RadHACApKParams p;
                     p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                     if (!mgr->BuildHMatrix(p)) throw std::runtime_error("analytic charge Gram H-matrix build failed");
                 }
                 // build=False leaves the H-matrix UNBUILT: only the geometry (outer quadrature, centroids,
                 // sizes) is set up in the ctor, so .entry() works as a cheap analytic ENTRY ORACLE while
                 // .matvec() is unavailable (would raise).  The Gauss near-correction uses this to sample
                 // exact analytic near entries WITHOUT paying the full O(N log N) analytic H-matrix build.
                 return mgr;
             }),
             py::arg("cell_verts"), py::arg("face_verts"), py::arg("n_el"),
             py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0, py::arg("near_factor") = 1e30,
             py::arg("image_masks") = std::vector<int>{}, py::arg("image_signs") = std::vector<double>{},
             py::arg("build") = true, py::arg("far_quad") = 0,
             "ANALYTIC mode (M2b): build the EXACT charge Gram as a HACApK H-matrix from per-charge "
             "geometry (cell_verts [n_el*12] tets, face_verts [n_bf*9] triangles). Entry = analytic "
             "PhiTet/TriPotential inner x outer quadrature (matches the independent analytic reference). "
             "near_factor (default 1e30 = all-analytic) gives the NEAR/FAR build speedup. far_quad selects the "
             "FAR evaluation: 0 = centroid-monopole (O((size/r)^2), slightly breaks symmetry); >0 = a low-order "
             "double-quadrature of 1/r (degree-2 4-pt tet / 3-pt tri, O((size/r)^4)) that reproduces the "
             "all-analytic Gram at ~monopole cost -- the precision-preserving build speedup (use near_factor~2 "
             "+ far_quad=4). build=False skips BuildHMatrix -> a geometry-only ENTRY ORACLE (.entry() only).")
        .def(py::init([](std::vector<double> cell_tris, std::vector<int> cell_troff,
                         std::vector<double> cell_cent, std::vector<double> cell_meas,
                         std::vector<double> face_tris, std::vector<int> face_troff,
                         std::vector<double> face_cent, std::vector<double> face_meas,
                         int n_el, double eps, int leaf, double eta, double near_factor,
                         std::vector<int> image_masks, std::vector<double> image_signs, int far_quad) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(
                     new RadHACApKChargeGram(std::move(cell_tris), std::move(cell_troff),
                                             std::move(cell_cent), std::move(cell_meas),
                                             std::move(face_tris), std::move(face_troff),
                                             std::move(face_cent), std::move(face_meas), n_el, near_factor,
                                             std::move(image_masks), std::move(image_signs), far_quad));
                 RadHACApKParams p;
                 p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                 if (!mgr->BuildHMatrix(p)) throw std::runtime_error("polytope charge Gram H-matrix build failed");
                 return mgr;
             }),
             py::arg("cell_tris"), py::arg("cell_troff"), py::arg("cell_cent"), py::arg("cell_meas"),
             py::arg("face_tris"), py::arg("face_troff"), py::arg("face_cent"), py::arg("face_meas"),
             py::arg("n_el"), py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0,
             py::arg("near_factor") = 1e30,
             py::arg("image_masks") = std::vector<int>{}, py::arg("image_signs") = std::vector<double>{},
             py::arg("far_quad") = 0,
             "POLYTOPE mode (hex/wedge cells + quad faces): the EXACT analytic charge Gram for any "
             "flat-faced convex cell, the triangulation supplied from Python.  cell_tris is a flat "
             "triangle soup (9 doubles/tri) of all cells' convex-hull triangles, cell_troff [n_el+1] the "
             "CSR offsets (in triangles), cell_cent [n_el*3] the vertex-mean centroid (fan apex + outward "
             "normal ref), cell_meas [n_el] the cell volume; face_tris/face_troff/face_cent/face_meas the "
             "boundary faces' sub-triangles (quad->2).  Entry = divergence-theorem polytope potential x "
             "centroid-fan / Dunavant outer quadrature (matches the independent analytic polytope reference). "
             "near_factor (default 1e30 = all-analytic); pass ~2 for the NEAR/FAR build speedup. far_quad>0 "
             "uses the precision-preserving low-order double-quad far (degree-2 on the sub-tets/sub-tris).")
        .def(py::init([](std::vector<double> cell_curved_nodes, std::vector<int> cell_subtet_off,
                         std::vector<double> cell_cent, std::vector<double> cell_meas,
                         std::vector<double> face_curved_nodes, std::vector<int> face_subtri_off,
                         std::vector<double> face_cent, std::vector<double> face_meas,
                         std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
                         std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
                         std::vector<double> curve_gl, std::vector<double> curve_gw, int n_el,
                         double eps, int leaf, double eta) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(
                     new RadHACApKChargeGram(std::move(cell_curved_nodes), std::move(cell_subtet_off),
                                             std::move(cell_cent), std::move(cell_meas),
                                             std::move(face_curved_nodes), std::move(face_subtri_off),
                                             std::move(face_cent), std::move(face_meas),
                                             std::move(ref_tet_pts), std::move(ref_tet_w),
                                             std::move(ref_tri_pts), std::move(ref_tri_w),
                                             std::move(curve_gl), std::move(curve_gw), n_el));
                 RadHACApKParams p;
                 p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                 if (!mgr->BuildHMatrix(p)) throw std::runtime_error("curved polytope charge Gram H-matrix build failed");
                 return mgr;
             }),
             py::arg("cell_curved_nodes"), py::arg("cell_subtet_off"), py::arg("cell_cent"), py::arg("cell_meas"),
             py::arg("face_curved_nodes"), py::arg("face_subtri_off"), py::arg("face_cent"), py::arg("face_meas"),
             py::arg("ref_tet_pts"), py::arg("ref_tet_w"), py::arg("ref_tri_pts"), py::arg("ref_tri_w"),
             py::arg("curve_gl"), py::arg("curve_gw"), py::arg("n_el"),
             py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0,
             "CURVED POLYTOPE mode (curved hex/wedge): FULLY curved -- curved CELL volume charge (sub-tets, "
             "CurvedTetPotential) + curved FACE surface charge (sub-tris, CurvedTriPotential). The cell volume "
             "charge is DOMINANT for the demag (the lowest-order curved charge cannot represent uniform M exactly -> div M != 0). "
             "cell_curved_nodes [n_cell_subtet*30] + cell_subtet_off [n_cell+1]; face_curved_nodes "
             "[n_bf_subtri*18] + face_subtri_off [n_bf+1]; ref_tet_pts/w + ref_tri_pts/w the outer quads; "
             "curve_gl/gw the inner Duffy rule. Both reuse the golden curved-tet/tri kernels.")
        .def(py::init([](std::vector<double> cell_verts, std::vector<double> face_verts, int n_el,
                         std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
                         std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
                         std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
                         std::vector<double> ref_tet_pts_lo, std::vector<double> ref_tet_w_lo,
                         std::vector<double> ref_tri_pts_lo, std::vector<double> ref_tri_w_lo,
                         double ho_far_factor,
                         std::vector<double> ref_tet_pts_in, std::vector<double> ref_tet_w_in,
                         std::vector<double> ref_tri_pts_in, std::vector<double> ref_tri_w_in,
                         std::vector<int> image_masks, std::vector<double> image_signs,
                         double eps, int leaf, double eta, bool build) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(
                     new RadHACApKChargeGram(std::move(cell_verts), std::move(face_verts), n_el,
                                             std::move(charge_host), std::move(charge_kind),
                                             std::move(charge_expo), std::move(ref_tet_pts),
                                             std::move(ref_tet_w), std::move(ref_tri_pts), std::move(ref_tri_w),
                                             std::move(ref_tet_pts_lo), std::move(ref_tet_w_lo),
                                             std::move(ref_tri_pts_lo), std::move(ref_tri_w_lo), ho_far_factor,
                                             std::move(ref_tet_pts_in), std::move(ref_tet_w_in),
                                             std::move(ref_tri_pts_in), std::move(ref_tri_w_in),
                                             std::move(image_masks), std::move(image_signs)));
                 if (build) {
                     RadHACApKParams p;
                     p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                     if (!mgr->BuildHMatrix(p)) throw std::runtime_error("high-order charge Gram H-matrix build failed");
                 }
                 // build=False -> geometry-only ENTRY ORACLE (.entry() only): the high-order Gauss near
                 // correction samples exact analytic high-order entries WITHOUT the full H-matrix build.
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
             py::arg("image_masks") = std::vector<int>{}, py::arg("image_signs") = std::vector<double>{},
             py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0, py::arg("build") = true,
             "HIGH-ORDER (order-p) mode: POLYNOMIAL charges (monomial basis per host). charge_host[c]/"
             "charge_kind[c] (0=cell,1=face)/charge_expo[3c] define each charge; ref_tet_pts[nqt*3]/ref_tet_w "
             "(sum 1/6) + ref_tri_pts[nqr*2]/ref_tri_w (sum 1/2) are the reference Gauss-Duffy rules. Entry = "
             "monomial-weighted outer quad x the subtraction inner potential (matches the independent high-order analytic reference). "
             "ref_*_lo + ho_far_factor (<inf) enable the accuracy-preserving NEAR/FAR adaptive quadrature: far "
             "pairs (|c_a-c_b| > ho_far_factor*(size_a+size_b)) use the cheap LOW-quad plain double-Gauss.")
        .def(py::init([](std::vector<double> cell_nodes, std::vector<double> face_nodes, int n_el, int curve_order,
                         std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
                         std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
                         std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
                         std::vector<double> curve_gl, std::vector<double> curve_gw,
                         double eps, int leaf, double eta, bool build) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(
                     new RadHACApKChargeGram(std::move(cell_nodes), std::move(face_nodes), n_el, curve_order,
                                             std::move(charge_host), std::move(charge_kind), std::move(charge_expo),
                                             std::move(ref_tet_pts), std::move(ref_tet_w),
                                             std::move(ref_tri_pts), std::move(ref_tri_w),
                                             std::move(curve_gl), std::move(curve_gw)));
                 if (build) {
                     RadHACApKParams p;
                     p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                     if (!mgr->BuildHMatrix(p)) throw std::runtime_error("curved charge Gram H-matrix build failed");
                 }
                 return mgr;
             }),
             py::arg("cell_nodes"), py::arg("face_nodes"), py::arg("n_el"), py::arg("curve_order"),
             py::arg("charge_host"), py::arg("charge_kind"), py::arg("charge_expo"),
             py::arg("ref_tet_pts"), py::arg("ref_tet_w"), py::arg("ref_tri_pts"), py::arg("ref_tri_w"),
             py::arg("curve_gl"), py::arg("curve_gw"),
             py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0, py::arg("build") = true,
             "CURVED HIGH-ORDER (isoparametric P2) mode: monomial charges on a mesh.Curve(2) geometry. "
             "cell_nodes [n_el*30] (10 P2 nodes/tet), face_nodes [n_bf*18] (6 P2 nodes/tri); curve_order=2. "
             "Outer quad = curved P2 map + curved measure; inner = the curved Duffy (curve_gl/gw = nq-pt "
             "Gauss-Legendre on [0,1]). Curved helps near-surface FIELD/flux accuracy, NOT the demag factor.")
        .def(py::init([](std::vector<double> hex_cell_nodes, std::vector<double> quad_face_nodes,
                         int n_el, int n_bf,
                         std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
                         std::vector<double> sym_tet_pts, std::vector<double> sym_tet_w,
                         std::vector<double> sym_tri_pts, std::vector<double> sym_tri_w,
                         std::vector<double> gl_out, std::vector<double> gw_out,
                         std::vector<double> gl_in, std::vector<double> gw_in,
                         std::vector<double> far_tet_pts, std::vector<double> far_tet_w,
                         std::vector<double> far_tri_pts, std::vector<double> far_tri_w,
                         double near_grade, double far_inner_factor,
                         std::vector<int> image_masks, std::vector<double> image_signs,
                         double eps, int leaf, double eta, bool build) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(
                     new RadHACApKChargeGram(std::move(hex_cell_nodes), std::move(quad_face_nodes), n_el, n_bf,
                                             std::move(charge_host), std::move(charge_kind), std::move(charge_expo),
                                             std::move(sym_tet_pts), std::move(sym_tet_w),
                                             std::move(sym_tri_pts), std::move(sym_tri_w),
                                             std::move(gl_out), std::move(gw_out),
                                             std::move(gl_in), std::move(gw_in),
                                             std::move(far_tet_pts), std::move(far_tet_w),
                                             std::move(far_tri_pts), std::move(far_tri_w),
                                             near_grade, far_inner_factor,
                                             std::move(image_masks), std::move(image_signs)));
                 if (build) {
                     RadHACApKParams p;
                     p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                     if (!mgr->BuildHMatrix(p)) throw std::runtime_error("hex RT1 charge Gram H-matrix build failed");
                 }
                 return mgr;
             }),
             py::arg("hex_cell_nodes"), py::arg("quad_face_nodes"), py::arg("n_el"), py::arg("n_bf"),
             py::arg("charge_host"), py::arg("charge_kind"), py::arg("charge_expo"),
             py::arg("sym_tet_pts"), py::arg("sym_tet_w"), py::arg("sym_tri_pts"), py::arg("sym_tri_w"),
             py::arg("gl_out"), py::arg("gw_out"), py::arg("gl_in"), py::arg("gw_in"),
             py::arg("far_tet_pts"), py::arg("far_tet_w"), py::arg("far_tri_pts"), py::arg("far_tri_w"),
             py::arg("near_grade") = 1.5, py::arg("far_inner_factor") = 1.5,
             py::arg("image_masks") = std::vector<int>{}, py::arg("image_signs") = std::vector<double>{},
             py::arg("eps") = 1e-4, py::arg("leaf") = 32, py::arg("eta") = 2.0, py::arg("build") = true,
             "HEX RT1 mode: Q1 monomial charges (8/hex volume + 4/quad-face surface) on the DIRECT Q2 "
             "isoparametric geometry -- hex_cell_nodes [n_el*81] = 27-node lattice, quad_face_nodes "
             "[n_bf*27] = 9-node lattice, both from GetTrafo at the reference lattice, so ONE path covers "
             "flat AND curved (mesh.Curve(2)) hexes.  Quadrature = the numpy-validated eig<=1 scheme: "
             "near sub pairs -> both-domains-graded Duffy (gl_out/gl_in 1D rules); far -> the regular "
             "symmetric rules (sym_* = Keast-15/Dunavant-7) + cheap far inner (far_*); the radial "
             "near/self inner fires only within far_inner_factor*size of a source sub (per outer point).  "
             "The H-matrix build is SYMMETRIC-FILL: strictly-lower leaves are skipped (all applies route "
             "through matvec_sym; plain matvec/matvec_transpose are routed to it).")
        .def(py::init([](std::vector<double> wedge_cell_nodes, std::vector<double> face_nodes,
                         std::vector<int> face_type, int n_el, int n_bf,
                         std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
                         std::vector<double> sym_tet_pts, std::vector<double> sym_tet_w,
                         std::vector<double> sym_tri_pts, std::vector<double> sym_tri_w,
                         std::vector<double> gl_out, std::vector<double> gw_out,
                         std::vector<double> gl_in, std::vector<double> gw_in,
                         std::vector<double> far_tet_pts, std::vector<double> far_tet_w,
                         std::vector<double> far_tri_pts, std::vector<double> far_tri_w,
                         double near_grade, double far_inner_factor,
                         std::vector<int> image_masks, std::vector<double> image_signs,
                         double eps, int leaf, double eta, bool build) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(
                     new RadHACApKChargeGram(std::move(wedge_cell_nodes), std::move(face_nodes),
                                             std::move(face_type), n_el, n_bf,
                                             std::move(charge_host), std::move(charge_kind), std::move(charge_expo),
                                             std::move(sym_tet_pts), std::move(sym_tet_w),
                                             std::move(sym_tri_pts), std::move(sym_tri_w),
                                             std::move(gl_out), std::move(gw_out),
                                             std::move(gl_in), std::move(gw_in),
                                             std::move(far_tet_pts), std::move(far_tet_w),
                                             std::move(far_tri_pts), std::move(far_tri_w),
                                             near_grade, far_inner_factor,
                                             std::move(image_masks), std::move(image_signs)));
                 if (build) {
                     RadHACApKParams p;
                     p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                     if (!mgr->BuildHMatrix(p)) throw std::runtime_error("wedge RT1 charge Gram H-matrix build failed");
                 }
                 return mgr;
             }),
             py::arg("wedge_cell_nodes"), py::arg("face_nodes"), py::arg("face_type"),
             py::arg("n_el"), py::arg("n_bf"),
             py::arg("charge_host"), py::arg("charge_kind"), py::arg("charge_expo"),
             py::arg("sym_tet_pts"), py::arg("sym_tet_w"), py::arg("sym_tri_pts"), py::arg("sym_tri_w"),
             py::arg("gl_out"), py::arg("gw_out"), py::arg("gl_in"), py::arg("gw_in"),
             py::arg("far_tet_pts"), py::arg("far_tet_w"), py::arg("far_tri_pts"), py::arg("far_tri_w"),
             py::arg("near_grade") = 0.6, py::arg("far_inner_factor") = 1.5,
             py::arg("image_masks") = std::vector<int>{}, py::arg("image_signs") = std::vector<double>{},
             py::arg("eps") = 1e-12, py::arg("leaf") = 64, py::arg("eta") = 2.0, py::arg("build") = true,
             "WEDGE (PRISM) RT1 mode: L2(prism,order=1) volume charges (6/prism = {1,x,y,z,xz,yz}, the "
             "prism div-image = a subset of the hex's 8 Q1 monomials) + SurfaceL2 face charges (tri P1 3, "
             "quad Q1 4) on the direct 18-node tri-P2(x)z-P2 (wedge_cell_nodes [n_el*54]) + mixed-face "
             "(face_nodes [n_bf*27] 9-node slots, a tri fills the first 6; face_type [n_bf] 0=tri/1=quad) "
             "isoparametric geometry via GetTrafo -> ONE path flat + curved.  Quadrature = the hex-mode "
             "both-domains-graded Duffy scheme on the 3-sub-tet / mixed 1-2 sub-tri decomposition "
             "(numpy de-risk eig(M^-1 N) in [0,1]: 0.989 @ n=2, 0.997 @ n=3; demag_z ~ 1/3).  Shares the "
             "hex block memo / symmetric-fill build; the golden hex path is byte-for-byte untouched.")
        .def(py::init([](int dim2, std::vector<double> cell_nodes9, std::vector<int> cell_type,
                         std::vector<double> edge_nodes3, int n_el, int n_be,
                         std::vector<int> charge_host, std::vector<int> charge_kind,
                         std::vector<int> charge_expo,
                         std::vector<double> sym_tri_pts, std::vector<double> sym_tri_w,
                         std::vector<double> gl_edge, std::vector<double> gw_edge,
                         std::vector<double> gl_in, std::vector<double> gw_in,
                         std::vector<double> far_tri_pts, std::vector<double> far_tri_w,
                         double near_grade, double far_inner_factor,
                         double eps, int leaf, double eta, bool build) {
                 auto mgr = std::unique_ptr<RadHACApKChargeGram>(new RadHACApKChargeGram(
                     dim2, std::move(cell_nodes9), std::move(cell_type), std::move(edge_nodes3),
                     n_el, n_be,
                     std::move(charge_host), std::move(charge_kind), std::move(charge_expo),
                     std::move(sym_tri_pts), std::move(sym_tri_w),
                     std::move(gl_edge), std::move(gw_edge), std::move(gl_in), std::move(gw_in),
                     std::move(far_tri_pts), std::move(far_tri_w),
                     near_grade, far_inner_factor));
                 if (build) {
                     RadHACApKParams p;
                     p.aca_eps = eps; p.leaf_size = leaf; p.eta = eta; p.print_level = 0;
                     if (!mgr->BuildHMatrix(p)) throw std::runtime_error("2D charge Gram H-matrix build failed");
                 }
                 return mgr;
             }),
             py::arg("dim2"), py::arg("cell_nodes9"), py::arg("cell_type"), py::arg("edge_nodes3"),
             py::arg("n_el"), py::arg("n_be"),
             py::arg("charge_host"), py::arg("charge_kind"), py::arg("charge_expo"),
             py::arg("sym_tri_pts"), py::arg("sym_tri_w"),
             py::arg("gl_edge"), py::arg("gw_edge"), py::arg("gl_in"), py::arg("gw_in"),
             py::arg("far_tri_pts"), py::arg("far_tri_w"),
             py::arg("near_grade") = 0.6, py::arg("far_inner_factor") = 1.5,
             py::arg("eps") = 1e-12, py::arg("leaf") = 64, py::arg("eta") = 2.0, py::arg("build") = true,
             "2D PLANAR mode (motor cross-sections; memory hdiv-vim-tri-quad-motor): charges rho = -div M "
             "on tri/quad cells (tri P0, quad Q1 -- the 2D hex-gotcha twin) + sigma = M.n on boundary "
             "EDGES (P1), Piola-exact REF measures, kernel -ln(r)/(2pi) (the ln-scale shift is killed by "
             "the zero-total-charge dof columns).  Geometry = P2 lattices via GetTrafo (tri 6-node in "
             "9-node slots, quad 9-node, edge 3-node; cell_type 0=tri 1=quad) -> flat + curved one path.  "
             "Regular symmetric outer everywhere (the log kernel needs NO graded outer -- numpy-validated); "
             "radial-cone inner for near/self, cheap far cloud otherwise.  Gates: eig(M^-1 N) in [0,1]; "
             "disk demag == 1/2 exact; ellipse a:b -> b/(a+b); 2D Clausius-Mossotti chi/(1+chi/2).")
        .def("ndof", [](RadHACApKChargeGram& s) { return s.GetNDOF(); })
        .def("matvec", [](RadHACApKChargeGram& s, const std::vector<double>& x) {
                 std::vector<double> y((size_t)s.GetNDOF(), 0.0);
                 s.MatVec(x, y);
                 return y;
             }, py::arg("x"), "G q (the O(N log N) Gram H-matvec).")
        .def("matvec_transpose", [](RadHACApKChargeGram& s, const std::vector<double>& x) {
                 std::vector<double> y((size_t)s.GetNDOF(), 0.0);
                 s.MatVecTranspose(x, y);
                 return y;
             }, py::arg("x"), "G^T q (transpose H-matvec; for symmetry probes / 0.5*(G+G^T) apply).")
        .def("matvec_sym", [](RadHACApKChargeGram& s, const std::vector<double>& x) {
                 std::vector<double> y((size_t)s.GetNDOF(), 0.0);
                 s.MatVecSym(x, y);
                 return y;
             }, py::arg("x"),
             "G_sym q -- EXACTLY symmetric H-matvec (upper-triangular leaves define both triangles), "
             "so CG/MINRES on B^T G_sym B use a machine-symmetric operator (the ACA-asymmetry failure mode is removed).")
        .def("entry", &RadHACApKChargeGram::GetInteractionMatrixElement, py::arg("i"), py::arg("j"),
             "Charge-Gram entry G[i,j] from the analytic / polytope / high-order kernel.")
        .def("hex_state_check", [](RadHACApKChargeGram& s) {
                 py::dict d;
                 d["ctor"] = s.HexStateCtorChecksum();
                 d["now"] = s.HexStateChecksum();
                 return d;
             },
             "Heap-stomp canary (hex mode): ctor-time vs recomputed checksum of every member array the "
             "block computation reads.  ctor != now proves the instance data was overwritten after "
             "construction (0xc0000374-class corruption), not computed wrong.")
        .def("hex_stored_nodes", [](RadHACApKChargeGram& s) {
                 py::dict d;
                 d["cell_nodes"] = s.HexStoredCellNodes();
                 d["face_nodes"] = s.HexStoredFaceNodes();
                 return d;
             },
             "The Q2 lattice node arrays THIS instance was constructed from (flake forensics: compare "
             "two instances' inputs directly).")
        .def("hex_state_breakdown", [](RadHACApKChargeGram& s) {
                 py::dict d;
                 for (const auto& kv : s.HexStateBreakdown()) d[kv.first.c_str()] = kv.second;
                 return d;
             },
             "Per-array checksum breakdown (flake forensics: WHICH array differs between instances).")
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
        .def("solve_linear_material_auto_prec",
             [](RadHACApKChargeGram& s,
                std::vector<int> B_indptr, std::vector<int> B_indices, std::vector<double> B_data,
                int n_face, std::vector<int> mI, std::vector<int> mJ, std::vector<double> mV,
                double inv_chi, std::vector<double> rhs,
                double tol, int maxit) {
                 if ((int)rhs.size() != n_face)
                     throw std::runtime_error("solve_linear_material_auto_prec: rhs size mismatch");
                 if ((int)B_indptr.size() < 1)
                     throw std::runtime_error("solve_linear_material_auto_prec: empty B_indptr");
                 std::vector<double> mass_diag((size_t)n_face, 0.0);
                 for (size_t k = 0; k < mV.size(); ++k) {
                     if (mI[k] == mJ[k] && mI[k] >= 0 && mI[k] < n_face) mass_diag[(size_t)mI[k]] += mV[k];
                 }
                 std::vector<std::vector<int>> supp_id((size_t)n_face);
                 std::vector<std::vector<double>> supp_val((size_t)n_face);
                 const int n_charge = (int)B_indptr.size() - 1;
                 for (int a = 0; a < n_charge; ++a) {
                     for (int k = B_indptr[(size_t)a]; k < B_indptr[(size_t)a + 1]; ++k) {
                         int f = B_indices[(size_t)k];
                         if (f < 0 || f >= n_face) throw std::runtime_error("solve_linear_material_auto_prec: B face index out of range");
                         supp_id[(size_t)f].push_back(a);
                         supp_val[(size_t)f].push_back(B_data[(size_t)k]);
                     }
                 }
                 std::vector<double> prec((size_t)n_face, 0.0);
                 for (int f = 0; f < n_face; ++f) {
                     double ndiag = 0.0;
                     const auto& ids = supp_id[(size_t)f];
                     const auto& vals = supp_val[(size_t)f];
                     for (size_t p = 0; p < ids.size(); ++p)
                         for (size_t q = 0; q < ids.size(); ++q)
                             ndiag += vals[p] * vals[q] * s.GetInteractionMatrixElement(ids[p], ids[q]);
                     double v = inv_chi * mass_diag[(size_t)f] + ndiag;
                     if (!(v > 0.0) || !std::isfinite(v)) v = 1.0;
                     prec[(size_t)f] = v;
                 }
                 int iters = 0;
                 std::vector<double> m = s.SolveLinearMaterial(B_indptr, B_indices, B_data, n_face,
                                                               mI, mJ, mV, inv_chi, prec, rhs,
                                                               tol, maxit, iters);
                 double pmin = n_face ? prec[0] : 0.0, pmax = pmin;
                 for (double v : prec) { if (v < pmin) pmin = v; if (v > pmax) pmax = v; }
                 py::dict timings;
                 for (const auto& kv : s.LastSolveTimings()) timings[py::str(kv.first)] = kv.second;
                 py::dict d;
                 d["m"] = m; d["iters"] = iters; d["prec_min"] = pmin; d["prec_max"] = pmax;
                 d["timings"] = timings;
                 return d;
             },
             py::arg("B_indptr"), py::arg("B_indices"), py::arg("B_data"), py::arg("n_face"),
             py::arg("mI"), py::arg("mJ"), py::arg("mV"), py::arg("inv_chi"), py::arg("rhs"),
             py::arg("tol") = 1e-9, py::arg("maxit") = 5000,
             "M3 production helper: build the exact Jacobi diagonal of ((1/chi)M_mass + B^T G B) in C++ "
             "from sparse B + mass COO, then run SolveLinearMaterial. Returns {m, iters, prec_min, prec_max}.")
        .def("solve_linear_material_auto_prec_arrays",
             [](RadHACApKChargeGram& s,
                py::array_t<int, py::array::c_style | py::array::forcecast> B_indptr_a,
                py::array_t<int, py::array::c_style | py::array::forcecast> B_indices_a,
                py::array_t<double, py::array::c_style | py::array::forcecast> B_data_a,
                int n_face,
                py::array_t<int, py::array::c_style | py::array::forcecast> mI_a,
                py::array_t<int, py::array::c_style | py::array::forcecast> mJ_a,
                py::array_t<double, py::array::c_style | py::array::forcecast> mV_a,
                double inv_chi,
                py::array_t<double, py::array::c_style | py::array::forcecast> rhs_a,
                double tol, int maxit, py::object x0_obj) {
                 auto B_indptr = to_1d_vector<int>(B_indptr_a, "B_indptr");
                 auto B_indices = to_1d_vector<int>(B_indices_a, "B_indices");
                 auto B_data = to_1d_vector<double>(B_data_a, "B_data");
                 auto mI = to_1d_vector<int>(mI_a, "mI");
                 auto mJ = to_1d_vector<int>(mJ_a, "mJ");
                 auto mV = to_1d_vector<double>(mV_a, "mV");
                 auto rhs = to_1d_vector<double>(rhs_a, "rhs");
                 if ((int)rhs.size() != n_face)
                     throw std::runtime_error("solve_linear_material_auto_prec_arrays: rhs size mismatch");
                 if ((int)B_indptr.size() < 1)
                     throw std::runtime_error("solve_linear_material_auto_prec_arrays: empty B_indptr");
                 std::vector<double> x0;
                 const std::vector<double>* x0_ptr = nullptr;
                 if (!x0_obj.is_none()) {
                     auto x0_a = py::cast<py::array_t<double, py::array::c_style | py::array::forcecast>>(x0_obj);
                     x0 = to_1d_vector<double>(x0_a, "x0");
                     if ((int)x0.size() != n_face)
                         throw std::runtime_error("solve_linear_material_auto_prec_arrays: x0 size mismatch");
                     x0_ptr = &x0;
                 }
                 std::vector<double> mass_diag((size_t)n_face, 0.0);
                 for (size_t k = 0; k < mV.size(); ++k) {
                     if (mI[k] == mJ[k] && mI[k] >= 0 && mI[k] < n_face) mass_diag[(size_t)mI[k]] += mV[k];
                 }
                 std::vector<std::vector<int>> supp_id((size_t)n_face);
                 std::vector<std::vector<double>> supp_val((size_t)n_face);
                 const int n_charge = (int)B_indptr.size() - 1;
                 for (int a = 0; a < n_charge; ++a) {
                     for (int k = B_indptr[(size_t)a]; k < B_indptr[(size_t)a + 1]; ++k) {
                         int f = B_indices[(size_t)k];
                         if (f < 0 || f >= n_face) throw std::runtime_error("solve_linear_material_auto_prec_arrays: B face index out of range");
                         supp_id[(size_t)f].push_back(a);
                         supp_val[(size_t)f].push_back(B_data[(size_t)k]);
                     }
                 }
                 std::vector<double> prec((size_t)n_face, 0.0);
                 for (int f = 0; f < n_face; ++f) {
                     double ndiag = 0.0;
                     const auto& ids = supp_id[(size_t)f];
                     const auto& vals = supp_val[(size_t)f];
                     for (size_t p = 0; p < ids.size(); ++p)
                         for (size_t q = 0; q < ids.size(); ++q)
                             ndiag += vals[p] * vals[q] * s.GetInteractionMatrixElement(ids[p], ids[q]);
                     double v = inv_chi * mass_diag[(size_t)f] + ndiag;
                     if (!(v > 0.0) || !std::isfinite(v)) v = 1.0;
                     prec[(size_t)f] = v;
                 }
                 int iters = 0;
                 std::vector<double> m = s.SolveLinearMaterial(B_indptr, B_indices, B_data, n_face,
                                                               mI, mJ, mV, inv_chi, prec, rhs,
                                                               tol, maxit, iters,
                                                               /*mass_riesz=*/false, /*symmetric=*/true,
                                                               x0_ptr);
                 double pmin = n_face ? prec[0] : 0.0, pmax = pmin;
                 for (double v : prec) { if (v < pmin) pmin = v; if (v > pmax) pmax = v; }
                 py::dict d;
                 d["m"] = m; d["iters"] = iters; d["prec_min"] = pmin; d["prec_max"] = pmax;
                 d["timings"] = solve_timings_dict(s);
                 return d;
             },
             py::arg("B_indptr"), py::arg("B_indices"), py::arg("B_data"), py::arg("n_face"),
             py::arg("mI"), py::arg("mJ"), py::arg("mV"), py::arg("inv_chi"), py::arg("rhs"),
             py::arg("tol") = 1e-9, py::arg("maxit") = 5000, py::arg("x0") = py::none(),
             "Array-input variant of solve_linear_material_auto_prec; avoids Python list materialization.")
        .def("solve_linear_material_mass_riesz",
             [](RadHACApKChargeGram& s,
                std::vector<int> B_indptr, std::vector<int> B_indices, std::vector<double> B_data,
                int n_face, std::vector<int> mI, std::vector<int> mJ, std::vector<double> mV,
                double inv_chi, std::vector<double> rhs,
                double tol, int maxit, bool symmetric) {
                 if ((int)rhs.size() != n_face)
                     throw std::runtime_error("solve_linear_material_mass_riesz: rhs size mismatch");
                 int iters = 0;
                 std::vector<double> noprec;   // mass_riesz=true ignores the diagonal prec
                 std::vector<double> m = s.SolveLinearMaterial(B_indptr, B_indices, B_data, n_face,
                                                               mI, mJ, mV, inv_chi, noprec, rhs,
                                                               tol, maxit, iters, /*mass_riesz=*/true, symmetric);
                 py::dict timings;
                 for (const auto& kv : s.LastSolveTimings()) timings[py::str(kv.first)] = kv.second;
                 py::dict d; d["m"] = m; d["iters"] = iters; d["timings"] = timings; return d;
             },
             py::arg("B_indptr"), py::arg("B_indices"), py::arg("B_data"), py::arg("n_face"),
             py::arg("mI"), py::arg("mJ"), py::arg("mV"), py::arg("inv_chi"),
             py::arg("rhs"), py::arg("tol") = 1e-8, py::arg("maxit") = 5000, py::arg("symmetric") = true,
             "DEFAULT linear demag solve ENTIRELY in C++: SPD +N system ((1/chi)M_mass + B^T G B) m = rhs "
             "by CG preconditioned with a PARDISO SPD factor of the HDiv mass M_mass (the MASS RIESZ map). "
             "symmetric=true (default) applies G via the EXACTLY-symmetric H-matvec so CG sees a symmetric operator; "
             "symmetric=false uses the general (asymmetric ACA) matvec. Returns {m, iters}.")
        .def("solve_linear_material_mass_riesz_arrays",
             [](RadHACApKChargeGram& s,
                py::array_t<int, py::array::c_style | py::array::forcecast> B_indptr_a,
                py::array_t<int, py::array::c_style | py::array::forcecast> B_indices_a,
                py::array_t<double, py::array::c_style | py::array::forcecast> B_data_a,
                int n_face,
                py::array_t<int, py::array::c_style | py::array::forcecast> mI_a,
                py::array_t<int, py::array::c_style | py::array::forcecast> mJ_a,
                 py::array_t<double, py::array::c_style | py::array::forcecast> mV_a,
                 double inv_chi,
                 py::array_t<double, py::array::c_style | py::array::forcecast> rhs_a,
                 double tol, int maxit, bool symmetric, py::object x0_obj) {
                 auto B_indptr = to_1d_vector<int>(B_indptr_a, "B_indptr");
                 auto B_indices = to_1d_vector<int>(B_indices_a, "B_indices");
                 auto B_data = to_1d_vector<double>(B_data_a, "B_data");
                 auto mI = to_1d_vector<int>(mI_a, "mI");
                 auto mJ = to_1d_vector<int>(mJ_a, "mJ");
                 auto mV = to_1d_vector<double>(mV_a, "mV");
                 auto rhs = to_1d_vector<double>(rhs_a, "rhs");
                 if ((int)rhs.size() != n_face)
                     throw std::runtime_error("solve_linear_material_mass_riesz_arrays: rhs size mismatch");
                 std::vector<double> x0;
                 const std::vector<double>* x0_ptr = nullptr;
                 if (!x0_obj.is_none()) {
                     auto x0_a = py::cast<py::array_t<double, py::array::c_style | py::array::forcecast>>(x0_obj);
                     x0 = to_1d_vector<double>(x0_a, "x0");
                     if ((int)x0.size() != n_face)
                         throw std::runtime_error("solve_linear_material_mass_riesz_arrays: x0 size mismatch");
                     x0_ptr = &x0;
                 }
                 int iters = 0;
                 std::vector<double> noprec;
                 std::vector<double> m = s.SolveLinearMaterial(B_indptr, B_indices, B_data, n_face,
                                                               mI, mJ, mV, inv_chi, noprec, rhs,
                                                               tol, maxit, iters, /*mass_riesz=*/true, symmetric,
                                                               x0_ptr);
                 py::dict d; d["m"] = m; d["iters"] = iters; d["timings"] = solve_timings_dict(s); return d;
              },
              py::arg("B_indptr"), py::arg("B_indices"), py::arg("B_data"), py::arg("n_face"),
              py::arg("mI"), py::arg("mJ"), py::arg("mV"), py::arg("inv_chi"),
              py::arg("rhs"), py::arg("tol") = 1e-8, py::arg("maxit") = 5000, py::arg("symmetric") = true,
               py::arg("x0") = py::none(),
               "Array-input variant of solve_linear_material_mass_riesz; avoids Python list materialization.")
        .def("apply_demag_arrays",
             [](RadHACApKChargeGram& s,
                py::array_t<int, py::array::c_style | py::array::forcecast> B_indptr_a,
                py::array_t<int, py::array::c_style | py::array::forcecast> B_indices_a,
                py::array_t<double, py::array::c_style | py::array::forcecast> B_data_a,
                int n_face,
                py::array_t<double, py::array::c_style | py::array::forcecast> x_a,
                bool symmetric) {
                 return s.ApplyDemagOperator(
                     to_1d_vector<int>(B_indptr_a, "B_indptr"),
                     to_1d_vector<int>(B_indices_a, "B_indices"),
                     to_1d_vector<double>(B_data_a, "B_data"), n_face,
                     to_1d_vector<double>(x_a, "x"), symmetric);
             },
             py::arg("B_indptr"), py::arg("B_indices"), py::arg("B_data"),
             py::arg("n_face"), py::arg("x"), py::arg("symmetric") = true,
             "Apply the matrix-free HDiv demagnetizing operator B^T G B to x in C++.")
        .def("apply_mass_riesz_arrays",
             [](RadHACApKChargeGram& s,
                py::array_t<int, py::array::c_style | py::array::forcecast> mI_a,
                py::array_t<int, py::array::c_style | py::array::forcecast> mJ_a,
                py::array_t<double, py::array::c_style | py::array::forcecast> mV_a,
                int n_face,
                py::array_t<double, py::array::c_style | py::array::forcecast> rhs_a) {
                 return s.ApplyMassRiesz(
                     to_1d_vector<int>(mI_a, "mI"),
                     to_1d_vector<int>(mJ_a, "mJ"),
                     to_1d_vector<double>(mV_a, "mV"), n_face,
                     to_1d_vector<double>(rhs_a, "rhs"));
             },
             py::arg("mI"), py::arg("mJ"), py::arg("mV"),
             py::arg("n_face"), py::arg("rhs"),
             "Apply the persistent C++ PARDISO mass-Riesz map M_mass^{-1} rhs.")
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
             "entirely in C++ (each step a mass-Riesz SolveLinearMaterial + closed-form chi update; G via "
             "the analytic H-matvec). Returns {m, Mavg, chi, Dscal, iters}.")
        .def("stats", [](RadHACApKChargeGram& s) {
                 const RadHACApKStats& st = s.GetStats();
                 py::dict d;
                 d["n_dof"] = st.n_dof; d["n_leaves"] = st.n_leaves; d["n_lowrank"] = st.n_lowrank;
                 d["n_dense"] = st.n_dense; d["max_rank"] = st.max_rank;
                 d["compression"] = st.compression; d["build_time"] = st.build_time;
                 d["memory_mb"] = st.memory_mb; d["dense_memory_mb"] = st.dense_memory_mb;
                 for (const auto& kv : s.HexCacheStats()) d[kv.first.c_str()] = kv.second;
                 return d;
             }, "H-matrix stats dict.");

    // ========================================================================
    // Object Creation
    // ========================================================================

    // NOTE: The Python-facing ObjRecMag constructor has been RETIRED. The
    // canonical rectangular permanent magnet is now the fixed-M surface-charge
    // hex via radia.magnet_box(...) (Python shim radia.magnet.ObjRecMag also
    // forwards to magnet_box). The internal C++ surface-current rectangular
    // block kernel (radTRecMag) and the C API RadObjRecMag are KEPT.

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
              5 surface-face coefficients for fixed-magnetization field evaluation.

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

    m.def("ObjPyramid", &radia_objects::ObjPyramid,
          py::arg("vertices"), py::arg("magnetization"),
          R"pbdoc(
              Create pyramid element from 5 vertices.

              Square-base pyramid: 1 quadrilateral base + 4 triangular sides (5 faces total).
              Retained for fixed-magnetization field evaluation.  Mesh-backed
              magnetic-material solves route through HDiv-VIM.

              Vertex convention (matches netgen_mesh_import.PYRAMID_FACES):
                  v0..v3 = base quad, v4 = apex.
                  Base : (v0, v3, v2, v1) -- outward normal away from the apex.
                  Sides: (v0,v1,v4) (v1,v2,v4) (v2,v3,v4) (v3,v0,v4).

              Args:
                  vertices: List of 5 vertex coordinates [[x,y,z], ...] (base v0..v3, then apex v4)
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

    m.def("PlanarChargeField", &radia_planar_charges::PlanarChargeField,
          py::arg("Xq"), py::arg("Q"), py::arg("P"),
          R"pbdoc(
              2D planar exterior field H at points P from a point-charge cloud
              (kernel -ln(r)/(2 pi)).  Shared by planar HDiv-VIM and dense planar helpers.
              Args:  Xq (nq,2), Q (nq,), P (nP,2) float64 -> H (nP,2) float64
          )pbdoc");

    m.def("PlanarChargeAz", &radia_planar_charges::PlanarChargeAz,
          py::arg("Xq"), py::arg("Q"), py::arg("P"),
          R"pbdoc(
              2D planar out-of-plane vector potential A_z at points P from a
              point-charge cloud: A_z = mu0/(2 pi) sum_a Q_a atan2(dy, dx).
              Args:  Xq (nq,2), Q (nq,), P (nP,2) float64 -> Az (nP,) float64
          )pbdoc");

    m.def("PlanarMaxwellTorqueCircle", &radia_planar_charges::PlanarMaxwellTorqueCircle,
          py::arg("Xq"), py::arg("Q"), py::arg("Rc"),
          py::arg("cx") = 0.0, py::arg("cy") = 0.0, py::arg("n") = 1440,
          py::arg("hextx") = 0.0, py::arg("hexty") = 0.0,
          R"pbdoc(
              Maxwell-stress torque per unit length about (cx,cy) on a circle of
              radius Rc in air, from a 2D point-charge cloud + uniform applied
              field (hextx,hexty): T = mu0 Rc^2 oint H_r H_phi dphi (n points).
          )pbdoc");

    m.def("PlanarMaxwellForceCircle", &radia_planar_charges::PlanarMaxwellForceCircle,
          py::arg("Xq"), py::arg("Q"), py::arg("Rc"),
          py::arg("cx") = 0.0, py::arg("cy") = 0.0, py::arg("n") = 1440,
          py::arg("hextx") = 0.0, py::arg("hexty") = 0.0,
          R"pbdoc(
              Maxwell-stress FORCE per unit length on a circle of radius Rc in air,
              from a 2D point-charge cloud + uniform applied field: Fout (2,) =
              mu0 Rc oint [H_r H - 1/2 |H|^2 n] dphi.  A uniform field gives ~0 net
              force; concatenate multiple bodies' clouds for the maglev inter-body force.
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
        (cHACApK_hlu_decomp + cHACApK_hlu_solve_vec). Call right after an
        H-LU probe (e.g. an HDiv-VIM H-LU round-trip) to get the
        factorization and solve wall times.
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
    // The shared-DOF approach was fundamentally incompatible with independent face coefficients.

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

    m.def("MatHysForwardBatch", &radia_material_ext::MatHysForwardBatch,
          py::arg("mat"), py::arg("B"), py::arg("states"),
          R"pbdoc(
              Batched Forward: H[i] = nu_rev * B[i] + H_irr(B[i]; states[i]).

              Restores each row's COMMITTED state before evaluating, so the call
              is pure w.r.t. the states (nothing commits).  Serial C++ loop with
              the GIL released -- the batched replacement for a per-row Python
              MatHysRestoreState + MatHysIrreversible loop.

              Args:
                  mat: Material handle from MatPlayHysteresis()
                  B: (n, 3) flux densities in Tesla
                  states: (n, state_len) committed states from MatHysSaveState()

              Returns:
                  (n, 3) magnetic field H in A/m
          )pbdoc");

    m.def("MatHysCommitBatch", &radia_material_ext::MatHysCommitBatch,
          py::arg("mat"), py::arg("B"), py::arg("states"),
          R"pbdoc(
              Batched commit: for each row, restore states[i], play to B[i],
              commit, and return the new committed state row.

              Args:
                  mat: Material handle from MatPlayHysteresis()
                  B: (n, 3) converged flux densities in Tesla
                  states: (n, state_len) committed states to advance

              Returns:
                  (n, state_len) new committed states
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
                  b_input_newton (bool): Enable B-input hysteresis stepping
                      for supported 3-DOF dipole relaxation (default: False)
                  b_input_hantila (bool): Enable Hantila B-input stepping
                      for supported 3-DOF dipole relaxation (default: False)
                  hantila_alpha (float): Hantila polarization parameter (0=auto, default: 0)
                  hantila_relax (float): Hantila under-relaxation (0=full step, default: 0)

              Example:
                  rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
                  rad.SolverConfig(bicgstab_tol=1e-6, relax_param=0.3)
                  rad.SolverConfig(newton_method=True, newton_damping=True)
                  rad.SolverConfig(b_input_newton=True)  # B-input hysteresis stepping
          )pbdoc");

    m.def("GetSolverConfig", &radia_solver_ext::GetSolverConfig,
          R"pbdoc(
              Get current solver configuration.

              Returns:
                  Dictionary with all solver parameters:
                  - bicgstab_tol
                  - relax_param, newton_method, b_input_newton, b_input_hantila
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
             "Pre-cache field values at given points for fast gf.Set() (direct rad.Fld per point).")
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
    // Resolves the former Phase 2 66% near-field undershoot from
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
            undershoot seen before the dyadic-kernel fix.

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
    // fixed-magnetization field kernels, ...).  ACA+ is HACApK's cHACApK_acaplus.
    // ========================================================================
    m.def("_stream_aca_tsvd",
        [](int M, int N, std::function<double(int, int)> entry,
           int modes, int kmax, double aca_eps)
        {
            if (M <= 0 || N <= 0)
                throw std::invalid_argument("M and N must be positive");
            if (!entry)
                throw std::invalid_argument("entry callback must be callable");
            radia::stream_function::TSVDResult r =
                radia::stream_function::ACATSVD(
                    M, N, entry, modes, kmax, aca_eps);

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
        py::arg("aca_eps") = 1.0e-4,
        R"pbdoc(
            (ACA+)+TSVD recompressed truncated SVD of an M x N matrix A whose
            entries are supplied on demand by the callback entry(i, j) ->
            A(i,j) (0-based row i in [0,M), col j in [0,N)).  Kernel-agnostic:
            the callback may use any Radia field computation (coil Biot-Savart,
            fixed-magnetization field kernels, ...).  ACA+ is delegated to HACApK
            (cHACApK_acaplus).  Returns (U, S, V, k_aca):
              U:     (M, modes) row-major
              S:     (modes,)
              V:     (N, modes) row-major
              k_aca: ACA+ rank found before truncation
            Recompression is the standard "SVD of a low-rank product": QR each
            tall-skinny ACA factor (C, D), then ONE small kt x kt SVD (peer
            review JIAM-2026-36).  The legacy manuscript Method 2/3 were removed.
        )pbdoc");
}
