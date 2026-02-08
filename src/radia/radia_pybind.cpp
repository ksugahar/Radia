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

// Radia core headers
#include "radentry.h"
#include "rad_constants.h"

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
 * @brief Compute field at a single point
 *
 * @param obj Object handle
 * @param field_type Field type: "b", "h", "a", "m", "bx", "by", "bz", etc.
 * @param point Evaluation point [x, y, z]
 * @return Field value (scalar or vector depending on field_type)
 */
py::object Fld(int obj, const std::string& field_type, py::array_t<double> point) {
    auto p = point.unchecked<1>();
    if (p.size() != 3) {
        throw std::runtime_error("point must have 3 coordinates");
    }

    double coords[3] = {p(0), p(1), p(2)};
    double result[6] = {0};  // Max size for any field type
    int nResult = 0;

    char* id = const_cast<char*>(field_type.c_str());
    int err = RadFld(result, &nResult, obj, id, coords, 1);
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
 * @brief Batch field computation at multiple points
 *
 * @param obj Object handle
 * @param points Numpy array of shape (N, 3)
 * @param method Computation method (0=direct, 1=FMM)
 * @return Dictionary with 'B' and 'H' arrays of shape (N, 3)
 */
py::dict FldBatch(int obj, py::array_t<double> points, int method = 0) {
    py::buffer_info buf = points.request();

    if (buf.ndim != 2 || buf.shape[1] != 3) {
        throw std::runtime_error("points must be array of shape (N, 3)");
    }

    int n_points = static_cast<int>(buf.shape[0]);
    double* pts = static_cast<double*>(buf.ptr);

    // Allocate output arrays
    std::vector<double> B_out(n_points * 3);
    std::vector<double> H_out(n_points * 3);

    // Release GIL for long computation
    {
        py::gil_scoped_release release;
        int err = RadFldBatch(B_out.data(), H_out.data(), n_points, pts, obj, method);
        check_error(err);
    }

    // Create numpy arrays
    py::array_t<double> B_arr({n_points, 3});
    py::array_t<double> H_arr({n_points, 3});

    auto B_buf = B_arr.mutable_unchecked<2>();
    auto H_buf = H_arr.mutable_unchecked<2>();

    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < 3; j++) {
            B_buf(i, j) = B_out[i * 3 + j];
            H_buf(i, j) = H_out[i * 3 + j];
        }
    }

    py::dict result;
    result["B"] = B_arr;
    result["H"] = H_arr;
    return result;
}

/**
 * @brief Compute vector potential A at multiple points
 *
 * @param obj Object handle
 * @param points Numpy array of shape (N, 3)
 * @return Numpy array of shape (N, 3) with A vectors
 */
py::array_t<double> FldA(int obj, py::array_t<double> points) {
    py::buffer_info buf = points.request();

    if (buf.ndim != 2 || buf.shape[1] != 3) {
        throw std::runtime_error("points must be array of shape (N, 3)");
    }

    int n_points = static_cast<int>(buf.shape[0]);
    double* pts = static_cast<double*>(buf.ptr);

    std::vector<double> A_out(n_points * 3);

    {
        py::gil_scoped_release release;
        int err = RadFldA(A_out.data(), n_points, pts, obj);
        check_error(err);
    }

    py::array_t<double> result({n_points, 3});
    auto r = result.mutable_unchecked<2>();

    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < 3; j++) {
            r(i, j) = A_out[i * 3 + j];
        }
    }

    return result;
}

/**
 * @brief Compute scalar potential Phi at multiple points
 *
 * @param obj Object handle
 * @param points Numpy array of shape (N, 3)
 * @return Numpy array of shape (N,) with Phi values
 */
py::array_t<double> FldPhi(int obj, py::array_t<double> points) {
    py::buffer_info buf = points.request();

    if (buf.ndim != 2 || buf.shape[1] != 3) {
        throw std::runtime_error("points must be array of shape (N, 3)");
    }

    int n_points = static_cast<int>(buf.shape[0]);
    double* pts = static_cast<double*>(buf.ptr);

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
std::string FldVTS(int obj, const std::string& filename,
                   py::array_t<double> x_range,
                   py::array_t<double> y_range,
                   py::array_t<double> z_range,
                   int nx, int ny, int nz,
                   int include_B, int include_H,
                   double unit_scale) {

    auto xr = x_range.unchecked<1>();
    auto yr = y_range.unchecked<1>();
    auto zr = z_range.unchecked<1>();

    if (xr.size() != 2 || yr.size() != 2 || zr.size() != 2) {
        throw std::runtime_error("ranges must have 2 elements [min, max]");
    }

    {
        py::gil_scoped_release release;
        int err = RadFldVTS(obj, filename.c_str(),
                           xr(0), xr(1), nx,
                           yr(0), yr(1), ny,
                           zr(0), zr(1), nz,
                           include_B, include_H, unit_scale);
        check_error(err);
    }

    return filename;
}

/**
 * @brief Set physical units
 *
 * @param unit_str Unit string: "mm" or "m"
 */
void FldUnits(const std::string& unit_str) {
    if (!unit_str.empty()) {
        int err = RadFldUnitsSet(unit_str.c_str());
        check_error(err);
    }
}

/**
 * @brief Get current units as string
 * @return Unit information string
 */
std::string FldUnitsGet() {
    char buf[1024] = {0};
    int err = RadFldUnits(buf);
    check_error(err);
    return std::string(buf);
}

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

/**
 * @brief Create SIBC material
 *
 * @param sigma Conductivity [S/m]
 * @param mu_r Relative permeability
 * @return Material handle
 */
int MatSIBC(double sigma, double mu_r) {
    int handle = 0;
    int err = RadMatSIBC(&handle, sigma, mu_r);
    check_error(err);
    return handle;
}

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
        return py::none().cast<py::dict>();
    }

    py::dict result;
    if (n >= 4) {
        result["t_matrix_build"] = stats[0];
        result["t_linear_solve"] = stats[1];
        result["linear_iterations"] = static_cast<int>(stats[2]);
        result["nonl_iterations"] = static_cast<int>(stats[3]);
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
 * @brief Get object geometry limits
 * @param obj Object handle
 * @return [xmin, xmax, ymin, ymax, zmin, zmax]
 */
py::array_t<double> ObjGeoLim(int obj) {
    double L[6] = {0};
    int err = RadObjGeoLim(L, obj);
    check_error(err);

    py::array_t<double> result(6);
    auto r = result.mutable_unchecked<1>();
    for (int i = 0; i < 6; i++) {
        r(i) = L[i];
    }
    return result;
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
 */
py::dict ObjM(int obj) {
    double M[6] = {0};  // [x, y, z, Mx, My, Mz]
    int arMesh[10] = {0};
    int err = RadObjM(M, arMesh, obj);
    check_error(err);

    py::dict result;
    result["center"] = py::make_tuple(M[0], M[1], M[2]);
    result["magnetization"] = py::make_tuple(M[3], M[4], M[5]);
    return result;
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
    result["compression"] = dOut[3];
    result["build_time"] = dOut[4];
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

// SetIMASymmetry, BuildIMAMatrix REMOVED (2026-01-31)
// Use BuildMatrix(obj, image="+x-z") or Solve(obj, ..., image="+x-z") instead

} // namespace radia_solver_ext


// ============================================================================
// Additional Field Functions
// ============================================================================

namespace radia_field_ext {

double FldEnr(int obj_dst, int obj_src) {
    double energy = 0;
    int SbdPar[3] = {0, 0, 0};
    int err = RadFldEnr(&energy, obj_dst, obj_src, SbdPar);
    check_error(err);
    return energy;
}

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

double FldFocPot(int obj, py::array_t<double> p1, py::array_t<double> p2, int np) {
    auto pt1 = p1.unchecked<1>();
    auto pt2 = p2.unchecked<1>();
    double P1[3] = {pt1(0), pt1(1), pt1(2)};
    double P2[3] = {pt2(0), pt2(1), pt2(2)};

    double d = 0;
    int err = RadFldFocPot(&d, obj, P1, P2, np);
    check_error(err);
    return d;
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


// ============================================================================
// Conductor PEEC Functions
// ============================================================================

namespace radia_conductor {

int CndRecBlock(py::array_t<double> center, py::array_t<double> dimensions, double sigma, int num_panels) {
    auto c = center.unchecked<1>();
    auto d = dimensions.unchecked<1>();
    double P[3] = {c(0), c(1), c(2)};
    double L[3] = {d(0), d(1), d(2)};

    int handle = 0;
    int err = RadCndRecBlock(&handle, P, L, sigma, num_panels);
    check_error(err);
    return handle;
}

int CndLoop(py::array_t<double> center, double radius, py::array_t<double> normal,
            const std::string& cross_section, double width, double height,
            double sigma, int num_around, int num_loop) {
    auto c = center.unchecked<1>();
    auto n = normal.unchecked<1>();
    double P[3] = {c(0), c(1), c(2)};
    double N[3] = {n(0), n(1), n(2)};
    char cs = cross_section[0];

    int handle = 0;
    int err = RadCndLoop(&handle, P, radius, N, cs, width, height, sigma, num_around, num_loop);
    check_error(err);
    return handle;
}

void CndSetFrequency(int cond, double freq) {
    int err = RadCndSetFrequency(cond, freq);
    check_error(err);
}

void CndSetMuR(int cond, double mu_r) {
    int err = RadCndSetMuR(cond, mu_r);
    check_error(err);
}

double CndGetSkinDepth(int cond) {
    double delta = 0;
    int err = RadCndGetSkinDepth(&delta, cond);
    check_error(err);
    return delta;
}

std::complex<double> CndGetSurfaceImpedance(int cond) {
    double Z_real = 0, Z_imag = 0;
    int err = RadCndGetSurfaceImpedance(&Z_real, &Z_imag, cond);
    check_error(err);
    return std::complex<double>(Z_real, Z_imag);
}

void CndSetVoltage(int cond, double V_real, double V_imag) {
    int err = RadCndSetVoltage(cond, V_real, V_imag);
    check_error(err);
}

void CndSetCurrent(int cond, double I_real, double I_imag) {
    int err = RadCndSetCurrent(cond, I_real, I_imag);
    check_error(err);
}

std::complex<double> CndGetTotalCurrent(int cond) {
    double I_real = 0, I_imag = 0;
    int err = RadCndGetTotalCurrent(&I_real, &I_imag, cond);
    check_error(err);
    return std::complex<double>(I_real, I_imag);
}

void CndSolve(int cond) {
    int err = RadCndSolve(cond);
    check_error(err);
}

std::complex<double> CndGetImpedance(int cond) {
    double Z_real = 0, Z_imag = 0;
    int err = RadCndGetImpedance(&Z_real, &Z_imag, cond);
    check_error(err);
    return std::complex<double>(Z_real, Z_imag);
}

void CndDefinePortAuto(int cond) {
    int err = RadCndDefinePortAuto(cond);
    check_error(err);
}

int CndNumPanels(int cond) {
    int n = 0;
    int err = RadCndNumPanels(&n, cond);
    check_error(err);
    return n;
}

void CndFmmSetEnabled(bool enabled) {
    int err = RadCndFmmSetEnabled(enabled ? 1 : 0);
    check_error(err);
}

bool CndFmmGetEnabled() {
    int enabled = 0;
    int err = RadCndFmmGetEnabled(&enabled);
    check_error(err);
    return enabled != 0;
}

void CndFmmSetParameters(int p, int ncrit, int threshold) {
    int err = RadCndFmmSetParameters(p, ncrit, threshold);
    check_error(err);
}

py::tuple CndFmmGetParameters() {
    int p = 0, ncrit = 0, threshold = 0;
    int err = RadCndFmmGetParameters(&p, &ncrit, &threshold);
    check_error(err);
    return py::make_tuple(p, ncrit, threshold);
}

} // namespace radia_conductor


// ============================================================================
// Coupled PEEC+MMM Solver Functions
// ============================================================================

namespace radia_coupled {

int CplMagCreate(int conductor, int magnet) {
    double pOut[10] = {0};
    int nOut = 0;
    int err = RadCplMagCreate(pOut, &nOut, conductor, magnet);
    check_error(err);
    return static_cast<int>(pOut[0]);  // Return solver handle
}

void CplMagSetFrequency(int solver, double freq) {
    int err = RadCplMagSetFrequency(solver, freq);
    check_error(err);
}

void CplMagSetVoltage(int solver, double V_real, double V_imag) {
    int err = RadCplMagSetVoltage(solver, V_real, V_imag);
    check_error(err);
}

void CplMagSetCurrent(int solver, double I_real, double I_imag) {
    int err = RadCplMagSetCurrent(solver, I_real, I_imag);
    check_error(err);
}

void CplMagSetExtField(int solver, double Hx, double Hy, double Hz) {
    int err = RadCplMagSetExtField(solver, Hx, Hy, Hz);
    check_error(err);
}

void CplMagSetMu(int solver, double mu_r_real, double mu_r_imag) {
    int err = RadCplMagSetMu(solver, mu_r_real, mu_r_imag);
    check_error(err);
}

void CplMagSetSymmetric(int solver, bool use_symmetric) {
    int err = RadCplMagSetSymmetric(solver, use_symmetric ? 1 : 0);
    check_error(err);
}

py::dict CplMagSolve(int solver) {
    double pOut[20] = {0};
    int nOut = 0;

    {
        py::gil_scoped_release release;
        int err = RadCplMagSolve(pOut, &nOut, solver);
        check_error(err);
    }

    py::dict result;
    result["Z_real"] = pOut[0];
    result["Z_imag"] = pOut[1];
    result["P_conductor"] = pOut[2];
    result["P_magnet"] = pOut[3];
    result["iterations"] = static_cast<int>(pOut[4]);
    return result;
}

std::complex<double> CplMagImpedance(int solver) {
    double pOut[2] = {0};
    int nOut = 0;
    int err = RadCplMagImpedance(pOut, &nOut, solver);
    check_error(err);
    return std::complex<double>(pOut[0], pOut[1]);
}

py::tuple CplMagPower(int solver) {
    double pOut[2] = {0};
    int nOut = 0;
    int err = RadCplMagPower(pOut, &nOut, solver);
    check_error(err);
    return py::make_tuple(pOut[0], pOut[1]);
}

py::dict CplMagFld(int solver, py::array_t<double> point) {
    auto p = point.unchecked<1>();
    double P[3] = {p(0), p(1), p(2)};

    double pOut[6] = {0};
    int nOut = 0;
    int err = RadCplMagFld(pOut, &nOut, solver, P);
    check_error(err);

    py::dict result;
    result["B_real"] = py::make_tuple(pOut[0], pOut[1], pOut[2]);
    result["B_imag"] = py::make_tuple(pOut[3], pOut[4], pOut[5]);
    return result;
}

void CplMagDelete(int solver) {
    int err = RadCplMagDelete(solver);
    check_error(err);
}

} // namespace radia_coupled


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
// Module Definition
// ============================================================================

PYBIND11_MODULE(_radia_pybind, m) {
    m.doc() = R"pbdoc(
        Radia - 3D Magnetostatics Library (pybind11 bindings)

        This module provides Python bindings for Radia using pybind11.
        It follows NGSolve design patterns for clean, efficient bindings.

        Example:
            import radia as rad
            rad.FldUnits('m')

            # Create rectangular magnet
            magnet = rad.ObjRecMag([0,0,0], [0.04, 0.04, 0.02], [0, 0, 954930])

            # Compute field
            B = rad.Fld(magnet, 'b', [0.05, 0, 0])
    )pbdoc";

    // Version info
    m.attr("__version__") = "1.4.0";

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
          py::arg("obj"), py::arg("field_type"), py::arg("point"),
          R"pbdoc(
              Compute field at a single point.

              Args:
                  obj: Object handle
                  field_type: "b", "h", "a", "m", "bx", "by", "bz", etc.
                  point: Evaluation point [x, y, z]

              Returns:
                  Field value (scalar or vector)
          )pbdoc");

    m.def("FldBatch", &radia_field::FldBatch,
          py::arg("obj"), py::arg("points"), py::arg("method") = 0,
          R"pbdoc(
              Batch field computation at multiple points.

              More efficient than calling Fld() in a loop.

              Args:
                  obj: Object handle
                  points: Numpy array of shape (N, 3)
                  method: Computation method (0=direct, 1=FMM)

              Returns:
                  Dictionary with 'B' and 'H' arrays of shape (N, 3)
          )pbdoc");

    m.def("FldA", &radia_field::FldA,
          py::arg("obj"), py::arg("points"),
          py::call_guard<py::gil_scoped_release>(),
          R"pbdoc(
              Compute vector potential A at multiple points.

              Args:
                  obj: Object handle
                  points: Numpy array of shape (N, 3)

              Returns:
                  Numpy array of shape (N, 3) with A vectors in T*m
          )pbdoc");

    m.def("FldPhi", &radia_field::FldPhi,
          py::arg("obj"), py::arg("points"),
          py::call_guard<py::gil_scoped_release>(),
          R"pbdoc(
              Compute scalar potential Phi at multiple points.

              Args:
                  obj: Object handle
                  points: Numpy array of shape (N, 3)

              Returns:
                  Numpy array of shape (N,) with Phi values in A
          )pbdoc");

    m.def("FldVTS", &radia_field::FldVTS,
          py::arg("obj"), py::arg("filename"),
          py::arg("x_range"), py::arg("y_range"), py::arg("z_range"),
          py::arg("nx") = 21, py::arg("ny") = 21, py::arg("nz") = 21,
          py::arg("include_B") = 1, py::arg("include_H") = 0,
          py::arg("unit_scale") = 1.0,
          py::call_guard<py::gil_scoped_release>(),
          R"pbdoc(
              Export field to VTS file.

              Args:
                  obj: Object handle
                  filename: Output filename
                  x_range: [xmin, xmax]
                  y_range: [ymin, ymax]
                  z_range: [zmin, zmax]
                  nx, ny, nz: Number of grid points
                  include_B: Include B field (default: 1)
                  include_H: Include H field (default: 0)
                  unit_scale: Coordinate scale factor (default: 1.0)

              Returns:
                  Filename
          )pbdoc");

    m.def("FldUnits", &radia_field::FldUnits,
          py::arg("unit_str") = "",
          R"pbdoc(
              Set or get physical units.

              Args:
                  unit_str: "mm" or "m" (empty to get current units)
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

    m.def("MatSIBC", &radia_material::MatSIBC,
          py::arg("sigma"), py::arg("mu_r"),
          R"pbdoc(
              Create SIBC (Surface Impedance Boundary Condition) material.

              Args:
                  sigma: Conductivity [S/m]
                  mu_r: Relative permeability

              Returns:
                  Material handle
          )pbdoc");

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

    m.def("ObjGeoLim", &radia_utility::ObjGeoLim,
          py::arg("obj"),
          "Get object geometry limits [xmin, xmax, ymin, ymax, zmin, zmax].");

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
                  Dictionary with 'center' and 'magnetization' tuples
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

    m.def("SetHACApKParams", &radia_solver_ext::SetHACApKParams,
          py::arg("eps"), py::arg("leaf_size"), py::arg("eta"),
          R"pbdoc(
              Set H-matrix (HACApK) parameters.

              Args:
                  eps: ACA tolerance (default: 1e-4)
                  leaf_size: Minimum cluster size (default: 10)
                  eta: Admissibility parameter (default: 2.0)
          )pbdoc");

    m.def("SetHMatrixEpsilon", &radia_solver_ext::SetHMatrixEpsilon,
          py::arg("eps"),
          "Set H-matrix ACA tolerance.");

    m.def("GetHACApKStats", &radia_solver_ext::GetHACApKStats,
          R"pbdoc(
              Get H-matrix statistics from last solve.

              Returns:
                  Dictionary with n_lowrank, n_dense, max_rank, compression, build_time
          )pbdoc");

    m.def("SetBiCGSTABTol", &radia_solver_ext::SetBiCGSTABTol,
          py::arg("tol"),
          "Set BiCGSTAB convergence tolerance.");

    m.def("GetBiCGSTABTol", &radia_solver_ext::GetBiCGSTABTol,
          "Get BiCGSTAB convergence tolerance.");

    m.def("SetRelaxParam", &radia_solver_ext::SetRelaxParam,
          py::arg("relax"),
          "Set under-relaxation parameter (0=full step, <1=damped).");

    m.def("GetRelaxParam", &radia_solver_ext::GetRelaxParam,
          "Get under-relaxation parameter.");

    // Image symmetry functions REMOVED (2026-01-31)
    // SetIMASymmetry, BuildIMAMatrix, etc. are replaced by the unified API:
    //   rad.Solve(obj, prec, maxiter, method, image='+x-z')
    //   rad.BuildMatrix(obj, image='+x-z')
    // The 'image' parameter specifies mirror symmetry: "+x", "-z", "+x-z", etc.

    // ========================================================================
    // Extended Field Functions
    // ========================================================================

    m.def("FldEnr", &radia_field_ext::FldEnr,
          py::arg("obj_dst"), py::arg("obj_src"),
          "Compute interaction energy between objects [J].");

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

    m.def("FldFocPot", &radia_field_ext::FldFocPot,
          py::arg("obj"), py::arg("p1"), py::arg("p2"), py::arg("np"),
          R"pbdoc(
              Compute focusing potential (integral of By*dz).

              Args:
                  obj: Object handle
                  p1: Start point [x, y, z]
                  p2: End point [x, y, z]
                  np: Number of integration points

              Returns:
                  Focusing potential value
          )pbdoc");

    m.def("FldCmpCrt", &radia_field_ext::FldCmpCrt,
          py::arg("prcB"), py::arg("prcA"), py::arg("prcBInt"),
          py::arg("prcFrc"), py::arg("prcTrjCrd"), py::arg("prcTrjAng"),
          "Set field computation precision criteria.");

    m.def("FldLenRndSw", &radia_field_ext::FldLenRndSw,
          py::arg("on_off"),
          "Switch field lens/rendering mode ('on' or 'off').");

    // ========================================================================
    // Conductor PEEC Functions
    // ========================================================================

    m.def("CndRecBlock", &radia_conductor::CndRecBlock,
          py::arg("center"), py::arg("dimensions"), py::arg("sigma"), py::arg("num_panels"),
          R"pbdoc(
              Create rectangular conductor block.

              Args:
                  center: Center point [x, y, z]
                  dimensions: Block dimensions [Lx, Ly, Lz]
                  sigma: Conductivity [S/m]
                  num_panels: Number of surface panels

              Returns:
                  Conductor handle
          )pbdoc");

    m.def("CndLoop", &radia_conductor::CndLoop,
          py::arg("center"), py::arg("radius"), py::arg("normal"),
          py::arg("cross_section"), py::arg("width"), py::arg("height"),
          py::arg("sigma"), py::arg("num_around"), py::arg("num_loop"),
          R"pbdoc(
              Create loop conductor (coil).

              Args:
                  center: Center point [x, y, z]
                  radius: Loop radius
                  normal: Normal vector [nx, ny, nz]
                  cross_section: "r" (rectangular) or "c" (circular)
                  width: Cross section width
                  height: Cross section height
                  sigma: Conductivity [S/m]
                  num_around: Panels around cross section
                  num_loop: Panels around loop

              Returns:
                  Conductor handle
          )pbdoc");

    m.def("CndSetFrequency", &radia_conductor::CndSetFrequency,
          py::arg("cond"), py::arg("freq"),
          "Set analysis frequency [Hz].");

    m.def("CndSetMuR", &radia_conductor::CndSetMuR,
          py::arg("cond"), py::arg("mu_r"),
          "Set relative permeability of conductor.");

    m.def("CndGetSkinDepth", &radia_conductor::CndGetSkinDepth,
          py::arg("cond"),
          "Get skin depth [m].");

    m.def("CndGetSurfaceImpedance", &radia_conductor::CndGetSurfaceImpedance,
          py::arg("cond"),
          "Get surface impedance [Ohm].");

    m.def("CndSetVoltage", &radia_conductor::CndSetVoltage,
          py::arg("cond"), py::arg("V_real"), py::arg("V_imag"),
          "Set port voltage [V] (complex).");

    m.def("CndSetCurrent", &radia_conductor::CndSetCurrent,
          py::arg("cond"), py::arg("I_real"), py::arg("I_imag"),
          "Set port current [A] (complex).");

    m.def("CndGetTotalCurrent", &radia_conductor::CndGetTotalCurrent,
          py::arg("cond"),
          "Get total current [A] (complex).");

    m.def("CndSolve", &radia_conductor::CndSolve,
          py::arg("cond"),
          "Solve PEEC system.");

    m.def("CndGetImpedance", &radia_conductor::CndGetImpedance,
          py::arg("cond"),
          "Get impedance [Ohm] (complex).");

    m.def("CndDefinePortAuto", &radia_conductor::CndDefinePortAuto,
          py::arg("cond"),
          "Automatically define port from conductor geometry.");

    m.def("CndNumPanels", &radia_conductor::CndNumPanels,
          py::arg("cond"),
          "Get number of surface panels.");

    m.def("CndFmmSetEnabled", &radia_conductor::CndFmmSetEnabled,
          py::arg("enabled"),
          "Enable/disable FMM acceleration.");

    m.def("CndFmmGetEnabled", &radia_conductor::CndFmmGetEnabled,
          "Check if FMM is enabled.");

    m.def("CndFmmSetParameters", &radia_conductor::CndFmmSetParameters,
          py::arg("p"), py::arg("ncrit"), py::arg("threshold"),
          "Set FMM parameters (expansion order, critical size, threshold).");

    m.def("CndFmmGetParameters", &radia_conductor::CndFmmGetParameters,
          "Get FMM parameters as (p, ncrit, threshold).");

    // ========================================================================
    // Coupled PEEC+MMM Solver
    // ========================================================================

    m.def("CplMagCreate", &radia_coupled::CplMagCreate,
          py::arg("conductor"), py::arg("magnet"),
          R"pbdoc(
              Create coupled PEEC+MMM solver.

              Args:
                  conductor: Conductor object handle
                  magnet: Magnetic object handle (container)

              Returns:
                  Solver handle
          )pbdoc");

    m.def("CplMagSetFrequency", &radia_coupled::CplMagSetFrequency,
          py::arg("solver"), py::arg("freq"),
          "Set analysis frequency [Hz].");

    m.def("CplMagSetVoltage", &radia_coupled::CplMagSetVoltage,
          py::arg("solver"), py::arg("V_real"), py::arg("V_imag"),
          "Set excitation voltage [V] (complex).");

    m.def("CplMagSetCurrent", &radia_coupled::CplMagSetCurrent,
          py::arg("solver"), py::arg("I_real"), py::arg("I_imag"),
          "Set excitation current [A] (complex).");

    m.def("CplMagSetExtField", &radia_coupled::CplMagSetExtField,
          py::arg("solver"), py::arg("Hx"), py::arg("Hy"), py::arg("Hz"),
          "Set external H field [A/m].");

    m.def("CplMagSetMu", &radia_coupled::CplMagSetMu,
          py::arg("solver"), py::arg("mu_r_real"), py::arg("mu_r_imag"),
          "Set complex relative permeability of magnetic material.");

    m.def("CplMagSetSymmetric", &radia_coupled::CplMagSetSymmetric,
          py::arg("solver"), py::arg("use_symmetric"),
          "Enable/disable symmetric matrix assumption.");

    m.def("CplMagSolve", &radia_coupled::CplMagSolve,
          py::arg("solver"),
          py::call_guard<py::gil_scoped_release>(),
          R"pbdoc(
              Solve coupled PEEC+MMM system.

              Returns:
                  Dictionary with Z_real, Z_imag, P_conductor, P_magnet, iterations
          )pbdoc");

    m.def("CplMagImpedance", &radia_coupled::CplMagImpedance,
          py::arg("solver"),
          "Get port impedance [Ohm] (complex).");

    m.def("CplMagPower", &radia_coupled::CplMagPower,
          py::arg("solver"),
          "Get power dissipation as (P_conductor, P_magnet) [W].");

    m.def("CplMagFld", &radia_coupled::CplMagFld,
          py::arg("solver"), py::arg("point"),
          R"pbdoc(
              Compute field at point from coupled solution.

              Args:
                  solver: Solver handle
                  point: Evaluation point [x, y, z]

              Returns:
                  Dictionary with B_real and B_imag tuples
          )pbdoc");

    m.def("CplMagDelete", &radia_coupled::CplMagDelete,
          py::arg("solver"),
          "Delete coupled solver and free resources.");

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

    // FldPtcTrj - Particle trajectory
    m.def("FldPtcTrj", [](int obj, double energy,
                          const py::list& init_cond,
                          const py::list& long_lim,
                          int np) -> py::list {
        auto ic = to_vector(init_cond.cast<py::object>());
        auto ll = to_vector(long_lim.cast<py::object>());
        if (ic.size() != 4) {
            throw std::runtime_error("Initial conditions must be [x0, dxdy0, z0, dzdy0]");
        }
        if (ll.size() != 2) {
            throw std::runtime_error("Longitudinal limits must be [y0, y1]");
        }

        std::vector<double> result(np * 5);  // y, x, dxdy, z, dzdy
        int nf = 0;

        int err = RadFldPtcTrj(result.data(), &nf, obj, energy, ic.data(), ll.data(), np);
        check_error(err);

        py::list out;
        for (int i = 0; i < np; ++i) {
            py::list row;
            for (int j = 0; j < 5; ++j) {
                row.append(result[i * 5 + j]);
            }
            out.append(row);
        }
        return out;
    },
    py::arg("obj"), py::arg("energy"), py::arg("init_cond"),
    py::arg("long_lim"), py::arg("np"),
    R"pbdoc(
        Compute relativistic particle trajectory.

        Args:
            obj: Object handle
            energy: Particle energy [GeV]
            init_cond: [x0, dxdy0, z0, dzdy0] in mm and radians
            long_lim: [y0, y1] longitudinal limits in mm
            np: Number of integration steps

        Returns:
            List of [y, x, dxdy, z, dzdy] at each step
    )pbdoc");

    // FldEnrFrc - Energy-based force
    m.def("FldEnrFrc", [](int obj_dst, int obj_src,
                          const std::string& component = "",
                          const py::object& subdiv = py::none()) -> py::object {
        double f[3] = {0, 0, 0};
        int nf = 0;
        int sbdPar[3] = {0, 0, 0};
        int* pSbdPar = nullptr;

        if (!subdiv.is_none()) {
            auto sd = to_vector(subdiv);
            if (sd.size() >= 3) {
                sbdPar[0] = static_cast<int>(sd[0]);
                sbdPar[1] = static_cast<int>(sd[1]);
                sbdPar[2] = static_cast<int>(sd[2]);
                pSbdPar = sbdPar;
            }
        }

        char id_cstr[16];
        strncpy(id_cstr, component.c_str(), 15);
        id_cstr[15] = '\0';

        int err = RadFldEnrFrc(f, &nf, obj_dst, obj_src, id_cstr, pSbdPar);
        check_error(err);

        if (nf == 1) {
            return py::cast(f[0]);
        } else {
            return py::make_tuple(f[0], f[1], f[2]);
        }
    },
    py::arg("obj_dst"), py::arg("obj_src"),
    py::arg("component") = "", py::arg("subdiv") = py::none(),
    R"pbdoc(
        Compute force on object from field of another object.

        Args:
            obj_dst: Destination object handle
            obj_src: Source object handle
            component: "fx", "fy", "fz", or "" for all
            subdiv: Optional subdivision [kx, ky, kz]

        Returns:
            Force in Newton (single value or tuple)
    )pbdoc");

    // FldEnrTrq - Energy-based torque
    m.def("FldEnrTrq", [](int obj_dst, int obj_src,
                          const std::string& component,
                          const py::list& point,
                          const py::object& subdiv = py::none()) -> py::object {
        auto p = to_vector(point.cast<py::object>());
        if (p.size() != 3) {
            throw std::runtime_error("Point must have 3 coordinates");
        }

        double f[3] = {0, 0, 0};
        int nf = 0;
        int sbdPar[3] = {0, 0, 0};
        int* pSbdPar = nullptr;

        if (!subdiv.is_none()) {
            auto sd = to_vector(subdiv);
            if (sd.size() >= 3) {
                sbdPar[0] = static_cast<int>(sd[0]);
                sbdPar[1] = static_cast<int>(sd[1]);
                sbdPar[2] = static_cast<int>(sd[2]);
                pSbdPar = sbdPar;
            }
        }

        char id_cstr[16];
        strncpy(id_cstr, component.c_str(), 15);
        id_cstr[15] = '\0';

        int err = RadFldEnrTrq(f, &nf, obj_dst, obj_src, id_cstr, p.data(), pSbdPar);
        check_error(err);

        if (nf == 1) {
            return py::cast(f[0]);
        } else {
            return py::make_tuple(f[0], f[1], f[2]);
        }
    },
    py::arg("obj_dst"), py::arg("obj_src"), py::arg("component"),
    py::arg("point"), py::arg("subdiv") = py::none(),
    R"pbdoc(
        Compute torque on object from field of another object.

        Args:
            obj_dst: Destination object handle
            obj_src: Source object handle
            component: "tx", "ty", "tz", or "" for all
            point: Torque reference point [x, y, z]
            subdiv: Optional subdivision [kx, ky, kz]

        Returns:
            Torque in Newton*mm (single value or tuple)
    )pbdoc");

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

    // ObjMltExtTri - Triangulated extruded polygon
    m.def("ObjMltExtTri", [](double xc, double lx,
                              const py::list& vertices,
                              const py::list& subdiv,
                              const std::string& axis = "x",
                              const py::list& magnetization = py::list(),
                              const std::string& opt = "") -> int {
        int nv = static_cast<int>(py::len(vertices));

        std::vector<double> flatVert;
        for (const auto& pt : vertices) {
            auto coord = to_vector(pt.cast<py::object>());
            if (coord.size() != 2) {
                throw std::runtime_error("Each 2D point must have 2 coordinates");
            }
            flatVert.push_back(coord[0]);
            flatVert.push_back(coord[1]);
        }

        std::vector<double> flatSubd;
        for (const auto& sd : subdiv) {
            auto params = to_vector(sd.cast<py::object>());
            if (params.size() != 2) {
                throw std::runtime_error("Each subdiv entry must be [k, q]");
            }
            flatSubd.push_back(params[0]);
            flatSubd.push_back(params[1]);
        }

        double M[3] = {0, 0, 0};
        if (py::len(magnetization) >= 3) {
            auto m = to_vector(magnetization.cast<py::object>());
            M[0] = m[0]; M[1] = m[1]; M[2] = m[2];
        }

        char a = axis.empty() ? 'x' : axis[0];
        char opt_cstr[256];
        strncpy(opt_cstr, opt.c_str(), 255);
        opt_cstr[255] = '\0';

        int n = 0;
        int err = RadObjMltExtTri(&n, xc, lx, flatVert.data(), flatSubd.data(),
                                   nv, a, M, opt_cstr);
        check_error(err);
        return n;
    },
    py::arg("xc"), py::arg("lx"), py::arg("vertices"), py::arg("subdiv"),
    py::arg("axis") = "x", py::arg("magnetization") = py::list(),
    py::arg("opt") = "",
    R"pbdoc(
        Create triangulated extruded polygon.

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

    // ========================================================================
    // Extended Conductor Functions
    // ========================================================================

    // CndHexahedron - Hexahedral conductor
    m.def("CndHexahedron", [](const py::list& vertices,
                               double sigma = 5.8e7,
                               int num_panels = 4) -> int {
        auto [flat, nv] = to_vertex_array(vertices);
        if (nv != 8) {
            throw std::runtime_error("Hexahedron requires exactly 8 vertices");
        }

        int ind = 0;
        int err = RadCndHexahedron(&ind, flat.data(), sigma, num_panels);
        check_error(err);
        return ind;
    },
    py::arg("vertices"), py::arg("sigma") = 5.8e7, py::arg("num_panels") = 4,
    R"pbdoc(
        Create hexahedral conductor.

        Args:
            vertices: 8 vertices [[x1,y1,z1], ..., [x8,y8,z8]]
            sigma: Conductivity [S/m] (default: copper)
            num_panels: Panels per face

        Returns:
            Conductor handle
    )pbdoc");

    // CndWire - Wire conductor
    m.def("CndWire", [](const py::list& points,
                        const std::string& cross_section,
                        double width, double height,
                        double sigma = 5.8e7,
                        int num_around = 8) -> int {
        auto [flat, np] = to_vertex_array(points);

        char cs = cross_section.empty() ? 'r' : cross_section[0];

        int ind = 0;
        int err = RadCndWire(&ind, flat.data(), np, cs, width, height, sigma, num_around);
        check_error(err);
        return ind;
    },
    py::arg("points"), py::arg("cross_section"),
    py::arg("width"), py::arg("height"),
    py::arg("sigma") = 5.8e7, py::arg("num_around") = 8,
    R"pbdoc(
        Create wire conductor along path.

        Args:
            points: Path points [[x,y,z], ...]
            cross_section: "r" (rectangular) or "c" (circular)
            width: Cross-section width
            height: Cross-section height
            sigma: Conductivity [S/m]
            num_around: Panels around cross-section

        Returns:
            Conductor handle
    )pbdoc");

    // CndSpiral - Spiral conductor
    m.def("CndSpiral", [](const py::list& center,
                          double inner_radius, double outer_radius,
                          double pitch, int num_turns,
                          const py::list& axis,
                          const std::string& cross_section,
                          double width, double height,
                          double sigma = 5.8e7,
                          int num_around = 8) -> int {
        auto c = to_vector(center.cast<py::object>());
        auto a = to_vector(axis.cast<py::object>());
        if (c.size() != 3 || a.size() != 3) {
            throw std::runtime_error("Center and axis must have 3 components");
        }

        char cs = cross_section.empty() ? 'r' : cross_section[0];

        int ind = 0;
        int err = RadCndSpiral(&ind, c.data(), inner_radius, outer_radius,
                               pitch, num_turns, a.data(), cs, width, height,
                               sigma, num_around);
        check_error(err);
        return ind;
    },
    py::arg("center"), py::arg("inner_radius"), py::arg("outer_radius"),
    py::arg("pitch"), py::arg("num_turns"), py::arg("axis"),
    py::arg("cross_section"), py::arg("width"), py::arg("height"),
    py::arg("sigma") = 5.8e7, py::arg("num_around") = 8,
    R"pbdoc(
        Create spiral conductor.

        Args:
            center: Center point [x, y, z]
            inner_radius: Inner radius
            outer_radius: Outer radius
            pitch: Height per turn
            num_turns: Number of turns
            axis: Spiral axis [ax, ay, az]
            cross_section: "r" or "c"
            width, height: Cross-section dimensions
            sigma: Conductivity [S/m]
            num_around: Panels around cross-section

        Returns:
            Conductor handle
    )pbdoc");

    // CndDefinePort - Define port terminals
    m.def("CndDefinePort", [](int cond,
                               const py::list& t1_panels,
                               const py::list& t2_panels) {
        std::vector<int> t1, t2;
        for (const auto& p : t1_panels) t1.push_back(p.cast<int>());
        for (const auto& p : t2_panels) t2.push_back(p.cast<int>());

        int err = RadCndDefinePort(cond, t1.data(), static_cast<int>(t1.size()),
                                    t2.data(), static_cast<int>(t2.size()));
        check_error(err);
    },
    py::arg("cond"), py::arg("t1_panels"), py::arg("t2_panels"),
    R"pbdoc(
        Define port terminals for impedance calculation.

        Args:
            cond: Conductor handle
            t1_panels: Panel indices for terminal 1
            t2_panels: Panel indices for terminal 2
    )pbdoc");

    // CndImpedanceSweep - Frequency sweep
    m.def("CndImpedanceSweep", [](int cond,
                                   const py::list& frequencies) -> py::list {
        std::vector<double> freqs;
        for (const auto& f : frequencies) freqs.push_back(f.cast<double>());

        int nf = static_cast<int>(freqs.size());
        std::vector<double> Z_real(nf), Z_imag(nf);

        int err = RadCndImpedanceSweep(Z_real.data(), Z_imag.data(),
                                        cond, freqs.data(), nf);
        check_error(err);

        py::list out;
        for (int i = 0; i < nf; ++i) {
            out.append(std::complex<double>(Z_real[i], Z_imag[i]));
        }
        return out;
    },
    py::arg("cond"), py::arg("frequencies"),
    R"pbdoc(
        Compute impedance at multiple frequencies.

        Args:
            cond: Conductor handle
            frequencies: List of frequencies [Hz]

        Returns:
            List of complex impedances [Ohm]
    )pbdoc");

    // ========================================================================
    // Extended Coupled Solver Functions
    // ========================================================================

    // CplMagSetConductor - Set conductor for existing solver
    m.def("CplMagSetConductor", [](int solver, int conductor) {
        int err = RadCplMagSetConductor(solver, conductor);
        check_error(err);
    },
    py::arg("solver"), py::arg("conductor"),
    "Set conductor for coupled solver.");

    // CplMagSweep - Frequency sweep for coupled solver
    m.def("CplMagSweep", [](int solver,
                             const py::list& frequencies) -> py::list {
        std::vector<double> freqs;
        for (const auto& f : frequencies) freqs.push_back(f.cast<double>());

        int nf = static_cast<int>(freqs.size());
        std::vector<double> pOut(nf * 2);  // Real and imag interleaved
        int nOut = 0;

        int err = RadCplMagSweep(pOut.data(), &nOut, solver, freqs.data(), nf);
        check_error(err);

        py::list out;
        for (int i = 0; i < nf; ++i) {
            out.append(std::complex<double>(pOut[2*i], pOut[2*i+1]));
        }
        return out;
    },
    py::arg("solver"), py::arg("frequencies"),
    R"pbdoc(
        Compute coupled impedance at multiple frequencies.

        Args:
            solver: Coupled solver handle
            frequencies: List of frequencies [Hz]

        Returns:
            List of complex impedances [Ohm]
    )pbdoc");

    // ========================================================================
    // Extended Utility Functions
    // ========================================================================

    m.def("UtiVer", &radia_utility_ext::UtiVer,
          "Get Radia library version number.");

    // UtiDmp - Dump object(s) to string
    m.def("UtiDmp", [](const py::object& elem,
                       const std::string& format = "asc") -> py::bytes {
        std::vector<int> objs;

        if (py::isinstance<py::int_>(elem)) {
            objs.push_back(elem.cast<int>());
        } else if (py::isinstance<py::list>(elem)) {
            for (const auto& e : elem.cast<py::list>()) {
                objs.push_back(e.cast<int>());
            }
        } else {
            throw std::runtime_error("elem must be int or list of ints");
        }

        char fmt[4];
        strncpy(fmt, format.c_str(), 3);
        fmt[3] = '\0';

        // First call to get size
        int size = 0;
        int err = RadUtiDmp(nullptr, &size, objs.data(),
                            static_cast<int>(objs.size()), fmt);
        check_error(err);

        // Allocate and get data
        std::vector<char> buffer(size + 1);
        err = RadUtiDmp(buffer.data(), &size, objs.data(),
                        static_cast<int>(objs.size()), fmt);
        check_error(err);

        return py::bytes(buffer.data(), size);
    },
    py::arg("elem"), py::arg("format") = "asc",
    R"pbdoc(
        Dump object(s) to byte string.

        Args:
            elem: Object handle or list of handles
            format: "asc" (ASCII) or "bin" (binary)

        Returns:
            Byte string with object data
    )pbdoc");

    // UtiDmpPrs - Parse dumped data
    m.def("UtiDmpPrs", [](const py::bytes& data) -> py::object {
        std::string s = data;
        std::vector<int> elems(1000);  // Max elements
        int nElem = 0;

        int err = RadUtiDmpPrs(elems.data(), &nElem,
                               reinterpret_cast<unsigned char*>(const_cast<char*>(s.data())),
                               static_cast<int>(s.size()));
        check_error(err);

        if (nElem == 1) {
            return py::cast(elems[0]);
        } else {
            py::list out;
            for (int i = 0; i < nElem; ++i) {
                out.append(elems[i]);
            }
            return out;
        }
    },
    py::arg("data"),
    R"pbdoc(
        Parse dumped byte string and recreate object(s).

        Args:
            data: Byte string from UtiDmp

        Returns:
            Object handle or list of handles
    )pbdoc");
}
