/*
 * rad_mmm_matrices_api.cpp
 *
 * pybind11 Python bindings for MMM (Magnetic Moment Method) matrix construction
 *
 * Module: mmm_core
 *
 * Designed for NGSolve mesh integration.
 * Style follows rad_peec_matrices_api.cpp and rad_cln_api.cpp patterns.
 *
 * Part of Radia project
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "../core/rad_mmm_matrices.h"

namespace py = pybind11;

using namespace radia::mmm;

//=========================================================================
// Python Wrapper Classes
//=========================================================================

/**
 * Python-friendly MMM Matrix Builder
 *
 * Wraps MMMBuilder with numpy array interface.
 */
class PyMMMBuilder {
public:
    PyMMMBuilder() : builder_(std::make_unique<MMMBuilder>()) {}

    /**
     * Add tetrahedron element
     * @param vertices: (4, 3) array of vertex coordinates [m]
     */
    void add_tetrahedron(py::array_t<double, py::array::c_style | py::array::forcecast> vertices) {
        auto v = vertices.unchecked<2>();
        if (v.shape(0) != 4 || v.shape(1) != 3) {
            throw std::runtime_error("Tetrahedron vertices must be shape (4, 3)");
        }

        std::vector<Vec3d> verts(4);
        for (int i = 0; i < 4; i++) {
            verts[i] = Vec3d(v(i, 0), v(i, 1), v(i, 2));
        }
        builder_->AddTetrahedron(verts);
    }

    /**
     * Add hexahedron element (MSC 6 DOF)
     * @param vertices: (8, 3) array of vertex coordinates [m]
     */
    void add_hexahedron(py::array_t<double, py::array::c_style | py::array::forcecast> vertices) {
        auto v = vertices.unchecked<2>();
        if (v.shape(0) != 8 || v.shape(1) != 3) {
            throw std::runtime_error("Hexahedron vertices must be shape (8, 3)");
        }

        std::vector<Vec3d> verts(8);
        for (int i = 0; i < 8; i++) {
            verts[i] = Vec3d(v(i, 0), v(i, 1), v(i, 2));
        }
        builder_->AddHexahedron(verts);
    }

    /**
     * Add multiple tetrahedra from NGSolve mesh
     * @param vertices: (n_verts, 3) vertex coordinates
     * @param elements: (n_elem, 4) vertex indices (0-based)
     */
    void add_tetrahedra_from_mesh(
            py::array_t<double, py::array::c_style | py::array::forcecast> vertices,
            py::array_t<int, py::array::c_style | py::array::forcecast> elements) {
        auto verts = vertices.unchecked<2>();
        auto elems = elements.unchecked<2>();

        if (verts.shape(1) != 3) {
            throw std::runtime_error("Vertices must have shape (n, 3)");
        }
        if (elems.shape(1) != 4) {
            throw std::runtime_error("Tetrahedron elements must have shape (n, 4)");
        }

        int n_elem = static_cast<int>(elems.shape(0));
        for (int e = 0; e < n_elem; e++) {
            std::vector<Vec3d> tet_verts(4);
            for (int i = 0; i < 4; i++) {
                int vi = elems(e, i);
                tet_verts[i] = Vec3d(verts(vi, 0), verts(vi, 1), verts(vi, 2));
            }
            builder_->AddTetrahedron(tet_verts);
        }
    }

    /**
     * Add multiple hexahedra from mesh
     * @param vertices: (n_verts, 3) vertex coordinates
     * @param elements: (n_elem, 8) vertex indices (0-based)
     */
    void add_hexahedra_from_mesh(
            py::array_t<double, py::array::c_style | py::array::forcecast> vertices,
            py::array_t<int, py::array::c_style | py::array::forcecast> elements) {
        auto verts = vertices.unchecked<2>();
        auto elems = elements.unchecked<2>();

        if (verts.shape(1) != 3) {
            throw std::runtime_error("Vertices must have shape (n, 3)");
        }
        if (elems.shape(1) != 8) {
            throw std::runtime_error("Hexahedron elements must have shape (n, 8)");
        }

        int n_elem = static_cast<int>(elems.shape(0));
        for (int e = 0; e < n_elem; e++) {
            std::vector<Vec3d> hex_verts(8);
            for (int i = 0; i < 8; i++) {
                int vi = elems(e, i);
                hex_verts[i] = Vec3d(verts(vi, 0), verts(vi, 1), verts(vi, 2));
            }
            builder_->AddHexahedron(hex_verts);
        }
    }

    /**
     * Add permanent magnet tetrahedron (fixed magnetization, 0 DOF)
     * @param vertices: (4, 3) array of vertex coordinates [m]
     * @param Mr: (3,) remanent magnetization [A/m]
     */
    void add_pm_tetrahedron(
            py::array_t<double, py::array::c_style | py::array::forcecast> vertices,
            py::array_t<double, py::array::c_style | py::array::forcecast> Mr) {
        auto v = vertices.unchecked<2>();
        auto m = Mr.unchecked<1>();

        if (v.shape(0) != 4 || v.shape(1) != 3) {
            throw std::runtime_error("Tetrahedron vertices must be shape (4, 3)");
        }
        if (m.shape(0) != 3) {
            throw std::runtime_error("Mr must be shape (3,)");
        }

        std::vector<Vec3d> verts(4);
        for (int i = 0; i < 4; i++) {
            verts[i] = Vec3d(v(i, 0), v(i, 1), v(i, 2));
        }
        Vec3d mr(m(0), m(1), m(2));
        builder_->AddPermanentMagnetTetrahedron(verts, mr);
    }

    /**
     * Add permanent magnet hexahedron (fixed magnetization, 0 DOF)
     * @param vertices: (8, 3) array of vertex coordinates [m]
     * @param Mr: (3,) remanent magnetization [A/m]
     */
    void add_pm_hexahedron(
            py::array_t<double, py::array::c_style | py::array::forcecast> vertices,
            py::array_t<double, py::array::c_style | py::array::forcecast> Mr) {
        auto v = vertices.unchecked<2>();
        auto m = Mr.unchecked<1>();

        if (v.shape(0) != 8 || v.shape(1) != 3) {
            throw std::runtime_error("Hexahedron vertices must be shape (8, 3)");
        }
        if (m.shape(0) != 3) {
            throw std::runtime_error("Mr must be shape (3,)");
        }

        std::vector<Vec3d> verts(8);
        for (int i = 0; i < 8; i++) {
            verts[i] = Vec3d(v(i, 0), v(i, 1), v(i, 2));
        }
        Vec3d mr(m(0), m(1), m(2));
        builder_->AddPermanentMagnetHexahedron(verts, mr);
    }

    /**
     * Add multiple PM tetrahedra from mesh
     * @param vertices: (n_verts, 3) vertex coordinates
     * @param elements: (n_elem, 4) vertex indices (0-based)
     * @param Mr: (3,) uniform remanent magnetization [A/m]
     */
    void add_pm_tetrahedra_from_mesh(
            py::array_t<double, py::array::c_style | py::array::forcecast> vertices,
            py::array_t<int, py::array::c_style | py::array::forcecast> elements,
            py::array_t<double, py::array::c_style | py::array::forcecast> Mr) {
        auto verts = vertices.unchecked<2>();
        auto elems = elements.unchecked<2>();
        auto m = Mr.unchecked<1>();

        if (verts.shape(1) != 3) {
            throw std::runtime_error("Vertices must have shape (n, 3)");
        }
        if (elems.shape(1) != 4) {
            throw std::runtime_error("Tetrahedron elements must have shape (n, 4)");
        }
        if (m.shape(0) != 3) {
            throw std::runtime_error("Mr must be shape (3,)");
        }

        Vec3d mr(m(0), m(1), m(2));
        int n_elem = static_cast<int>(elems.shape(0));
        for (int e = 0; e < n_elem; e++) {
            std::vector<Vec3d> tet_verts(4);
            for (int i = 0; i < 4; i++) {
                int vi = elems(e, i);
                tet_verts[i] = Vec3d(verts(vi, 0), verts(vi, 1), verts(vi, 2));
            }
            builder_->AddPermanentMagnetTetrahedron(tet_verts, mr);
        }
    }

    /**
     * Add multiple PM hexahedra from mesh
     * @param vertices: (n_verts, 3) vertex coordinates
     * @param elements: (n_elem, 8) vertex indices (0-based)
     * @param Mr: (3,) uniform remanent magnetization [A/m]
     */
    void add_pm_hexahedra_from_mesh(
            py::array_t<double, py::array::c_style | py::array::forcecast> vertices,
            py::array_t<int, py::array::c_style | py::array::forcecast> elements,
            py::array_t<double, py::array::c_style | py::array::forcecast> Mr) {
        auto verts = vertices.unchecked<2>();
        auto elems = elements.unchecked<2>();
        auto m = Mr.unchecked<1>();

        if (verts.shape(1) != 3) {
            throw std::runtime_error("Vertices must have shape (n, 3)");
        }
        if (elems.shape(1) != 8) {
            throw std::runtime_error("Hexahedron elements must have shape (n, 8)");
        }
        if (m.shape(0) != 3) {
            throw std::runtime_error("Mr must be shape (3,)");
        }

        Vec3d mr(m(0), m(1), m(2));
        int n_elem = static_cast<int>(elems.shape(0));
        for (int e = 0; e < n_elem; e++) {
            std::vector<Vec3d> hex_verts(8);
            for (int i = 0; i < 8; i++) {
                int vi = elems(e, i);
                hex_verts[i] = Vec3d(verts(vi, 0), verts(vi, 1), verts(vi, 2));
            }
            builder_->AddPermanentMagnetHexahedron(hex_verts, mr);
        }
    }

    /**
     * Build interaction matrix
     * @return Tuple (N, dof_offset)
     *   N: (total_dof, total_dof) interaction matrix
     *   dof_offset: (n_elem + 1,) DOF offset array
     */
    py::tuple build() {
        MMMMatrices matrices = builder_->Build();

        int n = matrices.total_dof;
        int n_elem = matrices.n_elements;

        // N matrix (convert from column-major to row-major for NumPy)
        py::array_t<double> N({n, n});
        auto N_buf = N.mutable_unchecked<2>();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                // Column-major in C++: matrices.N[j * n + i]
                N_buf(i, j) = matrices.N[static_cast<size_t>(j) * n + i];
            }
        }

        // DOF offset array
        py::array_t<int> dof_offset(n_elem + 1);
        auto dof_buf = dof_offset.mutable_unchecked<1>();
        for (int i = 0; i <= n_elem; i++) {
            dof_buf(i) = matrices.dof_offset[i];
        }

        return py::make_tuple(N, dof_offset);
    }

    // Properties
    int num_elements() const { return builder_->NumElements(); }
    int total_dof() const { return builder_->TotalDOF(); }

    void clear() { builder_->Clear(); }

    // Accessor for elements (for MMMFieldComputer)
    const std::vector<MMMElement>& get_elements() const { return builder_->GetElements(); }

private:
    std::unique_ptr<MMMBuilder> builder_;
};

/**
 * Python-friendly MMM Solver
 */
class PyMMMSolver {
public:
    PyMMMSolver() : solver_(std::make_unique<MMMSolver>()) {}

    /**
     * Set interaction matrix from numpy arrays
     */
    void set_matrix(py::array_t<double, py::array::c_style | py::array::forcecast> N,
                    py::array_t<int, py::array::c_style | py::array::forcecast> dof_offset) {
        auto N_info = N.request();
        auto dof_info = dof_offset.request();

        if (N_info.ndim != 2 || N_info.shape[0] != N_info.shape[1]) {
            throw std::runtime_error("N must be a square 2D array");
        }

        int total_dof = static_cast<int>(N_info.shape[0]);
        int n_elem = static_cast<int>(dof_info.shape[0]) - 1;

        // Convert from row-major (NumPy) to column-major (LAPACK)
        std::vector<double> N_col(total_dof * total_dof);
        auto N_ptr = static_cast<const double*>(N_info.ptr);
        for (int i = 0; i < total_dof; i++) {
            for (int j = 0; j < total_dof; j++) {
                // Row-major: N_ptr[i * total_dof + j]
                // Column-major: N_col[j * total_dof + i]
                N_col[static_cast<size_t>(j) * total_dof + i] = N_ptr[static_cast<size_t>(i) * total_dof + j];
            }
        }

        solver_->SetMatrix(N_col.data(),
                           static_cast<const int*>(dof_info.ptr),
                           total_dof, n_elem);
    }

    /**
     * Solve using LU decomposition
     * @param inv_chi: 1/chi values (total_dof,) or (n_elem,)
     * @param H_ext: External field (total_dof,)
     * @param chi_per_element: If True, inv_chi is (n_elem,) and broadcast
     * @return Magnetization M (total_dof,)
     */
    py::array_t<double> solve_lu(
            py::array_t<double, py::array::c_style | py::array::forcecast> inv_chi,
            py::array_t<double, py::array::c_style | py::array::forcecast> H_ext,
            bool chi_per_element = false) {

        auto inv_chi_info = inv_chi.request();
        auto H_ext_info = H_ext.request();

        int n = solver_->TotalDOF();
        int n_elem = solver_->NumElements();

        // Validate sizes
        if (chi_per_element) {
            if (inv_chi_info.size != n_elem) {
                throw std::runtime_error("inv_chi size must equal n_elements when chi_per_element=True");
            }
        } else {
            if (inv_chi_info.size != n) {
                throw std::runtime_error("inv_chi size must equal total_dof");
            }
        }

        if (H_ext_info.size != n) {
            throw std::runtime_error("H_ext size must equal total_dof");
        }

        std::vector<double> inv_chi_vec(static_cast<const double*>(inv_chi_info.ptr),
                                         static_cast<const double*>(inv_chi_info.ptr) + inv_chi_info.size);
        std::vector<double> H_ext_vec(static_cast<const double*>(H_ext_info.ptr),
                                       static_cast<const double*>(H_ext_info.ptr) + H_ext_info.size);
        std::vector<double> M;

        int info = solver_->SolveLU(inv_chi_vec, H_ext_vec, M, chi_per_element);
        if (info != 0) {
            throw std::runtime_error("LU solve failed with info=" + std::to_string(info));
        }

        return py::array_t<double>(M.size(), M.data());
    }

    /**
     * Solve using BiCGSTAB
     * @return Tuple (M, n_iter) - Magnetization and iteration count
     */
    py::tuple solve_bicgstab(
            py::array_t<double, py::array::c_style | py::array::forcecast> inv_chi,
            py::array_t<double, py::array::c_style | py::array::forcecast> H_ext,
            double tol = 1e-8,
            int max_iter = 1000,
            bool chi_per_element = false) {

        auto inv_chi_info = inv_chi.request();
        auto H_ext_info = H_ext.request();

        int n = solver_->TotalDOF();
        int n_elem = solver_->NumElements();

        if (chi_per_element) {
            if (inv_chi_info.size != n_elem) {
                throw std::runtime_error("inv_chi size must equal n_elements when chi_per_element=True");
            }
        } else {
            if (inv_chi_info.size != n) {
                throw std::runtime_error("inv_chi size must equal total_dof");
            }
        }

        if (H_ext_info.size != n) {
            throw std::runtime_error("H_ext size must equal total_dof");
        }

        std::vector<double> inv_chi_vec(static_cast<const double*>(inv_chi_info.ptr),
                                         static_cast<const double*>(inv_chi_info.ptr) + inv_chi_info.size);
        std::vector<double> H_ext_vec(static_cast<const double*>(H_ext_info.ptr),
                                       static_cast<const double*>(H_ext_info.ptr) + H_ext_info.size);
        std::vector<double> M;

        int n_iter = solver_->SolveBiCGSTAB(inv_chi_vec, H_ext_vec, M, tol, max_iter, chi_per_element);

        return py::make_tuple(
            py::array_t<double>(M.size(), M.data()),
            n_iter
        );
    }

    /**
     * Matrix-vector product: y = A * x where A = -diag(inv_chi) - N
     * For use in external iterative solvers (e.g., Schur complement coupling)
     */
    py::array_t<double> matvec(
            py::array_t<double, py::array::c_style | py::array::forcecast> inv_chi,
            py::array_t<double, py::array::c_style | py::array::forcecast> x,
            bool chi_per_element = false) {

        auto inv_chi_info = inv_chi.request();
        auto x_info = x.request();

        int n = solver_->TotalDOF();
        int n_elem = solver_->NumElements();

        if (chi_per_element) {
            if (inv_chi_info.size != n_elem) {
                throw std::runtime_error("inv_chi size must equal n_elements when chi_per_element=True");
            }
        } else {
            if (inv_chi_info.size != n) {
                throw std::runtime_error("inv_chi size must equal total_dof");
            }
        }

        if (x_info.size != n) {
            throw std::runtime_error("x size must equal total_dof");
        }

        std::vector<double> inv_chi_vec(static_cast<const double*>(inv_chi_info.ptr),
                                         static_cast<const double*>(inv_chi_info.ptr) + inv_chi_info.size);
        std::vector<double> x_vec(static_cast<const double*>(x_info.ptr),
                                   static_cast<const double*>(x_info.ptr) + x_info.size);
        std::vector<double> y;

        solver_->ApplySystemMatrix(inv_chi_vec, x_vec, y, chi_per_element);

        return py::array_t<double>(y.size(), y.data());
    }

    // Properties
    int total_dof() const { return solver_->TotalDOF(); }
    int num_elements() const { return solver_->NumElements(); }
    bool is_ready() const { return solver_->IsReady(); }

private:
    std::unique_ptr<MMMSolver> solver_;
};

/**
 * Python-friendly MMM Field Computer
 *
 * Computes B/H fields at observation points from solved magnetization.
 */
class PyMMMFieldComputer {
public:
    PyMMMFieldComputer() : computer_(std::make_unique<MMMFieldComputer>()) {}

    /**
     * Set elements from PyMMMBuilder
     */
    void set_elements_from_builder(const PyMMMBuilder& builder) {
        computer_->SetElements(builder.get_elements());
    }

    /**
     * Compute B field at observation points
     * @param M: Magnetization vector (total_dof,)
     * @param obs_points: Observation points (n_points, 3) [m]
     * @return B field (n_points, 3) [T]
     */
    py::array_t<double> compute_b_field(
            py::array_t<double, py::array::c_style | py::array::forcecast> M,
            py::array_t<double, py::array::c_style | py::array::forcecast> obs_points) {

        auto M_info = M.request();
        auto obs_info = obs_points.request();

        if (obs_info.ndim != 2 || obs_info.shape[1] != 3) {
            throw std::runtime_error("obs_points must have shape (n, 3)");
        }

        int n_points = static_cast<int>(obs_info.shape[0]);

        std::vector<double> M_vec(static_cast<const double*>(M_info.ptr),
                                   static_cast<const double*>(M_info.ptr) + M_info.size);
        std::vector<double> obs_vec(static_cast<const double*>(obs_info.ptr),
                                     static_cast<const double*>(obs_info.ptr) + obs_info.size);

        std::vector<double> B = computer_->ComputeBField(M_vec, obs_vec, n_points);

        // Reshape to (n_points, 3)
        py::array_t<double> result({n_points, 3});
        auto buf = result.mutable_unchecked<2>();
        for (int i = 0; i < n_points; i++) {
            buf(i, 0) = B[i * 3];
            buf(i, 1) = B[i * 3 + 1];
            buf(i, 2) = B[i * 3 + 2];
        }
        return result;
    }

    /**
     * Compute H field at observation points
     * @param M: Magnetization vector (total_dof,)
     * @param obs_points: Observation points (n_points, 3) [m]
     * @return H field (n_points, 3) [A/m]
     */
    py::array_t<double> compute_h_field(
            py::array_t<double, py::array::c_style | py::array::forcecast> M,
            py::array_t<double, py::array::c_style | py::array::forcecast> obs_points) {

        auto M_info = M.request();
        auto obs_info = obs_points.request();

        if (obs_info.ndim != 2 || obs_info.shape[1] != 3) {
            throw std::runtime_error("obs_points must have shape (n, 3)");
        }

        int n_points = static_cast<int>(obs_info.shape[0]);

        std::vector<double> M_vec(static_cast<const double*>(M_info.ptr),
                                   static_cast<const double*>(M_info.ptr) + M_info.size);
        std::vector<double> obs_vec(static_cast<const double*>(obs_info.ptr),
                                     static_cast<const double*>(obs_info.ptr) + obs_info.size);

        std::vector<double> H = computer_->ComputeHField(M_vec, obs_vec, n_points);

        // Reshape to (n_points, 3)
        py::array_t<double> result({n_points, 3});
        auto buf = result.mutable_unchecked<2>();
        for (int i = 0; i < n_points; i++) {
            buf(i, 0) = H[i * 3];
            buf(i, 1) = H[i * 3 + 1];
            buf(i, 2) = H[i * 3 + 2];
        }
        return result;
    }

    /**
     * Compute vector potential A at observation points
     * @param M: Magnetization vector (total_dof,)
     * @param obs_points: Observation points (n_points, 3) [m]
     * @return A field (n_points, 3) [T*m]
     *
     * Uses magnetic dipole formula. For hexahedra, reconstructs M
     * from surface charges via least-squares on face normals.
     */
    py::array_t<double> compute_a_field(
            py::array_t<double, py::array::c_style | py::array::forcecast> M,
            py::array_t<double, py::array::c_style | py::array::forcecast> obs_points) {

        auto M_info = M.request();
        auto obs_info = obs_points.request();

        if (obs_info.ndim != 2 || obs_info.shape[1] != 3) {
            throw std::runtime_error("obs_points must have shape (n, 3)");
        }

        int n_points = static_cast<int>(obs_info.shape[0]);

        std::vector<double> M_vec(static_cast<const double*>(M_info.ptr),
                                   static_cast<const double*>(M_info.ptr) + M_info.size);
        std::vector<double> obs_vec(static_cast<const double*>(obs_info.ptr),
                                     static_cast<const double*>(obs_info.ptr) + obs_info.size);

        std::vector<double> A = computer_->ComputeAField(M_vec, obs_vec, n_points);

        py::array_t<double> result({n_points, 3});
        auto buf = result.mutable_unchecked<2>();
        for (int i = 0; i < n_points; i++) {
            buf(i, 0) = A[i * 3];
            buf(i, 1) = A[i * 3 + 1];
            buf(i, 2) = A[i * 3 + 2];
        }
        return result;
    }

    /**
     * Get element centers
     * @return Element centers (n_elements, 3) [m]
     */
    py::array_t<double> get_element_centers() {
        std::vector<double> centers = computer_->GetElementCenters();
        int n_elem = computer_->NumElements();

        py::array_t<double> result({n_elem, 3});
        auto buf = result.mutable_unchecked<2>();
        for (int i = 0; i < n_elem; i++) {
            buf(i, 0) = centers[i * 3];
            buf(i, 1) = centers[i * 3 + 1];
            buf(i, 2) = centers[i * 3 + 2];
        }
        return result;
    }

    /**
     * Get element volumes
     * @return Element volumes (n_elements,) [m^3]
     */
    py::array_t<double> get_element_volumes() {
        std::vector<double> volumes = computer_->GetElementVolumes();
        return py::array_t<double>(volumes.size(), volumes.data());
    }

    // Properties
    int num_elements() const { return computer_->NumElements(); }
    int total_dof() const { return computer_->TotalDOF(); }

private:
    std::unique_ptr<MMMFieldComputer> computer_;
};

/**
 * Python-friendly Nonlinear MMM Solver
 *
 * Newton-Raphson iteration for nonlinear B-H curve materials.
 */
class PyMMMNonlinearSolver {
public:
    PyMMMNonlinearSolver() : solver_(std::make_unique<MMMNonlinearSolver>()) {}

    /**
     * Set interaction matrix from numpy arrays
     */
    void set_matrix(py::array_t<double, py::array::c_style | py::array::forcecast> N,
                    py::array_t<int, py::array::c_style | py::array::forcecast> dof_offset) {
        auto N_info = N.request();
        auto dof_info = dof_offset.request();

        if (N_info.ndim != 2 || N_info.shape[0] != N_info.shape[1]) {
            throw std::runtime_error("N must be a square 2D array");
        }

        int total_dof = static_cast<int>(N_info.shape[0]);
        int n_elem = static_cast<int>(dof_info.shape[0]) - 1;

        // Convert from row-major (NumPy) to column-major (LAPACK)
        std::vector<double> N_col(total_dof * total_dof);
        auto N_ptr = static_cast<const double*>(N_info.ptr);
        for (int i = 0; i < total_dof; i++) {
            for (int j = 0; j < total_dof; j++) {
                N_col[static_cast<size_t>(j) * total_dof + i] = N_ptr[static_cast<size_t>(i) * total_dof + j];
            }
        }

        solver_->SetMatrix(N_col.data(),
                           static_cast<const int*>(dof_info.ptr),
                           total_dof, n_elem);
    }

    /**
     * Set B-H curve for nonlinear material
     * @param bh_curve: (n_points, 2) array [[H0,B0], [H1,B1], ...]
     */
    void set_bh_curve(py::array_t<double, py::array::c_style | py::array::forcecast> bh_curve) {
        auto bh_info = bh_curve.request();

        if (bh_info.ndim != 2 || bh_info.shape[1] != 2) {
            throw std::runtime_error("bh_curve must have shape (n, 2)");
        }

        int n_pts = static_cast<int>(bh_info.shape[0]);
        auto bh_ptr = static_cast<const double*>(bh_info.ptr);

        std::vector<std::array<double, 2>> bh_vec(n_pts);
        for (int i = 0; i < n_pts; i++) {
            bh_vec[i][0] = bh_ptr[i * 2];      // H
            bh_vec[i][1] = bh_ptr[i * 2 + 1];  // B
        }

        solver_->SetBHCurve(bh_vec);
    }

    /**
     * Set solver parameters
     */
    void set_parameters(double tol = 1e-4, int max_iter = 100,
                        double relax = 1.0, int linear_solver = 0) {
        solver_->SetParameters(tol, max_iter, relax, linear_solver);
    }

    /**
     * Set initial chi guess (optional)
     */
    void set_initial_chi(py::array_t<double, py::array::c_style | py::array::forcecast> chi_init) {
        auto chi_info = chi_init.request();
        std::vector<double> chi_vec(static_cast<const double*>(chi_info.ptr),
                                     static_cast<const double*>(chi_info.ptr) + chi_info.size);
        solver_->SetInitialChi(chi_vec);
    }

    /**
     * Solve nonlinear system
     * @param H_ext: External field (3,) for uniform or (total_dof,)
     * @return dict with M, chi, iterations, final_error, converged
     */
    py::dict solve(py::array_t<double, py::array::c_style | py::array::forcecast> H_ext) {
        auto H_info = H_ext.request();
        std::vector<double> H_vec(static_cast<const double*>(H_info.ptr),
                                   static_cast<const double*>(H_info.ptr) + H_info.size);

        MMMNonlinearResult result = solver_->Solve(H_vec);

        py::dict ret;
        ret["M"] = py::array_t<double>(result.M.size(), result.M.data());
        ret["chi"] = py::array_t<double>(result.chi.size(), result.chi.data());
        ret["iterations"] = result.iterations;
        ret["final_error"] = result.final_error;
        ret["converged"] = result.converged;

        return ret;
    }

    /**
     * Get local H field after solve
     */
    py::array_t<double> get_local_h() {
        std::vector<double> H = solver_->GetLocalH();
        int n_elem = solver_->NumElements();

        py::array_t<double> result({n_elem, 3});
        auto buf = result.mutable_unchecked<2>();
        for (int i = 0; i < n_elem; i++) {
            buf(i, 0) = H[i * 3];
            buf(i, 1) = H[i * 3 + 1];
            buf(i, 2) = H[i * 3 + 2];
        }
        return result;
    }

    /**
     * Get local B field after solve
     */
    py::array_t<double> get_local_b() {
        std::vector<double> B = solver_->GetLocalB();
        int n_elem = solver_->NumElements();

        py::array_t<double> result({n_elem, 3});
        auto buf = result.mutable_unchecked<2>();
        for (int i = 0; i < n_elem; i++) {
            buf(i, 0) = B[i * 3];
            buf(i, 1) = B[i * 3 + 1];
            buf(i, 2) = B[i * 3 + 2];
        }
        return result;
    }

    // Properties
    int total_dof() const { return solver_->TotalDOF(); }
    int num_elements() const { return solver_->NumElements(); }

private:
    std::unique_ptr<MMMNonlinearSolver> solver_;
};

//=========================================================================
// pybind11 Module Definition
//=========================================================================

PYBIND11_MODULE(mmm_core, m) {
    m.doc() = R"doc(
MMM (Magnetic Moment Method) Core Module

Standalone matrix construction for magnetostatic analysis.
Designed for NGSolve mesh integration.

This module provides Radia-object-independent matrix construction,
allowing direct use with NGSolve meshes and external solvers.

Features:
- Tetrahedron elements (3 DOF: Mx, My, Mz)
- Hexahedron elements (6 DOF MSC: sigma per face)
- Mixed mesh support
- LU and BiCGSTAB solvers
- Nonlinear material helpers (B-H curve)

Example:
    from mmm_core import MMMBuilder, MMMSolver
    import numpy as np

    # Create builder
    builder = MMMBuilder()

    # Add tetrahedron (4 vertices)
    vertices = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0.5, 0.866, 0],
        [0.5, 0.289, 0.816]
    ])
    builder.add_tetrahedron(vertices)

    # Build interaction matrix
    N, dof_offset = builder.build()

    # Solve
    solver = MMMSolver()
    solver.set_matrix(N, dof_offset)

    mu_r = 1000
    inv_chi = np.full(3, 1.0 / (mu_r - 1))
    H_ext = np.array([0, 0, 1e5])

    M = solver.solve_lu(inv_chi, H_ext)
)doc";

    // MMMBuilder class
    py::class_<PyMMMBuilder>(m, "MMMBuilder",
        R"doc(
        MMM Interaction Matrix Builder

        Constructs interaction matrix from element geometry.
        Supports tetrahedra (3 DOF) and hexahedra (6 DOF MSC).
        )doc")
        .def(py::init<>(), "Create a new MMM matrix builder")

        .def("add_tetrahedron", &PyMMMBuilder::add_tetrahedron,
             py::arg("vertices"),
             R"doc(
             Add tetrahedron element.

             Args:
                 vertices: (4, 3) array of vertex coordinates [m]
                           Standard tetrahedron vertex ordering
             )doc")

        .def("add_hexahedron", &PyMMMBuilder::add_hexahedron,
             py::arg("vertices"),
             R"doc(
             Add hexahedron element (MSC 6 DOF).

             Args:
                 vertices: (8, 3) array of vertex coordinates [m]
                           Standard hex ordering: [0-3: bottom CCW, 4-7: top CCW]
             )doc")

        .def("add_tetrahedra_from_mesh", &PyMMMBuilder::add_tetrahedra_from_mesh,
             py::arg("vertices"), py::arg("elements"),
             R"doc(
             Add multiple tetrahedra from NGSolve mesh.

             Args:
                 vertices: (n_verts, 3) vertex coordinates [m]
                 elements: (n_elem, 4) vertex indices (0-based)
             )doc")

        .def("add_hexahedra_from_mesh", &PyMMMBuilder::add_hexahedra_from_mesh,
             py::arg("vertices"), py::arg("elements"),
             R"doc(
             Add multiple hexahedra from mesh.

             Args:
                 vertices: (n_verts, 3) vertex coordinates [m]
                 elements: (n_elem, 8) vertex indices (0-based)
             )doc")

        .def("add_pm_tetrahedron", &PyMMMBuilder::add_pm_tetrahedron,
             py::arg("vertices"), py::arg("Mr"),
             R"doc(
             Add permanent magnet tetrahedron (fixed magnetization).

             PM elements have 0 DOF in the solve - their magnetization
             is fixed at Mr and they contribute to the field directly.

             Args:
                 vertices: (4, 3) array of vertex coordinates [m]
                 Mr: (3,) remanent magnetization [A/m]
             )doc")

        .def("add_pm_hexahedron", &PyMMMBuilder::add_pm_hexahedron,
             py::arg("vertices"), py::arg("Mr"),
             R"doc(
             Add permanent magnet hexahedron (fixed magnetization).

             PM elements have 0 DOF in the solve - their magnetization
             is fixed at Mr and they contribute to the field directly.

             Args:
                 vertices: (8, 3) array of vertex coordinates [m]
                 Mr: (3,) remanent magnetization [A/m]
             )doc")

        .def("add_pm_tetrahedra_from_mesh", &PyMMMBuilder::add_pm_tetrahedra_from_mesh,
             py::arg("vertices"), py::arg("elements"), py::arg("Mr"),
             R"doc(
             Add multiple PM tetrahedra from mesh.

             All elements share the same remanent magnetization Mr.

             Args:
                 vertices: (n_verts, 3) vertex coordinates [m]
                 elements: (n_elem, 4) vertex indices (0-based)
                 Mr: (3,) remanent magnetization [A/m]
             )doc")

        .def("add_pm_hexahedra_from_mesh", &PyMMMBuilder::add_pm_hexahedra_from_mesh,
             py::arg("vertices"), py::arg("elements"), py::arg("Mr"),
             R"doc(
             Add multiple PM hexahedra from mesh.

             All elements share the same remanent magnetization Mr.

             Args:
                 vertices: (n_verts, 3) vertex coordinates [m]
                 elements: (n_elem, 8) vertex indices (0-based)
                 Mr: (3,) remanent magnetization [A/m]
             )doc")

        .def("build", &PyMMMBuilder::build,
             R"doc(
             Build interaction matrix.

             Returns:
                 Tuple (N, dof_offset):
                   N: (total_dof, total_dof) interaction matrix
                   dof_offset: (n_elem + 1,) DOF offset per element
             )doc")

        .def("clear", &PyMMMBuilder::clear, "Clear all elements")

        .def_property_readonly("num_elements", &PyMMMBuilder::num_elements,
            "Number of elements")
        .def_property_readonly("total_dof", &PyMMMBuilder::total_dof,
            "Total degrees of freedom");

    // MMMSolver class
    py::class_<PyMMMSolver>(m, "MMMSolver",
        R"doc(
        MMM Linear System Solver

        Solves (I/chi - N) * M = H_ext
        )doc")
        .def(py::init<>(), "Create a new MMM solver")

        .def("set_matrix", &PyMMMSolver::set_matrix,
             py::arg("N"), py::arg("dof_offset"),
             R"doc(
             Set interaction matrix.

             Args:
                 N: (total_dof, total_dof) interaction matrix
                 dof_offset: (n_elem + 1,) DOF offset array
             )doc")

        .def("solve_lu", &PyMMMSolver::solve_lu,
             py::arg("inv_chi"), py::arg("H_ext"),
             py::arg("chi_per_element") = false,
             R"doc(
             Solve using LU decomposition (LAPACK dgesv).

             Args:
                 inv_chi: 1/chi values (total_dof,) or (n_elem,)
                 H_ext: External field (total_dof,)
                 chi_per_element: If True, inv_chi is (n_elem,) and broadcast

             Returns:
                 M: Magnetization (total_dof,)
             )doc")

        .def("solve_bicgstab", &PyMMMSolver::solve_bicgstab,
             py::arg("inv_chi"), py::arg("H_ext"),
             py::arg("tol") = 1e-8, py::arg("max_iter") = 1000,
             py::arg("chi_per_element") = false,
             R"doc(
             Solve using BiCGSTAB iteration.

             Args:
                 inv_chi: 1/chi values
                 H_ext: External field
                 tol: Convergence tolerance
                 max_iter: Maximum iterations
                 chi_per_element: If True, inv_chi is per-element

             Returns:
                 Tuple (M, n_iter):
                   M: Magnetization (total_dof,)
                   n_iter: Number of iterations
             )doc")

        .def("matvec", &PyMMMSolver::matvec,
             py::arg("inv_chi"), py::arg("x"),
             py::arg("chi_per_element") = false,
             R"doc(
             Matrix-vector product: y = A * x where A = -diag(inv_chi) - N

             For use in external iterative solvers (e.g., Schur complement coupling).

             Args:
                 inv_chi: 1/chi values (total_dof,) or (n_elem,)
                 x: Input vector (total_dof,)
                 chi_per_element: If True, inv_chi is per-element

             Returns:
                 y: Output vector (total_dof,)
             )doc")

        .def_property_readonly("total_dof", &PyMMMSolver::total_dof,
            "Total degrees of freedom")
        .def_property_readonly("num_elements", &PyMMMSolver::num_elements,
            "Number of elements")
        .def_property_readonly("is_ready", &PyMMMSolver::is_ready,
            "True if matrix is set");

    // MMMFieldComputer class
    py::class_<PyMMMFieldComputer>(m, "MMMFieldComputer",
        R"doc(
        MMM Field Computer

        Computes B/H fields at observation points from solved magnetization.

        Example:
            builder = MMMBuilder()
            # ... add elements ...
            N, dof_offset = builder.build()

            solver = MMMSolver()
            solver.set_matrix(N, dof_offset)
            M = solver.solve_lu(inv_chi, H_ext)

            # Compute field
            field_computer = MMMFieldComputer()
            field_computer.set_elements_from_builder(builder)
            B = field_computer.compute_b_field(M, obs_points)
        )doc")
        .def(py::init<>(), "Create a new MMM field computer")

        .def("set_elements_from_builder", &PyMMMFieldComputer::set_elements_from_builder,
             py::arg("builder"),
             R"doc(
             Set element geometry from MMMBuilder.

             Args:
                 builder: MMMBuilder instance with elements added
             )doc")

        .def("compute_b_field", &PyMMMFieldComputer::compute_b_field,
             py::arg("M"), py::arg("obs_points"),
             R"doc(
             Compute B field at observation points.

             Args:
                 M: Magnetization vector (total_dof,) [A/m]
                 obs_points: Observation points (n_points, 3) [m]

             Returns:
                 B: Magnetic flux density (n_points, 3) [T]
             )doc")

        .def("compute_h_field", &PyMMMFieldComputer::compute_h_field,
             py::arg("M"), py::arg("obs_points"),
             R"doc(
             Compute H field at observation points.

             Args:
                 M: Magnetization vector (total_dof,) [A/m]
                 obs_points: Observation points (n_points, 3) [m]

             Returns:
                 H: Magnetic field intensity (n_points, 3) [A/m]
             )doc")

        .def("compute_a_field", &PyMMMFieldComputer::compute_a_field,
             py::arg("M"), py::arg("obs_points"),
             R"doc(
             Compute vector potential A at observation points.

             Uses magnetic dipole formula. For hexahedra (6 DOF), reconstructs
             equivalent magnetization M from surface charges via least-squares
             on face normals, then computes A = (mu_0/4pi) * m x r / |r|^3.

             Args:
                 M: Magnetization vector (total_dof,) [A/m or surface charge]
                 obs_points: Observation points (n_points, 3) [m]

             Returns:
                 A: Vector potential (n_points, 3) [T*m]
             )doc")

        .def("get_element_centers", &PyMMMFieldComputer::get_element_centers,
             R"doc(
             Get element center coordinates.

             Returns:
                 centers: Element centers (n_elements, 3) [m]
             )doc")

        .def("get_element_volumes", &PyMMMFieldComputer::get_element_volumes,
             R"doc(
             Get element volumes.

             Returns:
                 volumes: Element volumes (n_elements,) [m^3]
             )doc")

        .def_property_readonly("num_elements", &PyMMMFieldComputer::num_elements,
            "Number of elements")
        .def_property_readonly("total_dof", &PyMMMFieldComputer::total_dof,
            "Total degrees of freedom");

    // MMMNonlinearSolver class
    py::class_<PyMMMNonlinearSolver>(m, "MMMNonlinearSolver",
        R"doc(
        MMM Nonlinear Solver

        Newton-Raphson iteration for nonlinear B-H curve materials.

        Example:
            builder = MMMBuilder()
            # ... add elements ...
            N, dof_offset = builder.build()

            # B-H curve (soft iron)
            bh_curve = np.array([
                [0.0, 0.0],
                [100.0, 0.5],
                [1000.0, 1.5],
                [10000.0, 1.8],
            ])

            solver = MMMNonlinearSolver()
            solver.set_matrix(N, dof_offset)
            solver.set_bh_curve(bh_curve)
            solver.set_parameters(tol=1e-4, max_iter=100)

            result = solver.solve(H_ext)
            M = result['M']
            converged = result['converged']
        )doc")
        .def(py::init<>(), "Create a new nonlinear MMM solver")

        .def("set_matrix", &PyMMMNonlinearSolver::set_matrix,
             py::arg("N"), py::arg("dof_offset"),
             R"doc(
             Set interaction matrix.

             Args:
                 N: (total_dof, total_dof) interaction matrix
                 dof_offset: (n_elem + 1,) DOF offset array
             )doc")

        .def("set_bh_curve", &PyMMMNonlinearSolver::set_bh_curve,
             py::arg("bh_curve"),
             R"doc(
             Set B-H curve for nonlinear material.

             Args:
                 bh_curve: (n_points, 2) array [[H0,B0], [H1,B1], ...]
                           H in A/m, B in Tesla, sorted by H ascending
             )doc")

        .def("set_parameters", &PyMMMNonlinearSolver::set_parameters,
             py::arg("tol") = 1e-4, py::arg("max_iter") = 100,
             py::arg("relax") = 1.0, py::arg("linear_solver") = 0,
             R"doc(
             Set solver parameters.

             Args:
                 tol: Convergence tolerance (relative B change)
                 max_iter: Maximum Newton-Raphson iterations
                 relax: Relaxation factor (0 < relax <= 1, 1=no relaxation)
                 linear_solver: 0=LU, 1=BiCGSTAB
             )doc")

        .def("set_initial_chi", &PyMMMNonlinearSolver::set_initial_chi,
             py::arg("chi_init"),
             R"doc(
             Set initial chi guess (optional).

             Args:
                 chi_init: Initial chi per element (n_elem,)
             )doc")

        .def("solve", &PyMMMNonlinearSolver::solve,
             py::arg("H_ext"),
             R"doc(
             Solve nonlinear system.

             Args:
                 H_ext: External field (3,) for uniform or (total_dof,)

             Returns:
                 dict with:
                   M: Magnetization (total_dof,)
                   chi: Final chi per element (n_elem,)
                   iterations: Number of iterations used
                   final_error: Final convergence error
                   converged: True if converged
             )doc")

        .def("get_local_h", &PyMMMNonlinearSolver::get_local_h,
             R"doc(
             Get local H field at element centers after solve.

             Returns:
                 H_local: (n_elem, 3) array [A/m]
             )doc")

        .def("get_local_b", &PyMMMNonlinearSolver::get_local_b,
             R"doc(
             Get local B field at element centers after solve.

             Returns:
                 B_local: (n_elem, 3) array [T]
             )doc")

        .def_property_readonly("total_dof", &PyMMMNonlinearSolver::total_dof,
            "Total degrees of freedom")
        .def_property_readonly("num_elements", &PyMMMNonlinearSolver::num_elements,
            "Number of elements");

    // compute_chi_from_bh function
    m.def("compute_chi_from_bh",
          [](py::array_t<double, py::array::c_style | py::array::forcecast> M_flat,
             py::array_t<double, py::array::c_style | py::array::forcecast> H_flat,
             py::array_t<double, py::array::c_style | py::array::forcecast> bh_curve) {

              auto M_info = M_flat.request();
              auto H_info = H_flat.request();
              auto bh_info = bh_curve.request();

              if (M_info.size != H_info.size) {
                  throw std::runtime_error("M and H must have same size");
              }
              if (bh_info.ndim != 2 || bh_info.shape[1] != 2) {
                  throw std::runtime_error("bh_curve must have shape (n, 2)");
              }

              int n_elem = static_cast<int>(M_info.size / 3);
              int n_pts = static_cast<int>(bh_info.shape[0]);

              std::vector<double> M_vec(static_cast<const double*>(M_info.ptr),
                                        static_cast<const double*>(M_info.ptr) + M_info.size);
              std::vector<double> H_vec(static_cast<const double*>(H_info.ptr),
                                        static_cast<const double*>(H_info.ptr) + H_info.size);

              std::vector<std::array<double, 2>> bh_vec(n_pts);
              auto bh_ptr = static_cast<const double*>(bh_info.ptr);
              for (int i = 0; i < n_pts; i++) {
                  bh_vec[i][0] = bh_ptr[i * 2];
                  bh_vec[i][1] = bh_ptr[i * 2 + 1];
              }

              auto chi = ComputeChiFromBH(M_vec, H_vec, bh_vec, n_elem);
              return py::array_t<double>(chi.size(), chi.data());
          },
          py::arg("M"), py::arg("H"), py::arg("bh_curve"),
          R"doc(
          Compute chi from B-H curve.

          Args:
              M: Magnetization (n_elem * 3,)
              H: Local H field (n_elem * 3,)
              bh_curve: B-H data [[H0,B0], [H1,B1], ...] sorted by H

          Returns:
              chi: Susceptibility (n_elem,)
          )doc");

    // check_convergence function
    m.def("check_convergence",
          [](py::array_t<double, py::array::c_style | py::array::forcecast> B_new,
             py::array_t<double, py::array::c_style | py::array::forcecast> B_old,
             double B_sat) {

              auto B_new_info = B_new.request();
              auto B_old_info = B_old.request();

              if (B_new_info.size != B_old_info.size) {
                  throw std::runtime_error("B_new and B_old must have same size");
              }

              int n_elem = static_cast<int>(B_new_info.size / 3);

              std::vector<double> B_new_vec(static_cast<const double*>(B_new_info.ptr),
                                            static_cast<const double*>(B_new_info.ptr) + B_new_info.size);
              std::vector<double> B_old_vec(static_cast<const double*>(B_old_info.ptr),
                                            static_cast<const double*>(B_old_info.ptr) + B_old_info.size);

              return CheckConvergence(B_new_vec, B_old_vec, B_sat, n_elem);
          },
          py::arg("B_new"), py::arg("B_old"), py::arg("B_sat"),
          R"doc(
          Check convergence based on B-field change.

          Args:
              B_new: New B field (n_elem * 3,)
              B_old: Old B field (n_elem * 3,)
              B_sat: Saturation B for relative error

          Returns:
              max_rel_change: Maximum relative B change
          )doc");

    // Physical constants
    m.attr("MU_0") = constants::MU_0;
    m.attr("INV_FOUR_PI") = constants::INV_FOUR_PI;
}
