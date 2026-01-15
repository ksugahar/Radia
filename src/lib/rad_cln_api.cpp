/**
 * @file rad_cln_api.cpp
 * @brief pybind11 wrapper for CLN module
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "../core/rad_cln.h"

// Intel MKL for BLAS operations in create_compressed_cln
#include <mkl.h>

namespace py = pybind11;

/**
 * @brief Python wrapper for LanczosResult
 */
class PyLanczosResult {
public:
    radia::cln::LanczosResult result_;

    PyLanczosResult(radia::cln::LanczosResult&& r) : result_(std::move(r)) {}

    py::array_t<double> get_L() const {
        return py::array_t<double>(result_.L.size(), result_.L.data());
    }

    py::array_t<double> get_R() const {
        return py::array_t<double>(result_.R.size(), result_.R.data());
    }

    py::array_t<double> get_Q() const {
        // Return as 2D array (n_input x n_output)
        return py::array_t<double>(
            {result_.n_input, result_.n_output},
            {result_.n_output * sizeof(double), sizeof(double)},
            result_.Q.data()
        );
    }

    py::array_t<double> get_R_diag() const {
        return py::array_t<double>(
            {result_.n_output, result_.n_output},
            {result_.n_output * sizeof(double), sizeof(double)},
            result_.R_diag.data()
        );
    }

    py::array_t<double> get_L_tridiag() const {
        return py::array_t<double>(
            {result_.n_output, result_.n_output},
            {result_.n_output * sizeof(double), sizeof(double)},
            result_.L_tridiag.data()
        );
    }

    int get_n_input() const { return result_.n_input; }
    int get_n_output() const { return result_.n_output; }
    bool get_converged() const { return result_.converged; }

    // Backward compatibility aliases
    py::array_t<double> get_alpha() const { return get_L(); }
    py::array_t<double> get_beta() const { return get_R(); }
    py::array_t<double> get_LL() const { return get_R_diag(); }
    py::array_t<double> get_RR() const { return get_L_tridiag(); }
};


PYBIND11_MODULE(cln_core, m) {
    m.doc() = R"doc(
CLN (Cauer Ladder Network) core module

Provides Lanczos algorithm for PEEC model order reduction.

CLN I form reduction:
    lanczos(K=R, N=L) -> R_diag (diagonal), L_tridiag (tridiagonal)

Example:
    import numpy as np
    from cln_core import lanczos, compute_cln_impedance_sweep

    # PEEC matrices (R: diagonal, L: dense with mutual inductance)
    R = np.diag([0.001, 0.001, 0.001])  # 1 mOhm each
    L = np.array([[10e-9, 5e-9, 2e-9],
                  [5e-9, 10e-9, 5e-9],
                  [2e-9, 5e-9, 10e-9]])  # nH

    # Reduce to 2 stages
    result = lanczos(R, L, n_iter=2)

    # Compute impedance sweep
    freqs = np.logspace(0, 8, 100)
    Z = compute_cln_impedance_sweep(result.R_diag, result.L_tridiag, freqs)
)doc";

    // LanczosResult class
    py::class_<PyLanczosResult>(m, "LanczosResult", R"doc(
Result of Lanczos transformation for CLN I model reduction.

Attributes:
    L : ndarray
        Series inductance (alpha), size n_output
    R : ndarray
        Parallel resistance (beta), size n_output
    Q : ndarray
        Transformation matrix (n_input x n_output)
        K-orthogonal: Q^T * K * Q = R_diag
    R_diag : ndarray
        Diagonalized resistance matrix (n_output x n_output)
    L_tridiag : ndarray
        Tridiagonalized inductance matrix (n_output x n_output)
    n_input : int
        Original dimension
    n_output : int
        Reduced dimension
    converged : bool
        True if converged without breakdown

Backward compatibility aliases:
    alpha = L, beta = R, LL = R_diag, RR = L_tridiag
)doc")
        .def_property_readonly("L", &PyLanczosResult::get_L)
        .def_property_readonly("R", &PyLanczosResult::get_R)
        .def_property_readonly("Q", &PyLanczosResult::get_Q)
        .def_property_readonly("R_diag", &PyLanczosResult::get_R_diag)
        .def_property_readonly("L_tridiag", &PyLanczosResult::get_L_tridiag)
        .def_property_readonly("n_input", &PyLanczosResult::get_n_input)
        .def_property_readonly("n_output", &PyLanczosResult::get_n_output)
        .def_property_readonly("converged", &PyLanczosResult::get_converged)
        // Backward compatibility
        .def_property_readonly("alpha", &PyLanczosResult::get_alpha)
        .def_property_readonly("beta", &PyLanczosResult::get_beta)
        .def_property_readonly("LL", &PyLanczosResult::get_LL)
        .def_property_readonly("RR", &PyLanczosResult::get_RR);

    // lanczos function
    m.def("lanczos", [](py::array_t<double, py::array::c_style | py::array::forcecast> K,
                        py::array_t<double, py::array::c_style | py::array::forcecast> N,
                        int n_iter, double tol) -> PyLanczosResult {
        auto K_info = K.request();
        auto N_info = N.request();

        if (K_info.ndim != 2 || N_info.ndim != 2) {
            throw std::runtime_error("K and N must be 2D arrays");
        }

        int n = static_cast<int>(K_info.shape[0]);
        if (K_info.shape[1] != n || N_info.shape[0] != n || N_info.shape[1] != n) {
            throw std::runtime_error("K and N must be square matrices of same size");
        }

        auto result = radia::cln::lanczos(
            static_cast<const double*>(K_info.ptr),
            static_cast<const double*>(N_info.ptr),
            n, n_iter, tol
        );

        return PyLanczosResult(std::move(result));
    },
    py::arg("K"), py::arg("N"), py::arg("n_iter") = -1, py::arg("tol") = 1e-30,
    R"doc(
Lanczos algorithm for CLN I model reduction.

Transforms two Hermitian matrices (K, N) simultaneously:
    K -> R_diag (diagonal)
    N -> L_tridiag (tridiagonal)

For DC PEEC: lanczos(K=R, N=L) produces CLN I form.

Parameters:
    K : ndarray (n, n)
        Input matrix (resistance R)
    N : ndarray (n, n)
        Input matrix (inductance L)
    n_iter : int, optional
        Number of iterations (reduced dimension). -1 for full transform.
    tol : float, optional
        Convergence tolerance

Returns:
    LanczosResult with reduced CLN I model
)doc");

    // build_tridiagonal function
    m.def("build_tridiagonal", [](py::array_t<double, py::array::c_style | py::array::forcecast> diag) {
        auto diag_info = diag.request();
        if (diag_info.ndim != 1) {
            throw std::runtime_error("diag must be 1D array");
        }

        int n = static_cast<int>(diag_info.size);
        auto T = radia::cln::build_tridiagonal(
            static_cast<const double*>(diag_info.ptr), n
        );

        return py::array_t<double>(
            {n, n},
            {n * sizeof(double), sizeof(double)},
            T.data()
        );
    },
    py::arg("diag"),
    R"doc(
Build tridiagonal matrix from diagonal values.

Constructs CLN I inductance matrix structure:
    L[i,i] = diag[i] + diag[i+1]
    L[i,i+1] = L[i+1,i] = -diag[i+1]

Parameters:
    diag : ndarray (n,)
        Diagonal values

Returns:
    ndarray (n, n) : Tridiagonal matrix
)doc");

    // compute_cln_impedance function
    m.def("compute_cln_impedance", [](py::array_t<double, py::array::c_style | py::array::forcecast> R_diag,
                                       py::array_t<double, py::array::c_style | py::array::forcecast> L_tridiag,
                                       double freq) {
        auto R_info = R_diag.request();
        auto L_info = L_tridiag.request();

        if (R_info.ndim != 2 || L_info.ndim != 2) {
            throw std::runtime_error("R_diag and L_tridiag must be 2D arrays");
        }

        int n = static_cast<int>(R_info.shape[0]);

        return radia::cln::compute_cln_impedance(
            static_cast<const double*>(R_info.ptr),
            static_cast<const double*>(L_info.ptr),
            n, freq
        );
    },
    py::arg("R_diag"), py::arg("L_tridiag"), py::arg("freq"),
    R"doc(
Compute CLN I impedance at a single frequency.

Z(s) = 1 / I[0] where (R_diag + s*L_tridiag) * I = V, V[0]=1

Parameters:
    R_diag : ndarray (n, n)
        Diagonal resistance matrix
    L_tridiag : ndarray (n, n)
        Tridiagonal inductance matrix
    freq : float
        Frequency in Hz

Returns:
    complex : Impedance at given frequency
)doc");

    // compute_cln_impedance_sweep function
    m.def("compute_cln_impedance_sweep", [](py::array_t<double, py::array::c_style | py::array::forcecast> R_diag,
                                             py::array_t<double, py::array::c_style | py::array::forcecast> L_tridiag,
                                             py::array_t<double, py::array::c_style | py::array::forcecast> freqs) {
        auto R_info = R_diag.request();
        auto L_info = L_tridiag.request();
        auto freqs_info = freqs.request();

        if (R_info.ndim != 2 || L_info.ndim != 2) {
            throw std::runtime_error("R_diag and L_tridiag must be 2D arrays");
        }
        if (freqs_info.ndim != 1) {
            throw std::runtime_error("freqs must be 1D array");
        }

        int n = static_cast<int>(R_info.shape[0]);
        int n_freqs = static_cast<int>(freqs_info.size);

        // Allocate output
        auto Z_out = py::array_t<std::complex<double>>(n_freqs);
        auto Z_info = Z_out.request();

        radia::cln::compute_cln_impedance_sweep(
            static_cast<const double*>(R_info.ptr),
            static_cast<const double*>(L_info.ptr),
            n,
            static_cast<const double*>(freqs_info.ptr),
            n_freqs,
            static_cast<std::complex<double>*>(Z_info.ptr)
        );

        return Z_out;
    },
    py::arg("R_diag"), py::arg("L_tridiag"), py::arg("freqs"),
    R"doc(
Compute CLN I impedance over frequency sweep.

Parameters:
    R_diag : ndarray (n, n)
        Diagonal resistance matrix
    L_tridiag : ndarray (n, n)
        Tridiagonal inductance matrix
    freqs : ndarray (n_freqs,)
        Frequency array in Hz

Returns:
    ndarray (n_freqs,) complex : Impedance at each frequency
)doc");

    // transform_coupling function
    m.def("transform_coupling", [](py::array_t<double, py::array::c_style | py::array::forcecast> Q,
                                    py::array_t<double, py::array::c_style | py::array::forcecast> M_LS) {
        auto Q_info = Q.request();
        auto M_LS_info = M_LS.request();

        if (Q_info.ndim != 2 || M_LS_info.ndim != 2) {
            throw std::runtime_error("Q and M_LS must be 2D arrays");
        }

        int n_loop = static_cast<int>(Q_info.shape[0]);
        int n_reduced = static_cast<int>(Q_info.shape[1]);
        int n_star = static_cast<int>(M_LS_info.shape[1]);

        if (M_LS_info.shape[0] != n_loop) {
            throw std::runtime_error("Q and M_LS must have same number of rows (n_loop)");
        }

        auto result = radia::cln::transform_coupling(
            static_cast<const double*>(Q_info.ptr),
            static_cast<const double*>(M_LS_info.ptr),
            n_loop, n_reduced, n_star
        );

        return py::array_t<double>(
            {n_reduced, n_star},
            {n_star * sizeof(double), sizeof(double)},
            result.data()
        );
    },
    py::arg("Q"), py::arg("M_LS"),
    R"doc(
Transform Loop-Star coupling matrix with Lanczos Q matrix.

Computes: M_LS_reduced = Q^T * M_LS

Parameters:
    Q : ndarray (n_loop, n_reduced)
        Transformation matrix from Lanczos
    M_LS : ndarray (n_loop, n_star)
        Loop-Star coupling matrix

Returns:
    ndarray (n_reduced, n_star) : Transformed coupling matrix
)doc");

    // transform_port_vector function
    m.def("transform_port_vector", [](py::array_t<double, py::array::c_style | py::array::forcecast> Q,
                                       py::array_t<double, py::array::c_style | py::array::forcecast> v) {
        auto Q_info = Q.request();
        auto v_info = v.request();

        if (Q_info.ndim != 2) {
            throw std::runtime_error("Q must be 2D array");
        }
        if (v_info.ndim != 1) {
            throw std::runtime_error("v must be 1D array");
        }

        int n_loop = static_cast<int>(Q_info.shape[0]);
        int n_reduced = static_cast<int>(Q_info.shape[1]);

        if (v_info.size != n_loop) {
            throw std::runtime_error("v size must match Q rows (n_loop)");
        }

        auto result = radia::cln::transform_port_vector(
            static_cast<const double*>(Q_info.ptr),
            static_cast<const double*>(v_info.ptr),
            n_loop, n_reduced
        );

        return py::array_t<double>(result.size(), result.data());
    },
    py::arg("Q"), py::arg("v"),
    R"doc(
Transform port vector with Lanczos Q matrix.

Computes: v_reduced = Q^T * v

Parameters:
    Q : ndarray (n_loop, n_reduced)
        Transformation matrix from Lanczos
    v : ndarray (n_loop,)
        Port excitation vector

Returns:
    ndarray (n_reduced,) : Transformed port vector
)doc");

    // LoopStarCLN class
    py::class_<radia::cln::LoopStarCLN>(m, "LoopStarCLN", R"doc(
Full Loop-Star CLN system for PEEC analysis.

Represents the complete reduced PEEC system:
    [R_diag + s*L_tridiag   s*M_LS_reduced] [I_L_reduced]   [V_L_reduced]
    [s*M_LS_reduced^T       P/s           ] [I_S        ] = [V_S        ]

Attributes:
    R_diag : ndarray (n_reduced, n_reduced)
        Diagonal resistance matrix (CLN I)
    L_tridiag : ndarray (n_reduced, n_reduced)
        Tridiagonal inductance matrix (CLN I)
    M_LS_reduced : ndarray (n_reduced, n_star)
        Transformed Loop-Star coupling
    n_reduced : int
        Reduced loop dimension
    n_star : int
        Star dimension (unchanged)
    n_loop : int
        Original loop dimension
)doc")
        .def(py::init([](py::array_t<double, py::array::c_style | py::array::forcecast> R_diag,
                         py::array_t<double, py::array::c_style | py::array::forcecast> L_tridiag,
                         py::array_t<double, py::array::c_style | py::array::forcecast> M_LS_reduced,
                         py::array_t<double, py::array::c_style | py::array::forcecast> P,
                         py::array_t<double, py::array::c_style | py::array::forcecast> Q) {
            auto R_info = R_diag.request();
            auto L_info = L_tridiag.request();
            auto M_info = M_LS_reduced.request();
            auto P_info = P.request();
            auto Q_info = Q.request();

            radia::cln::LoopStarCLN cln;
            cln.n_reduced = static_cast<int>(R_info.shape[0]);
            cln.n_star = static_cast<int>(M_info.shape[1]);
            cln.n_loop = static_cast<int>(Q_info.shape[0]);

            // Copy data
            cln.R_diag.assign(static_cast<const double*>(R_info.ptr),
                              static_cast<const double*>(R_info.ptr) + R_info.size);
            cln.L_tridiag.assign(static_cast<const double*>(L_info.ptr),
                                 static_cast<const double*>(L_info.ptr) + L_info.size);
            cln.M_LS_reduced.assign(static_cast<const double*>(M_info.ptr),
                                    static_cast<const double*>(M_info.ptr) + M_info.size);
            cln.Q.assign(static_cast<const double*>(Q_info.ptr),
                         static_cast<const double*>(Q_info.ptr) + Q_info.size);

            // Store P pointer (NOTE: Python must keep P alive!)
            cln.P = static_cast<const double*>(P_info.ptr);

            return cln;
        }),
        py::arg("R_diag"), py::arg("L_tridiag"), py::arg("M_LS_reduced"), py::arg("P"), py::arg("Q"),
        "Construct LoopStarCLN from matrices. NOTE: P array must be kept alive in Python!")
        .def_readonly("n_reduced", &radia::cln::LoopStarCLN::n_reduced)
        .def_readonly("n_star", &radia::cln::LoopStarCLN::n_star)
        .def_readonly("n_loop", &radia::cln::LoopStarCLN::n_loop);

    // compute_loop_star_impedance function
    m.def("compute_loop_star_impedance", [](const radia::cln::LoopStarCLN& cln,
                                             double freq,
                                             py::array_t<double, py::array::c_style | py::array::forcecast> port_loop) {
        auto port_info = port_loop.request();
        if (port_info.ndim != 1 || port_info.size != cln.n_loop) {
            throw std::runtime_error("port_loop must be 1D array with n_loop elements");
        }

        return radia::cln::compute_loop_star_impedance(
            cln, freq, static_cast<const double*>(port_info.ptr)
        );
    },
    py::arg("cln"), py::arg("freq"), py::arg("port_loop"),
    R"doc(
Compute full Loop-Star CLN impedance at a frequency.

Solves the complete PEEC system including capacitive effects.

Parameters:
    cln : LoopStarCLN
        The reduced Loop-Star CLN system
    freq : float
        Frequency in Hz
    port_loop : ndarray (n_loop,)
        Port vector in original loop space

Returns:
    complex : Port impedance
)doc");

    // create_loop_star_cln helper function
    m.def("create_loop_star_cln", [](PyLanczosResult& lanczos_result,
                                      py::array_t<double, py::array::c_style | py::array::forcecast> M_LS,
                                      py::array_t<double, py::array::c_style | py::array::forcecast> P) {
        auto M_LS_info = M_LS.request();
        auto P_info = P.request();

        const auto& lr = lanczos_result.result_;

        // Transform M_LS
        int n_star = static_cast<int>(M_LS_info.shape[1]);
        auto M_LS_reduced = radia::cln::transform_coupling(
            lr.Q.data(),
            static_cast<const double*>(M_LS_info.ptr),
            lr.n_input, lr.n_output, n_star
        );

        radia::cln::LoopStarCLN cln;
        cln.n_reduced = lr.n_output;
        cln.n_star = n_star;
        cln.n_loop = lr.n_input;
        cln.R_diag = lr.R_diag;
        cln.L_tridiag = lr.L_tridiag;
        cln.M_LS_reduced = std::move(M_LS_reduced);
        cln.Q = lr.Q;
        cln.P = static_cast<const double*>(P_info.ptr);

        return cln;
    },
    py::arg("lanczos_result"), py::arg("M_LS"), py::arg("P"),
    R"doc(
Create LoopStarCLN from Lanczos result and PEEC matrices.

Convenience function that transforms M_LS using the Lanczos Q matrix.

Parameters:
    lanczos_result : LanczosResult
        Result from lanczos() function
    M_LS : ndarray (n_loop, n_star)
        Loop-Star coupling matrix
    P : ndarray (n_star, n_star)
        Potential coefficient matrix (kept by reference!)

Returns:
    LoopStarCLN : The complete reduced system
)doc");

    // ACAResult class
    py::class_<radia::cln::ACAResult>(m, "ACAResult", R"doc(
Result of ACA+ low-rank approximation.

P = U * V^T where U: (n x k), V: (n x k)

Attributes:
    U : ndarray (n, k)
        Left factor matrix
    V : ndarray (n, k)
        Right factor matrix
    n : int
        Original matrix dimension
    k : int
        Rank after compression
    compression_ratio : float
        k / n
    converged : bool
        True if ACA+ converged
)doc")
        .def_readonly("n", &radia::cln::ACAResult::n)
        .def_readonly("k", &radia::cln::ACAResult::k)
        .def_readonly("compression_ratio", &radia::cln::ACAResult::compression_ratio)
        .def_readonly("converged", &radia::cln::ACAResult::converged)
        .def_property_readonly("U", [](const radia::cln::ACAResult& r) {
            return py::array_t<double>(
                {r.n, r.k},
                {r.k * sizeof(double), sizeof(double)},
                r.U.data()
            );
        })
        .def_property_readonly("V", [](const radia::cln::ACAResult& r) {
            return py::array_t<double>(
                {r.n, r.k},
                {r.k * sizeof(double), sizeof(double)},
                r.V.data()
            );
        })
        .def("reconstruct", [](const radia::cln::ACAResult& r) {
            auto P_approx = r.reconstruct();
            return py::array_t<double>(
                {r.n, r.n},
                {r.n * sizeof(double), sizeof(double)},
                P_approx.data()
            );
        }, "Reconstruct full matrix P = U * V^T");

    // aca_compress function
    m.def("aca_compress", [](py::array_t<double, py::array::c_style | py::array::forcecast> P,
                              double eps, int kmax) {
        auto P_info = P.request();
        if (P_info.ndim != 2 || P_info.shape[0] != P_info.shape[1]) {
            throw std::runtime_error("P must be a square 2D array");
        }

        int n = static_cast<int>(P_info.shape[0]);
        return radia::cln::aca_compress(
            static_cast<const double*>(P_info.ptr),
            n, eps, kmax
        );
    },
    py::arg("P"), py::arg("eps") = 1e-4, py::arg("kmax") = -1,
    R"doc(
ACA+ low-rank approximation.

Computes P = U * V^T using ACA+ algorithm (via HACApK) or SVD fallback.

Parameters:
    P : ndarray (n, n)
        Input matrix
    eps : float
        ACA+ tolerance (default: 1e-4)
    kmax : int
        Maximum rank, -1 for automatic (default: -1)

Returns:
    ACAResult with U, V factors
)doc");

    // LoopStarCLNCompressed class
    py::class_<radia::cln::LoopStarCLNCompressed>(m, "LoopStarCLNCompressed", R"doc(
Full Loop-Star CLN with ACA+ compressed Star space.

Double reduction:
  - Loop: n_loop -> n_reduced (Lanczos)
  - Star: n_star -> k_star (ACA+)

Attributes:
    n_reduced : int
        Reduced loop dimension
    n_star : int
        Original star dimension
    k_star : int
        Compressed star dimension
    n_loop : int
        Original loop dimension
)doc")
        .def_readonly("n_reduced", &radia::cln::LoopStarCLNCompressed::n_reduced)
        .def_readonly("n_star", &radia::cln::LoopStarCLNCompressed::n_star)
        .def_readonly("k_star", &radia::cln::LoopStarCLNCompressed::k_star)
        .def_readonly("n_loop", &radia::cln::LoopStarCLNCompressed::n_loop);

    // create_compressed_cln helper function
    m.def("create_compressed_cln", [](PyLanczosResult& lanczos_result,
                                       py::array_t<double, py::array::c_style | py::array::forcecast> M_LS,
                                       py::array_t<double, py::array::c_style | py::array::forcecast> P,
                                       double eps, int kmax) {
        auto M_LS_info = M_LS.request();
        auto P_info = P.request();

        const auto& lr = lanczos_result.result_;

        int n_star = static_cast<int>(M_LS_info.shape[1]);

        // Step 1: ACA+ compress P -> U, V
        auto aca_result = radia::cln::aca_compress(
            static_cast<const double*>(P_info.ptr),
            n_star, eps, kmax
        );

        int k_star = aca_result.k;

        // Step 2: Transform M_LS to reduced loop space
        // M_LS_reduced = Q^T * M_LS (n_reduced x n_star)
        auto M_LS_reduced = radia::cln::transform_coupling(
            lr.Q.data(),
            static_cast<const double*>(M_LS_info.ptr),
            lr.n_input, lr.n_output, n_star
        );

        // Step 3: Transform M_LS_reduced to mode space
        // M_LS_mode = M_LS_reduced * V (n_reduced x k_star)
        std::vector<double> M_LS_mode(lr.n_output * k_star);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    lr.n_output, k_star, n_star,
                    1.0, M_LS_reduced.data(), n_star,
                    aca_result.V.data(), k_star,
                    0.0, M_LS_mode.data(), k_star);

        // Step 4: Compute P_mode = U^T * P * V (k_star x k_star)
        // First: temp = P * V (n_star x k_star)
        std::vector<double> PV(n_star * k_star);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    n_star, k_star, n_star,
                    1.0, static_cast<const double*>(P_info.ptr), n_star,
                    aca_result.V.data(), k_star,
                    0.0, PV.data(), k_star);

        // Then: P_mode = U^T * PV (k_star x k_star)
        std::vector<double> P_mode(k_star * k_star);
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                    k_star, k_star, n_star,
                    1.0, aca_result.U.data(), k_star,
                    PV.data(), k_star,
                    0.0, P_mode.data(), k_star);

        // Build result
        radia::cln::LoopStarCLNCompressed cln;
        cln.n_reduced = lr.n_output;
        cln.n_star = n_star;
        cln.k_star = k_star;
        cln.n_loop = lr.n_input;

        cln.R_diag = lr.R_diag;
        cln.L_tridiag = lr.L_tridiag;
        cln.M_LS_mode = std::move(M_LS_mode);
        cln.P_mode = std::move(P_mode);
        cln.U = aca_result.U;
        cln.V = aca_result.V;
        cln.Q = lr.Q;

        return cln;
    },
    py::arg("lanczos_result"), py::arg("M_LS"), py::arg("P"),
    py::arg("eps") = 1e-4, py::arg("kmax") = -1,
    R"doc(
Create LoopStarCLNCompressed with double reduction.

Applies:
  - Lanczos reduction to Loop space (from lanczos_result)
  - ACA+ compression to Star space

Parameters:
    lanczos_result : LanczosResult
        Result from lanczos() function
    M_LS : ndarray (n_loop, n_star)
        Loop-Star coupling matrix
    P : ndarray (n_star, n_star)
        Potential coefficient matrix
    eps : float
        ACA+ tolerance (default: 1e-4)
    kmax : int
        Maximum rank for ACA+, -1 for automatic

Returns:
    LoopStarCLNCompressed : Doubly-reduced system
)doc");

    // compute_compressed_impedance function
    m.def("compute_compressed_impedance", [](const radia::cln::LoopStarCLNCompressed& cln,
                                              double freq,
                                              py::array_t<double, py::array::c_style | py::array::forcecast> port_loop) {
        auto port_info = port_loop.request();
        if (port_info.ndim != 1 || port_info.size != cln.n_loop) {
            throw std::runtime_error("port_loop must be 1D array with n_loop elements");
        }

        return radia::cln::compute_compressed_impedance(
            cln, freq, static_cast<const double*>(port_info.ptr)
        );
    },
    py::arg("cln"), py::arg("freq"), py::arg("port_loop"),
    R"doc(
Compute compressed Loop-Star CLN impedance.

Uses doubly-reduced system (Lanczos + ACA+).

Parameters:
    cln : LoopStarCLNCompressed
        The compressed system
    freq : float
        Frequency in Hz
    port_loop : ndarray (n_loop,)
        Port vector in original loop space

Returns:
    complex : Port impedance
)doc");
}
