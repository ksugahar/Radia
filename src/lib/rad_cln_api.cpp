/**
 * @file rad_cln_api.cpp
 * @brief pybind11 wrapper for CLN module
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "../core/rad_cln.h"

// Historical native dependency; generic core helpers remain available.
#include <mkl.h>

namespace py = pybind11;

PYBIND11_MODULE(cln_core, m) {
    m.doc() = "Generic reduced-matrix, impedance, coupling and ACA helpers.";




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
Transform Loop-Star coupling matrix with a supplied projection Q matrix.

Computes: M_LS_reduced = Q^T * M_LS

Parameters:
    Q : ndarray (n_loop, n_reduced)
        Supplied transformation matrix
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
Transform port vector with a supplied projection Q matrix.

Computes: v_reduced = Q^T * v

Parameters:
    Q : ndarray (n_loop, n_reduced)
        Supplied transformation matrix
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

    // Historical compressed-result type: no Python constructor remains after
    // retiring its CLN-specific producer. Core C++ matrix helpers are retained.
    // LoopStarCLNCompressed class
    py::class_<radia::cln::LoopStarCLNCompressed>(m, "LoopStarCLNCompressed", R"doc(
Full Loop-Star CLN with ACA+ compressed Star space.

Double reduction:
  - Loop: n_loop -> n_reduced (projection)
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

Uses doubly-reduced system (projection + ACA+).

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
