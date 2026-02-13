/*
 * rad_peec_matrices_api.cpp
 *
 * pybind11 Python bindings for PEEC matrix construction
 *
 * Style follows radia_ngsolve.cpp pattern for consistency.
 *
 * Part of Radia project
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include "rad_peec_matrices.h"

namespace py = pybind11;

using namespace radia;

//=========================================================================
// Python Wrapper Classes
//=========================================================================

/**
 * Python-friendly PEEC Matrix Builder
 *
 * Wraps PEECMatrixBuilder with numpy array returns.
 */
class PyPEECBuilder {
public:
    PyPEECBuilder() : builder_(std::make_unique<PEECMatrixBuilder>()),
                      matrices_built_(false) {}

    /**
     * Add a conductor segment
     *
     * Args:
     *   center: Segment center [x, y, z] in meters
     *   direction: Unit direction vector [dx, dy, dz]
     *   length: Segment length [m]
     *   width: Cross-section width [m]
     *   height: Cross-section height [m]
     *   sigma: Conductivity [S/m] (default: 5.8e7 for copper)
     *   cross_section_type: 0 = rectangular (default), 1 = circular
     */
    void add_segment(py::array_t<double> center,
                     py::array_t<double> direction,
                     double length,
                     double width,
                     double height,
                     double sigma = 5.8e7,
                     int cross_section_type = 0) {
        auto c = center.unchecked<1>();
        auto d = direction.unchecked<1>();

        CrossSectionType type = (cross_section_type == 1)
                                    ? CrossSectionType::CIRCULAR
                                    : CrossSectionType::RECTANGULAR;

        PEECSegment seg;
        seg.center = TVector3d(c(0), c(1), c(2));
        seg.direction = TVector3d(d(0), d(1), d(2));
        seg.length = length;
        seg.width = width;
        seg.height = height;
        seg.sigma = sigma;
        seg.cross_section_type = type;

        builder_->AddSegment(seg);
        matrices_built_ = false;
    }

    /**
     * Add a charge node
     */
    void add_node(py::array_t<double> position, double area) {
        auto p = position.unchecked<1>();

        PEECNode node;
        node.position = TVector3d(p(0), p(1), p(2));
        node.area = area;

        builder_->AddNode(node);
        matrices_built_ = false;
    }

    /**
     * Add a 2D surface panel (triangle or quadrilateral)
     *
     * Args:
     *   vertices: List of 3 or 4 vertices [[x,y,z], ...] in meters
     */
    void add_panel(py::list vertices) {
        std::vector<TVector3d> verts;
        for (auto item : vertices) {
            auto v = item.cast<py::array_t<double>>().unchecked<1>();
            verts.push_back(TVector3d(v(0), v(1), v(2)));
        }

        if (verts.size() != 3 && verts.size() != 4) {
            throw std::invalid_argument("Panel must have 3 or 4 vertices");
        }

        PEECPanel panel(verts);
        builder_->AddPanel(panel);
        matrices_built_ = false;
    }

    /**
     * Create wire segments along a straight path
     *
     * Returns number of segments created.
     */
    int create_wire(py::array_t<double> start,
                    py::array_t<double> end,
                    double width,
                    double height,
                    int n_segments,
                    double sigma = 5.8e7) {
        auto s = start.unchecked<1>();
        auto e = end.unchecked<1>();

        TVector3d p1(s(0), s(1), s(2));
        TVector3d p2(e(0), e(1), e(2));

        auto segments = CreateWireSegments(p1, p2, width, height, n_segments, sigma);
        builder_->AddSegments(segments);
        matrices_built_ = false;

        return static_cast<int>(segments.size());
    }

    /**
     * Create circular loop segments
     *
     * Returns number of segments created.
     */
    int create_loop(py::array_t<double> center,
                    double radius,
                    py::array_t<double> normal,
                    double width,
                    double height,
                    int n_segments,
                    double sigma = 5.8e7) {
        auto c = center.unchecked<1>();
        auto n = normal.unchecked<1>();

        TVector3d ctr(c(0), c(1), c(2));
        TVector3d nrm(n(0), n(1), n(2));

        auto segments = CreateLoopSegments(ctr, radius, nrm, width, height, n_segments, sigma);
        builder_->AddSegments(segments);
        matrices_built_ = false;

        return static_cast<int>(segments.size());
    }

    /**
     * Build PEEC matrices (frequency-independent: L, R_dc, P, M_LS)
     *
     * Returns tuple: (L, R, P, M_LS)
     *   L: Inductance matrix [H] (n_loop x n_loop)
     *   R: Resistance vector [Ohm] (n_loop,) diagonal only
     *   P: Potential coefficient [1/F] (n_star x n_star) or None
     *   M_LS: Loop-Star coupling (n_loop x n_star) or None
     */
    py::tuple build(bool include_star = true) {
        // Auto-generate nodes if not provided
        if (builder_->NumNodes() == 0) {
            builder_->AutoGenerateNodes();
        }

        matrices_ = builder_->Build(include_star);
        matrices_built_ = true;

        int n_loop = matrices_.n_loop;
        int n_star = matrices_.n_star;

        // L matrix
        py::array_t<double> L({n_loop, n_loop});
        auto L_buf = L.mutable_unchecked<2>();
        for (int i = 0; i < n_loop; ++i) {
            for (int j = 0; j < n_loop; ++j) {
                L_buf(i, j) = matrices_.L[i * n_loop + j];
            }
        }

        // R vector (diagonal)
        py::array_t<double> R(n_loop);
        auto R_buf = R.mutable_unchecked<1>();
        for (int i = 0; i < n_loop; ++i) {
            R_buf(i) = matrices_.R[i];
        }

        if (include_star && n_star > 0) {
            // P matrix
            py::array_t<double> P({n_star, n_star});
            auto P_buf = P.mutable_unchecked<2>();
            for (int i = 0; i < n_star; ++i) {
                for (int j = 0; j < n_star; ++j) {
                    P_buf(i, j) = matrices_.P[i * n_star + j];
                }
            }

            // M_LS matrix
            py::array_t<double> M_LS({n_loop, n_star});
            auto MLS_buf = M_LS.mutable_unchecked<2>();
            for (int i = 0; i < n_loop; ++i) {
                for (int j = 0; j < n_star; ++j) {
                    MLS_buf(i, j) = matrices_.M_LS[i * n_star + j];
                }
            }

            return py::make_tuple(L, R, P, M_LS);
        } else {
            return py::make_tuple(L, R, py::none(), py::none());
        }
    }

    /**
     * Compute port impedance at given frequency
     *
     * Args:
     *   freq_hz: Frequency in Hz
     *   port_vector: Excitation vector (n_loop,)
     *
     * Returns:
     *   Complex impedance [Ohm]
     */
    std::complex<double> compute_impedance(double freq_hz,
                                            py::array_t<double> port_vector) {
        if (!matrices_built_) {
            build(true);
        }

        PEECSolver solver;
        solver.SetMatrices(matrices_);
        solver.SetFrequency(freq_hz);

        auto pv = port_vector.unchecked<1>();
        std::vector<double> pv_vec(pv.size());
        for (py::ssize_t i = 0; i < pv.size(); ++i) {
            pv_vec[i] = pv(i);
        }

        return solver.ComputePortImpedance(pv_vec);
    }

    // Properties
    int n_loop() const { return builder_->NumSegments(); }
    int n_star() const { return builder_->NumNodes(); }
    int num_segments() const { return builder_->NumSegments(); }
    int num_nodes() const { return builder_->NumNodes(); }
    int num_panels() const { return builder_->NumPanels(); }

    /**
     * Clear all geometry
     */
    void clear() {
        builder_->Clear();
        matrices_built_ = false;
    }

private:
    std::unique_ptr<PEECMatrixBuilder> builder_;
    PEECMatrices matrices_;
    bool matrices_built_;
};

//=========================================================================
// pybind11 Module Definition
//=========================================================================

PYBIND11_MODULE(peec_matrices, m) {
    m.doc() = R"doc(
PEEC Matrix Construction Module

Loop-Star decomposition for quasi-static electromagnetic analysis.
Uses Darwin approximation: G(r) = 1/(4*pi*r)

Matrix Types:
  L: Inductance matrix (Loop-Loop) - Neumann formula [H]
  P: Potential coefficient matrix (Star-Star) [1/F], C = P^{-1}
  R: Resistance matrix (diagonal) [Ohm]
  M_LS: Loop-Star coupling matrix

PEEC System Equation:
  [Z_LL   Z_LS] [I_L]   [V_L]
  [Z_SL   Z_SS] [I_S] = [V_S]

  Z_LL = R_dc + jw*L + Z_s  (Z_s computed in Python: Bessel/Dowell/ESIM)
  Z_SS = P / jw
  Z_LS = jw*M_LS
  Z_SL = Z_LS^T (reciprocity)

Note: build() returns frequency-independent matrices (L, R_dc, P, M_LS).
      Frequency-dependent surface impedance Z_s is computed in Python
      using scipy.special.jv (circular) or Dowell formula (rectangular).

Valid Frequency Range: DC to ~100 MHz (Darwin approximation)

Example:
    from peec_matrices import PEECBuilder
    import numpy as np
    from scipy.special import jv

    builder = PEECBuilder()
    builder.create_wire([0,0,0], [0.1,0,0], 1e-3, 1e-3, 10, 5.8e7)
    L, R_dc, P, M_LS = builder.build()

    # AC impedance with Bessel SIBC (Python-side)
    freq = 1e6
    omega = 2 * np.pi * freq
    Z = np.sum(R_dc) + 1j * omega * np.sum(L)  # + Z_s from Bessel
)doc";

    py::class_<PyPEECBuilder>(m, "PEECBuilder",
        R"doc(
        PEEC Matrix Builder

        Constructs L, P, R, M_LS matrices from conductor segments.
        )doc")
        .def(py::init<>(), "Create a new PEEC matrix builder")

        .def("add_segment", &PyPEECBuilder::add_segment,
             py::arg("center"),
             py::arg("direction"),
             py::arg("length"),
             py::arg("width"),
             py::arg("height"),
             py::arg("sigma") = 5.8e7,
             py::arg("cross_section_type") = 0,
             R"doc(
             Add a conductor segment.

             Args:
                 center: Segment center [x, y, z] in meters
                 direction: Unit direction vector [dx, dy, dz]
                 length: Segment length [m]
                 width: Cross-section width [m]
                 height: Cross-section height [m]
                 sigma: Conductivity [S/m] (default: copper 5.8e7)
                 cross_section_type: 0 = rectangular (default), 1 = circular
             )doc")

        .def("add_node", &PyPEECBuilder::add_node,
             py::arg("position"),
             py::arg("area"),
             R"doc(
             Add a charge node for capacitive effects.

             Args:
                 position: Node position [x, y, z] in meters
                 area: Associated area [m^2] for self-potential
             )doc")

        .def("add_panel", &PyPEECBuilder::add_panel,
             py::arg("vertices"),
             R"doc(
             Add a 2D surface panel (triangle or quadrilateral).

             Args:
                 vertices: List of 3 or 4 vertices [[x,y,z], ...] in meters
                           3 vertices: Triangle
                           4 vertices: Quadrilateral

             Example:
                 # Triangle panel
                 builder.add_panel([[0,0,0], [0.01,0,0], [0,0.01,0]])

                 # Quadrilateral panel
                 builder.add_panel([[0,0,0], [0.01,0,0], [0.01,0.01,0], [0,0.01,0]])
             )doc")

        .def("create_wire", &PyPEECBuilder::create_wire,
             py::arg("start"),
             py::arg("end"),
             py::arg("width"),
             py::arg("height"),
             py::arg("n_segments"),
             py::arg("sigma") = 5.8e7,
             R"doc(
             Create wire segments along a straight path.

             Args:
                 start: Start point [x, y, z] in meters
                 end: End point [x, y, z] in meters
                 width: Cross-section width [m]
                 height: Cross-section height [m]
                 n_segments: Number of segments
                 sigma: Conductivity [S/m]

             Returns:
                 Number of segments created
             )doc")

        .def("create_loop", &PyPEECBuilder::create_loop,
             py::arg("center"),
             py::arg("radius"),
             py::arg("normal"),
             py::arg("width"),
             py::arg("height"),
             py::arg("n_segments"),
             py::arg("sigma") = 5.8e7,
             R"doc(
             Create circular loop segments.

             Args:
                 center: Loop center [x, y, z] in meters
                 radius: Loop radius [m]
                 normal: Normal vector [nx, ny, nz]
                 width: Wire cross-section width [m]
                 height: Wire cross-section height [m]
                 n_segments: Number of segments
                 sigma: Conductivity [S/m]

             Returns:
                 Number of segments created
             )doc")

        .def("build", &PyPEECBuilder::build,
             py::arg("include_star") = true,
             R"doc(
             Build PEEC matrices.

             Args:
                 include_star: If True, compute P and M_LS matrices

             Returns:
                 Tuple (L, R_dc, P, M_LS):
                   L: Inductance matrix [H] (n_loop x n_loop)
                   R_dc: DC resistance vector [Ohm] (n_loop,)
                   P: Potential coefficient [1/F] (n_star x n_star) or None
                   M_LS: Loop-Star coupling (n_loop x n_star) or None

             Note:
                 All matrices are frequency-independent.
                 For AC analysis, compute Z_s in Python (scipy Bessel/Dowell/ESIM)
                 and add to R_dc.
             )doc")

        .def("compute_impedance", &PyPEECBuilder::compute_impedance,
             py::arg("freq_hz"),
             py::arg("port_vector"),
             R"doc(
             Compute port impedance at given frequency.

             Args:
                 freq_hz: Frequency in Hz
                 port_vector: Excitation vector (n_loop,)

             Returns:
                 Complex port impedance [Ohm]
             )doc")

        .def("clear", &PyPEECBuilder::clear, "Clear all geometry")

        .def_property_readonly("n_loop", &PyPEECBuilder::n_loop,
            "Number of Loop elements (segments)")
        .def_property_readonly("n_star", &PyPEECBuilder::n_star,
            "Number of Star elements (nodes)")
        .def_property_readonly("num_segments", &PyPEECBuilder::num_segments,
            "Number of segments")
        .def_property_readonly("num_nodes", &PyPEECBuilder::num_nodes,
            "Number of nodes")
        .def_property_readonly("num_panels", &PyPEECBuilder::num_panels,
            "Number of panels");

    // Convenience function for wire PEEC
    m.def("create_wire_peec",
          [](py::array_t<double> start,
             py::array_t<double> end,
             double width,
             double height,
             int n_segments,
             double sigma,
             bool include_star) {
              PyPEECBuilder builder;
              builder.create_wire(start, end, width, height, n_segments, sigma);
              return builder.build(include_star);
          },
          py::arg("start"),
          py::arg("end"),
          py::arg("width"),
          py::arg("height"),
          py::arg("n_segments"),
          py::arg("sigma") = 5.8e7,
          py::arg("include_star") = true,
          R"doc(
          Create PEEC matrices for a straight wire.

          Convenience function that creates a builder, adds wire segments,
          and returns the matrices.

          Args:
              start: Start point [x, y, z] in meters
              end: End point [x, y, z] in meters
              width: Cross-section width [m]
              height: Cross-section height [m]
              n_segments: Number of segments
              sigma: Conductivity [S/m]
              include_star: Include P and M_LS matrices

          Returns:
              Tuple (L, R, P, M_LS)
          )doc");

    // Convenience function for loop PEEC
    m.def("create_loop_peec",
          [](py::array_t<double> center,
             double radius,
             py::array_t<double> normal,
             double width,
             double height,
             int n_segments,
             double sigma,
             bool include_star) {
              PyPEECBuilder builder;
              builder.create_loop(center, radius, normal, width, height, n_segments, sigma);
              return builder.build(include_star);
          },
          py::arg("center"),
          py::arg("radius"),
          py::arg("normal"),
          py::arg("width"),
          py::arg("height"),
          py::arg("n_segments"),
          py::arg("sigma") = 5.8e7,
          py::arg("include_star") = true,
          R"doc(
          Create PEEC matrices for a circular loop.

          Convenience function that creates a builder, adds loop segments,
          and returns the matrices.

          Args:
              center: Loop center [x, y, z] in meters
              radius: Loop radius [m]
              normal: Normal vector [nx, ny, nz]
              width: Wire cross-section width [m]
              height: Wire cross-section height [m]
              n_segments: Number of segments
              sigma: Conductivity [S/m]
              include_star: Include P and M_LS matrices

          Returns:
              Tuple (L, R, P, M_LS)
          )doc");
}
