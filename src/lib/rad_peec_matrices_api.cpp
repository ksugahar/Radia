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
     * Generate face panels from segment box surfaces
     *
     * Automatically creates quad panels from each segment's rectangular
     * cross-section extruded along its length.
     *
     * Args:
     *   mode: Panel generation mode
     *     0 = 6 faces (all box surfaces)
     *     1 = 2 faces (top + bottom) — inter-layer capacitance
     *     2 = 4 faces (top + bottom + left + right) — inter-turn + inter-layer
     *   eps_r: Relative permittivity for P matrix scaling (e.g., 4.4 for FR4)
     *
     * Returns number of panels generated.
     */
    int generate_face_panels(int mode = 2, double eps_r = 1.0) {
        builder_->GenerateFacePanels(mode, eps_r);
        matrices_built_ = false;
        return builder_->NumPanels();
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

    // ========== Topology-aware methods ==========

    /**
     * Add a node at a given position
     *
     * Returns node ID (auto-incrementing).
     */
    int add_node_at(double x, double y, double z, double area = 0.0) {
        TVector3d pos(x, y, z);
        int id = builder_->AddNodeAt(pos, area);
        matrices_built_ = false;
        return id;
    }

    /**
     * Add a segment connecting two nodes
     */
    void add_connected_segment(int node_from, int node_to,
                               double width, double height,
                               double sigma = 5.8e7,
                               int cross_section_type = 0,
                               int nwinc = 1, int nhinc = 1) {
        CrossSectionType type = (cross_section_type == 1)
                                    ? CrossSectionType::CIRCULAR
                                    : CrossSectionType::RECTANGULAR;
        builder_->AddConnectedSegment(node_from, node_to, width, height, sigma, type, nwinc, nhinc);
        matrices_built_ = false;
    }

    /**
     * Add a port between two nodes
     *
     * Returns port ID.
     */
    int add_port(int node_positive, int node_negative) {
        int id = builder_->AddPort(node_positive, node_negative);
        matrices_built_ = false;
        return id;
    }

    /**
     * Build PEEC matrices with topology information
     *
     * Returns dict with L, R, P (optional), incidence matrix (CSR), ports.
     */
    py::dict build_topology(bool include_star = false, double eps_eff = 1.0) {
        matrices_ = builder_->Build(include_star, eps_eff);
        matrices_built_ = true;

        int n_loop = matrices_.n_loop;
        int n_star = matrices_.n_star;

        py::dict result;

        // L matrix
        py::array_t<double> L({n_loop, n_loop});
        auto L_buf = L.mutable_unchecked<2>();
        for (int i = 0; i < n_loop; ++i) {
            for (int j = 0; j < n_loop; ++j) {
                L_buf(i, j) = matrices_.L[i * n_loop + j];
            }
        }
        result["L"] = L;

        // R vector (diagonal)
        py::array_t<double> R(n_loop);
        auto R_buf = R.mutable_unchecked<1>();
        for (int i = 0; i < n_loop; ++i) {
            R_buf(i) = matrices_.R[i];
        }
        result["R"] = R;

        // P matrix (optional)
        if (include_star && n_star > 0) {
            py::array_t<double> P({n_star, n_star});
            auto P_buf = P.mutable_unchecked<2>();
            for (int i = 0; i < n_star; ++i) {
                for (int j = 0; j < n_star; ++j) {
                    P_buf(i, j) = matrices_.P[i * n_star + j];
                }
            }
            result["P"] = P;

            // M_LS matrix (Loop-Star coupling)
            py::array_t<double> M_LS({n_loop, n_star});
            auto MLS_buf = M_LS.mutable_unchecked<2>();
            for (int i = 0; i < n_loop; ++i) {
                for (int j = 0; j < n_star; ++j) {
                    MLS_buf(i, j) = matrices_.M_LS[i * n_star + j];
                }
            }
            result["M_LS"] = M_LS;
        } else {
            result["P"] = py::none();
            result["M_LS"] = py::none();
        }

        // Incidence matrix (CSR format)
        int n_junction = matrices_.n_junction;
        result["n_junction"] = n_junction;

        {
            int indptr_size = static_cast<int>(matrices_.incidence_indptr.size());
            py::array_t<int> indptr(indptr_size);
            auto ip_buf = indptr.mutable_unchecked<1>();
            for (int i = 0; i < indptr_size; ++i) {
                ip_buf(i) = matrices_.incidence_indptr[i];
            }
            result["incidence_indptr"] = indptr;
        }

        {
            int indices_size = static_cast<int>(matrices_.incidence_indices.size());
            py::array_t<int> indices(indices_size);
            auto idx_buf = indices.mutable_unchecked<1>();
            for (int i = 0; i < indices_size; ++i) {
                idx_buf(i) = matrices_.incidence_indices[i];
            }
            result["incidence_indices"] = indices;
        }

        {
            int data_size = static_cast<int>(matrices_.incidence_data.size());
            py::array_t<double> data(data_size);
            auto d_buf = data.mutable_unchecked<1>();
            for (int i = 0; i < data_size; ++i) {
                d_buf(i) = matrices_.incidence_data[i];
            }
            result["incidence_data"] = data;
        }

        // Port definitions
        py::list port_list;
        for (const auto& port : matrices_.ports) {
            port_list.append(py::make_tuple(port.node_positive, port.node_negative, port.port_id));
        }
        result["ports"] = port_list;

        result["n_loop"] = n_loop;
        result["n_star"] = n_star;

        // Segment connectivity and geometry
        {
            const auto& segs = builder_->GetSegments();
            int n_seg = static_cast<int>(segs.size());

            // Connectivity (node_from, node_to)
            py::array_t<int> seg_nodes({n_seg, 2});
            auto sn_buf = seg_nodes.mutable_unchecked<2>();

            // Geometry arrays
            py::array_t<double> seg_centers({n_seg, 3});
            py::array_t<double> seg_directions({n_seg, 3});
            py::array_t<double> seg_lengths(n_seg);
            py::array_t<double> seg_widths(n_seg);
            py::array_t<double> seg_heights(n_seg);
            py::array_t<double> seg_sigmas(n_seg);

            auto sc_buf = seg_centers.mutable_unchecked<2>();
            auto sd_buf = seg_directions.mutable_unchecked<2>();
            auto sl_buf = seg_lengths.mutable_unchecked<1>();
            auto sw_buf = seg_widths.mutable_unchecked<1>();
            auto sh_buf = seg_heights.mutable_unchecked<1>();
            auto ss_buf = seg_sigmas.mutable_unchecked<1>();

            for (int i = 0; i < n_seg; ++i) {
                sn_buf(i, 0) = segs[i].node_from;
                sn_buf(i, 1) = segs[i].node_to;

                sc_buf(i, 0) = segs[i].center.x;
                sc_buf(i, 1) = segs[i].center.y;
                sc_buf(i, 2) = segs[i].center.z;

                sd_buf(i, 0) = segs[i].direction.x;
                sd_buf(i, 1) = segs[i].direction.y;
                sd_buf(i, 2) = segs[i].direction.z;

                sl_buf(i) = segs[i].length;
                sw_buf(i) = segs[i].width;
                sh_buf(i) = segs[i].height;
                ss_buf(i) = segs[i].sigma;
            }

            result["segment_nodes"] = seg_nodes;
            result["segment_centers"] = seg_centers;
            result["segment_directions"] = seg_directions;
            result["segment_lengths"] = seg_lengths;
            result["segment_widths"] = seg_widths;
            result["segment_heights"] = seg_heights;
            result["segment_sigmas"] = seg_sigmas;
        }

        // Total number of nodes
        result["n_nodes"] = static_cast<int>(builder_->GetNodes().size());

        // Gathered capacitance C_eff = G × (P/eps_eff)⁻¹ × G^T
        if (!matrices_.C_gathered.empty() && matrices_.n_nodes_gathered > 0) {
            int ng = matrices_.n_nodes_gathered;
            py::array_t<double> C_np({ng, ng});
            auto C_buf = C_np.mutable_unchecked<2>();
            for (int i = 0; i < ng; ++i) {
                for (int j = 0; j < ng; ++j) {
                    C_buf(i, j) = matrices_.C_gathered[i * ng + j];
                }
            }
            result["C_gathered"] = C_np;
            // Override n_nodes from C_gathered if it provides a better count
            if (ng > static_cast<int>(builder_->GetNodes().size())) {
                result["n_nodes"] = ng;
            }
        } else {
            result["C_gathered"] = py::none();
        }

        return result;
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

    // ========== MNA Multi-Port Solver (C++ LAPACK) ==========

    /**
     * Compute multi-port Z-parameter matrix at a single frequency.
     *
     * Uses LAPACK zgetrf_/zgetrs_.
     *
     * Args:
     *   freq: Frequency in Hz
     *   Zs: Optional surface impedance array (n_loop,) complex.
     *        If None, only DC resistance is used.
     *
     * Returns:
     *   Z_matrix: (n_ports, n_ports) complex numpy array
     */
    py::array_t<std::complex<double>> solve_z_matrix(double freq,
                                                      py::object Zs_obj = py::none()) {
        ensure_solver_ready();

        int n_ports = solver_->NumPorts();
        if (n_ports == 0) {
            return py::array_t<std::complex<double>>({0, 0});
        }

        // Extract Zs from numpy or None
        const std::complex<double>* Zs_ptr = nullptr;
        int n_Zs = 0;
        py::array_t<std::complex<double>> Zs_arr;
        if (!Zs_obj.is_none()) {
            Zs_arr = Zs_obj.cast<py::array_t<std::complex<double>>>();
            auto Zs_buf = Zs_arr.unchecked<1>();
            Zs_ptr = Zs_buf.data(0);
            n_Zs = static_cast<int>(Zs_buf.size());
        }

        // Allocate output
        py::array_t<std::complex<double>> Z_out({n_ports, n_ports});
        auto Z_buf = Z_out.mutable_unchecked<2>();

        // Solve
        std::vector<std::complex<double>> Z_result(n_ports * n_ports);
        solver_->ComputeZMatrix(freq, Zs_ptr, n_Zs, Z_result.data(), n_ports);

        // Copy to numpy
        for (int i = 0; i < n_ports; ++i) {
            for (int j = 0; j < n_ports; ++j) {
                Z_buf(i, j) = Z_result[i * n_ports + j];
            }
        }

        return Z_out;
    }

    /**
     * Compute coupling coefficient from Z-parameters.
     *
     * Args:
     *   freq: Frequency in Hz
     *   Zs: Optional surface impedance array (n_loop,) complex
     *
     * Returns:
     *   dict with: Z_matrix, L_matrix, k_matrix, k, L1, L2, M
     */
    py::dict solve_coupling(double freq, py::object Zs_obj = py::none()) {
        ensure_solver_ready();

        int n_ports = solver_->NumPorts();
        py::dict result;

        if (n_ports < 2) {
            result["Z_matrix"] = py::array_t<std::complex<double>>({n_ports, n_ports});
            result["L_matrix"] = py::array_t<double>({n_ports, n_ports});
            result["k_matrix"] = py::array_t<double>({n_ports, n_ports});
            result["k"] = 0.0;
            result["L1"] = 0.0;
            result["L2"] = 0.0;
            result["M"] = 0.0;
            return result;
        }

        // Extract Zs
        const std::complex<double>* Zs_ptr = nullptr;
        int n_Zs = 0;
        py::array_t<std::complex<double>> Zs_arr;
        if (!Zs_obj.is_none()) {
            Zs_arr = Zs_obj.cast<py::array_t<std::complex<double>>>();
            auto Zs_buf = Zs_arr.unchecked<1>();
            Zs_ptr = Zs_buf.data(0);
            n_Zs = static_cast<int>(Zs_buf.size());
        }

        // Compute Z matrix
        std::vector<std::complex<double>> Z_result(n_ports * n_ports);
        solver_->ComputeZMatrix(freq, Zs_ptr, n_Zs, Z_result.data(), n_ports);

        // Compute L and k
        std::vector<double> L_result(n_ports * n_ports, 0.0);
        std::vector<double> k_result(n_ports * n_ports, 0.0);
        double omega = 2.0 * 3.14159265358979323846 * freq;

        if (omega > 1e-10) {
            for (int i = 0; i < n_ports; ++i) {
                for (int j = 0; j < n_ports; ++j) {
                    L_result[i * n_ports + j] = Z_result[i * n_ports + j].imag() / omega;
                }
            }
            for (int i = 0; i < n_ports; ++i) {
                for (int j = 0; j < n_ports; ++j) {
                    double Li = L_result[i * n_ports + i];
                    double Lj = L_result[j * n_ports + j];
                    if (Li > 0 && Lj > 0) {
                        k_result[i * n_ports + j] = L_result[i * n_ports + j] / std::sqrt(Li * Lj);
                    }
                }
            }
        }

        // Build numpy arrays
        py::array_t<std::complex<double>> Z_np({n_ports, n_ports});
        py::array_t<double> L_np({n_ports, n_ports});
        py::array_t<double> k_np({n_ports, n_ports});

        auto Zb = Z_np.mutable_unchecked<2>();
        auto Lb = L_np.mutable_unchecked<2>();
        auto kb = k_np.mutable_unchecked<2>();

        for (int i = 0; i < n_ports; ++i) {
            for (int j = 0; j < n_ports; ++j) {
                Zb(i, j) = Z_result[i * n_ports + j];
                Lb(i, j) = L_result[i * n_ports + j];
                kb(i, j) = k_result[i * n_ports + j];
            }
        }

        result["Z_matrix"] = Z_np;
        result["L_matrix"] = L_np;
        result["k_matrix"] = k_np;
        result["k"] = (n_ports >= 2) ? k_result[0 * n_ports + 1] : 0.0;
        result["L1"] = L_result[0 * n_ports + 0];
        result["L2"] = (n_ports >= 2) ? L_result[1 * n_ports + 1] : 0.0;
        result["M"] = (n_ports >= 2) ? L_result[0 * n_ports + 1] : 0.0;

        return result;
    }

    /**
     * Batch frequency sweep for multi-port Z-parameters.
     *
     * Args:
     *   freqs: Frequency array [Hz]
     *   Zs_all: Optional (n_freq, n_loop) complex array of per-frequency
     *           surface impedance. If None, only DC resistance is used.
     *
     * Returns:
     *   dict with: freqs, Z_matrix, L_matrix, k, L1, L2, M arrays
     */
    py::dict solve_multiport_sweep(py::array_t<double> freqs,
                                    py::object Zs_all_obj = py::none()) {
        ensure_solver_ready();

        auto f_buf = freqs.unchecked<1>();
        int n_freq = static_cast<int>(f_buf.size());
        int n_ports = solver_->NumPorts();
        int n_loop = matrices_.n_loop;

        // Extract Zs_all (optional)
        const std::complex<double>* Zs_all_ptr = nullptr;
        py::array_t<std::complex<double>> Zs_all_arr;
        if (!Zs_all_obj.is_none()) {
            Zs_all_arr = Zs_all_obj.cast<py::array_t<std::complex<double>>>();
            Zs_all_ptr = Zs_all_arr.data();
        }

        // Allocate results
        std::vector<std::complex<double>> Z_all(n_freq * n_ports * n_ports);
        solver_->FrequencySweep(f_buf.data(0), n_freq, Zs_all_ptr,
                                Z_all.data(), n_ports);

        // Build output arrays
        py::array_t<std::complex<double>> Z_np({n_freq, n_ports, n_ports});
        py::array_t<double> L_np({n_freq, n_ports, n_ports});
        py::array_t<double> k_arr(n_freq);
        py::array_t<double> L1_arr(n_freq);
        py::array_t<double> L2_arr(n_freq);
        py::array_t<double> M_arr(n_freq);

        auto Zb = Z_np.mutable_unchecked<3>();
        auto Lb = L_np.mutable_unchecked<3>();
        auto kb = k_arr.mutable_unchecked<1>();
        auto L1b = L1_arr.mutable_unchecked<1>();
        auto L2b = L2_arr.mutable_unchecked<1>();
        auto Mb = M_arr.mutable_unchecked<1>();

        for (int f = 0; f < n_freq; ++f) {
            double omega = 2.0 * 3.14159265358979323846 * f_buf(f);
            for (int i = 0; i < n_ports; ++i) {
                for (int j = 0; j < n_ports; ++j) {
                    int idx = f * n_ports * n_ports + i * n_ports + j;
                    Zb(f, i, j) = Z_all[idx];
                    Lb(f, i, j) = (omega > 1e-10) ? Z_all[idx].imag() / omega : 0.0;
                }
            }

            double L1_val = Lb(f, 0, 0);
            double L2_val = (n_ports >= 2) ? Lb(f, 1, 1) : 0.0;
            double M_val = (n_ports >= 2) ? Lb(f, 0, 1) : 0.0;

            L1b(f) = L1_val;
            L2b(f) = L2_val;
            Mb(f) = M_val;
            kb(f) = (L1_val > 0 && L2_val > 0) ? M_val / std::sqrt(L1_val * L2_val) : 0.0;
        }

        py::dict result;
        result["freqs"] = freqs;
        result["Z_matrix"] = Z_np;
        result["L_matrix"] = L_np;
        result["k"] = k_arr;
        result["L1"] = L1_arr;
        result["L2"] = L2_arr;
        result["M"] = M_arr;

        return result;
    }

private:
    std::unique_ptr<PEECMatrixBuilder> builder_;
    PEECMatrices matrices_;
    bool matrices_built_;
    std::unique_ptr<PEECSolver> solver_;

    /**
     * Ensure solver is initialized from current matrices and topology.
     */
    void ensure_solver_ready() {
        if (!matrices_built_) {
            matrices_ = builder_->Build(true);
            matrices_built_ = true;
        }

        if (!solver_) {
            solver_ = std::make_unique<PEECSolver>();
        }

        solver_->SetMatrices(matrices_);

        // Set topology (segment_nodes, n_nodes, ports)
        const auto& segs = builder_->GetSegments();
        int n_seg = static_cast<int>(segs.size());
        std::vector<std::pair<int,int>> seg_nodes(n_seg);
        for (int i = 0; i < n_seg; ++i) {
            seg_nodes[i] = {segs[i].node_from, segs[i].node_to};
        }
        solver_->SetSegmentNodes(seg_nodes, builder_->NumNodes());
        solver_->SetPorts(matrices_.ports);
    }
};

/**
 * Standalone MNA Solver (accepts raw numpy arrays, no PEECBuilder dependency)
 *
 * This class enables C++ LAPACK MNA solve for ALL callers:
 * - PEECBuilder topology path
 * - ngbem bridge (ngbem_interface.py)
 * - Shielded and reduced PEEC workflows
 * - Shielded PEEC (peec_shielded.py)
 *
 * Uses LAPACK zgesv_/zgetrf_/zgetrs_ + MKL cblas_zgemm
 * (shared MKL/LAPACK infrastructure).
 */
class PyMNASolver {
public:
    /**
     * Initialize from raw numpy arrays (topology_dict fields).
     *
     * Args:
     *   L: Inductance matrix (n_loop x n_loop) float64
     *   R: DC resistance vector (n_loop,) float64
     *   segment_nodes: Connectivity (n_seg x 2) int32, [node_from, node_to]
     *   n_nodes: Total number of nodes
     *   ports: List of (node_positive, node_negative, port_id) tuples
     *   P: Optional potential coefficient matrix (n_star x n_star) float64
     *   M_LS: Optional Loop-Star coupling (n_loop x n_star) float64
     */
    PyMNASolver(py::array_t<double> L,
                py::array_t<double> R,
                py::array_t<int> segment_nodes,
                int n_nodes,
                py::list ports,
                py::object P_obj = py::none(),
                py::object M_LS_obj = py::none(),
                py::object C_gathered_obj = py::none()) {

        auto L_buf = L.unchecked<2>();
        auto R_buf = R.unchecked<1>();
        auto sn_buf = segment_nodes.unchecked<2>();

        int n_loop = static_cast<int>(L_buf.shape(0));
        int n_seg = static_cast<int>(sn_buf.shape(0));

        // Build PEECMatrices
        matrices_.n_loop = n_loop;
        matrices_.L.resize(n_loop * n_loop);
        for (int i = 0; i < n_loop; ++i) {
            for (int j = 0; j < n_loop; ++j) {
                matrices_.L[i * n_loop + j] = L_buf(i, j);
            }
        }

        matrices_.R.resize(n_loop);
        for (int i = 0; i < n_loop; ++i) {
            matrices_.R[i] = R_buf(i);
        }

        // Optional P matrix
        if (!P_obj.is_none()) {
            auto P_arr = P_obj.cast<py::array_t<double>>();
            auto P_buf = P_arr.unchecked<2>();
            int n_star = static_cast<int>(P_buf.shape(0));
            matrices_.n_star = n_star;
            matrices_.P.resize(n_star * n_star);
            for (int i = 0; i < n_star; ++i) {
                for (int j = 0; j < n_star; ++j) {
                    matrices_.P[i * n_star + j] = P_buf(i, j);
                }
            }

            // Optional M_LS matrix
            if (!M_LS_obj.is_none()) {
                auto MLS_arr = M_LS_obj.cast<py::array_t<double>>();
                auto MLS_buf = MLS_arr.unchecked<2>();
                matrices_.M_LS.resize(n_loop * n_star);
                for (int i = 0; i < n_loop; ++i) {
                    for (int j = 0; j < n_star; ++j) {
                        matrices_.M_LS[i * n_star + j] = MLS_buf(i, j);
                    }
                }
            }
        } else {
            matrices_.n_star = 0;
        }

        // Segment connectivity
        std::vector<std::pair<int,int>> seg_nodes_vec(n_seg);
        for (int i = 0; i < n_seg; ++i) {
            seg_nodes_vec[i] = {sn_buf(i, 0), sn_buf(i, 1)};
        }

        // Ports
        std::vector<PEECPort> port_vec;
        for (auto item : ports) {
            auto tup = item.cast<py::tuple>();
            int pos = tup[0].cast<int>();
            int neg = tup[1].cast<int>();
            int pid = tup[2].cast<int>();
            port_vec.push_back(PEECPort(pos, neg, pid));
        }
        matrices_.ports = port_vec;

        // Initialize solver
        solver_ = std::make_unique<PEECSolver>();
        solver_->SetMatrices(matrices_);
        solver_->SetSegmentNodes(seg_nodes_vec, n_nodes);
        solver_->SetPorts(port_vec);

        // Set gathered capacitance (enables proper MNA with resonance)
        if (!C_gathered_obj.is_none()) {
            auto C_arr = C_gathered_obj.cast<py::array_t<double>>();
            auto C_buf = C_arr.unchecked<2>();
            int ng = static_cast<int>(C_buf.shape(0));
            std::vector<double> C_vec(ng * ng);
            for (int i = 0; i < ng; ++i) {
                for (int j = 0; j < ng; ++j) {
                    C_vec[i * ng + j] = C_buf(i, j);
                }
            }
            solver_->SetGatheredCapacitance(C_vec, ng);
        }

        n_ports_ = static_cast<int>(port_vec.size());
    }

    /**
     * Compute multi-port Z-parameter matrix at a single frequency.
     */
    py::array_t<std::complex<double>> solve_z_matrix(double freq,
                                                      py::object Zs_obj = py::none()) {
        if (n_ports_ == 0) {
            return py::array_t<std::complex<double>>({0, 0});
        }

        const std::complex<double>* Zs_ptr = nullptr;
        int n_Zs = 0;
        py::array_t<std::complex<double>> Zs_arr;
        if (!Zs_obj.is_none()) {
            Zs_arr = Zs_obj.cast<py::array_t<std::complex<double>>>();
            auto Zs_buf = Zs_arr.unchecked<1>();
            Zs_ptr = Zs_buf.data(0);
            n_Zs = static_cast<int>(Zs_buf.size());
        }

        py::array_t<std::complex<double>> Z_out({n_ports_, n_ports_});
        auto Z_buf = Z_out.mutable_unchecked<2>();

        std::vector<std::complex<double>> Z_result(n_ports_ * n_ports_);
        solver_->ComputeZMatrix(freq, Zs_ptr, n_Zs, Z_result.data(), n_ports_);

        for (int i = 0; i < n_ports_; ++i) {
            for (int j = 0; j < n_ports_; ++j) {
                Z_buf(i, j) = Z_result[i * n_ports_ + j];
            }
        }

        return Z_out;
    }

    /**
     * Batch frequency sweep for multi-port Z-parameters.
     */
    py::dict frequency_sweep(py::array_t<double> freqs,
                              py::object Zs_all_obj = py::none()) {
        auto f_buf = freqs.unchecked<1>();
        int n_freq = static_cast<int>(f_buf.size());

        const std::complex<double>* Zs_all_ptr = nullptr;
        py::array_t<std::complex<double>> Zs_all_arr;
        if (!Zs_all_obj.is_none()) {
            Zs_all_arr = Zs_all_obj.cast<py::array_t<std::complex<double>>>();
            Zs_all_ptr = Zs_all_arr.data();
        }

        std::vector<std::complex<double>> Z_all(n_freq * n_ports_ * n_ports_);
        solver_->FrequencySweep(f_buf.data(0), n_freq, Zs_all_ptr,
                                Z_all.data(), n_ports_);

        // Build output
        py::array_t<std::complex<double>> Z_np({n_freq, n_ports_, n_ports_});
        auto Zb = Z_np.mutable_unchecked<3>();

        for (int f = 0; f < n_freq; ++f) {
            for (int i = 0; i < n_ports_; ++i) {
                for (int j = 0; j < n_ports_; ++j) {
                    Zb(f, i, j) = Z_all[f * n_ports_ * n_ports_ + i * n_ports_ + j];
                }
            }
        }

        py::dict result;
        result["freqs"] = freqs;
        result["Z_matrix"] = Z_np;
        return result;
    }

    /**
     * Solve for per-filament complex branch currents at a given frequency,
     * given user-specified port current injection.
     *
     * I_branch[f] = current flowing from node_from[f] to node_to[f] (complex).
     * Useful for the "N parallel filaments between two ports, I_total injected"
     * workflow: pass port_currents=[I_total] and read off the AC current
     * distribution including skin + proximity coupling.
     */
    py::array_t<std::complex<double>> solve_branch_currents(
        double freq,
        py::array_t<std::complex<double>> port_currents,
        py::object Zs_obj = py::none()) {

        auto pc_buf = port_currents.unchecked<1>();
        int n_pc = static_cast<int>(pc_buf.size());

        const std::complex<double>* Zs_ptr = nullptr;
        int n_Zs = 0;
        py::array_t<std::complex<double>> Zs_arr;
        if (!Zs_obj.is_none()) {
            Zs_arr = Zs_obj.cast<py::array_t<std::complex<double>>>();
            auto Zs_buf = Zs_arr.unchecked<1>();
            Zs_ptr = Zs_buf.data(0);
            n_Zs = static_cast<int>(Zs_buf.size());
        }

        int n_loop = matrices_.n_loop;
        py::array_t<std::complex<double>> I_out(n_loop);
        auto I_buf = I_out.mutable_unchecked<1>();

        std::vector<std::complex<double>> I_branch(n_loop);
        solver_->ComputeBranchCurrents(freq,
                                        pc_buf.data(0), n_pc,
                                        Zs_ptr, n_Zs,
                                        I_branch.data(), n_loop);

        for (int i = 0; i < n_loop; ++i) I_buf(i) = I_branch[i];
        return I_out;
    }

    int num_ports() const { return n_ports_; }

    void set_solver_method(int method) { solver_->SetSolverMethod(method); }
    int get_solver_method() const { return solver_->GetSolverMethod(); }
    void set_bicgstab_params(double tol, int max_iter) {
        solver_->SetBiCGSTABParams(tol, max_iter);
    }

    /**
     * Set gathered capacitance for proper MNA (enables resonance).
     * C_gathered = G × (P/eps_eff)⁻¹ × G^T
     */
    void set_gathered_capacitance(py::array_t<double> C_gathered) {
        auto C_buf = C_gathered.unchecked<2>();
        int ng = static_cast<int>(C_buf.shape(0));
        std::vector<double> C_vec(ng * ng);
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < ng; ++j) {
                C_vec[i * ng + j] = C_buf(i, j);
            }
        }
        solver_->SetGatheredCapacitance(C_vec, ng);
    }

    bool has_gathered_capacitance() const { return solver_->HasGatheredCapacitance(); }

private:
    PEECMatrices matrices_;
    std::unique_ptr<PEECSolver> solver_;
    int n_ports_;
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

        .def("generate_face_panels", &PyPEECBuilder::generate_face_panels,
             py::arg("mode") = 2,
             py::arg("eps_r") = 1.0,
             R"doc(
             Generate face panels from segment box surfaces.

             Automatically creates quad panels from each segment's rectangular
             cross-section extruded along its length. Must be called after
             adding all segments and before build/build_topology.

             Args:
                 mode: Panel generation mode
                     0 = 6 faces (all box surfaces)
                     1 = 2 faces (top + bottom) — inter-layer capacitance
                     2 = 4 faces (top + bottom + left + right) — inter-turn + inter-layer
                     3 = 2 faces (left + right) — inter-turn only (single-layer PCB)
                     4 = 1 face (top only) — simplified single-side
                 eps_r: Relative permittivity for P matrix scaling (e.g., 4.4 for FR4)

             Returns:
                 Number of panels generated.

             Example:
                 builder.generate_face_panels(mode=1, eps_r=4.4)  # FR4 PCB
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

        // Topology-aware methods
        .def("add_node_at", &PyPEECBuilder::add_node_at,
             py::arg("x"), py::arg("y"), py::arg("z"),
             py::arg("area") = 0.0,
             R"doc(
             Add a node at a given position.

             Args:
                 x, y, z: Node position in meters
                 area: Associated area [m^2] (for capacitive effects)

             Returns:
                 Node ID (auto-incrementing, starting from 0)
             )doc")

        .def("add_connected_segment", &PyPEECBuilder::add_connected_segment,
             py::arg("node_from"), py::arg("node_to"),
             py::arg("width"), py::arg("height"),
             py::arg("sigma") = 5.8e7,
             py::arg("cross_section_type") = 0,
             py::arg("nwinc") = 1, py::arg("nhinc") = 1,
             R"doc(
             Add a segment connecting two nodes.

             Computes center, direction, and length from node positions.
             Parallel segments share the same pair of nodes.

             If nwinc*nhinc > 1, the cross-section is subdivided into
             nwinc*nhinc parallel sub-filaments during build(). Each
             sub-filament has width/nwinc x height/nhinc cross-section
             and shares the same node_from/node_to (parallel connection).
             The mutual inductance between sub-filaments naturally
             captures skin and proximity effects.

             Args:
                 node_from: Source node ID
                 node_to: Destination node ID
                 width: Cross-section width [m]
                 height: Cross-section height [m]
                 sigma: Conductivity [S/m] (default: copper 5.8e7)
                 cross_section_type: 0 = rectangular, 1 = circular
                 nwinc: Width subdivisions (default: 1, no subdivision)
                 nhinc: Height subdivisions (default: 1, no subdivision)

             Example:
                 n1 = builder.add_node_at(0, 0, 0)
                 n2 = builder.add_node_at(0.1, 0, 0)
                 # Single filament
                 builder.add_connected_segment(n1, n2, 1e-3, 1e-3)
                 # 3x3 = 9 sub-filaments (solid conductor)
                 builder.add_connected_segment(n1, n2, 3e-3, 3e-3, nwinc=3, nhinc=3)
             )doc")

        .def("add_port", &PyPEECBuilder::add_port,
             py::arg("node_positive"), py::arg("node_negative"),
             R"doc(
             Add a port between two nodes.

             Args:
                 node_positive: Positive terminal node ID
                 node_negative: Negative terminal node ID

             Returns:
                 Port ID (auto-incrementing, starting from 0)
             )doc")

        .def("build_topology", &PyPEECBuilder::build_topology,
             py::arg("include_star") = false,
             py::arg("eps_eff") = 1.0,
             R"doc(
             Build PEEC matrices with topology information.

             Args:
                 include_star: If True, compute P and M_LS matrices
                 eps_eff: Effective permittivity for gathered capacitance.
                     For PCB half-space model: eps_eff = (1 + eps_r) / 2.
                     If > 0 and panels present, also computes C_gathered.

             Returns a dict containing:
                 L: Inductance matrix [H] (n_loop x n_loop)
                 R: DC resistance vector [Ohm] (n_loop,)
                 P: Potential coefficient [1/F] or None
                 C_gathered: Gathered capacitance (n_nodes x n_nodes) or None
                 incidence_indptr: CSR row pointers for incidence matrix
                 incidence_indices: CSR column indices
                 incidence_data: CSR values (+1 or -1)
                 n_junction: Number of internal junction nodes
                 ports: List of (node_positive, node_negative, port_id)
                 n_loop: Number of filaments
                 n_star: Number of star elements
             )doc")

        .def_property_readonly("n_loop", &PyPEECBuilder::n_loop,
            "Number of Loop elements (segments)")
        .def_property_readonly("n_star", &PyPEECBuilder::n_star,
            "Number of Star elements (nodes)")
        .def_property_readonly("num_segments", &PyPEECBuilder::num_segments,
            "Number of segments")
        .def_property_readonly("num_nodes", &PyPEECBuilder::num_nodes,
            "Number of nodes")
        .def_property_readonly("num_panels", &PyPEECBuilder::num_panels,
            "Number of panels")

        // MNA Multi-Port Solver (C++ LAPACK zgesv_/zgetrf_/zgetrs_)
        .def("solve_z_matrix", &PyPEECBuilder::solve_z_matrix,
             py::arg("freq"),
             py::arg("Zs") = py::none(),
             R"doc(
             Compute multi-port Z-parameter matrix at a single frequency.

             Uses LAPACK zgetrf_/zgetrs_.

             Args:
                 freq: Frequency in Hz
                 Zs: Optional surface impedance array (n_loop,) complex.
                      If None, only DC resistance is used.

             Returns:
                 Z_matrix: (n_ports, n_ports) complex numpy array

             Example:
                 builder = PEECBuilder()
                 n1 = builder.add_node_at(0, 0, 0)
                 n2 = builder.add_node_at(0.1, 0, 0)
                 builder.add_connected_segment(n1, n2, 1e-3, 1e-3)
                 builder.add_port(n1, n2)
                 Z = builder.solve_z_matrix(1e6)  # (1, 1) complex array
             )doc")

        .def("solve_coupling", &PyPEECBuilder::solve_coupling,
             py::arg("freq"),
             py::arg("Zs") = py::none(),
             R"doc(
             Compute coupling coefficient from Z-parameters.

             Convenience method that extracts L, M, and k from the Z-matrix.

             Args:
                 freq: Frequency in Hz
                 Zs: Optional surface impedance array (n_loop,) complex

             Returns:
                 dict with keys:
                   Z_matrix: (n_ports, n_ports) complex numpy array
                   L_matrix: (n_ports, n_ports) inductance [H]
                   k_matrix: (n_ports, n_ports) coupling coefficients
                   k: Scalar coupling coefficient k12
                   L1: Self-inductance of port 1 [H]
                   L2: Self-inductance of port 2 [H]
                   M: Mutual inductance M12 [H]
             )doc")

        .def("solve_multiport_sweep", &PyPEECBuilder::solve_multiport_sweep,
             py::arg("freqs"),
             py::arg("Zs_all") = py::none(),
             R"doc(
             Batch frequency sweep for multi-port Z-parameters.

             Computes Z-matrix at each frequency in the sweep array.

             Args:
                 freqs: Frequency array [Hz], shape (n_freq,)
                 Zs_all: Optional (n_freq, n_loop) complex array of per-frequency
                          surface impedance. If None, only DC resistance is used.

             Returns:
                 dict with keys:
                   freqs: Input frequency array
                   Z_matrix: (n_freq, n_ports, n_ports) complex array
                   L_matrix: (n_freq, n_ports, n_ports) inductance [H]
                   k: (n_freq,) coupling coefficient k12
                   L1: (n_freq,) self-inductance of port 1 [H]
                   L2: (n_freq,) self-inductance of port 2 [H]
                   M: (n_freq,) mutual inductance M12 [H]
             )doc");

    // ========== Standalone MNA Solver (no PEECBuilder dependency) ==========

    py::class_<PyMNASolver>(m, "MNASolver",
        R"doc(
        Standalone MNA (Modified Nodal Analysis) Solver.

        Accepts raw L, R, segment_nodes, ports arrays directly (no PEECBuilder).
        Uses LAPACK zgesv_/zgetrf_/zgetrs_.

        Used internally by PEECCircuitSolver for all MNA linear algebra.
        Also usable directly by ngbem bridge and other callers.
        )doc")
        .def(py::init<py::array_t<double>, py::array_t<double>,
                       py::array_t<int>, int, py::list,
                       py::object, py::object, py::object>(),
             py::arg("L"),
             py::arg("R"),
             py::arg("segment_nodes"),
             py::arg("n_nodes"),
             py::arg("ports"),
             py::arg("P") = py::none(),
             py::arg("M_LS") = py::none(),
             py::arg("C_gathered") = py::none(),
             R"doc(
             Create MNA solver from raw matrices.

             Args:
                 L: Inductance matrix (n_loop x n_loop) float64
                 R: DC resistance vector (n_loop,) float64
                 segment_nodes: Connectivity (n_seg x 2) int32
                 n_nodes: Total number of nodes
                 ports: List of (node_positive, node_negative, port_id)
                 P: Optional potential coefficient (n_star x n_star)
                 M_LS: Optional Loop-Star coupling (n_loop x n_star)
                 C_gathered: Optional gathered capacitance (n_nodes x n_nodes).
                     Enables proper MNA with resonance.
             )doc")

        .def("solve_z_matrix", &PyMNASolver::solve_z_matrix,
             py::arg("freq"),
             py::arg("Zs") = py::none(),
             R"doc(
             Compute multi-port Z-parameter matrix at a single frequency.

             Args:
                 freq: Frequency in Hz
                 Zs: Optional surface impedance (n_loop,) complex

             Returns:
                 Z_matrix: (n_ports, n_ports) complex numpy array
             )doc")

        .def("frequency_sweep", &PyMNASolver::frequency_sweep,
             py::arg("freqs"),
             py::arg("Zs_all") = py::none(),
             R"doc(
             Batch frequency sweep for multi-port Z-parameters.

             Args:
                 freqs: Frequency array [Hz]
                 Zs_all: Optional (n_freq, n_loop) complex array

             Returns:
                 dict with 'freqs' and 'Z_matrix' (n_freq, n_ports, n_ports)
             )doc")

        .def("solve_branch_currents", &PyMNASolver::solve_branch_currents,
             py::arg("freq"),
             py::arg("port_currents"),
             py::arg("Zs") = py::none(),
             R"doc(
             Solve for per-filament complex branch currents under user-specified
             port current injection.

             Args:
                 freq: Frequency in Hz
                 port_currents: Complex currents injected at each port (n_ports,)
                 Zs: Optional per-filament surface impedance (n_loop,) complex

             Returns:
                 I_branch: (n_loop,) complex numpy array. I_branch[f] > 0 means
                     current flows from node_from[f] to node_to[f]. For a bundle
                     of parallel filaments between port+ and port-, sum(I_branch)
                     equals port_currents[0].
             )doc")

        .def_property_readonly("num_ports", &PyMNASolver::num_ports,
            "Number of ports")

        .def("set_solver_method", &PyMNASolver::set_solver_method,
             py::arg("method"),
             R"doc(
             Set solver method for Z_branch inversion.

             Args:
                 method: 0 = LU (default, LAPACK zgesv_), 1 = BiCGSTAB
             )doc")

        .def("get_solver_method", &PyMNASolver::get_solver_method,
             "Get current solver method (0=LU, 1=BiCGSTAB)")

        .def("set_bicgstab_params", &PyMNASolver::set_bicgstab_params,
             py::arg("tol") = 1e-10,
             py::arg("max_iter") = 1000,
             R"doc(
             Set BiCGSTAB solver parameters.

             Args:
                 tol: Convergence tolerance (default 1e-10)
                 max_iter: Maximum iterations (default 1000)
             )doc")

        .def("set_gathered_capacitance", &PyMNASolver::set_gathered_capacitance,
             py::arg("C_gathered"),
             R"doc(
             Set gathered capacitance for proper MNA.

             When set, MNA uses: Y_node = A * Y_branch^{-1} * A^T + jw * C_gathered
             instead of the Schur complement, enabling resonance (SRF) prediction.

             Args:
                 C_gathered: Gathered capacitance matrix (n_nodes x n_nodes) float64.
                     C_gathered = G * inv(P / eps_eff) * G^T
             )doc")

        .def_property_readonly("has_gathered_capacitance",
             &PyMNASolver::has_gathered_capacitance,
             "Whether gathered capacitance is set (proper MNA mode)");

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
