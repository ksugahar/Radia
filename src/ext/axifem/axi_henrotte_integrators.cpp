// axi_henrotte_integrators.cpp — Closed-form BilinearFormIntegrators.
//
// Element-matrix formulas ported from the validated Python fixtures under
// tests/axifem/_reference_python/ (research copies are retained under
// validation_test/axifem/research/prototypes/):
//   axifem_quad.py
//     _element_matrices_quad_closed_form  (Q1 stiffness)
//     element_sigma_mass_quad             (Q1 sigma mass)
// P1/P2 triangles and curved Q2 quads use the s-isoparametric quadrature of
// axi_henrotte_siso.hpp.
//
// DOF convention (V-DOF): u_j = A_phi at vertex j. This matches the Python
// references (validation_test/axifem/research/verification/disk_axisym_reference.py).
// - Q1 closed form: M_V = T M_phi T with T_jj = 2 pi r_j; axis-vertex rows/
//   cols become zero (caller MUST Dirichlet axis edges).
// - Triangles / curved quads: V-DOF integrals with N_A_j = r_j psi_j / r over
//   the (s, z) element map.
//
// The user MUST apply Dirichlet on every axis edge ("axis|...") for axis-
// touching elements; otherwise the V-DOF zero rows/cols make the system
// singular.

#include <comp.hpp>
#include <python_comp.hpp>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cmath>
#include "axi_henrotte_integrators.hpp"
#include "axi_henrotte_fe.hpp"
#include "axi_henrotte_numeric.hpp"
#include "q2_henrotte_generated.hpp"

namespace axifem {

using namespace ngfem;
using std::log;
using std::abs;

namespace {

constexpr double PI = 3.14159265358979323846;

// ---------------------------------------------------------------------------
// Sample a CoefficientFunction at the element centroid (constant per element).
// ---------------------------------------------------------------------------
double SampleAtCentroid(const CoefficientFunction & cf,
                        const ElementTransformation & eltrans,
                        ELEMENT_TYPE et, LocalHeap & lh)
{
    IntegrationPoint ip_center;
    if (et == ET_QUAD)      ip_center = IntegrationPoint(0.5, 0.5, 0.0);
    else /* ET_TRIG */      ip_center = IntegrationPoint(1.0/3.0, 1.0/3.0, 0.0);

    HeapReset hr(lh);
    auto & mip = eltrans(ip_center, lh);
    return cf.Evaluate(mip);
}

void CopyMatrix4(const numeric::Matrix4 & source, FlatMatrix<double> target)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            target(i, j) = source[4 * i + j];
}

void CopyMatrix9(const numeric::Matrix9 & source, FlatMatrix<double> target)
{
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            target(i, j) = source[9 * i + j];
}

void Q1InverseVandermonde(double ra, double rb, double za, double zb,
                          Mat<4,4> & inverse)
{
    const auto values =
        numeric::ComputeQ1InverseVandermonde(ra, rb, za, zb);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            inverse(i, j) = values[4 * i + j];
}

// ---------------------------------------------------------------------------
// Q2 element matrix (interior or axis-touching, 9 DOFs in V-DOF basis).
//
//   M_phi (monomial)  =  KPhiGeneral / MSigmaPhiGeneral / KPhiAxis / ...
//                        (delivered by q2_henrotte_generated.hpp)
//   M_node (Lagrange) =  Vinv^T * M_phi * Vinv
//   M_V (psi-DOF)     =  T * M_node * T,    T_jj = 2 pi r_node_j
//
// For axis elements the axis nodes (local indices 0, 3, 7) get zero rows and
// columns — caller MUST Dirichlet those global DOFs.
// ---------------------------------------------------------------------------

// Returns r-coord at each of the 9 local nodes (sm convention).
inline void Q2NodeRCoords(double ra, double rb, double r_node[9])
{
    double sa = ra * ra, sb = rb * rb;
    double rm = std::sqrt(0.5 * (sa + sb));   // r at s-midpoint = sqrt((ra^2+rb^2)/2)
    r_node[0] = ra;
    r_node[1] = rb;
    r_node[2] = rb;
    r_node[3] = ra;
    r_node[4] = rm;
    r_node[5] = rb;
    r_node[6] = rm;
    r_node[7] = ra;
    r_node[8] = rm;
}

// Generic V-DOF transform driver.
//   For interior:  9x9 monomial matrix M_phi[81] -> 9x9 V-DOF matrix elmat.
//   For axis    :  6x6 monomial matrix M_phi[36] embedded into 9x9 elmat.
//                  Uses fe.Vinv[k=0..5][j_local] for the 6-monomial axis basis.
template <int N_PHI>
void Q2MonomialToVDof(const AxiHenrotteFE_Q2_AxisAligned & fe,
                      const double M_phi_flat[],
                      FlatMatrix<double> elmat)
{
    static_assert(N_PHI == 9 || N_PHI == 6, "Q2MonomialToVDof: N_PHI must be 6 or 9");

    elmat = 0.0;

    double r_node[9];
    Q2NodeRCoords(fe.r_a, fe.r_b, r_node);
    double T[9];
    for (int j = 0; j < 9; ++j) T[j] = 2.0 * PI * r_node[j];

    // M_phi flat -> M_phi[N_PHI][N_PHI]
    double Mphi[N_PHI][N_PHI];
    for (int i = 0; i < N_PHI; ++i)
        for (int j = 0; j < N_PHI; ++j)
            Mphi[i][j] = M_phi_flat[i * N_PHI + j];

    // Step 1: tmp[k, j_local] = sum_l M_phi[k, l] * Vinv[l, j_local]
    // (Vinv stored as Vinv[mono_index, local_node_index]; for axis case only
    // the first 6 rows of Vinv are non-zero, indices into nz_idx.)
    int n_active = fe.n_nz;
    double tmp[N_PHI][9];
    for (int k = 0; k < N_PHI; ++k) {
        for (int q = 0; q < n_active; ++q) {
            int j = fe.nz_idx[q];
            double s = 0.0;
            for (int l = 0; l < N_PHI; ++l) s += Mphi[k][l] * fe.Vinv[l][j];
            tmp[k][j] = s;
        }
    }

    // Step 2: M_node[i_local, j_local] = sum_k Vinv[k, i_local] * tmp[k, j_local]
    // Then M_V[i, j] = T[i] * M_node[i, j] * T[j].
    for (int p = 0; p < n_active; ++p) {
        int i = fe.nz_idx[p];
        for (int q = 0; q < n_active; ++q) {
            int j = fe.nz_idx[q];
            double s = 0.0;
            for (int k = 0; k < N_PHI; ++k) s += fe.Vinv[k][i] * tmp[k][j];
            elmat(i, j) = T[i] * s * T[j];
        }
    }

    // Symmetrise (numerical noise).
    for (int i = 0; i < 9; ++i)
        for (int j = i + 1; j < 9; ++j) {
            double avg = 0.5 * (elmat(i, j) + elmat(j, i));
            elmat(i, j) = avg;
            elmat(j, i) = avg;
        }
}

void Q2StiffnessElement(const AxiHenrotteFE_Q2_AxisAligned & fe, double mu,
                        FlatMatrix<double> elmat)
{
    CopyMatrix9(numeric::ComputeQ2MagneticElementMatrices(
        fe.r_a, fe.r_b, fe.z_a, fe.z_b, mu, 0.0).stiffness, elmat);
}

void Q2SigmaMassElement(const AxiHenrotteFE_Q2_AxisAligned & fe, double sigma,
                        FlatMatrix<double> elmat)
{
    CopyMatrix9(numeric::ComputeQ2MagneticElementMatrices(
        fe.r_a, fe.r_b, fe.z_a, fe.z_b, 1.0, sigma).sigma_mass, elmat);
}

}  // anonymous namespace

// ---------------------------------------------------------------------------
// AxiHenrotteStiffnessBFI
// ---------------------------------------------------------------------------

AxiHenrotteStiffnessBFI::AxiHenrotteStiffnessBFI(
    shared_ptr<CoefficientFunction> amu_cf)
  : mu_cf(amu_cf)
{ SetName("AxiHenrotteStiffnessBFI"); }

void AxiHenrotteStiffnessBFI::CalcElementMatrix(
    const FiniteElement & fel,
    const ElementTransformation & eltrans,
    FlatMatrix<double> elmat,
    LocalHeap & lh) const
{
    elmat = 0.0;
    ELEMENT_TYPE et = fel.ElementType();
    double mu = SampleAtCentroid(*mu_cf, eltrans, et, lh);

    if (auto * q1 = dynamic_cast<const AxiHenrotteFE_Q1_AxisAligned*>(&fel)) {
        CopyMatrix4(numeric::ComputeQ1MagneticStiffness(
                        q1->r_a, q1->r_b, q1->z_a, q1->z_b, mu),
                    elmat);
        return;
    }
    if (auto * q2 = dynamic_cast<const AxiHenrotteFE_Q2_AxisAligned*>(&fel)) {
        Q2StiffnessElement(*q2, mu, elmat);
        return;
    }
    if (auto * q2c = dynamic_cast<const AxiHenrotteFE_Q2_Curved*>(&fel)) {
        q2c->geo.Stiffness(mu, elmat);
        return;
    }
    if (auto * p1 = dynamic_cast<const AxiHenrotteFE_P1_Triangle*>(&fel)) {
        p1->geo.Stiffness(mu, elmat);
        return;
    }
    if (auto * p2 = dynamic_cast<const AxiHenrotteFE_P2_Triangle*>(&fel)) {
        p2->geo.Stiffness(mu, elmat);
        return;
    }
    throw Exception(string("AxiHenrotteStiffnessBFI: unsupported FE type ")
                    + typeid(fel).name());
}

// ---------------------------------------------------------------------------
// AxiHenrotteSigmaMassBFI
// ---------------------------------------------------------------------------

AxiHenrotteSigmaMassBFI::AxiHenrotteSigmaMassBFI(
    shared_ptr<CoefficientFunction> asigma_cf)
  : sigma_cf(asigma_cf)
{ SetName("AxiHenrotteSigmaMassBFI"); }

void AxiHenrotteSigmaMassBFI::CalcElementMatrix(
    const FiniteElement & fel,
    const ElementTransformation & eltrans,
    FlatMatrix<double> elmat,
    LocalHeap & lh) const
{
    elmat = 0.0;
    ELEMENT_TYPE et = fel.ElementType();
    double sigma = SampleAtCentroid(*sigma_cf, eltrans, et, lh);
    if (sigma == 0.0) return;  // air / non-conductor: zero contribution

    if (auto * q1 = dynamic_cast<const AxiHenrotteFE_Q1_AxisAligned*>(&fel)) {
        CopyMatrix4(numeric::ComputeQ1SigmaMass(
                        q1->r_a, q1->r_b, q1->z_a, q1->z_b, sigma),
                    elmat);
        return;
    }
    if (auto * q2 = dynamic_cast<const AxiHenrotteFE_Q2_AxisAligned*>(&fel)) {
        Q2SigmaMassElement(*q2, sigma, elmat);
        return;
    }
    if (auto * q2c = dynamic_cast<const AxiHenrotteFE_Q2_Curved*>(&fel)) {
        q2c->geo.SigmaMass(sigma, elmat);
        return;
    }
    if (auto * p1 = dynamic_cast<const AxiHenrotteFE_P1_Triangle*>(&fel)) {
        p1->geo.SigmaMass(sigma, elmat);
        return;
    }
    if (auto * p2 = dynamic_cast<const AxiHenrotteFE_P2_Triangle*>(&fel)) {
        p2->geo.SigmaMass(sigma, elmat);
        return;
    }
    throw Exception(string("AxiHenrotteSigmaMassBFI: unsupported FE type ")
                    + typeid(fel).name());
}

// ---------------------------------------------------------------------------
// Python bindings
// ---------------------------------------------------------------------------

void ExportAxiHenrotteIntegrators(pybind11::module & m)
{
    namespace py = pybind11;

    m.def(
        "q1_magnetic_element_matrices",
        [](double ra, double rb, double za, double zb,
           double mu, double sigma) {
            const auto matrices = numeric::ComputeQ1MagneticElementMatrices(
                ra, rb, za, zb, mu, sigma);
            auto as_numpy = [](const numeric::Matrix4 & values) {
                py::array_t<double> result({
                    static_cast<py::ssize_t>(4), static_cast<py::ssize_t>(4)});
                auto view = result.mutable_unchecked<2>();
                for (py::ssize_t i = 0; i < 4; ++i)
                    for (py::ssize_t j = 0; j < 4; ++j)
                        view(i, j) = values[4 * i + j];
                return result;
            };
            py::dict result;
            result["stiffness"] = as_numpy(matrices.stiffness);
            result["sigma_mass"] = as_numpy(matrices.sigma_mass);
            result["backend"] = "native-pybind";
            result["dof_convention"] = "nodal A_phi (V-DOF)";
            result["node_order"] =
                "(ra,za),(rb,za),(rb,zb),(ra,zb)";
            return result;
        },
        py::arg("ra"), py::arg("rb"), py::arg("za"), py::arg("zb"),
        py::arg("mu"), py::arg("sigma"),
        "Return shared-native Q1 Henrotte stiffness and sigma-mass matrices.");

    m.def(
        "q2_magnetic_element_matrices",
        [](double ra, double rb, double za, double zb,
           double mu, double sigma) {
            const auto matrices = numeric::ComputeQ2MagneticElementMatrices(
                ra, rb, za, zb, mu, sigma);
            auto as_numpy = [](const numeric::Matrix9 & values) {
                py::array_t<double> result({
                    static_cast<py::ssize_t>(9), static_cast<py::ssize_t>(9)});
                auto view = result.mutable_unchecked<2>();
                for (py::ssize_t i = 0; i < 9; ++i)
                    for (py::ssize_t j = 0; j < 9; ++j)
                        view(i, j) = values[9 * i + j];
                return result;
            };
            py::dict result;
            result["stiffness"] = as_numpy(matrices.stiffness);
            result["sigma_mass"] = as_numpy(matrices.sigma_mass);
            result["backend"] = "native-pybind";
            result["dof_convention"] = "nodal A_phi (V-DOF)";
            result["node_order"] =
                "(ra,za),(rb,za),(rb,zb),(ra,zb),(rm,za),(rb,zm),(rm,zb),(ra,zm),(rm,zm)";
            result["axis_touching"] = matrices.axis_touching;
            return result;
        },
        py::arg("ra"), py::arg("rb"), py::arg("za"), py::arg("zb"),
        py::arg("mu"), py::arg("sigma"),
        "Return shared-native Q2 Henrotte stiffness and sigma-mass matrices.");

    py::class_<AxiHenrotteStiffnessBFI, BilinearFormIntegrator,
               shared_ptr<AxiHenrotteStiffnessBFI>>(
        m, "AxiHenrotteStiffnessBFI",
        "Closed-form axisymmetric stiffness integrator (Henrotte basis).\n"
        "Bypasses NGSolve quadrature; uses Mathematica-derived element matrices\n"
        "exact for the 1/r-weighted axisymmetric weak form.")
        .def(py::init<shared_ptr<CoefficientFunction>>(),
             py::arg("mu"));

    py::class_<AxiHenrotteSigmaMassBFI, BilinearFormIntegrator,
               shared_ptr<AxiHenrotteSigmaMassBFI>>(
        m, "AxiHenrotteSigmaMassBFI",
        "Closed-form axisymmetric sigma-mass integrator (Henrotte basis).\n"
        "Axis-aligned Q1/Q2 quads use closed forms; P1/P2 triangles and general\n"
        "Q2 quads use the s-isoparametric quadrature (axi_henrotte_siso.hpp).")
        .def(py::init<shared_ptr<CoefficientFunction>>(),
             py::arg("sigma"));
}

}  // namespace axifem
