// axi_henrotte_integrators.hpp — Closed-form BilinearFormIntegrators that
// bypass NGSolve's Gauss quadrature for the FEMM/Henrotte axisymmetric weak
// form. The 1/r weight in the axisymmetric integrals would otherwise lose
// 5-10% accuracy when integrated by Gauss quadrature near the symmetry axis;
// here we use Mathematica-derived closed forms instead (port of
// W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/axifemm/axifemm_quad.py and
// axifemm_core.py).
//
// Two integrators are exposed:
//
//   AxiHenrotteStiffnessBFI(mu_cf)
//       elmat_ij = (1/(2 pi mu)) * Integrate( grad(shape_i) . grad(shape_j) / r dx )
//
//   AxiHenrotteSigmaMassBFI(sigma_cf)
//       elmat_ij = (sigma/(2 pi))  * Integrate( shape_i shape_j / r dx )
//
// Conventions match Phase 1b Python (DOF u_i = psi_i = 2 pi r_i A_phi(r_i, z_i)),
// so the eigenvalue problem K v = lambda M v gives the same tau_n as the
// validated reference (axifemm/disk_hiruma_quad.py).
//
// Material coefficients (mu, sigma) are sampled at the element centroid and
// treated as element-constant; the closed forms assume per-element mu, sigma.
// For curved or graded materials use SymbolicBFI instead.

#ifndef RADIA_AXIFEMM_AXI_HENROTTE_INTEGRATORS_HPP
#define RADIA_AXIFEMM_AXI_HENROTTE_INTEGRATORS_HPP

#include <fem.hpp>
#include <integrator.hpp>

namespace axifem {

using namespace ngfem;

// ---------------------------------------------------------------------------
// AxiHenrotteStiffnessBFI
// ---------------------------------------------------------------------------

class AxiHenrotteStiffnessBFI : public BilinearFormIntegrator {
public:
    shared_ptr<CoefficientFunction> mu_cf;

    AxiHenrotteStiffnessBFI(shared_ptr<CoefficientFunction> amu_cf);

    xbool IsSymmetric() const override { return true; }
    VorB VB() const override { return VOL; }
    int DimElement() const override { return 2; }
    int DimSpace()   const override { return 2; }
    string Name() const override { return "AxiHenrotteStiffnessBFI"; }

    void CalcElementMatrix(
        const FiniteElement & fel,
        const ElementTransformation & eltrans,
        FlatMatrix<double> elmat,
        LocalHeap & lh) const override;
};

// ---------------------------------------------------------------------------
// AxiHenrotteSigmaMassBFI
// ---------------------------------------------------------------------------

class AxiHenrotteSigmaMassBFI : public BilinearFormIntegrator {
public:
    shared_ptr<CoefficientFunction> sigma_cf;

    AxiHenrotteSigmaMassBFI(shared_ptr<CoefficientFunction> asigma_cf);

    xbool IsSymmetric() const override { return true; }
    VorB VB() const override { return VOL; }
    int DimElement() const override { return 2; }
    int DimSpace()   const override { return 2; }
    string Name() const override { return "AxiHenrotteSigmaMassBFI"; }

    void CalcElementMatrix(
        const FiniteElement & fel,
        const ElementTransformation & eltrans,
        FlatMatrix<double> elmat,
        LocalHeap & lh) const override;
};

// ---------------------------------------------------------------------------
// AxiHenrotteHeatStiffnessBFI -- axisymmetric heat conduction stiffness
// ---------------------------------------------------------------------------
//
// Closed-form heat stiffness for the axisymmetric Henrotte basis.  The
// weak form (after substitution s = r^2, ds = 2r dr, 2 pi r dr dz =
// pi ds dz) is
//
//   a(T, v) = pi * Integrate( k * [ 4 s d_s T d_s v + d_z T d_z v ] ds dz )
//
// Both terms are POLYNOMIAL in (s, z) -- no 1/s factor (unlike the
// magnetic case).  Element matrices come from sympy-derived closed forms
// in q_heat_henrotte_generated.hpp.
//
// DOF semantics: nodal temperature T_j (NOT the magnetic flux psi_j =
// 2 pi r_j A_phi).  No T = diag(2 pi r_node) wrap is applied -- T is
// regular at the axis and stays in nodal form.  This is the key
// difference from AxiHenrotteStiffnessBFI.
class AxiHenrotteHeatStiffnessBFI : public BilinearFormIntegrator {
public:
    shared_ptr<CoefficientFunction> k_cf;

    AxiHenrotteHeatStiffnessBFI(shared_ptr<CoefficientFunction> ak_cf);

    xbool IsSymmetric() const override { return true; }
    VorB VB() const override { return VOL; }
    int DimElement() const override { return 2; }
    int DimSpace()   const override { return 2; }
    string Name() const override { return "AxiHenrotteHeatStiffnessBFI"; }

    void CalcElementMatrix(
        const FiniteElement & fel,
        const ElementTransformation & eltrans,
        FlatMatrix<double> elmat,
        LocalHeap & lh) const override;
};

// ---------------------------------------------------------------------------
// AxiHenrotteHeatMassBFI -- axisymmetric heat capacity (transient term)
// ---------------------------------------------------------------------------
//
//   m(T, v) = pi * Integrate( rho_c * T * v ) ds dz
//
// rho_c = rho * c_p [J/(m^3 K)]; pass via CoefficientFunction.  Same
// nodal-T DOF semantics as AxiHenrotteHeatStiffnessBFI.
class AxiHenrotteHeatMassBFI : public BilinearFormIntegrator {
public:
    shared_ptr<CoefficientFunction> rho_c_cf;

    AxiHenrotteHeatMassBFI(shared_ptr<CoefficientFunction> arho_c_cf);

    xbool IsSymmetric() const override { return true; }
    VorB VB() const override { return VOL; }
    int DimElement() const override { return 2; }
    int DimSpace()   const override { return 2; }
    string Name() const override { return "AxiHenrotteHeatMassBFI"; }

    void CalcElementMatrix(
        const FiniteElement & fel,
        const ElementTransformation & eltrans,
        FlatMatrix<double> elmat,
        LocalHeap & lh) const override;
};

}  // namespace axifem

#endif
