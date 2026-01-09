/*
 * rad_conductor.h
 *
 * Conductor element classes for electromagnetic analysis
 * Supports DC (analytical), MQS, and Full-wave formulations
 *
 * Compatible with Radia geometry (ObjRecMag, ObjHexahedron, etc.)
 *
 * References:
 * [1] Z. Zhu et al., "Algorithms in FastImp", IEEE TCAD, 2005
 * [2] S. Bilicz et al., "Nonlocal SIBC", ISEM 2023
 *
 * Part of Radia project
 */

#ifndef RAD_CONDUCTOR_H
#define RAD_CONDUCTOR_H

#include "gmvect.h"
#include "rad_pfft.h"
#include "rad_constants.h"
#include "rad_rwg_basis.h"
#include <vector>
#include <array>
#include <complex>
#include <memory>
#include <functional>

namespace radia {

// Physical constants for conductor analysis
constexpr double MU_0 = RadConst::MU_0;               // Permeability of free space [H/m]
constexpr double EPS_0 = 8.854187817e-12;             // Permittivity of free space [F/m]
constexpr double INV_FOUR_PI = RadConst::INV_FOUR_PI; // 1/(4*pi)

// Forward declarations
class radTg3d;

/**
 * @brief Formulation type for conductor analysis
 */
enum class ConductorFormulation {
    DC,         // DC: Use existing Radia analytical formulas (Biot-Savart)
    MQS,        // Magneto-Quasi-Static: G(r) = 1/(4*pi*r), no displacement current
    EMQS,       // Electro-Magneto-Quasi-Static: Low frequency with capacitance
    FullWave    // Full-wave: G(r) = exp(-jkr)/(4*pi*r), includes resonance
};

/**
 * @brief Surface panel for integral equation discretization
 */
struct SurfacePanel {
    TVector3d center;         // Panel center
    TVector3d normal;         // Outward normal
    double area;              // Panel area
    std::vector<TVector3d> vertices;  // Panel vertices (3 or 4)

    // Panel type
    enum Type { Triangle, Quadrilateral } type;

    SurfacePanel() : area(0), type(Quadrilateral) {
        center = TVector3d(0, 0, 0);
        normal = TVector3d(0, 0, 1);
    }
};

/**
 * @brief Conductor element for IE formulation
 *
 * This class represents a conductor region discretized into surface panels.
 * It can be created from existing Radia objects (ObjRecMag, ObjHexahedron, etc.)
 * or from external mesh data.
 */
class radTConductor {
public:
    radTConductor();
    ~radTConductor();

    /**
     * @brief Set conductor material properties
     * @param conductivity Electrical conductivity [S/m]
     */
    void SetConductivity(double conductivity);

    /**
     * @brief Get conductivity
     */
    double GetConductivity() const { return conductivity_; }

    /**
     * @brief Set analysis formulation
     */
    void SetFormulation(ConductorFormulation form) { formulation_ = form; }
    ConductorFormulation GetFormulation() const { return formulation_; }

    /**
     * @brief Set frequency for AC analysis
     * @param frequency Frequency in Hz (0 for DC)
     */
    void SetFrequency(double frequency);
    double GetFrequency() const { return frequency_; }

    // ========== Geometry creation (Radia-compatible) ==========

    /**
     * @brief Create conductor from rectangular block
     * Similar to Radia ObjRecMag but for conductor
     * @param center Block center [x, y, z]
     * @param dimensions Block dimensions [Lx, Ly, Lz]
     * @param numPanelsPerFace Number of panels per face direction
     */
    void CreateFromRecBlock(const TVector3d& center,
                            const TVector3d& dimensions,
                            int numPanelsPerFace = 4);

    /**
     * @brief Create conductor from hexahedron vertices
     * Similar to Radia ObjHexahedron
     * @param vertices 8 vertices of hexahedron
     * @param numPanelsPerFace Number of panels per face direction
     */
    void CreateFromHexahedron(const std::vector<TVector3d>& vertices,
                              int numPanelsPerFace = 4);

    /**
     * @brief Create conductor from existing Radia object
     * @param radiaObj Pointer to existing Radia geometry object
     * @param numPanelsPerFace Surface discretization density
     */
    void CreateFromRadiaObject(radTg3d* radiaObj, int numPanelsPerFace = 4);

    /**
     * @brief Create conductor from surface mesh
     * @param vertices Mesh vertices
     * @param triangles Triangle connectivity (3 indices per triangle)
     */
    void CreateFromSurfaceMesh(const std::vector<TVector3d>& vertices,
                               const std::vector<std::array<int, 3>>& triangles);

    /**
     * @brief Create conductor from quadrilateral panels
     * @param panels List of surface panels
     */
    void CreateFromPanels(const std::vector<SurfacePanel>& panels);

    // ========== Wire/Coil creation ==========

    /**
     * @brief Create wire conductor along a path
     * @param path Center line of wire
     * @param crossSection Cross-section type ("circular" or "rectangular")
     * @param width Width (or diameter for circular)
     * @param height Height (ignored for circular)
     * @param numPanelsAround Number of panels around circumference
     * @param numPanelsAlongLength Panels along wire length (0 = auto)
     */
    void CreateWire(const std::vector<TVector3d>& path,
                    const std::string& crossSection,
                    double width,
                    double height = 0,
                    int numPanelsAround = 8,
                    int numPanelsAlongLength = 0);

    /**
     * @brief Create circular loop coil
     * @param center Loop center
     * @param radius Loop radius
     * @param normal Loop normal direction
     * @param crossSection Wire cross-section
     * @param wireWidth Wire width
     * @param wireHeight Wire height
     * @param numPanelsAround Panels around wire
     * @param numPanelsLoop Panels along loop
     */
    void CreateLoop(const TVector3d& center,
                    double radius,
                    const TVector3d& normal,
                    const std::string& crossSection,
                    double wireWidth,
                    double wireHeight,
                    int numPanelsAround = 8,
                    int numPanelsLoop = 30);

    /**
     * @brief Create spiral coil
     */
    void CreateSpiral(const TVector3d& center,
                      double innerRadius,
                      double outerRadius,
                      double pitch,
                      int numTurns,
                      const TVector3d& axis,
                      const std::string& crossSection,
                      double wireWidth,
                      double wireHeight,
                      int numPanelsAround = 8);

    // ========== Access to discretization ==========

    /**
     * @brief Get surface panels
     */
    const std::vector<SurfacePanel>& GetPanels() const { return panels_; }

    /**
     * @brief Get number of panels
     */
    int NumPanels() const { return static_cast<int>(panels_.size()); }

    /**
     * @brief Get panel centers
     */
    std::vector<TVector3d> GetPanelCenters() const;

    /**
     * @brief Get panel areas
     */
    std::vector<double> GetPanelAreas() const;

    // ========== Solution data ==========

    /**
     * @brief Surface current density K [A/m]
     * Complex for AC analysis
     */
    std::vector<std::complex<double>>& SurfaceCurrent() { return surfaceCurrent_; }
    const std::vector<std::complex<double>>& SurfaceCurrent() const { return surfaceCurrent_; }

    /**
     * @brief Surface charge density sigma [C/m^2]
     * Complex for AC analysis
     */
    std::vector<std::complex<double>>& SurfaceCharge() { return surfaceCharge_; }
    const std::vector<std::complex<double>>& SurfaceCharge() const { return surfaceCharge_; }

    // ========== Field computation ==========

    /**
     * @brief Compute magnetic field B at a point
     * Uses appropriate formulation (DC/MQS/Full-wave)
     * @param point Evaluation point
     * @return B field (complex for AC)
     */
    std::complex<double> ComputeB(const TVector3d& point, int component) const;

    /**
     * @brief Compute all 3 components of B field
     */
    void ComputeB(const TVector3d& point,
                  std::complex<double>& Bx,
                  std::complex<double>& By,
                  std::complex<double>& Bz) const;

    /**
     * @brief Compute electric field E at a point
     */
    void ComputeE(const TVector3d& point,
                  std::complex<double>& Ex,
                  std::complex<double>& Ey,
                  std::complex<double>& Ez) const;

    /**
     * @brief Compute vector potential A at a point
     * A = mu_0/(4*pi) * integral{ K / |r - r'| } dA'
     * @param point Evaluation point
     * @param Ax, Ay, Az Output vector potential components [T*m]
     */
    void ComputeA(const TVector3d& point,
                  std::complex<double>& Ax,
                  std::complex<double>& Ay,
                  std::complex<double>& Az) const;

    /**
     * @brief Compute scalar potential Phi at a point
     * Phi = 1/(4*pi*eps_0) * integral{ sigma / |r - r'| } dA'
     * @param point Evaluation point
     * @return Scalar electric potential [V]
     */
    std::complex<double> ComputePhi(const TVector3d& point) const;

    // ========== Port definition (for impedance calculation) ==========

    /**
     * @brief Define input/output port
     * @param terminal1 First terminal panel indices
     * @param terminal2 Second terminal panel indices
     */
    void DefinePort(const std::vector<int>& terminal1,
                    const std::vector<int>& terminal2);

    /**
     * @brief Get port impedance after solve
     */
    std::complex<double> GetPortImpedance() const { return portImpedance_; }

    /**
     * @brief Set port impedance (used by solver)
     */
    void SetPortImpedance(std::complex<double> Z) { portImpedance_ = Z; }

    // ========== Excitation settings ==========

    /**
     * @brief Set voltage excitation at port
     * @param V_real Real part of voltage [V]
     * @param V_imag Imaginary part of voltage [V]
     */
    void SetVoltageExcitation(double V_real, double V_imag);

    /**
     * @brief Set current excitation at port
     * @param I_real Real part of current [A]
     * @param I_imag Imaginary part of current [A]
     */
    void SetCurrentExcitation(double I_real, double I_imag);

    /**
     * @brief Get excitation type
     * @return 0=none, 1=voltage, 2=current
     */
    int GetExcitationType() const { return excitationType_; }

    /**
     * @brief Get excitation value (voltage or current)
     */
    std::complex<double> GetExcitationValue() const { return excitationValue_; }

    /**
     * @brief Get total current through conductor (after solve)
     */
    std::complex<double> GetTotalCurrent() const { return totalCurrent_; }

    /**
     * @brief Set total current (used by solver)
     */
    void SetTotalCurrent(std::complex<double> I) { totalCurrent_ = I; }

private:
    // Material properties
    double conductivity_;     // [S/m]
    double frequency_;        // [Hz]

    // Formulation
    ConductorFormulation formulation_;

    // Surface discretization
    std::vector<SurfacePanel> panels_;

    // Solution vectors
    std::vector<std::complex<double>> surfaceCurrent_;   // K [A/m], 3 components per panel
    std::vector<std::complex<double>> surfaceCharge_;    // sigma [C/m^2]

    // Port definition
    std::vector<int> portTerminal1_;
    std::vector<int> portTerminal2_;
    std::complex<double> portImpedance_;

    // Excitation settings
    int excitationType_;                    // 0=none, 1=voltage, 2=current
    std::complex<double> excitationValue_;  // V or I depending on type
    std::complex<double> totalCurrent_;     // Total current after solve

    // Wire centerline path (for Biot-Savart field computation)
    // Stored when creating wire/loop conductors
    std::vector<TVector3d> wireCenterPath_;

    // Helper methods
    void GenerateRectangularFacePanels(const TVector3d& corner,
                                       const TVector3d& edge1,
                                       const TVector3d& edge2,
                                       int n1, int n2);

    void ComputePanelNormals();

    // Green's function based on formulation
    std::complex<double> GreenFunction(double r) const;
    std::complex<double> GreenFunctionGrad(double r, double dr_dx) const;
};

/**
 * @brief Frequency solver for conductor analysis
 */
class radTConductorSolver {
public:
    radTConductorSolver();
    ~radTConductorSolver();

    /**
     * @brief Add conductor to solver
     */
    void AddConductor(std::shared_ptr<radTConductor> conductor);

    /**
     * @brief Clear all conductors
     */
    void Clear();

    /**
     * @brief Set analysis frequency
     */
    void SetFrequency(double frequency);

    /**
     * @brief Set formulation for all conductors
     */
    void SetFormulation(ConductorFormulation form);

    /**
     * @brief Solve the system at current frequency
     */
    void Solve();

    /**
     * @brief Solve and return port impedance
     */
    std::complex<double> SolveImpedance(double frequency);

    /**
     * @brief Frequency sweep
     */
    std::vector<std::complex<double>> ImpedanceSweep(
        const std::vector<double>& frequencies);

    /**
     * @brief Use pFFT acceleration
     */
    void EnablePfft(bool enable = true) { usePfft_ = enable; }

    /**
     * @brief Compute B field at a point (after solve)
     */
    void ComputeB(const TVector3d& point,
                  std::complex<double>& Bx,
                  std::complex<double>& By,
                  std::complex<double>& Bz) const;

    /**
     * @brief Get port impedance after solve
     */
    std::complex<double> GetPortImpedance() const { return portImpedance_; }

    /**
     * @brief Get total current after solve
     */
    std::complex<double> GetTotalCurrent() const { return totalCurrent_; }

private:
    std::vector<std::shared_ptr<radTConductor>> conductors_;
    double frequency_;
    ConductorFormulation formulation_;
    bool usePfft_;
    bool directImpedanceMode_;  // True when impedance computed directly

    // pFFT accelerator
    std::unique_ptr<radTPfft> pfft_;

    // RWG EFIE solver
    std::shared_ptr<RWGMesh> rwgMesh_;
    std::unique_ptr<RWGEFIESolver> rwgSolver_;
    std::complex<double> portImpedance_;
    std::complex<double> totalCurrent_;

    // System matrix and RHS (for backward compatibility)
    std::vector<std::complex<double>> systemMatrix_;
    std::vector<std::complex<double>> rhs_;
    std::vector<std::complex<double>> solution_;

    void BuildSystemMatrix();
    void SolveLinearSystem();
    void ExtractSolution();
    void ComputePortImpedance();
    void BuildRWGMesh();
};

} // namespace radia

#endif // RAD_CONDUCTOR_H
