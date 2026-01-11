/*
 * rad_peec_mmm_api.cpp
 *
 * Python API implementation for Coupled solver (PEEC + MMM).
 *
 * This provides the interface between Radia Python bindings and
 * the PEECMMMCoupledSolver class (accessed via CplMag* API).
 *
 * Part of Radia project
 */

#include "rad_peec_mmm_coupled.h"
#include "rad_io_buffer.h"
#include "rad_application.h"
#include "rad_interaction.h"
#include <map>
#include <memory>
#include <cstring>

using namespace radia;

// External I/O buffer (defined in radentry.cpp)
extern radTIOBuffer ioBuffer;

// External Radia application instance (defined in radentry.cpp)
extern radTApplication rad;

// External functions from rad_conductor_api.cpp
extern bool IsConductorHandle(int handle);
extern std::shared_ptr<radTConductor> GetConductorByHandle(int handle);

//=========================================================================
// Helper: Create radTInteraction from magnet handle
//=========================================================================

/**
 * @brief Create radTInteraction from a magnet object handle
 *
 * This creates the interaction matrix needed for MMM solver from a
 * Radia magnet object (ObjHexahedron, ObjRecMag, etc.)
 *
 * @param magnetHandle Handle to magnet object
 * @return Pointer to radTInteraction or nullptr on error
 */
static radTInteraction* CreateInteractionFromMagnetHandle(int magnetHandle)
{
    // Get object from handle using PreRelax approach
    // PreRelax creates the interaction and returns its key
    int interactKey = rad.PreRelax(magnetHandle, 0, false);
    if (interactKey <= 0) {
        return nullptr;  // Failed to create interaction
    }

    // Retrieve the interaction using public API method
    radTInteraction* interactPtr = rad.GetInteractionByKey(interactKey);
    if (interactPtr == nullptr) {
        return nullptr;
    }

    // Note: The interaction is owned by the global handle table
    // We return the pointer but don't own it
    return interactPtr;
}

//=========================================================================
// Coupled Solver Handle Management
//=========================================================================

// Global solver storage
static std::map<int, std::shared_ptr<PEECMMMCoupledSolver>> gCplMagSolverMap;
static int gNextCoupledHandle = 20000;  // Start from 20000 to avoid conflicts

// Store a solver and return its handle
static int StoreCplMagSolver(std::shared_ptr<PEECMMMCoupledSolver> solver)
{
    int handle = gNextCoupledHandle++;
    gCplMagSolverMap[handle] = solver;
    return handle;
}

// Get a solver by handle
static std::shared_ptr<PEECMMMCoupledSolver> GetCplMagSolver(int handle)
{
    auto it = gCplMagSolverMap.find(handle);
    if (it == gCplMagSolverMap.end()) {
        return nullptr;
    }
    return it->second;
}

//=========================================================================
// External API Functions (called from radpy_pyapi.cpp)
//=========================================================================

/**
 * @brief Create coupled solver
 *
 * Python: solver = rad.CoupledCreate(conductor_handle, magnet_handle)
 *
 * @param pOut Output array (solver handle)
 * @param pNout Number of outputs (1)
 * @param conductor Conductor handle (from CndRecBlock, CndLoop, etc.)
 * @param magnet Magnet handle (from ObjRecMag, ObjHexahedron, etc.)
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagCreate(double* pOut, int* pNout, int conductor, int magnet)
{
    try {
        auto solver = std::make_shared<PEECMMMCoupledSolver>();

        // Set conductor from handle
        auto cond = GetConductorByHandle(conductor);
        if (!cond) {
            ioBuffer.StoreErrorMessage("CplMagCreate: Invalid conductor handle");
            return -1;
        }
        solver->SetConductor(cond);

        // Set magnet handle (resolved during solve)
        solver->SetMagnet(magnet);

        // Store and return handle
        int handle = StoreCplMagSolver(solver);
        pOut[0] = static_cast<double>(handle);
        *pNout = 1;

        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Set coupled solver frequency
 *
 * Python: rad.CplMagSetFrequency(solver, frequency)
 *
 * @param solver Solver handle
 * @param frequency Analysis frequency [Hz]
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagSetFrequency(int solver, double frequency)
{
    try {
        auto pSolver = GetCplMagSolver(solver);
        if (!pSolver) {
            ioBuffer.StoreErrorMessage("CplMagSetFrequency: Invalid solver handle");
            return -1;
        }

        pSolver->SetFrequency(frequency);
        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Set coupled solver voltage excitation
 *
 * Python: rad.CplMagSetVoltage(solver, V_real, V_imag)
 *
 * @param solver Solver handle
 * @param V_real Real part of voltage [V]
 * @param V_imag Imaginary part of voltage [V]
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagSetVoltage(int solver, double V_real, double V_imag)
{
    try {
        auto pSolver = GetCplMagSolver(solver);
        if (!pSolver) {
            ioBuffer.StoreErrorMessage("CplMagSetVoltage: Invalid solver handle");
            return -1;
        }

        pSolver->SetVoltageExcitation(std::complex<double>(V_real, V_imag));
        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Set coupled solver current excitation
 *
 * Python: rad.CplMagSetCurrent(solver, I_real, I_imag)
 *
 * @param solver Solver handle
 * @param I_real Real part of current [A]
 * @param I_imag Imaginary part of current [A]
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagSetCurrent(int solver, double I_real, double I_imag)
{
    try {
        auto pSolver = GetCplMagSolver(solver);
        if (!pSolver) {
            ioBuffer.StoreErrorMessage("CplMagSetCurrent: Invalid solver handle");
            return -1;
        }

        pSolver->SetCurrentExcitation(std::complex<double>(I_real, I_imag));
        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Set external field for coupled solver
 *
 * Python: rad.CplMagSetExtField(solver, Hx, Hy, Hz)
 *
 * @param solver Solver handle
 * @param Hx, Hy, Hz External H field [A/m]
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagSetExtField(int solver, double Hx, double Hy, double Hz)
{
    try {
        auto pSolver = GetCplMagSolver(solver);
        if (!pSolver) {
            ioBuffer.StoreErrorMessage("CplMagSetExtField: Invalid solver handle");
            return -1;
        }

        pSolver->SetExternalField(Hx, Hy, Hz);
        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Solve coupled system
 *
 * Python: result = rad.CplMagSolve(solver)
 *
 * @param pOut Output array [Z_real, Z_imag, P_cond, P_mag, n_iter]
 * @param pNout Number of outputs (5)
 * @param solver Solver handle
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagSolve(double* pOut, int* pNout, int solver)
{
    try {
        auto pSolver = GetCplMagSolver(solver);
        if (!pSolver) {
            ioBuffer.StoreErrorMessage("CplMagSolve: Invalid solver handle");
            return -1;
        }

        // Create radTInteraction from magnet handle if not already set
        int magnetHandle = pSolver->GetMagnetHandle();
        if (magnetHandle > 0 && !pSolver->HasMagnetInteraction()) {
            radTInteraction* interaction = CreateInteractionFromMagnetHandle(magnetHandle);
            if (interaction) {
                pSolver->SetMagnetInteraction(interaction);
                // Interaction is owned by Radia global table, not the solver
                pSolver->SetOwnsMagnetInteraction(false);
            } else {
                ioBuffer.StoreErrorMessage("CplMagSolve: Failed to create interaction from magnet handle");
                return -1;
            }
        }

        PEECMMMResult result = pSolver->Solve();

        pOut[0] = result.impedance.real();
        pOut[1] = result.impedance.imag();
        pOut[2] = result.conductorPower;
        pOut[3] = result.magnetPower;
        pOut[4] = static_cast<double>(result.iterations);
        *pNout = 5;

        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Get impedance from coupled solver
 *
 * Python: Z = rad.CplMagImpedance(solver)
 *
 * @param pOut Output array [Z_real, Z_imag]
 * @param pNout Number of outputs (2)
 * @param solver Solver handle
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagImpedance(double* pOut, int* pNout, int solver)
{
    try {
        auto pSolver = GetCplMagSolver(solver);
        if (!pSolver) {
            ioBuffer.StoreErrorMessage("CplMagImpedance: Invalid solver handle");
            return -1;
        }

        std::complex<double> Z = pSolver->GetImpedance();
        pOut[0] = Z.real();
        pOut[1] = Z.imag();
        *pNout = 2;

        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Compute B field from coupled system
 *
 * Python: B = rad.CplMagFld(solver, [x, y, z])
 *
 * @param pOut Output array [Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]
 * @param pNout Number of outputs (6)
 * @param solver Solver handle
 * @param pCoord Evaluation point [x, y, z]
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagFld(double* pOut, int* pNout, int solver, double* pCoord)
{
    try {
        auto pSolver = GetCplMagSolver(solver);
        if (!pSolver) {
            ioBuffer.StoreErrorMessage("CplMagFld: Invalid solver handle");
            return -1;
        }

        TVector3d point(pCoord[0], pCoord[1], pCoord[2]);
        std::complex<double> Bx, By, Bz;
        pSolver->ComputeB(point, Bx, By, Bz);

        pOut[0] = Bx.real();
        pOut[1] = By.real();
        pOut[2] = Bz.real();
        pOut[3] = Bx.imag();
        pOut[4] = By.imag();
        pOut[5] = Bz.imag();
        *pNout = 6;

        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Delete coupled solver
 *
 * Python: rad.CplMagDelete(solver)
 *
 * @param solver Solver handle
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagDelete(int solver)
{
    try {
        auto it = gCplMagSolverMap.find(solver);
        if (it != gCplMagSolverMap.end()) {
            gCplMagSolverMap.erase(it);
        }
        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Set complex permeability for coupled solver
 *
 * Python: rad.CplMagSetMu(solver, mu_r_real, mu_r_imag)
 *
 * Complex permeability: mu = mu' - j*mu"
 * - mu_r_real: Real part of relative permeability (mu'_r)
 * - mu_r_imag: Imaginary part (mu"_r, positive for loss)
 *
 * @param solver Solver handle
 * @param mu_r_real Real part of relative permeability
 * @param mu_r_imag Imaginary part of relative permeability
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagSetMu(int solver, double mu_r_real, double mu_r_imag)
{
    try {
        auto pSolver = GetCplMagSolver(solver);
        if (!pSolver) {
            ioBuffer.StoreErrorMessage("CplMagSetMu: Invalid solver handle");
            return -1;
        }

        pSolver->SetComplexPermeability(mu_r_real, mu_r_imag);
        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Set conductor for coupled solver
 *
 * Python: rad.CplMagSetConductor(solver, conductor_handle)
 *
 * @param solver Solver handle
 * @param conductor Conductor handle (from CndLoop, CndRecBlock, etc.)
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagSetConductor(int solver, int conductor)
{
    try {
        auto pSolver = GetCplMagSolver(solver);
        if (!pSolver) {
            ioBuffer.StoreErrorMessage("CplMagSetConductor: Invalid solver handle");
            return -1;
        }

        // Get conductor by handle
        auto cond = GetConductorByHandle(conductor);
        if (!cond) {
            ioBuffer.StoreErrorMessage("CplMagSetConductor: Invalid conductor handle");
            return -1;
        }

        pSolver->SetConductor(cond);
        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Get power losses from coupled solver
 *
 * Python: [P_cond, P_mag] = rad.CplMagPower(solver)
 *
 * @param pOut Output array [P_conductor, P_magnet]
 * @param pNout Number of outputs (2)
 * @param solver Solver handle
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagPower(double* pOut, int* pNout, int solver)
{
    try {
        auto pSolver = GetCplMagSolver(solver);
        if (!pSolver) {
            ioBuffer.StoreErrorMessage("CplMagPower: Invalid solver handle");
            return -1;
        }

        pOut[0] = pSolver->GetConductorPower();
        pOut[1] = pSolver->GetMagnetPower();
        *pNout = 2;

        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}

/**
 * @brief Frequency sweep for coupled solver
 *
 * Python: results = rad.CplMagSweep(solver, [f1, f2, ...])
 *
 * @param pOut Output array [Z_re1, Z_im1, Z_re2, Z_im2, ...]
 * @param pNout Number of outputs (2 * n_freq)
 * @param solver Solver handle
 * @param pFreq Array of frequencies [Hz]
 * @param nFreq Number of frequencies
 * @return 0 on success, error code otherwise
 */
extern "C" int RadCplMagSweep(double* pOut, int* pNout, int solver, double* pFreq, int nFreq)
{
    try {
        auto pSolver = GetCplMagSolver(solver);
        if (!pSolver) {
            ioBuffer.StoreErrorMessage("CplMagSweep: Invalid solver handle");
            return -1;
        }

        std::vector<double> frequencies(pFreq, pFreq + nFreq);
        std::vector<PEECMMMResult> results = pSolver->FrequencySweep(frequencies);

        int idx = 0;
        for (const auto& r : results) {
            pOut[idx++] = r.impedance.real();
            pOut[idx++] = r.impedance.imag();
        }
        *pNout = idx;

        return 0;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        return -1;
    }
}
