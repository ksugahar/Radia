/*
 * rad_conductor_api.cpp
 *
 * Internal API implementation for conductor analysis functions.
 * These functions are called by radentry.cpp and use rad_conductor.h classes.
 *
 * Part of Radia project - FastImp integration
 */

#include "rad_conductor.h"
#include "rad_io_buffer.h"
#include <map>
#include <memory>
#include <cstring>

using namespace radia;

// External I/O buffer (defined in radentry.cpp)
extern radTIOBuffer ioBuffer;

//=========================================================================
// Conductor Handle Management
//=========================================================================

// Global conductor storage (similar to Radia's magnetic object storage)
static std::map<int, std::shared_ptr<radTConductor>> gConductorMap;
static std::map<int, std::shared_ptr<radTConductorSolver>> gConductorSolverMap;
static int gNextConductorHandle = 10000;  // Start from 10000 to avoid conflicts with mag objects

// Store a conductor and return its handle
static int StoreConductor(std::shared_ptr<radTConductor> cond)
{
    int handle = gNextConductorHandle++;
    gConductorMap[handle] = cond;
    return handle;
}

// Get a conductor by handle
static std::shared_ptr<radTConductor> GetConductor(int handle)
{
    auto it = gConductorMap.find(handle);
    if (it == gConductorMap.end()) {
        return nullptr;
    }
    return it->second;
}

// Container storage for conductor groups
static std::map<int, std::vector<int>> gConductorContainerMap;

//=========================================================================
// Conductor Creation Functions
//=========================================================================

void CndRecBlock(double* P, double* L, double sigma, int num_panels)
{
    try {
        auto cond = std::make_shared<radTConductor>();

        TVector3d center(P[0], P[1], P[2]);
        TVector3d dims(L[0], L[1], L[2]);

        cond->SetConductivity(sigma);
        cond->CreateFromRecBlock(center, dims, num_panels > 0 ? num_panels : 4);

        int handle = StoreConductor(cond);
        ioBuffer.StoreInt(handle);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        ioBuffer.StoreInt(0);
    }
    catch (...) {
        ioBuffer.StoreErrorMessage("Unknown error in CndRecBlock");
        ioBuffer.StoreInt(0);
    }
}

//-------------------------------------------------------------------------

void CndHexahedron(double* FlatVert, double sigma, int num_panels)
{
    try {
        auto cond = std::make_shared<radTConductor>();

        std::vector<TVector3d> vertices(8);
        for (int i = 0; i < 8; i++) {
            vertices[i] = TVector3d(FlatVert[3*i], FlatVert[3*i+1], FlatVert[3*i+2]);
        }

        cond->SetConductivity(sigma);
        cond->CreateFromHexahedron(vertices, num_panels > 0 ? num_panels : 4);

        int handle = StoreConductor(cond);
        ioBuffer.StoreInt(handle);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        ioBuffer.StoreInt(0);
    }
    catch (...) {
        ioBuffer.StoreErrorMessage("Unknown error in CndHexahedron");
        ioBuffer.StoreInt(0);
    }
}

//-------------------------------------------------------------------------

void CndWire(double* FlatPath, int np, char cross_section, double width, double height,
             double sigma, int num_panels_around)
{
    try {
        auto cond = std::make_shared<radTConductor>();

        std::vector<TVector3d> path(np);
        for (int i = 0; i < np; i++) {
            path[i] = TVector3d(FlatPath[3*i], FlatPath[3*i+1], FlatPath[3*i+2]);
        }

        std::string cs_type = (cross_section == 'c' || cross_section == 'C') ? "circular" : "rectangular";

        cond->SetConductivity(sigma);
        cond->CreateWire(path, cs_type, width, height, num_panels_around > 0 ? num_panels_around : 8, 0);

        int handle = StoreConductor(cond);
        ioBuffer.StoreInt(handle);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        ioBuffer.StoreInt(0);
    }
    catch (...) {
        ioBuffer.StoreErrorMessage("Unknown error in CndWire");
        ioBuffer.StoreInt(0);
    }
}

//-------------------------------------------------------------------------

void CndLoop(double* center, double radius, double* normal, char cross_section,
             double wire_width, double wire_height, double sigma,
             int num_panels_around, int num_panels_loop)
{
    try {
        auto cond = std::make_shared<radTConductor>();

        TVector3d ctr(center[0], center[1], center[2]);
        TVector3d nrm(normal[0], normal[1], normal[2]);
        std::string cs_type = (cross_section == 'c' || cross_section == 'C') ? "circular" : "rectangular";

        cond->SetConductivity(sigma);
        cond->CreateLoop(ctr, radius, nrm, cs_type, wire_width, wire_height,
                         num_panels_around > 0 ? num_panels_around : 8,
                         num_panels_loop > 0 ? num_panels_loop : 30);

        int handle = StoreConductor(cond);
        ioBuffer.StoreInt(handle);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        ioBuffer.StoreInt(0);
    }
    catch (...) {
        ioBuffer.StoreErrorMessage("Unknown error in CndLoop");
        ioBuffer.StoreInt(0);
    }
}

//-------------------------------------------------------------------------

void CndSpiral(double* center, double inner_radius, double outer_radius, double pitch,
               int num_turns, double* axis, char cross_section, double wire_width,
               double wire_height, double sigma, int num_panels_around)
{
    try {
        auto cond = std::make_shared<radTConductor>();

        TVector3d ctr(center[0], center[1], center[2]);
        TVector3d ax(axis[0], axis[1], axis[2]);
        std::string cs_type = (cross_section == 'c' || cross_section == 'C') ? "circular" : "rectangular";

        cond->SetConductivity(sigma);
        cond->CreateSpiral(ctr, inner_radius, outer_radius, pitch, num_turns, ax,
                           cs_type, wire_width, wire_height,
                           num_panels_around > 0 ? num_panels_around : 8);

        int handle = StoreConductor(cond);
        ioBuffer.StoreInt(handle);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        ioBuffer.StoreInt(0);
    }
    catch (...) {
        ioBuffer.StoreErrorMessage("Unknown error in CndSpiral");
        ioBuffer.StoreInt(0);
    }
}

//=========================================================================
// Container Functions
//=========================================================================

void CndCnt(int* Conds, int ncond)
{
    try {
        int container_handle = gNextConductorHandle++;
        std::vector<int> members;

        for (int i = 0; i < ncond; i++) {
            if (GetConductor(Conds[i])) {
                members.push_back(Conds[i]);
            }
        }

        gConductorContainerMap[container_handle] = members;
        ioBuffer.StoreInt(container_handle);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        ioBuffer.StoreInt(0);
    }
}

//-------------------------------------------------------------------------

void CndCntAdd(int cnt, int* Conds, int ncond)
{
    try {
        auto it = gConductorContainerMap.find(cnt);
        if (it == gConductorContainerMap.end()) {
            ioBuffer.StoreErrorMessage("Invalid conductor container handle");
            return;
        }

        for (int i = 0; i < ncond; i++) {
            if (GetConductor(Conds[i])) {
                it->second.push_back(Conds[i]);
            }
        }
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//=========================================================================
// Configuration Functions
//=========================================================================

void CndSetFormulation(int cond, const char* formulation)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            return;
        }

        ConductorFormulation form = ConductorFormulation::MQS;
        if (strcmp(formulation, "dc") == 0 || strcmp(formulation, "DC") == 0) {
            form = ConductorFormulation::DC;
        } else if (strcmp(formulation, "mqs") == 0 || strcmp(formulation, "MQS") == 0) {
            form = ConductorFormulation::MQS;
        } else if (strcmp(formulation, "emqs") == 0 || strcmp(formulation, "EMQS") == 0) {
            form = ConductorFormulation::EMQS;
        } else if (strcmp(formulation, "fullwave") == 0 || strcmp(formulation, "FULLWAVE") == 0) {
            form = ConductorFormulation::FullWave;
        }

        conductor->SetFormulation(form);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//-------------------------------------------------------------------------

void CndSetFrequency(int cond, double frequency)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            return;
        }

        conductor->SetFrequency(frequency);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//-------------------------------------------------------------------------

void CndSetPfft(int cond, int enable)
{
    try {
        // pFFT setting is typically on the solver, not individual conductor
        // For now, this is a no-op placeholder
        (void)cond;
        (void)enable;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//=========================================================================
// Excitation Functions
//=========================================================================

void CndSetVoltage(int cond, double V_real, double V_imag)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            return;
        }

        conductor->SetVoltageExcitation(V_real, V_imag);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//-------------------------------------------------------------------------

void CndSetCurrent(int cond, double I_real, double I_imag)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            return;
        }

        conductor->SetCurrentExcitation(I_real, I_imag);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//-------------------------------------------------------------------------

void CndGetTotalCurrent(double* I_real, double* I_imag, int cond)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            *I_real = 0;
            *I_imag = 0;
            return;
        }

        std::complex<double> I = conductor->GetTotalCurrent();
        *I_real = I.real();
        *I_imag = I.imag();
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        *I_real = 0;
        *I_imag = 0;
    }
}

//=========================================================================
// Port Definition Functions
//=========================================================================

void CndDefinePort(int cond, int* terminal1, int n1, int* terminal2, int n2)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            return;
        }

        std::vector<int> t1, t2;
        if (terminal1 && n1 > 0) {
            t1.assign(terminal1, terminal1 + n1);
        }
        if (terminal2 && n2 > 0) {
            t2.assign(terminal2, terminal2 + n2);
        }

        conductor->DefinePort(t1, t2);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//-------------------------------------------------------------------------

void CndDefinePortAuto(int cond)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            return;
        }

        // Auto-detect: use first and last panels
        int np = conductor->NumPanels();
        if (np >= 2) {
            std::vector<int> t1 = {0};
            std::vector<int> t2 = {np - 1};
            conductor->DefinePort(t1, t2);
        }
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//=========================================================================
// Solver Functions
//=========================================================================

void CndSolve(int cond)
{
    try {
        // Check if it's a container
        auto cnt_it = gConductorContainerMap.find(cond);
        if (cnt_it != gConductorContainerMap.end()) {
            // Solve container: create solver with all conductors
            auto solver = std::make_shared<radTConductorSolver>();
            for (int handle : cnt_it->second) {
                auto c = GetConductor(handle);
                if (c) {
                    solver->AddConductor(c);
                }
            }
            solver->Solve();
            gConductorSolverMap[cond] = solver;
        } else {
            // Single conductor
            auto conductor = GetConductor(cond);
            if (!conductor) {
                ioBuffer.StoreErrorMessage("Invalid conductor handle");
                return;
            }

            auto solver = std::make_shared<radTConductorSolver>();
            solver->AddConductor(conductor);
            solver->SetFrequency(conductor->GetFrequency());
            solver->SetFormulation(conductor->GetFormulation());
            solver->Solve();
            gConductorSolverMap[cond] = solver;
        }
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//-------------------------------------------------------------------------

void CndGetImpedance(double* Z_real, double* Z_imag, int cond)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            *Z_real = 0;
            *Z_imag = 0;
            return;
        }

        std::complex<double> Z = conductor->GetPortImpedance();
        *Z_real = Z.real();
        *Z_imag = Z.imag();
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        *Z_real = 0;
        *Z_imag = 0;
    }
}

//-------------------------------------------------------------------------

void CndImpedanceSweep(double* Z_real, double* Z_imag, int cond, double* freqs, int nfreq)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            return;
        }

        auto solver = std::make_shared<radTConductorSolver>();
        solver->AddConductor(conductor);

        std::vector<double> freq_vec(freqs, freqs + nfreq);
        std::vector<std::complex<double>> Z_vec = solver->ImpedanceSweep(freq_vec);

        for (int i = 0; i < nfreq; i++) {
            Z_real[i] = Z_vec[i].real();
            Z_imag[i] = Z_vec[i].imag();
        }
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//=========================================================================
// Field Computation Functions
//=========================================================================

void CndFldB(double* B_real, double* B_imag, int cond, double* point)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            B_real[0] = B_real[1] = B_real[2] = 0;
            B_imag[0] = B_imag[1] = B_imag[2] = 0;
            return;
        }

        TVector3d pt(point[0], point[1], point[2]);
        std::complex<double> Bx, By, Bz;
        conductor->ComputeB(pt, Bx, By, Bz);

        B_real[0] = Bx.real();
        B_real[1] = By.real();
        B_real[2] = Bz.real();
        B_imag[0] = Bx.imag();
        B_imag[1] = By.imag();
        B_imag[2] = Bz.imag();
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        B_real[0] = B_real[1] = B_real[2] = 0;
        B_imag[0] = B_imag[1] = B_imag[2] = 0;
    }
}

//-------------------------------------------------------------------------

void CndFldE(double* E_real, double* E_imag, int cond, double* point)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            E_real[0] = E_real[1] = E_real[2] = 0;
            E_imag[0] = E_imag[1] = E_imag[2] = 0;
            return;
        }

        TVector3d pt(point[0], point[1], point[2]);
        std::complex<double> Ex, Ey, Ez;
        conductor->ComputeE(pt, Ex, Ey, Ez);

        E_real[0] = Ex.real();
        E_real[1] = Ey.real();
        E_real[2] = Ez.real();
        E_imag[0] = Ex.imag();
        E_imag[1] = Ey.imag();
        E_imag[2] = Ez.imag();
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        E_real[0] = E_real[1] = E_real[2] = 0;
        E_imag[0] = E_imag[1] = E_imag[2] = 0;
    }
}

//-------------------------------------------------------------------------

void CndFldBBatch(double* B_real, double* B_imag, int cond, double* points, int npts)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            return;
        }

        for (int i = 0; i < npts; i++) {
            TVector3d pt(points[3*i], points[3*i+1], points[3*i+2]);
            std::complex<double> Bx, By, Bz;
            conductor->ComputeB(pt, Bx, By, Bz);

            B_real[3*i]   = Bx.real();
            B_real[3*i+1] = By.real();
            B_real[3*i+2] = Bz.real();
            B_imag[3*i]   = Bx.imag();
            B_imag[3*i+1] = By.imag();
            B_imag[3*i+2] = Bz.imag();
        }
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//=========================================================================
// Solution Access Functions
//=========================================================================

void CndGetSurfaceCurrent(double* K_real, double* K_imag, int* npanels, int cond)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            *npanels = 0;
            return;
        }

        const auto& K = conductor->SurfaceCurrent();
        int np = conductor->NumPanels();
        *npanels = np;

        // Surface current is stored as 3 complex components per panel
        for (int i = 0; i < np * 3; i++) {
            K_real[i] = K[i].real();
            K_imag[i] = K[i].imag();
        }
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        *npanels = 0;
    }
}

//-------------------------------------------------------------------------

void CndGetSurfaceCharge(double* sigma_real, double* sigma_imag, int* npanels, int cond)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            *npanels = 0;
            return;
        }

        const auto& sigma = conductor->SurfaceCharge();
        int np = conductor->NumPanels();
        *npanels = np;

        for (int i = 0; i < np; i++) {
            sigma_real[i] = sigma[i].real();
            sigma_imag[i] = sigma[i].imag();
        }
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        *npanels = 0;
    }
}

//-------------------------------------------------------------------------

void CndNumPanels(int* npanels, int cond)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            *npanels = 0;
            return;
        }

        *npanels = conductor->NumPanels();
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        *npanels = 0;
    }
}

//-------------------------------------------------------------------------

void CndGetPanelInfo(double* centers, double* areas, int cond)
{
    try {
        auto conductor = GetConductor(cond);
        if (!conductor) {
            ioBuffer.StoreErrorMessage("Invalid conductor handle");
            return;
        }

        std::vector<TVector3d> panel_centers = conductor->GetPanelCenters();
        std::vector<double> panel_areas = conductor->GetPanelAreas();

        int np = conductor->NumPanels();
        for (int i = 0; i < np; i++) {
            centers[3*i]   = panel_centers[i].x;
            centers[3*i+1] = panel_centers[i].y;
            centers[3*i+2] = panel_centers[i].z;
            areas[i] = panel_areas[i];
        }
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//=========================================================================
// Coupled Solver Functions
//=========================================================================

void CoupledSolve(int cond_cnt, int mag_cnt, double precision, int max_iter)
{
    try {
        // TODO: Implement coupled conductor-magnetic solver
        // This requires integration with existing Radia magnetic solver
        ioBuffer.StoreErrorMessage("CoupledSolve not yet implemented");
        (void)cond_cnt;
        (void)mag_cnt;
        (void)precision;
        (void)max_iter;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//=========================================================================
// SIBC Material Functions
//=========================================================================

// SIBC material storage
static std::map<int, std::pair<double, double>> gSIBCMaterialMap;  // handle -> (sigma, mu_r)
static int gNextSIBCHandle = 20000;

void MatSIBC(double sigma, double mu_r)
{
    try {
        int handle = gNextSIBCHandle++;
        gSIBCMaterialMap[handle] = std::make_pair(sigma, mu_r);
        ioBuffer.StoreInt(handle);
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
        ioBuffer.StoreInt(0);
    }
}

//-------------------------------------------------------------------------

void SIBCSetType(int mat, const char* sibc_type)
{
    try {
        auto it = gSIBCMaterialMap.find(mat);
        if (it == gSIBCMaterialMap.end()) {
            ioBuffer.StoreErrorMessage("Invalid SIBC material handle");
            return;
        }

        // TODO: Store SIBC type for the material
        (void)sibc_type;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}

//-------------------------------------------------------------------------

void SIBCSetCrossSection(int mat, const char* shape, double* params, int nparams)
{
    try {
        auto it = gSIBCMaterialMap.find(mat);
        if (it == gSIBCMaterialMap.end()) {
            ioBuffer.StoreErrorMessage("Invalid SIBC material handle");
            return;
        }

        // TODO: Store cross-section info for nonlocal SIBC
        (void)shape;
        (void)params;
        (void)nparams;
    }
    catch (const std::exception& e) {
        ioBuffer.StoreErrorMessage(e.what());
    }
}
