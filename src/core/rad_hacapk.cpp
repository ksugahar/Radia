/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk.cpp
*
* Project:        RADIA
*
* Description:    HACApK (H-matrix with ACA+) interface for BiCGSTAB solver
*                 Implementation of RadHACApKBase / RadHACApKMSCManager and callback functions
*
* First release:  2025
*
* Reference:      ppOpen-HPC project (MIT License)
*                 https://github.com/ppohHPC/ppOpen-HPC
*
-------------------------------------------------------------------------*/

#include "rad_hacapk.h"
#include "rad_interaction.h"
#include "rad_constants.h"
#include "rad_poly_analytical.h"
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <algorithm>
#include <cstring>
#include <mkl.h>   // cblas_dgemm + LAPACKE_dsyevd for the antisym-IMA plane-slab null space (SVD-free)
#include <cstdio>
#include <cstdlib>  // std::getenv / std::atoi (RADIA_LS_BLOCK_CELLS research toggle)
#include <iostream>
#include <chrono>
#include <atomic>
#include "rad_parallel.h"

// Include C++ compatible HACApK wrapper header
extern "C" {
#include "../ext/HACApK/cHACApK_cpp.h"
// H-LU preconditioner wrappers (defined in cHACApK_harith.c; declared in
// cHACApK_harith.h, which is not included here): factor A_SS once, apply many.
void* cHACApK_hlu_factor_leafmtxp(void* leafmtxp_void, void* control_void, int nffc);
int   cHACApK_hlu_apply(void* root_void, void* control_void,
                        const double* r, double* z, int nd);
void  cHACApK_hlu_free_factors(void* root_void);
void  cHACApK_hlu_set_trunc_tol(double tol);
}

//=========================================================================
// Constants for 6DOF MSC field computation
// Now using unified constants from rad_constants.h
//=========================================================================

//=========================================================================
// Global callback state (required by HACApK C interface)
//=========================================================================

namespace {
    // Global callback state shared across all OpenMP threads.
    // NOT thread_local: HACApK calls cHACApK_entry_ij from OpenMP worker threads,
    // which must see the same manager/invChi/interaction set by the main thread.
    // Thread safety for concurrent BuildHMatrix calls (multiple Python threads)
    // is ensured by Python's GIL; standalone C++ use would require a mutex.
    RadHACApKBase* g_currentManager = nullptr;
    std::vector<double> g_invChi;
    radTInteraction* g_interaction = nullptr;
    int g_nElem = 0;
    int g_nffc = 3;  // DOF per element (default 3 for standard elements)
    // Note: lod (DOF permutation array) is now accessed via C-side accessors:
    //   HACApK_get_current_lod() and HACApK_get_current_lod_size()
    // Set during H-matrix build by HACApK_build_hmatrix_varDOF_wrapper

    // FIX (2025-02-04): Generation counter for thread-local cache invalidation
    // Incremented each time a new HACApK manager is created or BuildHMatrix is called.
    // Thread-local caches check this counter to invalidate stale entries.
    // Using a counter instead of pointer comparison prevents issues with memory reuse.
    std::atomic<uint64_t> g_hacapk_generation{0};
}

namespace RadHACApKCallback {

void SetCurrentManager(RadHACApKBase* manager) {
    g_currentManager = manager;
}

RadHACApKBase* GetCurrentManager() {
    return g_currentManager;
}

void SetInvChi(const std::vector<double>& inv_chi) {
    g_invChi = inv_chi;
}

const std::vector<double>& GetInvChi() {
    return g_invChi;
}

void SetInteraction(radTInteraction* interaction, int n_elem, int nffc) {
    g_interaction = interaction;
    g_nElem = n_elem;
    g_nffc = nffc;
}

void SetLod(int* lod, int size) {
    // Deprecated: lod is now set/cleared on the C side by HACApK_build_hmatrix_varDOF_wrapper
    // This function is kept for API compatibility but does nothing
    (void)lod;
    (void)size;
}

void ClearLod() {
    // Deprecated: lod is now cleared on the C side by HACApK_build_hmatrix_varDOF_wrapper
    // This function is kept for API compatibility but does nothing
}

void ClearGlobalState() {
    // Clear all global callback state to prevent interference between solves
    g_currentManager = nullptr;
    g_interaction = nullptr;
    g_nElem = 0;
    g_nffc = 3;
    g_invChi.clear();
    // Note: g_hacapk_generation is NOT reset here - it continues incrementing
}

// FIX (2025-02-04): Get current generation for cache invalidation
uint64_t GetGeneration() {
    return g_hacapk_generation.load(std::memory_order_acquire);
}

// FIX (2025-02-04): Increment generation to invalidate all thread-local caches
void IncrementGeneration() {
    g_hacapk_generation.fetch_add(1, std::memory_order_release);
}

double ComputeEntry(int i, int j) {
    // i, j are 1-based ORIGINAL indices from HACApK.
    // The lod conversion is already done by cHACApK_fill_leafmtx_hyp:
    //   val = cHACApK_entry_ij(lodl[permuted_pos], lodt[permuted_pos], i_bemv)
    // so we receive original indices, NOT permuted indices.
    //
    // The kernel-specific system-matrix convention is delegated to
    // RadHACApKBase::ComputeSystemEntry so that each subclass (MSC,
    // PEEC, future BEM) can store exactly what HACApK needs.

    if (g_currentManager == nullptr) {
        std::cerr << "[HACApK] Error: g_currentManager is null in ComputeEntry" << std::endl;
        return 0.0;
    }

    int ndof = g_currentManager->GetNDOF();

    // Direct 0-based conversion (indices are already original, NOT permuted)
    int i0 = i - 1;
    int j0 = j - 1;

    if (i0 < 0 || i0 >= ndof || j0 < 0 || j0 >= ndof) {
        std::cerr << "[HACApK] Error: Invalid DOF indices: i0=" << i0
                  << " j0=" << j0 << " ndof=" << ndof << std::endl;
        return 0.0;
    }

    return g_currentManager->ComputeSystemEntry(i0, j0);
}

}  // namespace RadHACApKCallback

//=========================================================================
// C callback function for HACApK
// This is the function HACApK calls to get matrix elements
//=========================================================================

extern "C" {

// Optional per-call entry-function override (implements the HACApK_set_entry_func
// API declared in cHACApK_cpp.h, previously unimplemented).  When non-null,
// cHACApK_acaplus / fill routines fetch matrix entries from this kernel instead
// of the default MMM/MSC system matrix.  This lets callers (e.g. the
// stream-function (ACA+)+TSVD solver) factor an arbitrary rectangular kernel
// block with HACApK's ACA+ -- keeping ACA+ a single source of truth instead of
// re-porting it.  Default null => unchanged MMM behaviour.  Set/cleared
// synchronously around one factorization (GIL-serialized; no concurrent MMM
// build), so a plain (non-thread_local) pointer is sufficient.
static HACApK_entry_func g_entry_override = NULL;

void HACApK_set_entry_func(HACApK_entry_func func) { g_entry_override = func; }
void HACApK_clear_entry_func(void) { g_entry_override = NULL; }

double cHACApK_entry_ij(int i, int j, int i_bemv) {
    if (g_entry_override != NULL) return g_entry_override(i, j, i_bemv);
    (void)i_bemv;  // Unused in Radia
    return RadHACApKCallback::ComputeEntry(i, j);
}

}  // extern "C"

//=========================================================================
// RadHACApKBase Implementation (kernel-agnostic H-matrix lifecycle)
//=========================================================================

RadHACApKBase::RadHACApKBase()
    : m_leafmtxp(nullptr)
    , m_control(nullptr)
    , m_valid(false)
    , m_ndof(0)
    , m_n_elem(0)
    , m_diag_cached(false)
    , m_defl_alpha(0.0)
    , m_defl_nplaq(0)
{
}

RadHACApKBase::~RadHACApKBase() {
    FreeResources();
}

//=========================================================================
// RadHACApKMSCManager Implementation (MMM / MSC kernel)
//=========================================================================

RadHACApKMSCManager::RadHACApKMSCManager(radTInteraction* interaction)
    : RadHACApKBase()
    , m_interaction(interaction)
    , m_nffc(3)
    , m_is_6dof(false)
    , m_is_5dof(false)
    , m_is_mixed_dof(false)
    , m_geometry_3dof_ready(false)
    , m_flat_N_ready(false)
{
    // Hash-based cache is initialized automatically
}

void RadHACApKMSCManager::FreeStarHMatrix() {
    if (m_star_hlu_root) { cHACApK_hlu_free_factors(m_star_hlu_root); m_star_hlu_root = nullptr; }
    if (m_star_leafmtxp || m_star_control)
        HACApK_free_hmatrix_wrapper(m_star_leafmtxp, m_star_control);
    if (m_star_leafmtxp) { HACApK_free_leafmtxp(m_star_leafmtxp); m_star_leafmtxp = nullptr; }
    if (m_star_control)  { HACApK_free_lcontrol(m_star_control);  m_star_control  = nullptr; }
    m_star_hmat_built = false;
}

RadHACApKMSCManager::~RadHACApKMSCManager() { FreeStarHMatrix(); }

void RadHACApKMSCManager::BuildLoopBasis(double alpha) {
    // Loop (cycle) deflation basis from the element-adjacency graph. The bulk is
    // the SHORT CYCLES (3- and 4-cycles), which span the contractible part of
    // the null space. For multiply-connected bodies (first Betti number b1 > 0,
    // e.g. a closed core / ring) these are completed by b1 global "belt" loops
    // via a tree-cotree fundamental-cycle GF(2) pass (a belted tree) so the
    // basis spans the FULL null space for any topology. Adjacency + shared-face
    // DOF come from coincident face centroids (poly->FaceCenter[f]); DOF of
    // element e face f is dofOffset(e) + f. General for any conforming mesh;
    // O(N) for the short cycles (bounded element degree). Installs via SetDeflation.
    if (!m_interaction || m_n_elem <= 0) { SetDeflation({}, {}, {}, 0.0); return; }

    auto dofOffset = [&](int e) {
        return m_dof_offset.empty() ? e * GetUniformNFFC() : m_dof_offset[e];
    };

    // 1. Collect (centroid, dof, elem) for every element face + length scale.
    struct FaceRec { double c[3]; int dof; int elem; };
    std::vector<FaceRec> faces;
    faces.reserve((size_t)m_n_elem * 6);
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (int e = 0; e < m_n_elem; e++) {
        radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(m_interaction->g3dRelaxPtrVect[e]);
        if (!poly) continue;
        int nf = poly->AmOfFaces;
        for (int f = 0; f < nf; f++) {
            const TVector3d& fc = poly->FaceCenter[f];
            FaceRec r; r.c[0] = fc.x; r.c[1] = fc.y; r.c[2] = fc.z;
            r.dof = dofOffset(e) + f; r.elem = e;
            faces.push_back(r);
            for (int d = 0; d < 3; d++) { lo[d] = std::min(lo[d], r.c[d]); hi[d] = std::max(hi[d], r.c[d]); }
        }
    }
    double diag = std::sqrt((hi[0]-lo[0])*(hi[0]-lo[0]) + (hi[1]-lo[1])*(hi[1]-lo[1]) + (hi[2]-lo[2])*(hi[2]-lo[2]));
    double tol = 1.0e-7 * (diag > 0.0 ? diag : 1.0);
    auto key = [&](const double* c) {
        return std::make_tuple((long long)std::llround(c[0]/tol),
                               (long long)std::llround(c[1]/tol),
                               (long long)std::llround(c[2]/tol));
    };

    // 2. Group faces by centroid; internal faces (exactly 2) define adjacency.
    std::map<std::tuple<long long, long long, long long>, std::vector<int>> grp;
    for (int i = 0; i < (int)faces.size(); i++) grp[key(faces[i].c)].push_back(i);
    std::vector<std::vector<int>> nbr(m_n_elem);
    std::map<std::pair<int,int>, std::pair<int,int>> sdof;  // (ea,eb) -> (dofa,dofb)
    for (auto& kv : grp) {
        if (kv.second.size() != 2) continue;
        const FaceRec& A = faces[kv.second[0]];
        const FaceRec& B = faces[kv.second[1]];
        if (A.elem == B.elem) continue;
        nbr[A.elem].push_back(B.elem);
        nbr[B.elem].push_back(A.elem);
        sdof[{A.elem, B.elem}] = {A.dof, B.dof};
        sdof[{B.elem, A.elem}] = {B.dof, A.dof};
    }
    for (auto& v : nbr) { std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); }
    auto isAdj = [&](int a, int b) { return sdof.count({a, b}) > 0; };

    // 2b. Index undirected graph edges (for the GF(2) belt completion, step 4).
    std::map<std::pair<int,int>, int> eidx;
    for (auto& kv : sdof) {
        int a = kv.first.first, b = kv.first.second;
        if (a < b) eidx.emplace(std::make_pair(a, b), (int)eidx.size());
    }
    auto edgeId = [&](int a, int b) {
        return eidx[std::make_pair(std::min(a, b), std::max(a, b))];
    };

    // 3. Enumerate short cycles (3 + 4), dedup by sorted node set, build CSR.
    //    installCycle also records each added cycle's GF(2) edge-vector so the
    //    belt-loop completion (step 4) can test homological independence.
    std::set<std::vector<int>> seen;
    std::vector<int> offs; offs.push_back(0);
    std::vector<int> cdofs; std::vector<double> csign;
    std::vector<std::vector<int>> cycleEdges;   // GF(2) edge-vector per cycle
    auto installCycle = [&](const std::vector<int>& ordered) {
        int k = (int)ordered.size();
        std::vector<int> dbuf; std::vector<double> sbuf; std::vector<int> ebuf;
        for (int t = 0; t < k; t++) {
            auto it = sdof.find({ordered[t], ordered[(t+1) % k]});
            if (it == sdof.end()) return;        // not a closed walk: skip
            // +/-1 (unit density) on the two coincident DOFs of each shared face.
            // This is the discrete NULL mode of N (zero normal field at every
            // collocation point: N L ~ 0, verified field-exact in the linear
            // regime). NOTE: it is NOT charge-neutral (m_e = A_out - A_in != 0)
            // on distorted hexes, but charge-neutrality is the WRONG requirement
            // for deflation -- field-nullness (N L = 0) is what matters, and a
            // 1/area "charge-neutral" reweighting breaks N L = 0 (changes the
            // field by ~100%, verified). The residual nonlinear B_z drift comes
            // from a different mechanism (non-uniform chi makes span(L) not
            // A-invariant), addressed at the solver level, not the basis.
            dbuf.push_back(it->second.first);  sbuf.push_back(1.0);
            dbuf.push_back(it->second.second); sbuf.push_back(-1.0);
            ebuf.push_back(edgeId(ordered[t], ordered[(t+1) % k]));
        }
        for (size_t i = 0; i < dbuf.size(); i++) { cdofs.push_back(dbuf[i]); csign.push_back(sbuf[i]); }
        offs.push_back((int)cdofs.size());
        std::sort(ebuf.begin(), ebuf.end());
        cycleEdges.push_back(ebuf);
    };
    auto addCycle = [&](std::vector<int> ordered) {
        std::vector<int> kk = ordered; std::sort(kk.begin(), kk.end());
        if (seen.count(kk)) return; seen.insert(kk);
        installCycle(ordered);
    };
    for (int e = 0; e < m_n_elem; e++) {
        const std::vector<int>& Ne = nbr[e];
        for (size_t ai = 0; ai < Ne.size(); ai++) {
            int a = Ne[ai];
            for (size_t bi = ai + 1; bi < Ne.size(); bi++) {
                int b = Ne[bi];
                if (isAdj(a, b)) addCycle({e, a, b});                 // 3-cycle
                for (int d : nbr[a]) {                                // 4-cycle e-a-d-b
                    if (d == e || d == a || d == b) continue;
                    if (isAdj(d, b)) addCycle({e, a, d, b});
                }
            }
        }
    }

    // 4. Belt loops for multiply-connected bodies (first Betti number b1 > 0).
    //    Short cycles span only the contractible part (boundaries of plaquette
    //    2-cells). The full null space (= cycle space of the element graph) also
    //    contains b1 = (E - V + C) - rank(short) global "belt" loops around the
    //    holes. We complete it with a tree-cotree fundamental-cycle GF(2) pass
    //    (a belted tree). Each accepted belt loop is a genuine exact null mode
    //    (a divergence-free circulation). For simply-connected bodies the short
    //    cycles already span everything, so beltNeeded == 0 and step 4 is a no-op
    //    (only the cheap GF(2) rank check runs; the installed DOF vectors are
    //    unchanged from the short-cycle-only path).
    {
        // 4a. BFS spanning forest over the active sub-graph (deg >= 1):
        //     parent/depth (for fundamental-cycle paths) + V_active + #components.
        std::vector<int> parent(m_n_elem, -1), depth(m_n_elem, -1);
        int Vp = 0, C = 0;
        std::vector<int> q;
        for (int s = 0; s < m_n_elem; s++) {
            if (!nbr[s].empty()) Vp++;
            if (depth[s] >= 0 || nbr[s].empty()) continue;
            C++; depth[s] = 0; q.clear(); q.push_back(s);
            for (size_t h = 0; h < q.size(); h++) {
                int u = q[h];
                for (int v : nbr[u]) if (depth[v] < 0) {
                    depth[v] = depth[u] + 1; parent[v] = u; q.push_back(v);
                }
            }
        }
        int E = (int)eidx.size();
        int b1 = E - Vp + C;   // graph cycle rank = dim ker N

        // 4b. GF(2) basis seeded with short cycles (sparse rows keyed by pivot
        //     edge); rankShort = dimension of the contractible subspace.
        std::map<int, std::vector<int>> basis;
        auto symdiff = [](const std::vector<int>& A, const std::vector<int>& B) {
            std::vector<int> R; size_t i = 0, j = 0;
            while (i < A.size() && j < B.size()) {
                if (A[i] < B[j]) R.push_back(A[i++]);
                else if (B[j] < A[i]) R.push_back(B[j++]);
                else { i++; j++; }
            }
            while (i < A.size()) R.push_back(A[i++]);
            while (j < B.size()) R.push_back(B[j++]);
            return R;
        };
        auto reduceRow = [&](std::vector<int> row) {
            while (!row.empty()) {
                auto it = basis.find(row[0]);
                if (it == basis.end()) break;
                row = symdiff(row, it->second);
            }
            return row;
        };
        int rankShort = 0;
        for (auto& ev : cycleEdges) {
            if (rankShort >= b1) break;   // cycle space saturated: rest are dependent
            std::vector<int> r = reduceRow(ev);
            if (!r.empty()) { basis[r[0]] = r; rankShort++; }
        }
        int beltNeeded = b1 - rankShort;

        // 4c. complete with cotree fundamental cycles until beltNeeded added.
        if (beltNeeded > 0) {
            auto treePath = [&](int a, int b) {
                std::vector<int> pa, pb; int x = a, y = b;
                while (depth[x] > depth[y]) { pa.push_back(x); x = parent[x]; }
                while (depth[y] > depth[x]) { pb.push_back(y); y = parent[y]; }
                while (x != y) { pa.push_back(x); x = parent[x]; pb.push_back(y); y = parent[y]; }
                pa.push_back(x);   // LCA
                for (int i = (int)pb.size() - 1; i >= 0; i--) pa.push_back(pb[i]);
                return pa;         // ordered walk a -> ... -> b (closes via {a,b})
            };
            for (auto& kv : eidx) {
                if (beltNeeded <= 0) break;
                int a = kv.first.first, b = kv.first.second;
                if (depth[a] < 0 || depth[b] < 0) continue;
                if (parent[a] == b || parent[b] == a) continue;       // tree edge
                std::vector<int> walk = treePath(a, b);
                std::vector<int> ev;                                  // fundamental cycle edges
                for (size_t t = 0; t + 1 < walk.size(); t++) ev.push_back(edgeId(walk[t], walk[t+1]));
                ev.push_back(edgeId(a, b));
                std::sort(ev.begin(), ev.end());
                std::vector<int> r = reduceRow(ev);
                if (r.empty()) continue;                              // homologically trivial
                std::vector<int> kk = walk; std::sort(kk.begin(), kk.end());
                if (!seen.count(kk)) { seen.insert(kk); installCycle(walk); }
                basis[r[0]] = r; beltNeeded--;
            }
        }
    }

    // 5. Auto-scale the shift strength when alpha <= 0 (recommended default).
    //    The matrix-free shift moves a loop mode v by alpha * (its L L^T
    //    eigenvalue), so the WORST-shifted mode is lifted by alpha * lambda_max
    //    where lambda_max = spectral radius of L L^T. Capping that worst lift at
    //    the physical eigenvalue scale (~ mean |N_ii|, the self-demagnetization)
    //    gives  alpha = mean|N_ii| / lambda_max(L L^T).  This is MESH-ROBUST: a
    //    fixed alpha=1 (fine on a regular cube, where L L^T is near-uniform with
    //    small lambda_max) over-shoots on an irregular mesh (chamfers, varying
    //    connectivity) whose overcomplete cycles correlate, giving a large
    //    lambda_max -- throwing the most-overlapped loops far past the physical
    //    band and WORSENING conditioning. Using the diagonal d_max instead of
    //    the true lambda_max underestimates the lift (overlapping cycles make
    //    lambda_max >> d_max) and over-shoots ~10x; the power iteration below is
    //    the principled cap. It is also chi-independent and self-regulating:
    //    where the iron saturates (small chi, large 1/chi) the loop cluster is
    //    no longer the smallest-eigenvalue bottleneck and a shift of order
    //    mean|N_ii| << 1/chi barely moves it, so deflation does no harm.
    double alpha_use = alpha;
    if (alpha <= 0.0 && (int)offs.size() > 1) {
        int nplaq = (int)offs.size() - 1;
        // lambda_max(L L^T) via power iteration (L L^T x: per cycle c=sum sign*x,
        // then y[dof] += sign*c). Cheap: ~20 matvecs of O(nnz(L)) = O(N).
        std::vector<double> xv(m_ndof, 0.0), yv(m_ndof, 0.0);
        // Seed inside range(L) with a VARIED (pseudo-random) pattern: a uniform
        // seed is annihilated by the balanced +/-1 cycle signs (each cycle's
        // signs sum to zero, so L L^T * uniform = 0) and the power iteration
        // would stall at the init value. A deterministic LCG keeps it reproducible.
        unsigned int rng = 2463534242u;
        for (int d : cdofs) { rng = rng * 1664525u + 1013904223u;
                              xv[d] = (double)(rng >> 9) / 8388608.0 - 0.5; }
        double nrm = 0.0; for (double v : xv) nrm += v * v; nrm = std::sqrt(nrm);
        if (nrm > 0.0) for (double& v : xv) v /= nrm;
        double lam_max = 1.0;
        for (int it = 0; it < 20; it++) {
            std::fill(yv.begin(), yv.end(), 0.0);
            for (int p = 0; p < nplaq; p++) {
                double c = 0.0;
                for (int k = offs[p]; k < offs[p + 1]; k++) c += csign[k] * xv[cdofs[k]];
                for (int k = offs[p]; k < offs[p + 1]; k++) yv[cdofs[k]] += csign[k] * c;
            }
            double ny = 0.0; for (double v : yv) ny += v * v; ny = std::sqrt(ny);
            if (ny <= 0.0) break;
            lam_max = ny;                           // ||L L^T x|| with ||x||=1
            for (int i = 0; i < m_ndof; i++) xv[i] = yv[i] / ny;
        }
        double sumN = 0.0;
        for (int i = 0; i < m_ndof; i++) sumN += std::fabs(GetInteractionMatrixElement(i, i));
        double meanN = (m_ndof > 0) ? sumN / (double)m_ndof : 1.0;
        // CONSERVATIVE safety factor: keep alpha WELL BELOW the over-shift cliff.
        // meanN/lam_max alone lands ~0.07 on the C-type, which is AT the cliff
        // (non-reproducible: 3871 vs 9019 iters across runs at 6^3 mu=1e5) and
        // over the cliff at 10^3 (diverges). The cliff moves DOWN with mesh size
        // while lam_max(L L^T) stays ~constant, so a formula cannot fully track
        // it -- this factor is a STOPGAP that lands ~0.0035 (safe at 6^3 and
        // 10^3: 1.75x and 4x speedup, field-exact). NOTE: auto-alpha is a fragile
        // convenience and may be retired in favour of an alpha-free cure
        // (tree-cotree gauge) or default-off; do not treat it as load-bearing.
        const double DEFL_ALPHA_SAFETY = 0.05;
        alpha_use = (lam_max > 0.0) ? DEFL_ALPHA_SAFETY * meanN / lam_max : 0.0;
    }
    SetDeflation(offs, cdofs, csign, alpha_use);
}

//=========================================================================
// Loop-star gauge (ALPHA-FREE): sparse star basis + reduced BiCGSTAB sandwich.
// Build T (columns span the complement of the loop null space) and solve the
// reduced non-singular system T^T A T y = T^T b; the H-matrix is UNCHANGED.
//=========================================================================

/* Diagnostic counters for the star-basis build (set by BuildStarBasis, read via
 * RadGetStarBasisStats). Used to localize the low-mu_r loop-star field error to
 * star-basis imperfection on complex (non-cube) geometry. */
static int g_sb_ndof = 0, g_sb_nface = 0, g_sb_boundary = 0, g_sb_internal = 0,
           g_sb_skip_ge3 = 0, g_sb_skip_same = 0, g_sb_nstar = 0, g_sb_ncomp = 0;
extern "C" void RadGetStarBasisStats(int* out /* len 8 */) {
    if (!out) return;
    out[0]=g_sb_ndof; out[1]=g_sb_nface; out[2]=g_sb_boundary; out[3]=g_sb_internal;
    out[4]=g_sb_skip_ge3; out[5]=g_sb_skip_same; out[6]=g_sb_nstar; out[7]=g_sb_ncomp;
}

/* Diagnostic for the KEEP-LOOPS block Gauss-Seidel in SolveLoopStar:
 * [0]=n_loop, [1]=res0 (||b - A sigma_S|| after the initial star solve),
 * [2]=res_final_rel (||b - A sigma|| / ||b|| after GS), [3]=GS sweeps,
 * [4]=total CG iterations (loop blocks), [5]=res_final (absolute). */
static double g_kl_dbg[6] = {0,0,0,0,0,0};
extern "C" void RadGetKeepLoopStats(double* out /* len 6 */) {
    if (!out) return; for (int i=0;i<6;i++) out[i]=g_kl_dbg[i];
}

/* A_SS-as-H-matrix path. mode 0 = off (dense K-dense LU, the default);
 * mode 1 = build A_SS as a HACApK H-matrix in SolveLoopStar and VERIFY its
 * MatVec against the exact T^T A T (sandwich) operator (milestone: compression
 * accuracy). g_ass_dbg: [0]=n_star, [1]=H-matrix rel MatVec error vs exact,
 * [2]=BuildStarHMatrix return code, [3]=build time (s). */
static int g_loopstar_hmat_mode = 0;
extern "C" void RadSetLoopStarHMatMode(int m) { g_loopstar_hmat_mode = m; }
extern "C" int  RadGetLoopStarHMatMode(void) { return g_loopstar_hmat_mode; }
// H-ILU truncation tolerance for the A_SS H-LU factorization (mode 2). Larger =
// more aggressive rk truncation during the Schur updates (cheaper, "incomplete"
// -> H-ILU); smaller = more accurate (-> H-LU). Applied via
// cHACApK_hlu_set_trunc_tol immediately before the factor.
static double g_star_hlu_trunc_tol = 1.0e-3;
extern "C" void   RadSetStarHLUTruncTol(double t) { g_star_hlu_trunc_tol = (t > 0.0) ? t : 1.0e-3; }
extern "C" double RadGetStarHLUTruncTol(void)     { return g_star_hlu_trunc_tol; }
// g_ass_dbg: [0]=n_star, [1]=A_SS H-matrix rel MatVec error vs exact T^T A T,
// [2]=BuildStarHMatrix rc, [3]=A_SS H-build time (s), [4]=H-ILU-preconditioned
// BiCGSTAB rel err vs K-dense y, [5]=BiCGSTAB iters (-1 = factor failed),
// [6]=H-LU factor time (s), [7]=H-ILU BiCGSTAB solve time (s),
// [8]=preconditioner residual ||A_SS (M^-1 b) - b||/||b||, [9]=trunc tol used,
// [10]=A_SS nlf (total leaves), [11]=A_SS nlfkt (low-rank leaves),
// [12]=A_SS ktmax (max rank), [13]=A_SS compression (H-bytes / dense-bytes).
// Compression stats [10..13] are recorded by BuildStarHMatrix BEFORE the H-LU
// factor converts the leaves (probe for "does A_SS actually compress?").
static double g_ass_dbg[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
extern "C" void RadGetStarHMatStats(double* out /* len 14 */) {
    if (!out) return; for (int i=0;i<14;i++) out[i]=g_ass_dbg[i];
}
// Recorded by BuildStarHMatrix so the probe survives even when SolveLoopStar
// is not the entry point. Read into g_ass_dbg[10..13] in SolveLoopStar.
static double g_star_compress[4] = {0,0,0,0};   // nlf, nlfkt, ktmax, compression

int RadHACApKMSCManager::BuildStarBasis() {
    m_star_offsets.clear(); m_star_dofs.clear(); m_star_coeffs.clear(); m_n_star = 0;
    if (!m_interaction || m_n_elem <= 0) return 0;

    auto dofOffset = [&](int e) {
        return m_dof_offset.empty() ? e * GetUniformNFFC() : m_dof_offset[e];
    };

    // 1. faces (centroid, dof, elem) + length scale (same as BuildLoopBasis).
    struct FaceRec { double c[3]; int dof; int elem; };
    std::vector<FaceRec> faces; faces.reserve((size_t)m_n_elem * 6);
    double lo[3] = {1e300,1e300,1e300}, hi[3] = {-1e300,-1e300,-1e300};
    for (int e = 0; e < m_n_elem; e++) {
        radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(m_interaction->g3dRelaxPtrVect[e]);
        if (!poly) continue;
        for (int f = 0; f < poly->AmOfFaces; f++) {
            const TVector3d& fc = poly->FaceCenter[f];
            FaceRec r; r.c[0]=fc.x; r.c[1]=fc.y; r.c[2]=fc.z; r.dof=dofOffset(e)+f; r.elem=e;
            faces.push_back(r);
            for (int d=0; d<3; d++){ lo[d]=std::min(lo[d],r.c[d]); hi[d]=std::max(hi[d],r.c[d]); }
        }
    }
    double diag = std::sqrt((hi[0]-lo[0])*(hi[0]-lo[0])+(hi[1]-lo[1])*(hi[1]-lo[1])+(hi[2]-lo[2])*(hi[2]-lo[2]));
    double tol = 1.0e-7 * (diag > 0.0 ? diag : 1.0);
    auto key = [&](const double* c){ return std::make_tuple((long long)std::llround(c[0]/tol),(long long)std::llround(c[1]/tol),(long long)std::llround(c[2]/tol)); };

    // NOTE on IMA: symmetric planes (sign>0) need NO special handling here -- they
    // double the boundary-face charge (2*sigma, not null), preserving the topological
    // +-1 loop null structure, and the standard star is field-exact (verified +x:
    // dB/B ~ 2e-11). ANTISYMMETRIC planes: the topological loops STAY null under IMA
    // (verified ||N_ima @ loop|| ~ 1e-15), so this standard star is still correct for
    // the loop part; ker(N_ima) just gains a FEW extra LOCAL modes confined to the
    // antisymmetric-plane slab. Those are handled separately by
    // ComputePlaneSlabAddedModes() + a reduced-space deflation in SolveLoopStar (NOT
    // by changing the star here). So no antisym branch is needed in this function.

    // 2. group by centroid: size-2 = internal (adjacency + shared DOF pair),
    //    size-1 = boundary face.
    std::map<std::tuple<long long,long long,long long>, std::vector<int>> grp;
    for (int i=0; i<(int)faces.size(); i++) grp[key(faces[i].c)].push_back(i);
    std::vector<std::vector<int>> nbr(m_n_elem);
    std::map<std::pair<int,int>, std::pair<int,int>> sdof;
    std::vector<int> boundaryDofs;
    std::vector<std::pair<int,int>> internalPairs;
    g_sb_skip_ge3 = 0; g_sb_skip_same = 0;
    for (auto& kv : grp) {
        if (kv.second.size() == 1) { boundaryDofs.push_back(faces[kv.second[0]].dof); continue; }
        if (kv.second.size() != 2) { g_sb_skip_ge3++; continue; }
        const FaceRec& A = faces[kv.second[0]]; const FaceRec& B = faces[kv.second[1]];
        if (A.elem == B.elem) { g_sb_skip_same++; continue; }
        nbr[A.elem].push_back(B.elem); nbr[B.elem].push_back(A.elem);
        sdof[{A.elem,B.elem}] = {A.dof,B.dof}; sdof[{B.elem,A.elem}] = {B.dof,A.dof};
        internalPairs.push_back({A.dof, B.dof});
    }
    g_sb_nface = (int)faces.size();
    g_sb_boundary = (int)boundaryDofs.size();
    g_sb_internal = (int)internalPairs.size();
    g_sb_ndof = m_ndof;
    for (auto& v : nbr){ std::sort(v.begin(),v.end()); v.erase(std::unique(v.begin(),v.end()),v.end()); }

    // 3. BFS components -> drop exactly one divergence per connected component
    //    (the per-component sum of internal divergences vanishes: 1 redundancy).
    std::vector<int> comp(m_n_elem, -1); int ncomp = 0; std::vector<int> q;
    for (int s=0; s<m_n_elem; s++){
        if (comp[s] >= 0 || nbr[s].empty()) continue;
        q.clear(); q.push_back(s); comp[s]=ncomp;
        for (size_t h=0; h<q.size(); h++){ int u=q[h]; for(int v: nbr[u]) if(comp[v]<0){ comp[v]=ncomp; q.push_back(v);} }
        ncomp++;
    }

    // 4. assemble the sparse star CSR (columns): boundary unit vectors,
    //    symmetric internal (sigma_A+sigma_B), per-element divergence.
    m_star_offsets.push_back(0);
    auto addCol = [&](const std::vector<std::pair<int,double>>& entries){
        for (auto& e : entries){ m_star_dofs.push_back(e.first); m_star_coeffs.push_back(e.second); }
        m_star_offsets.push_back((int)m_star_dofs.size());
    };
    for (int d : boundaryDofs) addCol({{d, 1.0}});
    for (auto& p : internalPairs) addCol({{p.first, 1.0},{p.second, 1.0}});
    std::vector<char> compHasRoot(ncomp > 0 ? ncomp : 1, 0);
    for (int e=0; e<m_n_elem; e++){
        if (nbr[e].empty()) continue;
        if (comp[e] >= 0 && !compHasRoot[comp[e]]) { compHasRoot[comp[e]] = 1; continue; }  // drop 1 / comp
        std::vector<std::pair<int,double>> col;
        for (int e2 : nbr[e]){ auto it=sdof.find({e,e2}); if(it==sdof.end()) continue;
            col.push_back({it->second.first, 1.0}); col.push_back({it->second.second, -1.0}); }
        if (!col.empty()) addCol(col);
    }
    m_n_star = (int)m_star_offsets.size() - 1;
    g_sb_nstar = m_n_star; g_sb_ncomp = ncomp;

    // Star-column centroids (for the A_SS H-matrix cluster tree): each column's
    // geometric center = average of its support faces' centers. dofCoord maps a
    // sigma DOF -> its face center (faces still in scope).
    std::vector<double> dofCoord((size_t)m_ndof * 3, 0.0);
    for (const auto& fr : faces)
        if (fr.dof >= 0 && fr.dof < m_ndof) {
            dofCoord[3*(size_t)fr.dof+0] = fr.c[0];
            dofCoord[3*(size_t)fr.dof+1] = fr.c[1];
            dofCoord[3*(size_t)fr.dof+2] = fr.c[2];
        }
    m_star_coords.assign((size_t)m_n_star * 3, 0.0);
    for (int k = 0; k < m_n_star; k++) {
        double cx=0,cy=0,cz=0; int cnt=0;
        for (int j = m_star_offsets[k]; j < m_star_offsets[k+1]; j++) {
            int d = m_star_dofs[j];
            if (d>=0 && d<m_ndof) { cx+=dofCoord[3*(size_t)d]; cy+=dofCoord[3*(size_t)d+1]; cz+=dofCoord[3*(size_t)d+2]; cnt++; }
        }
        if (cnt>0) { m_star_coords[3*(size_t)k]=cx/cnt; m_star_coords[3*(size_t)k+1]=cy/cnt; m_star_coords[3*(size_t)k+2]=cz/cnt; }
    }
    return m_n_star;
}

//=========================================================================
// A_SS = S^T A S as a HACApK H-matrix ("HACApK-only star" path).
// The reduced star operator is built as its OWN HACApK leafmtxp via the ACA
// entry-function override (StarSSEntry = S_i^T A S_j) on a cluster tree over the
// star-column centroids m_star_coords. This replaces the dense K-dense LU (caps
// ~15^3 / 8 GB) so the star block can be H-LU/H-ILU preconditioned and solved
// with HACApK alone at 20^3+.  Tree-cotree loop-star H-matrix (Ida-endorsed).
//=========================================================================

// A_SS entry: small double sum over the sparse supports of star columns si, sj.
double RadHACApKMSCManager::StarSSEntry(int si, int sj) const {
    if (si < 0 || sj < 0 || si >= m_n_star || sj >= m_n_star) return 0.0;
    double acc = 0.0;
    for (int a = m_star_offsets[si]; a < m_star_offsets[si+1]; a++) {
        int dp = m_star_dofs[a]; double cp = m_star_coeffs[a];
        for (int b = m_star_offsets[sj]; b < m_star_offsets[sj+1]; b++)
            acc += cp * m_star_coeffs[b] * ComputeSystemEntry(dp, m_star_dofs[b]);
    }
    return acc;
}

// Manager bound for the A_SS ACA entry override (set only during BuildStarHMatrix).
// HACApK passes 1-based original indices (see ComputeEntry); convert to 0-based.
static RadHACApKMSCManager* g_star_entry_mgr = nullptr;
static double RadStarSSEntryOverride(int i, int j, int i_bemv) {
    (void)i_bemv;
    return g_star_entry_mgr ? g_star_entry_mgr->StarSSEntry(i - 1, j - 1) : 0.0;
}

int RadHACApKMSCManager::BuildStarHMatrix(double eps, int leaf, double eta) {
    if (m_n_star <= 0) { if (BuildStarBasis() <= 0) return -1; }
    if (m_star_hmat_built) return m_n_star;
    if ((int)m_star_coords.size() != m_n_star * 3) return -1;
    m_star_leafmtxp = HACApK_alloc_leafmtxp();
    m_star_control  = HACApK_alloc_lcontrol();
    if (!m_star_leafmtxp || !m_star_control) return -1;
    // Build A_SS via the entry override. nffc=1 (each star column = one DOF at its
    // centroid); ndim=3. ACA calls RadStarSSEntryOverride(i,j) = S_i^T A S_j.
    g_star_entry_mgr = this;
    HACApK_set_entry_func(RadStarSSEntryOverride);
    int rc = HACApK_build_hmatrix_wrapper(
        m_star_leafmtxp, m_star_control, m_star_coords.data(),
        m_n_star, /*nffc=*/1, /*ndim=*/3, eps, leaf, eta, /*print=*/0);
    HACApK_clear_entry_func();
    g_star_entry_mgr = nullptr;
    if (rc != 0) { return -1; }
    m_star_hmat_built = true;
    // Compression probe (read-only, BEFORE the H-LU factor converts the leaves):
    // does A_SS actually compress as an H-matrix?  nlfkt/nlf near 0 -> near-dense
    // -> H-LU degenerates to dense LU (the materialize wall).  H-bytes/dense-bytes
    // is the storage ratio (1.0 = no compression).
    {
        int nlf   = HACApK_leafmtxp_get_nlf(m_star_leafmtxp);
        int nlfkt = HACApK_leafmtxp_get_nlfkt(m_star_leafmtxp);
        int ktmax = HACApK_leafmtxp_get_ktmax(m_star_leafmtxp);
        int64_t hb = 0, db = 0;
        HACApK_get_memory_stats(m_star_leafmtxp, &hb, &db);
        g_star_compress[0] = (double)nlf;
        g_star_compress[1] = (double)nlfkt;
        g_star_compress[2] = (double)ktmax;
        g_star_compress[3] = (db > 0) ? (double)hb / (double)db : 0.0;
    }
    return m_n_star;
}

void RadHACApKMSCManager::StarHMatVec(const std::vector<double>& x, std::vector<double>& y) {
    y.assign(m_n_star, 0.0);
    if (!m_star_hmat_built || !m_star_leafmtxp || !m_star_control) return;
    int nd = HACApK_leafmtxp_get_nd(m_star_leafmtxp);
    HACApK_matvec_wrapper(m_star_leafmtxp, m_star_control,
                          const_cast<double*>(x.data()), y.data(), nd);
}

// z = A_SS^{-1} r (approx) via the cached A_SS H-LU factors. NOTE: the H-LU factor
// (cHACApK_hlu_factor_leafmtxp) converts the leafmtxp leaves to internal layout, so
// StarHMatVec must NOT be used after the H-LU is factored -- use the exact T-sandwich
// (redMatVec) as the operator and this as the preconditioner.
bool RadHACApKMSCManager::StarHLUApply(const std::vector<double>& r, std::vector<double>& z) {
    z.assign(m_n_star, 0.0);
    if (!m_star_hlu_root || !m_star_control || m_n_star <= 0) return false;
    return cHACApK_hlu_apply(m_star_hlu_root, m_star_control,
                             r.data(), z.data(), m_n_star) == 0;
}

// True iff the interaction has an ANTISYMMETRIC IMA plane (sign<0). These add a
// few LOCAL plane-coupled modes to ker(N_ima); ComputePlaneSlabAddedModes finds
// them and SolveLoopStar deflates them so the gauge is field-exact.
bool RadHACApKMSCManager::HasAntisymmetricIMAPlane() const {
    if (!m_interaction) return false;
    int sym = m_interaction->GetIMASymmetry();
    return ((sym & radTInteraction::IMA_X) && m_interaction->GetIMASignX() < 0)
        || ((sym & radTInteraction::IMA_Y) && m_interaction->GetIMASignY() < 0)
        || ((sym & radTInteraction::IMA_Z) && m_interaction->GetIMASignZ() < 0);
}

// Antisym-IMA O(N) gauge: the topological star gauges the loops (which stay null
// under IMA -- verified ||N_ima @ loop|| ~ 1e-15), but ker(N_ima) ALSO contains a
// few LOCAL "plane-coupled" modes confined to the antisymmetric-plane slab (count
// = nP-1, nP = #element-faces ON the plane). Find them, SVD-FREE, in two stages:
//   (1) densify N_ima on the slab (plane-touching elements) via
//       GetInteractionMatrixElement (IMA-folded), then ker(N_sl) via the Gram
//       eigendecomposition G = N_sl^T N_sl (LAPACKE_dsyevd, NOT dgesdd) -- the null
//       gap WIDENS under squaring (s~1e-11 -> lambda~1e-22) so detection is robust.
//   (2) isolate the PLANE-COUPLED subspace BASIS-INDEPENDENTLY: the intra-slab loops
//       have ZERO weight on the antisym-plane faces and are ALREADY gauged by the star
//       T, so deflate only their orthogonal complement within ker(N_sl). Build the
//       small plane-restriction Gram H = (Q_plane)^T (Q_plane) (m x m, Q_plane = the
//       plane-face rows of the null basis Q); the eigenvectors of H with NONZERO
//       eigenvalue are the combos of Q carrying plane weight = K_plane (count nP-1).
// (The earlier per-vector "plane-weight >= 5%" filter was basis-DEPENDENT -- with a
//  different but equally valid orthonormal null basis it could drop modes spanning
//  K_plane; the H-restriction above is basis-independent. Verified in
//  C:/temp/gram_vs_svd_indexing.py: option B recovers K_plane to angle 0 from BOTH the
//  SVD and Gram bases, while the per-vector filter misses it.)
// Cache the kept (orthonormal: Q orthonormal x hvec orthonormal) modes as a sigma-DOF
// CSR. SolveLoopStar deflates them via SetDeflation(alpha L L^T): robust because they
// are FEW + ORTHONORMAL, and the T^T(.)T sandwich annihilates any loop-overlap part
// (T excludes the loops). One-time (cached via m_added_built). Validated field-exact
// in examples/mmm_eigenvalue_study + ctype_loopstar_test.py ('+x-z').
int RadHACApKMSCManager::ComputePlaneSlabAddedModes() {
    if (m_added_built) return m_n_added;
    m_added_built = true;
    m_added_offsets.clear(); m_added_dofs.clear(); m_added_coeffs.clear(); m_n_added = 0;
    if (!m_interaction || m_n_elem <= 0 || !HasAntisymmetricIMAPlane()) return 0;

    int sym = m_interaction->GetIMASymmetry();
    bool antiX = (sym & radTInteraction::IMA_X) && m_interaction->GetIMASignX() < 0;
    bool antiY = (sym & radTInteraction::IMA_Y) && m_interaction->GetIMASignY() < 0;
    bool antiZ = (sym & radTInteraction::IMA_Z) && m_interaction->GetIMASignZ() < 0;
    auto dofOffset = [&](int e){ return m_dof_offset.empty() ? e * GetUniformNFFC() : m_dof_offset[e]; };

    // length scale for the plane tolerance
    double lo[3]={1e300,1e300,1e300}, hi[3]={-1e300,-1e300,-1e300};
    for (int e=0;e<m_n_elem;e++){
        radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(m_interaction->g3dRelaxPtrVect[e]);
        if(!poly) continue;
        for(int f=0; f<poly->AmOfFaces; f++){ const TVector3d& fc=poly->FaceCenter[f];
            double c[3]={fc.x,fc.y,fc.z}; for(int d=0;d<3;d++){lo[d]=std::min(lo[d],c[d]); hi[d]=std::max(hi[d],c[d]);} }
    }
    double diag=std::sqrt((hi[0]-lo[0])*(hi[0]-lo[0])+(hi[1]-lo[1])*(hi[1]-lo[1])+(hi[2]-lo[2])*(hi[2]-lo[2]));
    double ptol = 1.0e-6 * (diag>0.0?diag:1.0);

    // plane-touching elements + antisym-plane face DOFs (for the plane-weight filter)
    std::vector<char> isSlab(m_n_elem, 0);
    std::set<int> planeFaceDofs;
    for (int e=0;e<m_n_elem;e++){
        radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(m_interaction->g3dRelaxPtrVect[e]);
        if(!poly) continue;
        for(int f=0; f<poly->AmOfFaces; f++){ const TVector3d& fc=poly->FaceCenter[f];
            if((antiX && std::fabs(fc.x)<ptol) || (antiY && std::fabs(fc.y)<ptol) || (antiZ && std::fabs(fc.z)<ptol)){
                isSlab[e]=1; planeFaceDofs.insert(dofOffset(e)+f);
            }
        }
    }
    std::vector<int> slabDofs;
    for (int e=0;e<m_n_elem;e++){ if(!isSlab[e]) continue;
        radTPolyhedron* poly=dynamic_cast<radTPolyhedron*>(m_interaction->g3dRelaxPtrVect[e]);
        int nf = poly? poly->AmOfFaces : 0;
        for(int f=0;f<nf;f++) slabDofs.push_back(dofOffset(e)+f);
    }
    int nsl = (int)slabDofs.size();
    if (nsl < 1) return 0;

    // densify N_ima on the slab (IMA-folded), column-major for LAPACK
    std::vector<double> Nslab((size_t)nsl*nsl, 0.0);
    for (int a=0;a<nsl;a++)
        for (int b=0;b<nsl;b++)
            Nslab[(size_t)b*nsl + a] = GetInteractionMatrixElement(slabDofs[a], slabDofs[b]);

    // Stage 1: ker(N_sl) via the Gram eigendecomposition G = N_sl^T N_sl (SPD), SVD-free.
    // dsyevd returns orthonormal eigenVECTORS as COLUMNS, eigenvalues ASCENDING; the null
    // space = eigenvectors with lambda < (1e-6*smax)^2 = 1e-12*lmax (same threshold as the
    // old dgesdd snull=1e-6*smax, squared). Eigvec k component j is G[j + k*nsl] = G(j,k).
    std::vector<double> G((size_t)nsl*nsl, 0.0);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, nsl, nsl, nsl,
                1.0, Nslab.data(), nsl, Nslab.data(), nsl, 0.0, G.data(), nsl);
    std::vector<double> evals(nsl);
    if (LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'U', nsl, G.data(), nsl, evals.data()) != 0) return 0;
    double lmax = (nsl>0) ? evals[nsl-1] : 0.0;          // ascending -> last is largest
    double lnull = 1.0e-12 * (lmax>0.0 ? lmax : 1.0);
    std::vector<int> nullCols;                           // eigenvector cols spanning ker(N_sl)
    for (int k=0;k<nsl;k++) if (evals[k] < lnull) nullCols.push_back(k);
    int m = (int)nullCols.size();
    if (m < 1) return 0;                                 // Q(j,t) = G[j + nullCols[t]*nsl]

    // Stage 2: isolate K_plane (basis-independent). H = (Q_plane)^T (Q_plane) (m x m),
    // Q_plane = plane-face rows of Q; eigenvectors of H with NONZERO eigenvalue carry
    // plane weight = K_plane (intra-slab loops have zero plane weight -> H-kernel, and
    // are already gauged by the star T). dsyevd again (SVD-free, m x m tiny).
    std::vector<char> isPlaneRow(nsl, 0);
    for (int j=0;j<nsl;j++) if (planeFaceDofs.count(slabDofs[j])) isPlaneRow[j]=1;
    std::vector<double> H((size_t)m*m, 0.0);
    for (int j=0;j<nsl;j++){
        if (!isPlaneRow[j]) continue;
        for (int t2=0;t2<m;t2++){ double q2 = G[(size_t)j + (size_t)nullCols[t2]*nsl];
            for (int t1=0;t1<m;t1++) H[(size_t)t1 + (size_t)t2*m] += G[(size_t)j + (size_t)nullCols[t1]*nsl]*q2; } }
    std::vector<double> hvals(m);
    if (LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'U', m, H.data(), m, hvals.data()) != 0) return 0;
    double hmax = (m>0) ? hvals[m-1] : 0.0;
    double hkeep = 1.0e-9 * (hmax>0.0 ? hmax : 1.0);

    // assemble each kept plane-coupled mode mode_c[j] = sum_t Q(j,t) * hvec_c[t]
    // (orthonormal). hvec_c component t is H[t + c*m]. Store as a sigma-DOF CSR column.
    m_added_offsets.push_back(0);
    std::vector<double> mode(nsl);
    for (int c=0;c<m;c++){
        if (hvals[c] <= hkeep) continue;                 // pure intra-slab loop -> skip (gauged by T)
        for (int j=0;j<nsl;j++){ double acc=0.0;
            for (int t=0;t<m;t++) acc += G[(size_t)j + (size_t)nullCols[t]*nsl]*H[(size_t)t + (size_t)c*m];
            mode[j]=acc; }
        for (int j=0;j<nsl;j++) if (std::fabs(mode[j]) > 1e-14){ m_added_dofs.push_back(slabDofs[j]); m_added_coeffs.push_back(mode[j]); }
        m_added_offsets.push_back((int)m_added_dofs.size());
    }
    m_n_added = (int)m_added_offsets.size() - 1;
    return m_n_added;
}

int RadHACApKMSCManager::SolveLoopStar(const std::vector<double>& b, std::vector<double>& sigma,
                                       double tol, int max_iter,
                                       const std::vector<double>& blockInverse,
                                       const std::vector<int>& blockOffsets) {
    if (m_n_star <= 0) { if (BuildStarBasis() <= 0) return -1; }
    const int n = m_ndof, ns = m_n_star;

    // KEEP-LOOPS: cache the topological loop (cycle) basis L for the loop-mode
    // back-substitution at the end of this routine. The reduced star solve gives
    // only sigma_S = S y_S (the loop-free part); the physical solution also has a
    // loop component sigma_L = L y_L which sigma_S drops. Because L = ker(N)
    // (N L ~ 0) and the star is loop-orthogonal (S^T L ~ 0), the missing part lives
    // purely in span(L) and is recovered exactly (linear / uniform-chi) by the
    // block back-substitution below. Build via BuildLoopBasis (alpha=1 skips its
    // auto-alpha power iteration), copy the cycle CSR out, then CLEAR the deflation
    // so MatVec stays the clean A during the reduced solve.
    if (!m_loop_built) {
        BuildLoopBasis(1.0);
        m_loop_offsets = m_defl_offsets;
        m_loop_dofs    = m_defl_dofs;
        m_loop_coeffs  = m_defl_signs;
        SetDeflation({}, {}, {}, 0.0);
        m_n_loop = m_loop_offsets.empty() ? 0 : (int)m_loop_offsets.size() - 1;
        m_loop_built = true;
    }

    auto applyT = [&](const std::vector<double>& y, std::vector<double>& s){
        s.assign(n, 0.0);
        for (int k=0; k<ns; k++){ double yk=y[k];
            for(int j=m_star_offsets[k]; j<m_star_offsets[k+1]; j++) s[m_star_dofs[j]] += m_star_coeffs[j]*yk; }
    };
    auto applyTt = [&](const std::vector<double>& r, std::vector<double>& y){
        y.assign(ns, 0.0);
        for (int k=0; k<ns; k++){ double acc=0.0;
            for(int j=m_star_offsets[k]; j<m_star_offsets[k+1]; j++) acc += m_star_coeffs[j]*r[m_star_dofs[j]];
            y[k]=acc; }
    };
    std::vector<double> ts(n), tAs(n);
    auto redMatVec = [&](const std::vector<double>& z, std::vector<double>& out){
        applyT(z, ts); MatVec(ts, tAs); applyTt(tAs, out);   // out = T^T A T z
    };

    // A_SS-as-H-matrix milestone (mode 1): build A_SS as its own HACApK H-matrix
    // and VERIFY its MatVec vs the exact T^T A T sandwich. Deflation is not set
    // yet here, so MatVec is the clean A. (mode 2 -- H-LU precondition + solve --
    // is the next increment.)
    if (g_loopstar_hmat_mode >= 1 && ns > 0) {
        auto th0 = std::chrono::high_resolution_clock::now();
        int rc = BuildStarHMatrix(1.0e-4, 10, 2.0);
        auto th1 = std::chrono::high_resolution_clock::now();
        g_ass_dbg[0] = (double)ns;
        g_ass_dbg[2] = (double)rc;
        g_ass_dbg[3] = std::chrono::duration<double>(th1 - th0).count();
        // A_SS compression probe recorded by BuildStarHMatrix (before any factor).
        for (int z=0; z<4; z++) g_ass_dbg[10+z] = g_star_compress[z];
        // StarHMatVec is only valid BEFORE the H-LU factor converts the leafmtxp
        // leaves to internal layout (mode 2).  Skip the verify on re-entry once the
        // A_SS H-LU is already factored (m_star_hlu_root != null).
        if (rc > 0 && !m_star_hlu_root) {
            std::vector<double> yv(ns), hy(ns), ey(ns);
            unsigned int rng = 2246822519u;
            for (int k=0;k<ns;k++){ rng=rng*1664525u+1013904223u; yv[k]=(double)(rng>>9)/8388608.0-0.5; }
            StarHMatVec(yv, hy);                                  // A_SS_Hmat * y
            applyT(yv, ts); MatVec(ts, tAs); applyTt(tAs, ey);    // exact A_SS * y
            double en=0.0, dn=0.0;
            for (int k=0;k<ns;k++){ en+=ey[k]*ey[k]; double d=hy[k]-ey[k]; dn+=d*d; }
            g_ass_dbg[1] = (en>0.0) ? std::sqrt(dn/en) : 0.0;
        }
    }

    // Preconditioner M^-1 for the reduced system T^T A T:
    //  (a) BLOCK-JACOBI SANDWICH  M^-1 = T^T M_bj^-1 T  (when blockInverse/blockOffsets
    //      supplied). EMPIRICALLY MUCH WORSE than (b) -- linear 50->902, nonlinear
    //      3023->49141 iters -- because the FULL-A block-diag does NOT match the reduced
    //      operator T^T A T (see rad_relaxation_methods.cpp ~2972). The caller never
    //      passes block data; this code path is left for completeness/experimentation.
    //  (b) CONSERVATIVE COLUMN-SCALE Jacobi  M^-1 = 1/d_k, d_k = sum c_i^2 |A_{i,i}|.
    //      Robust to sign cancellations between 1/chi diagonal and -N cross-terms.
    //      (Two stronger variants were prototyped and removed: B6 -- exact bilinear
    //       d_k = sum_{i,j} c_i c_j A_{i,j} -- diverged on C-type mu_r=1e5 when |d|~0
    //       was inverted. B7 -- true reduced block-Jacobi grouped by element-owner with
    //       LAPACK-pivoted block inversion -- slowed convergence on every working case
    //       (cube 1.77x->1.05x, antisym 49->81) and did not crack the C-type cap. The
    //       per-element block does not align with the natural block structure of T^T A T,
    //       which mixes element DOFs through the boundary/sym-internal/divergence star
    //       columns.  Both rolled back 2026-05-30.)
    bool useBlock = (!blockOffsets.empty() && !blockInverse.empty() && m_interaction);
    std::vector<double> dinv(ns, 1.0);
    if (!useBlock) {
        bool haveN = (m_diag_cached && (int)m_diag_N.size() == n);
        std::vector<double> diagA(n);
        for (int i = 0; i < n; i++){
            double Nii = haveN ? m_diag_N[i] : GetInteractionMatrixElement(i, i);
            double ic  = (i < (int)m_inv_chi.size()) ? m_inv_chi[i] : 0.0;
            diagA[i] = -Nii + ic;
        }
        for (int k = 0; k < ns; k++){
            double d = 0.0;
            for (int j = m_star_offsets[k]; j < m_star_offsets[k+1]; j++)
                d += m_star_coeffs[j]*m_star_coeffs[j]*std::fabs(diagA[m_star_dofs[j]]);
            dinv[k] = (d > 1e-300) ? 1.0/d : 1.0;
        }
    }
    // (d) K-DENSE preconditioner state. Populated AFTER any SetDeflation (so the
    // reduced operator T^T A T includes the antisym slab deflation alpha L L^T).
    // Build path runs after the deflation setup below; applyMinv prefers K-dense
    // when ready, else falls back to (b) col-scale Jacobi.
    std::vector<double>     kdense_lu;     // ns x ns col-major, overwritten with LU
    std::vector<lapack_int> kdense_ipiv;
    int                     kdense_n = 0;
    bool                    kdense_ready = false;

    // (e) FEW-ELEMENT BLOCK-JACOBI reduced preconditioner (Ida's S-graph block,
    // 2026-06-03).  Toggle via env RADIA_LS_BLOCK_CELLS:
    //   > 0 : cells-per-element-edge for the spatial-cell blocks (3 = 3^3-element
    //         clusters).  The cube study (C:/temp/test_ass_*.py) showed PER-element
    //         blocks FAIL (1/chi-vs--N sign cancellation -> small blocks near-singular)
    //         but FEW-element blocks WORK (2-3x over col-scale Jacobi, near
    //         mesh-independent).  Scalable O(N) -- the alternative to the K-dense LU,
    //         which caps ~15^3 / 8 GB.  Each block = A_SS restricted to a cell's star
    //         columns (StarSSEntry = S_i^T A S_j, exact for no-IMA / symmetric IMA),
    //         LAPACK LU.
    //   = -1 : force the (b) col-scale Jacobi baseline (skip K-dense) -- for comparison.
    //   0/unset : default (K-dense LU if it fits, else col-scale Jacobi).
    int blockCells = 0;
    if (const char* ev = std::getenv("RADIA_LS_BLOCK_CELLS")) blockCells = std::atoi(ev);
    std::vector<std::vector<int>>        bjCols;
    std::vector<std::vector<double>>     bjLU;
    std::vector<std::vector<lapack_int>> bjIpiv;
    bool bjReady = false;

    std::vector<double> mf(n), mg(n);
    auto applyMinv = [&](const std::vector<double>& rstar, std::vector<double>& zstar){
        if (kdense_ready) {
            // M^-1 r = (T^T A T)^-1 r  via LAPACK dgetrs (forward + backward subst).
            std::copy(rstar.begin(), rstar.end(), zstar.begin());
            LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', kdense_n, 1,
                           kdense_lu.data(), kdense_n, kdense_ipiv.data(),
                           zstar.data(), kdense_n);
            return;
        }
        if (bjReady) {
            // (e) few-element block-Jacobi: per spatial cell, dgetrs on the cell's
            // reduced sub-vector (an empty LU = a singular cell -> identity).
            for (size_t c = 0; c < bjCols.size(); c++) {
                const std::vector<int>& cs = bjCols[c]; int m = (int)cs.size();
                if (bjLU[c].empty()) { for (int i=0;i<m;i++) zstar[cs[i]] = rstar[cs[i]]; continue; }
                std::vector<double> rb(m);
                for (int i=0;i<m;i++) rb[i] = rstar[cs[i]];
                LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', m, 1, bjLU[c].data(),
                               m, bjIpiv[c].data(), rb.data(), m);
                for (int i=0;i<m;i++) zstar[cs[i]] = rb[i];
            }
            return;
        }
        if (!useBlock) { for (int k=0;k<ns;k++) zstar[k] = dinv[k]*rstar[k]; return; }
        applyT(rstar, mf);                                   // full = T r
        int ne = m_interaction->AmOfMainElem;
        std::fill(mg.begin(), mg.end(), 0.0);
        for (int e=0; e<ne; e++){
            int dof = m_interaction->GetElementDOF(e);
            int off = m_interaction->GetElementDOFOffset(e);
            int boff = blockOffsets[e];
            for (int i=0;i<dof;i++){ double acc=0.0;
                for (int j=0;j<dof;j++) acc += blockInverse[boff + i*dof + j]*mf[off+j];
                mg[off+i]=acc; }
        }
        applyTt(mg, zstar);                                  // z = T^T M_bj^-1 T r
    };
    auto dot = [](const std::vector<double>& a, const std::vector<double>& c){ double s=0; for(size_t i=0;i<a.size();i++) s+=a[i]*c[i]; return s; };
    auto nrm = [&](const std::vector<double>& a){ return std::sqrt(dot(a,a)); };

    std::vector<double> bred(ns); applyTt(b, bred);
    double bnorm = nrm(bred);
    if (bnorm < 1e-300){ sigma.assign(n, 0.0); return 0; }

    // ANTISYM IMA: deflate the few LOCAL plane-coupled null modes (the topological
    // star already gauges the loops, which stay null under IMA). MatVec then applies
    // A + alpha L L^T; the T^T(.)T sandwich lifts these modes out of the reduced
    // null space (loop-overlap part is annihilated by T^T). alpha ~ mean |A_ii| over
    // the modes' support so they land in the physical band. Robust (FEW + orthonormal).
    bool antisymDefl = false;
    if (HasAntisymmetricIMAPlane() && ComputePlaneSlabAddedModes() > 0){
        double asum=0.0; int acnt=0;
        bool haveN = (m_diag_cached && (int)m_diag_N.size()==n);
        for (int k=0;k<m_n_added;k++)
            for (int j=m_added_offsets[k]; j<m_added_offsets[k+1]; j++){
                int d=m_added_dofs[j];
                double Nii = haveN ? m_diag_N[d] : GetInteractionMatrixElement(d,d);
                double ic  = (d<(int)m_inv_chi.size()) ? m_inv_chi[d] : 0.0;
                asum += std::fabs(-Nii+ic); acnt++;
            }
        double alpha = (acnt>0) ? asum/(double)acnt : 1.0;
        SetDeflation(m_added_offsets, m_added_dofs, m_added_coeffs, alpha);
        antisymDefl = true;
    }

    // (d) K-DENSE preconditioner BUILD: explicit dense T^T A T + LAPACK LU.
    //   Cost: ns HACApK matvecs (build) + O(ns^3/3) factor (dgetrf). Memory ns^2*8 B.
    //   Built AFTER SetDeflation so MatVec applies (A + alpha L L^T) -> the LU is
    //   on the FULL reduced operator including antisym slab deflation.
    //   Skipped (-> falls back to (b) col-scale Jacobi via applyMinv) when:
    //     - ns^2 * 8 > KDENSE_BYTES_LIMIT (currently 8 GB; rules out ~20^3)
    //     - allocation throws (std::bad_alloc)
    //     - dgetrf returns info != 0 (singular even after deflation, should not happen
    //       for the symmetric / no-IMA / antisym-deflated reduced operator)
    //   Optimal at moderate scale: cube 4^3 (n=384) and C-type 6^3 (n=7200) factor in
    //   under a few seconds; iteration count drops to ~2 (BiCGSTAB direct-solve regime).
    if (blockCells == 0) {
        const size_t KDENSE_BYTES_LIMIT = (size_t)8 * 1024 * 1024 * 1024;   // 8 GB
        size_t kdense_bytes = (size_t)ns * (size_t)ns * sizeof(double);
        if (ns > 0 && kdense_bytes <= KDENSE_BYTES_LIMIT) {
            bool alloc_ok = true;
            try {
                kdense_lu.assign((size_t)ns * (size_t)ns, 0.0);
                kdense_ipiv.assign((size_t)ns, 0);
            } catch (std::bad_alloc&) { alloc_ok = false; kdense_lu.clear(); kdense_ipiv.clear(); }
            if (alloc_ok) {
                std::vector<double> kd_v(n, 0.0), kd_u(n, 0.0);
                for (int k = 0; k < ns; k++) {
                    // kd_v = T[:,k] (sparse): scatter star col k into full-space vector
                    for (int j = m_star_offsets[k]; j < m_star_offsets[k+1]; j++)
                        kd_v[m_star_dofs[j]] = m_star_coeffs[j];
                    MatVec(kd_v, kd_u);                              // kd_u = (A + alpha L L^T) kd_v
                    // (T^T kd_u) gives column k of T^T A T (col-major: kdense_lu[k2 + k*ns])
                    for (int k2 = 0; k2 < ns; k2++) {
                        double s = 0.0;
                        for (int j = m_star_offsets[k2]; j < m_star_offsets[k2+1]; j++)
                            s += m_star_coeffs[j] * kd_u[m_star_dofs[j]];
                        kdense_lu[(size_t)k2 + (size_t)k * (size_t)ns] = s;
                    }
                    // un-scatter (restore kd_v to zero for next col -- cheaper than full fill)
                    for (int j = m_star_offsets[k]; j < m_star_offsets[k+1]; j++)
                        kd_v[m_star_dofs[j]] = 0.0;
                }
                // LAPACK dense LU factorization
                int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, ns, ns,
                                          kdense_lu.data(), ns, kdense_ipiv.data());
                if (info == 0) { kdense_n = ns; kdense_ready = true; }
                else { kdense_lu.clear(); kdense_ipiv.clear(); }
            }
        }
    }

    // (e) FEW-ELEMENT BLOCK-JACOBI BUILD (when RADIA_LS_BLOCK_CELLS > 0): sort the
    // star columns by a FINE (1-element) spatial cell, then chunk into BALANCED
    // fixed-size blocks (spatially coherent AND balanced -- a graded mesh's dense
    // region does NOT collapse into one huge block, which it does with uniform
    // cells).  Each block = A_SS restricted to its columns via StarSSEntry (exact
    // for no-IMA / symmetric IMA), LAPACK LU.  Scalable O(N): fixed block size.
    if (blockCells > 0 && ns > 0) {
        double mn[3] = {1e30,1e30,1e30}, mx[3] = {-1e30,-1e30,-1e30};
        for (int k=0;k<ns;k++) for(int d=0;d<3;d++){ double vv=m_star_coords[3*(size_t)k+d];
            if(vv<mn[d])mn[d]=vv; if(vv>mx[d])mx[d]=vv; }
        int nelem = m_ndof/6; if (nelem < 1) nelem = 1;
        double ex=std::max(1e-30,mx[0]-mn[0]), ey=std::max(1e-30,mx[1]-mn[1]), ez=std::max(1e-30,mx[2]-mn[2]);
        double h = std::cbrt(ex*ey*ez/(double)nelem);     // typical element size
        std::vector<std::tuple<long,long,long,int>> keyed; keyed.reserve(ns);
        for (int k=0;k<ns;k++)
            keyed.emplace_back((long)std::floor((m_star_coords[3*(size_t)k+0]-mn[0])/h),
                               (long)std::floor((m_star_coords[3*(size_t)k+1]-mn[1])/h),
                               (long)std::floor((m_star_coords[3*(size_t)k+2]-mn[2])/h), k);
        std::sort(keyed.begin(), keyed.end());
        int Bsz = std::max(32, blockCells*blockCells*blockCells*4);   // target block size (cols)
        for (int start=0; start<ns; start+=Bsz){
            int m = std::min(Bsz, ns-start);
            std::vector<int> cols(m);
            for (int i=0;i<m;i++) cols[i]=std::get<3>(keyed[start+i]);
            std::vector<double> Bm((size_t)m*(size_t)m);
            for (int j=0;j<m;j++) for(int i=0;i<m;i++)
                Bm[(size_t)i + (size_t)j*(size_t)m] = StarSSEntry(cols[i], cols[j]);  // exact A_SS sub-block
            std::vector<lapack_int> ip(m);
            int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, m, Bm.data(), m, ip.data());
            bjCols.push_back(cols);
            if (info==0){ bjLU.push_back(std::move(Bm)); bjIpiv.push_back(std::move(ip)); }
            else { bjLU.push_back({}); bjIpiv.push_back({}); }   // singular block -> identity
        }
        bjReady = true;
    }

    // reducedSolve: right-preconditioned BiCGSTAB on the reduced system A_SS y = rhs
    // (operator = redMatVec = the matrix-free T^T A T sandwich, preconditioner =
    // applyMinv).  Wrapped in a lambda so the keep-loops star-correction below uses
    // the SAME (converged) solve instead of a single applyMinv application -- critical
    // for the few-element block-J, where one applyMinv is a WEAK preconditioner step,
    // NOT the exact A_SS^-1 the keep-loops refinement assumes (a single application
    // there amplifies the residual -> divergence).  For K-dense (applyMinv = exact
    // A_SS^-1) this converges in 1 iter, so the production path is unchanged.
    auto reducedSolve = [&](const std::vector<double>& rhs, std::vector<double>& yout)->int{
        double bn = nrm(rhs);
        if (bn < 1e-300){ yout.assign(ns, 0.0); return 0; }
        std::vector<double> y(ns,0.0), r(rhs), rhat(rhs), v(ns,0.0), p(ns,0.0),
                            s(ns), t(ns), phat(ns), shat(ns);
        double rho=1.0, alpha=1.0, omega=1.0;
        int it = 0;
        for (; it < max_iter; it++){
            double rho_new = dot(rhat, r);
            if (std::fabs(rho_new) < 1e-300) break;
            double beta = (rho_new/rho)*(alpha/omega);
            for (int i=0;i<ns;i++) p[i] = r[i] + beta*(p[i]-omega*v[i]);
            applyMinv(p, phat);
            redMatVec(phat, v);
            double rv = dot(rhat, v); if (std::fabs(rv) < 1e-300) break;
            alpha = rho_new/rv;
            for (int i=0;i<ns;i++) s[i] = r[i] - alpha*v[i];
            if (nrm(s)/bn < tol){ for(int i=0;i<ns;i++) y[i]+=alpha*phat[i]; it++; break; }
            applyMinv(s, shat);
            redMatVec(shat, t);
            double tt = dot(t,t); if (tt < 1e-300) break;
            omega = dot(t,s)/tt;
            for (int i=0;i<ns;i++){ y[i] += alpha*phat[i] + omega*shat[i]; r[i] = s[i] - omega*t[i]; }
            if (nrm(r)/bn < tol){ it++; break; }
            rho = rho_new;
        }
        yout = y; return it;
    };
    std::vector<double> y(ns, 0.0);
    int it = reducedSolve(bred, y);
    if (antisymDefl) SetDeflation({}, {}, {}, 0.0);   // clear: deflation was per-solve

    // A_SS-as-H-matrix MILESTONE 2 (mode 2): solve the star block with the A_SS
    // H-matrix factored by H-ILU as the PRECONDITIONER (Ida's scalable path).  The
    // dense K-dense LU caps at ~15^3 / 8 GB; H-ILU (low-rank-truncated H-LU of the
    // A_SS H-matrix) is the only >15^3 option.  THE LAW: A_SS's ill-conditioning is
    // the long-range 1/r demag coupling (non-local) -- a scalar threshold ILU(eps)
    // on entries fails (entries are not individually small), but H-ILU thresholds the
    // block SINGULAR VALUES (rank), which DO decay for the smooth 1/r kernel.  Steps:
    //   1. factor the A_SS H-matrix in place -> H-ILU factors (trunc tol =
    //      g_star_hlu_trunc_tol; larger = more aggressive truncation = cheaper),
    //   2. right-preconditioned BiCGSTAB with operator = redMatVec (the EXACT T^T A T
    //      sandwich, so the solution is field-exact, NOT the ACA-approx StarHMatVec)
    //      and preconditioner = StarHLUApply (the H-ILU factors).
    // The H-LU factor CONVERTS the leafmtxp leaves to internal layout, so StarHMatVec
    // is invalid afterwards -- that is why the operator is redMatVec, not StarHMatVec.
    // FreeStarHMatrix frees the H-LU root then the (converted) leafmtxp in order.
    // Diagnostics: [4]=rel err vs K-dense y (field-exactness), [5]=BiCGSTAB iters
    // (mu_r headline; -1 = factor failed), [6]=factor time, [7]=solve time,
    // [8]=preconditioner residual ||A_SS (M^-1 b) - b||/||b||, [9]=trunc tol.
    // Measurement side-channel: the production solution still uses the K-dense y
    // (when available, <=15^3); >15^3 (no K-dense) is the H-ILU's production regime.
    if (g_loopstar_hmat_mode >= 2 && m_star_hmat_built && ns > 0) {
        g_ass_dbg[9] = g_star_hlu_trunc_tol;
        // (1) H-ILU factorization of the A_SS H-matrix (factor once; reuse on re-entry)
        cHACApK_hlu_set_trunc_tol(g_star_hlu_trunc_tol);
        if (!m_star_hlu_root) {
            auto tf0 = std::chrono::high_resolution_clock::now();
            m_star_hlu_root = cHACApK_hlu_factor_leafmtxp(
                m_star_leafmtxp, m_star_control, /*nffc=*/1);
            auto tf1 = std::chrono::high_resolution_clock::now();
            g_ass_dbg[6] = std::chrono::duration<double>(tf1 - tf0).count();
        }
        if (m_star_hlu_root) {
            // (8) preconditioner quality probe: z = M^-1 b, resid = ||A_SS z - b||/||b||
            {
                std::vector<double> zp(ns), azp(ns);
                StarHLUApply(bred, zp);
                redMatVec(zp, azp);
                double en=0.0, dn=0.0;
                for (int i=0;i<ns;i++){ en+=bred[i]*bred[i]; double d=azp[i]-bred[i]; dn+=d*d; }
                g_ass_dbg[8] = (en>0.0)? std::sqrt(dn/en):0.0;
            }
            // (2) right-preconditioned BiCGSTAB: operator=redMatVec (exact A_SS),
            //     preconditioner=StarHLUApply (H-ILU factors).
            auto ts0 = std::chrono::high_resolution_clock::now();
            std::vector<double> yh(ns,0.0), rr(bred), rh(bred), vv(ns,0.0), pp(ns,0.0),
                                ss(ns), tt2(ns), ph(ns), sh(ns);
            double rho2=1.0, al2=1.0, om2=1.0; int it2=0;
            for (; it2 < max_iter; it2++){
                double rn2 = dot(rh, rr);
                if (std::fabs(rn2) < 1e-300) break;
                double beta2 = (rn2/rho2)*(al2/om2);
                for(int i=0;i<ns;i++) pp[i]=rr[i]+beta2*(pp[i]-om2*vv[i]);
                StarHLUApply(pp, ph); redMatVec(ph, vv);      // M^-1 then exact A_SS
                double rv2=dot(rh,vv); if(std::fabs(rv2)<1e-300) break;
                al2=rn2/rv2;
                for(int i=0;i<ns;i++) ss[i]=rr[i]-al2*vv[i];
                if(nrm(ss)/bnorm<tol){ for(int i=0;i<ns;i++) yh[i]+=al2*ph[i]; it2++; break; }
                StarHLUApply(ss, sh); redMatVec(sh, tt2);
                double ttd=dot(tt2,tt2); if(ttd<1e-300) break;
                om2=dot(tt2,ss)/ttd;
                for(int i=0;i<ns;i++){ yh[i]+=al2*ph[i]+om2*sh[i]; rr[i]=ss[i]-om2*tt2[i]; }
                if(nrm(rr)/bnorm<tol){ it2++; break; }
                rho2=rn2;
            }
            auto ts1 = std::chrono::high_resolution_clock::now();
            g_ass_dbg[7] = std::chrono::duration<double>(ts1 - ts0).count();
            double en=0.0, dn=0.0;
            for(int i=0;i<ns;i++){ en+=y[i]*y[i]; double d=yh[i]-y[i]; dn+=d*d; }
            g_ass_dbg[4] = (en>0.0)? std::sqrt(dn/en):0.0;
            g_ass_dbg[5] = (double)it2;
        } else {
            g_ass_dbg[5] = -1.0;   // H-LU factorization failed
        }
    }

    applyT(y, sigma);                                 // sigma = sigma_S (loop-free star part)

    // KEEP-LOOPS via block Gauss-Seidel (star <-> loop). The reduced star solve
    // alone (sigma = S y_S) (i) drops the loop content and (ii) -- because the
    // sparse star S and the loop basis L are NOT exactly mutually orthogonal
    // (S^T L != 0, the per-element divergence star columns overlap the cycles) --
    // leaves the star coordinates slightly off (A_SL != 0 breaks the strict
    // block-triangular structure). Both show up as a ~0.5% external-field error on
    // the C-type at low/moderate mu_r. We fix BOTH by a few Gauss-Seidel sweeps on
    // the full system A sigma = b, alternating:
    //   (i)  STAR correction:  A_SS d_S = S^T r  via the same K-dense reduced solve
    //        (applyMinv); sigma += S d_S,
    //   (ii) LOOP correction:  A_LL y_L = L^T r  with A_LL = L^T diag(inv_chi) L
    //        (= (1/chi) L^T L for uniform chi, SPD), solved by CG on the matrix-free
    //        operator z -> L^T( inv_chi . (L z) ); sigma += L y_L.
    // Each sweep recomputes the true residual r = b - A sigma (one HACApK MatVec),
    // so convergence is to the DIRECT solution (field-exact, loops kept) regardless
    // of basis orthogonality. The coupling is weak (~0.5%) so a handful of sweeps
    // reach tol. This is iterative refinement with the star+loop blocks as a
    // (cheap, mu_r-robust) preconditioner -- far fewer MatVecs than plain BiCGSTAB.
    if (m_n_loop > 0) {
        const int nl = m_n_loop;
        auto applyL = [&](const std::vector<double>& yL, std::vector<double>& full){
            std::fill(full.begin(), full.end(), 0.0);
            for (int p=0; p<nl; p++){ double c=yL[p];
                for (int k=m_loop_offsets[p]; k<m_loop_offsets[p+1]; k++)
                    full[m_loop_dofs[k]] += m_loop_coeffs[k]*c; }
        };
        auto applyLt = [&](const std::vector<double>& full, std::vector<double>& yL){
            for (int p=0; p<nl; p++){ double acc=0.0;
                for (int k=m_loop_offsets[p]; k<m_loop_offsets[p+1]; k++)
                    acc += m_loop_coeffs[k]*full[m_loop_dofs[k]];
                yL[p]=acc; }
        };
        std::vector<double> lz(n);
        auto gL = [&](const std::vector<double>& z, std::vector<double>& out){
            applyL(z, lz);                                  // lz = L z
            for (int i=0;i<n;i++)
                lz[i] *= (i<(int)m_inv_chi.size() ? m_inv_chi[i] : 0.0);
            applyLt(lz, out);                               // out = L^T (inv_chi . L z)
        };
        std::vector<double> As(n), r(n), tmp(n), rstar(ns), dyS(ns);
        std::vector<double> rhsL(nl), yL(nl), rk(nl), pk(nl), Ap(nl);
        double bn = 0.0; for (int i=0;i<n;i++) bn += b[i]*b[i]; bn = std::sqrt(bn) + 1e-300;
        const int MAXSWEEP = 50;
        double res0 = 0.0, resf = 0.0; int sweeps = 0, cg_tot = 0;
        for (int sw=0; sw<MAXSWEEP; sw++){
            MatVec(sigma, As); for (int i=0;i<n;i++) r[i] = b[i]-As[i];   // full residual
            double rn=0.0; for (int i=0;i<n;i++) rn += r[i]*r[i]; rn = std::sqrt(rn);
            if (sw==0) res0 = rn;
            resf = rn; sweeps = sw;
            if (rn <= tol*bn) break;
            // (i) STAR correction d_S = A_SS^{-1} S^T r ; sigma += S d_S.
            // Use the CONVERGED reducedSolve (not a single applyMinv): for block-J,
            // applyMinv is a weak step and a single application would leave d_S far
            // from A_SS^{-1} S^T r, so the refinement amplifies (-> 1e88).  For
            // K-dense (applyMinv exact) reducedSolve converges in 1 iter (unchanged).
            applyTt(r, rstar); reducedSolve(rstar, dyS); applyT(dyS, tmp);
            for (int i=0;i<n;i++) sigma[i] += tmp[i];
            // (ii) LOOP correction y_L = A_LL^{-1} L^T r2 ; sigma += L y_L (CG)
            MatVec(sigma, As); for (int i=0;i<n;i++) r[i] = b[i]-As[i];
            applyLt(r, rhsL);
            std::fill(yL.begin(), yL.end(), 0.0); rk = rhsL; pk = rhsL;
            double rsold=0.0; for(int i=0;i<nl;i++) rsold += rk[i]*rk[i];
            double rhn = std::sqrt(rsold);
            if (rhn > 1e-300){
                int cgmax = 2*nl + 50;
                for (int cgit=0; cgit<cgmax; cgit++){
                    cg_tot++;
                    gL(pk, Ap);
                    double pAp=0.0; for(int i=0;i<nl;i++) pAp += pk[i]*Ap[i];
                    if (std::fabs(pAp) < 1e-300) break;
                    double a = rsold/pAp;
                    for(int i=0;i<nl;i++){ yL[i]+=a*pk[i]; rk[i]-=a*Ap[i]; }
                    double rsnew=0.0; for(int i=0;i<nl;i++) rsnew += rk[i]*rk[i];
                    if (std::sqrt(rsnew) <= tol*rhn) break;
                    double bcg = rsnew/rsold;
                    for(int i=0;i<nl;i++) pk[i]=rk[i]+bcg*pk[i];
                    rsold = rsnew;
                }
                applyL(yL, tmp); for (int i=0;i<n;i++) sigma[i] += tmp[i];
            }
        }
        g_kl_dbg[0]=(double)nl; g_kl_dbg[1]=res0; g_kl_dbg[2]=resf/bn;
        g_kl_dbg[3]=(double)sweeps; g_kl_dbg[4]=(double)cg_tot; g_kl_dbg[5]=resf;
    }
    return it;
}

//=========================================================================
// HELMHOLTZ-HODGE LOOP REMOVAL: subtract the non-physical loop (circulating
// surface-charge) component from a solved sigma.  This is NOT a solver/convergence
// trick (plain BiCGSTAB already converges); it gives a PHYSICAL answer with zero
// loop content.  The loop space is span(L), L = the topological cycle basis
// (m_loop_*, = ker(N), the circulating charges).  Hodge projection off span(L):
//   c solves (L^T L) c = L^T sigma   (CG; L^T L is the loop GRAM matrix)
//   sigma <- sigma - L c
// L^T L is geometric, sparse, mu_r-INDEPENDENT and WELL-CONDITIONED (cond ~ O(1),
// the cycle basis is near-orthonormal) -- unlike A_SS (cond ~ mu_r) -- so the CG
// converges in a handful of iterations: NO significant slowdown.  Because N L = 0,
// removing the loop part does NOT change N*sigma (the on-element field identity is
// preserved); only the non-physical circulation is removed.  Scalable: uses only
// the O(N) cycle CSR (NOT the +3 non-topological near-null modes, which grow with N
// and are NOT loops; NOT A_SL; NOT A).  Returns CG iters, or -1 if no loops.
// g_loopproj_dbg: [0]=n_loop, [1]=CG iters, [2]=||L^T sigma|| before,
// [3]=||L^T sigma|| after (loop residual), [4]=loop fraction ||L c||/||sigma||.
//=========================================================================
static double g_loopproj_dbg[5] = {0,0,0,0,0};
extern "C" void RadGetLoopProjStats(double* out /* len 5 */) {
    if (!out) return; for (int i=0;i<5;i++) out[i]=g_loopproj_dbg[i];
}

int RadHACApKMSCManager::ProjectOutLoops(std::vector<double>& sigma, double tol, int max_iter) {
    const int n = m_ndof;
    if ((int)sigma.size() != n) return -1;
    // Build/cache the topological loop cycle CSR (same pattern as SolveLoopStar).
    if (!m_loop_built) {
        BuildLoopBasis(1.0);
        m_loop_offsets = m_defl_offsets;
        m_loop_dofs    = m_defl_dofs;
        m_loop_coeffs  = m_defl_signs;
        SetDeflation({}, {}, {}, 0.0);   // clear: keep MatVec the clean A
        m_n_loop = m_loop_offsets.empty() ? 0 : (int)m_loop_offsets.size() - 1;
        m_loop_built = true;
    }
    const int nl = m_n_loop;
    if (nl <= 0) { g_loopproj_dbg[0]=0; return -1; }

    auto applyL = [&](const std::vector<double>& yL, std::vector<double>& full){
        std::fill(full.begin(), full.end(), 0.0);
        for (int p=0;p<nl;p++){ double c=yL[p];
            for (int k=m_loop_offsets[p]; k<m_loop_offsets[p+1]; k++)
                full[m_loop_dofs[k]] += m_loop_coeffs[k]*c; }
    };
    auto applyLt = [&](const std::vector<double>& full, std::vector<double>& yL){
        for (int p=0;p<nl;p++){ double acc=0.0;
            for (int k=m_loop_offsets[p]; k<m_loop_offsets[p+1]; k++)
                acc += m_loop_coeffs[k]*full[m_loop_dofs[k]];
            yL[p]=acc; }
    };
    auto dot=[](const std::vector<double>&a,const std::vector<double>&b){
        double s=0; for(size_t i=0;i<a.size();i++) s+=a[i]*b[i]; return s; };

    // rhs = L^T sigma ; solve (L^T L) c = rhs by CG on the matrix-free Gram operator
    // G z = L^T (L z).  G is SPD + well-conditioned (cycle basis near-orthonormal).
    std::vector<double> rhs(nl), c(nl,0.0), r(nl), p(nl), Gp(nl), Lz(n);
    applyLt(sigma, rhs);
    double loop_before = std::sqrt(dot(rhs,rhs));
    g_loopproj_dbg[0] = (double)nl;
    g_loopproj_dbg[2] = loop_before;
    if (loop_before < 1e-300) { g_loopproj_dbg[1]=0; g_loopproj_dbg[3]=0; g_loopproj_dbg[4]=0; return 0; }
    auto Gapply = [&](const std::vector<double>& z, std::vector<double>& out){
        applyL(z, Lz); applyLt(Lz, out);
    };
    r = rhs; p = rhs;     // c0 = 0 -> r0 = rhs
    double rs = dot(r,r); double rs0 = rs;
    int it = 0;
    for (; it < max_iter; it++){
        Gapply(p, Gp);
        double pGp = dot(p,Gp);
        if (std::fabs(pGp) < 1e-300) break;
        double a = rs/pGp;
        for (int i=0;i<nl;i++){ c[i]+=a*p[i]; r[i]-=a*Gp[i]; }
        double rs_new = dot(r,r);
        if (std::sqrt(rs_new) <= tol*std::sqrt(rs0)) { it++; break; }
        double beta = rs_new/rs;
        for (int i=0;i<nl;i++) p[i] = r[i] + beta*p[i];
        rs = rs_new;
    }
    // sigma <- sigma - L c
    applyL(c, Lz);
    double lc_norm = std::sqrt(dot(Lz,Lz));
    double sig_norm = std::sqrt(dot(sigma,sigma)) + 1e-300;
    for (int i=0;i<n;i++) sigma[i] -= Lz[i];
    // loop residual after
    std::vector<double> rhs_after(nl); applyLt(sigma, rhs_after);
    g_loopproj_dbg[1] = (double)it;
    g_loopproj_dbg[3] = std::sqrt(dot(rhs_after,rhs_after));
    g_loopproj_dbg[4] = lc_norm / sig_norm;
    return it;
}

//=========================================================================
// MSC system-matrix convention: A(i, j) = -N(i, j) + delta_ij / chi_i
// where N is the physical demagnetization tensor entry returned by
// GetInteractionMatrixElement (MSC/MMM convention: +N).
//=========================================================================

double RadHACApKMSCManager::ComputeSystemEntry(int dof_i, int dof_j) const {
    double N_val = GetInteractionMatrixElement(dof_i, dof_j);
    double A_val = -N_val;
    if (dof_i == dof_j && dof_i < (int)m_inv_chi.size()) {
        A_val += m_inv_chi[dof_i];
    }
    return A_val;
}

void RadHACApKBase::FreeResources() {
    // Clear global callback state first (prevents stale state in next solve)
    RadHACApKCallback::ClearGlobalState();

    if (m_leafmtxp || m_control) {
        HACApK_free_hmatrix_wrapper(m_leafmtxp, m_control);
    }
    if (m_leafmtxp) {
        HACApK_free_leafmtxp(m_leafmtxp);
        m_leafmtxp = nullptr;
    }
    if (m_control) {
        HACApK_free_lcontrol(m_control);
        m_control = nullptr;
    }

    // Reset HACApK global state (persistent matvec buffers, lod)
    HACApK_reset_global_state();

    m_valid = false;
}

//=========================================================================
// STANDALONE loop removal for the NON-HACApK solvers (method 0 = LU, method 1 =
// dense BiCGSTAB). Those solver classes do not own a RadHACApKMSCManager, but the
// Helmholtz-Hodge projection needs ONLY the cycle basis (BuildLoopBasis) + sparse
// matrix-free CG on L^T L -- NO H-matrix, NO MatVec. So build a TEMPORARY manager
// bound to the same interaction, run ExtractCoordinates() (cheap: just dof offsets)
// + ProjectOutLoops on the element sigma gathered from the objects, write back.
// This makes "loop is non-physical" consistent across ALL solver paths, not only
// HACApK -- so method 0/1/2 all return the same loop-free sigma. Returns CG iters
// (-1 if no loops / not applicable). Reads the loop-projection diagnostics into the
// same g_loopproj_dbg as the HACApK path.
int RadHACApKMSCManager::ProjectOutLoopsStandalone(radTInteraction* interaction,
                                                   double tol, int max_iter)
{
    if (!interaction || interaction->GetAmOfMainElem() <= 0) return -1;
    RadHACApKMSCManager mgr(interaction);
    mgr.ExtractCoordinates();                  // populate m_n_elem, m_dof_offset, m_ndof
    int n = mgr.GetNDOF();
    if (n <= 0) return -1;

    // Gather sigma from element objects into a flat vector (3DOF Magn / 6DOF Sigma[]).
    std::vector<double> sigma(n, 0.0);
    int ne = interaction->GetAmOfMainElem();
    for (int e = 0; e < ne; e++) {
        int dof = interaction->GetElementDOF(e);
        int off = interaction->GetElementDOFOffset(e);
        radTg3dRelax* g = interaction->GetElement(e);
        if (!g) continue;
        if (dof == 3) {
            sigma[off+0] = g->Magn.x; sigma[off+1] = g->Magn.y; sigma[off+2] = g->Magn.z;
        } else if (dof >= 5) {
            radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g);
            if (poly && poly->Use6DOF_MSC)
                for (int k = 0; k < dof; k++) sigma[off+k] = poly->Sigma[k];
        }
    }

    int it = mgr.ProjectOutLoops(sigma, tol, max_iter);
    if (it < 0) return -1;

    // Scatter the loop-free sigma back to the element objects.
    for (int e = 0; e < ne; e++) {
        int dof = interaction->GetElementDOF(e);
        int off = interaction->GetElementDOFOffset(e);
        radTg3dRelax* g = interaction->GetElement(e);
        if (!g) continue;
        if (dof == 3) {
            g->Magn.x = sigma[off+0]; g->Magn.y = sigma[off+1]; g->Magn.z = sigma[off+2];
        } else if (dof >= 5) {
            radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g);
            if (poly && poly->Use6DOF_MSC)
                for (int k = 0; k < dof; k++) poly->Sigma[k] = sigma[off+k];
        }
    }
    return it;
}

void RadHACApKMSCManager::ExtractCoordinates() {
    // Extract element center coordinates for clustering
    // Supports both 3DOF tetrahedra and 6DOF MSC hexahedra
    if (!m_interaction) return;

    m_n_elem = m_interaction->AmOfMainElem;
    m_coordinates.resize(m_n_elem * 3);
    m_dof_offset.resize(m_n_elem + 1);

    // Check if using variable DOF (6DOF hexahedra) or uniform 3DOF (tetrahedra)
    bool has_variable_dof = m_interaction->HasVariableDOF();

    int total_dof = 0;
    int n_3dof = 0;
    int n_6dof = 0;

    for (int i = 0; i < m_n_elem; i++) {
        m_dof_offset[i] = total_dof;
        int elem_dof = m_interaction->GetElementDOF(i);
        total_dof += elem_dof;

        if (elem_dof == 3) {
            n_3dof++;
        } else if (elem_dof >= 5) {
            n_6dof++;  // Count all MSC elements (5-DOF wedges and 6-DOF hexahedra)
        } else {
            std::cerr << "[HACApK] Error: Element " << i << " has " << elem_dof
                      << " DOF, expected 3, 5, or 6" << std::endl;
            m_ndof = 0;
            m_nffc = 0;
            m_is_6dof = false;
            return;
        }
    }

    m_dof_offset[m_n_elem] = total_dof;
    m_ndof = total_dof;

    // Support mixed DOF elements (hex + tetra + wedge)
    if (n_3dof > 0 && n_6dof > 0) {
        // Mixed mode: variable DOF per element (tetra + MSC)
        m_nffc = 0;  // Indicates variable DOF
        m_is_6dof = false;
        m_is_mixed_dof = true;
#ifdef HACAPK_RADIA_LOGGING
        std::cout << "[HACApK] Mixed DOF mode: " << n_3dof << " tetrahedra (3DOF) + "
                  << n_6dof << " MSC elements (5/6DOF), total " << total_dof << " DOF" << std::endl;
#endif
    } else if (n_6dof > 0) {
        // Check element types: pure hex, pure wedge, or mixed
        bool allHex = true, allWedge = true;
        for (int i = 0; i < m_n_elem; i++) {
            int dof = m_interaction->GetElementDOF(i);
            if (dof != 6) allHex = false;
            if (dof != 5) allWedge = false;
        }
        if (allHex) {
            m_nffc = 6;
            m_is_6dof = true;   // Pure hex (fast path with Compute6x6BlockFast)
            m_is_5dof = false;
            m_is_mixed_dof = false;
        } else if (allWedge) {
            m_nffc = 5;
            m_is_6dof = false;
            m_is_5dof = true;   // Pure wedge (fast path with Compute5x5BlockFast)
            m_is_mixed_dof = false;
        } else {
            // Mixed hex+wedge: use variable DOF mode (flat matrix)
            m_nffc = 0;
            m_is_6dof = false;
            m_is_5dof = false;
            m_is_mixed_dof = true;
        }
    } else {
        m_nffc = 3;
        m_is_6dof = false;
        m_is_5dof = false;
        m_is_mixed_dof = false;
    }

    // Get element centers from g3dRelaxPtrVect
    for (int i = 0; i < m_n_elem; i++) {
        radTg3dRelax* elem = m_interaction->g3dRelaxPtrVect[i];
        if (elem) {
            TVector3d center = elem->CentrPoint;
            m_coordinates[i * 3 + 0] = center.x;
            m_coordinates[i * 3 + 1] = center.y;
            m_coordinates[i * 3 + 2] = center.z;
        }
    }

    // Build O(1) DOF lookup table
    BuildDOFLookupTable();
}

//=========================================================================
// BuildDOFLookupTable: Create O(1) DOF-to-element lookup (ELF-style)
//=========================================================================

void RadHACApKMSCManager::BuildDOFLookupTable() {
    if (m_ndof == 0) return;

    m_dof_to_elem.resize(m_ndof);
    m_dof_to_local.resize(m_ndof);

    for (int e = 0; e < m_n_elem; e++) {
        int offset = m_dof_offset[e];
        int elem_dof = m_dof_offset[e + 1] - offset;
        for (int d = 0; d < elem_dof; d++) {
            m_dof_to_elem[offset + d] = e;
            m_dof_to_local[offset + d] = d;
        }
    }
}

//=========================================================================
// PrecomputeGeometry3DOF: ELF-style pre-computed geometry for 3DOF tetrahedra
// Extracts all tetrahedron vertices and face geometry into contiguous arrays
// This avoids calling B_comp() which has significant overhead during H-matrix build
//=========================================================================

void RadHACApKMSCManager::PrecomputeGeometry3DOF() {
    if (m_geometry_3dof_ready || m_n_elem == 0 || !m_interaction) return;

    // Allocate arrays for tetrahedra (4 triangular faces, 3 vertices each)
    m_tetra_centers.resize(m_n_elem * 3);
    m_tetra_face_vertices.resize(m_n_elem * 4 * 3 * 3);  // 4 faces, 3 vertices, 3 coords
    m_tetra_face_normals.resize(m_n_elem * 4 * 3);       // 4 faces, 3 coords (outward normals)
    m_tetra_face_areas.resize(m_n_elem * 4);              // 4 faces

    // Parallel per-element fill: each iteration writes to disjoint slices of
    // the pre-allocated flat arrays (m_tetra_centers, m_tetra_face_vertices,
    // m_tetra_face_normals, m_tetra_face_areas), so no synchronisation needed.
    ngcore::ParallelFor(ngcore::IntRange(m_n_elem), [&](size_t e_size) {
        int e = (int)e_size;
        radTg3dRelax* elem = m_interaction->g3dRelaxPtrVect[e];
        if (!elem) return;

        radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(elem);
        if (!poly || poly->AmOfFaces != 4) return;  // Skip non-tetrahedra

        // Store element center
        int cIdx = e * 3;
        m_tetra_centers[cIdx + 0] = poly->CentrPoint.x;
        m_tetra_centers[cIdx + 1] = poly->CentrPoint.y;
        m_tetra_centers[cIdx + 2] = poly->CentrPoint.z;

        // Store face data for each of the 4 triangular faces
        for (int f = 0; f < 4; f++) {
            const radTHandlePgnAndTrans& hpt = poly->VectHandlePgnAndTrans[f];
            radTPolygon* pgn = hpt.PgnHndl.rep;
            radTrans* tr = hpt.TransHndl.rep;

            // Get 3 vertices of this triangular face
            const radTVect2dVect& verts2d = pgn->EdgePointsVector;
            if (verts2d.size() < 3) continue;

            int fvIdx = (e * 4 + f) * 3 * 3;  // Starting index for this face's vertices
            TVector3d V[3];
            for (int v = 0; v < 3; v++) {
                V[v] = tr->TrPoint(TVector3d(verts2d[v].x, verts2d[v].y, pgn->CoordZ));
                m_tetra_face_vertices[fvIdx + v * 3 + 0] = V[v].x;
                m_tetra_face_vertices[fvIdx + v * 3 + 1] = V[v].y;
                m_tetra_face_vertices[fvIdx + v * 3 + 2] = V[v].z;
            }

            // Compute face normal (outward pointing)
            TVector3d e1 = {V[1].x - V[0].x, V[1].y - V[0].y, V[1].z - V[0].z};
            TVector3d e2 = {V[2].x - V[0].x, V[2].y - V[0].y, V[2].z - V[0].z};
            TVector3d n = {e1.y*e2.z - e1.z*e2.y, e1.z*e2.x - e1.x*e2.z, e1.x*e2.y - e1.y*e2.x};
            double nLen = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);

            // Face area = 0.5 * |cross product|
            m_tetra_face_areas[e * 4 + f] = 0.5 * nLen;

            // Normalize and check orientation (outward from centroid)
            if (nLen > 1e-20) {
                n.x /= nLen; n.y /= nLen; n.z /= nLen;

                // Face center
                TVector3d fc = {(V[0].x + V[1].x + V[2].x) / 3.0,
                                (V[0].y + V[1].y + V[2].y) / 3.0,
                                (V[0].z + V[1].z + V[2].z) / 3.0};
                // Vector from centroid to face center
                TVector3d toFace = {fc.x - poly->CentrPoint.x,
                                    fc.y - poly->CentrPoint.y,
                                    fc.z - poly->CentrPoint.z};
                // If normal points inward, flip it
                if (n.x*toFace.x + n.y*toFace.y + n.z*toFace.z < 0) {
                    n.x = -n.x; n.y = -n.y; n.z = -n.z;
                }
            }

            // Store normalized outward normal
            int fnIdx = (e * 4 + f) * 3;
            m_tetra_face_normals[fnIdx + 0] = n.x;
            m_tetra_face_normals[fnIdx + 1] = n.y;
            m_tetra_face_normals[fnIdx + 2] = n.z;
        }
    });

    m_geometry_3dof_ready = true;
}

//=========================================================================
// PrecomputeFlatInteractMatrix: Flatten InteractMatrix for 3DOF tetrahedra
// Converts 2D pointer array TMatrix3df** to contiguous double array
// This eliminates pointer chasing during matrix element access
//=========================================================================

void RadHACApKMSCManager::PrecomputeFlatInteractMatrix() {
    if (m_flat_N_ready || m_n_elem == 0 || !m_interaction) return;
    if (!m_interaction->InteractMatrix) {
        return;  // InteractMatrix not computed
    }

    // Allocate flat storage: n_elem * n_elem * 9 doubles (3x3 block per element pair)
    int64_t total_size = (int64_t)m_n_elem * m_n_elem * 9;
    m_flat_N_data.resize(total_size);

    // Copy from InteractMatrix[i][j] to flat array as +N (physical quantity).
    // ComputeEntry() handles the sign flip to -N for the system matrix.
    //
    // InteractMatrix[i][j] stores TMatrix3df where:
    //   Str0 = dH/dMx, Str1 = dH/dMy, Str2 = dH/dMz (COLUMN vectors)
    // We store row-major: N[row][col] = dH_row / dM_col
    int total_collapse = m_n_elem * m_n_elem;
    ngcore::ParallelFor(ngcore::IntRange(total_collapse), [&](size_t idx)
    {
        int i = (int)(idx / m_n_elem);
        int j = (int)(idx % m_n_elem);
            const TMatrix3df& M = m_interaction->InteractMatrix[i][j];
            int64_t base_idx = ((int64_t)i * m_n_elem + j) * 9;

            m_flat_N_data[base_idx + 0] = static_cast<double>(M.Str0.x);  // dHx/dMx
            m_flat_N_data[base_idx + 1] = static_cast<double>(M.Str1.x);  // dHx/dMy
            m_flat_N_data[base_idx + 2] = static_cast<double>(M.Str2.x);  // dHx/dMz
            m_flat_N_data[base_idx + 3] = static_cast<double>(M.Str0.y);  // dHy/dMx
            m_flat_N_data[base_idx + 4] = static_cast<double>(M.Str1.y);  // dHy/dMy
            m_flat_N_data[base_idx + 5] = static_cast<double>(M.Str2.y);  // dHy/dMz
            m_flat_N_data[base_idx + 6] = static_cast<double>(M.Str0.z);  // dHz/dMx
            m_flat_N_data[base_idx + 7] = static_cast<double>(M.Str1.z);  // dHz/dMy
            m_flat_N_data[base_idx + 8] = static_cast<double>(M.Str2.z);  // dHz/dMz
    });

    m_flat_N_ready = true;
}

//=========================================================================
// Compute6x6BlockFast: Delegates to radTInteraction for unified implementation
// 6DOF hexahedra use radTInteraction::Compute6x6BlockFast (shared with LU/BiCGSTAB)
//=========================================================================

void RadHACApKMSCManager::Compute6x6BlockFast(int elem_i, int elem_j, double* K_mat) const {
    // Use radTInteraction's unified implementation (shared with LU/BiCGSTAB)
    if (m_interaction && m_interaction->IsHexaTriDataReady()) {
        m_interaction->Compute6x6BlockFast(elem_i, elem_j, K_mat);
        return;
    }

    // Error: radTInteraction geometry not ready
    std::memset(K_mat, 0, 36 * sizeof(double));
}

bool RadHACApKBase::BuildHMatrix(const RadHACApKParams& params) {
    FreeResources();

    auto start_time = std::chrono::high_resolution_clock::now();

    ExtractCoordinates();

    if (m_n_elem == 0) {
        std::cerr << "[HACApK] Error: No elements" << std::endl;
        return false;
    }
    if (m_ndof == 0) {
        std::cerr << "[HACApK] Error: Invalid DOF configuration (ndof=0)" << std::endl;
        return false;
    }

    // Set global callback state (MUST happen before OnBeforeBuild / InitializeInvChi
    // so kernel-side hooks may register additional callback state if needed)
    RadHACApKCallback::SetCurrentManager(this);

    // FIX (2025-02-04): Increment generation counter to invalidate thread-local caches
    // so that matrix-element caches from previous solves are not reused.
    RadHACApKCallback::IncrementGeneration();

    // Kernel-specific precomputation (e.g. PrecomputeHexaGeometry for MSC hex)
    OnBeforeBuild();

    // Kernel-specific initial chi (for MSC: initial susceptibility from material state)
    InitializeInvChi();

    RadHACApKCallback::SetInvChi(m_inv_chi);

    // Allocate opaque structures
    m_leafmtxp = HACApK_alloc_leafmtxp();
    m_control = HACApK_alloc_lcontrol();

    if (!m_leafmtxp || !m_control) {
        std::cerr << "[HACApK] Error: Failed to allocate structures" << std::endl;
        return false;
    }

    // Build H-matrix via HACApK wrapper. Variable-DOF mode routes through the
    // dof_offset array; uniform-DOF mode uses GetUniformNFFC().
    int ndim = 3;

    auto t_hmatrix_start = std::chrono::high_resolution_clock::now();
    int result;

    if (IsVariableDOF()) {
        result = HACApK_build_hmatrix_varDOF_wrapper(
            m_leafmtxp,
            m_control,
            m_coordinates.data(),
            m_n_elem,
            m_dof_offset.data(),
            m_ndof,
            ndim,
            params.aca_eps,
            params.leaf_size,
            params.eta,
            params.print_level
        );
    } else {
        result = HACApK_build_hmatrix_wrapper(
            m_leafmtxp,
            m_control,
            m_coordinates.data(),
            m_n_elem,
            GetUniformNFFC(),
            ndim,
            params.aca_eps,
            params.leaf_size,
            params.eta,
            params.print_level
        );
    }
    auto t_hmatrix_end = std::chrono::high_resolution_clock::now();
    double t_hmatrix = std::chrono::duration<double>(t_hmatrix_end - t_hmatrix_start).count();
    if (params.print_level > 0) {
        std::cout << "[HACApK] HACApK_build_hmatrix_wrapper: " << t_hmatrix << " s" << std::endl;
    }

    if (result != 0) {
        std::cerr << "[HACApK] Error: H-matrix build failed with code " << result << std::endl;
        FreeResources();
        return false;
    }

    // Store permutation from control structure
    int* lod = HACApK_lcontrol_get_lod(m_control);
    if (lod) {
        m_permutation.resize(m_ndof);
        for (int i = 0; i < m_ndof; i++) {
            m_permutation[i] = lod[i + 1];
        }
    }

    // Get statistics from leafmtxp
    m_stats.n_dof = HACApK_leafmtxp_get_nd(m_leafmtxp);
    m_stats.n_leaves = HACApK_leafmtxp_get_nlf(m_leafmtxp);
    m_stats.n_lowrank = HACApK_leafmtxp_get_nlfkt(m_leafmtxp);
    m_stats.n_dense = m_stats.n_leaves - m_stats.n_lowrank;
    m_stats.max_rank = HACApK_leafmtxp_get_ktmax(m_leafmtxp);

    int64_t hmat_bytes = 0;
    int64_t dense_bytes = 0;
    HACApK_get_memory_stats(m_leafmtxp, &hmat_bytes, &dense_bytes);

    m_stats.memory_mb = (double)hmat_bytes / (1024.0 * 1024.0);
    m_stats.dense_memory_mb = (double)dense_bytes / (1024.0 * 1024.0);
    m_stats.compression = (dense_bytes > 0) ?
        (double)hmat_bytes / (double)dense_bytes : 1.0;

    auto end_time = std::chrono::high_resolution_clock::now();
    m_stats.build_time = std::chrono::duration<double>(end_time - start_time).count();

    if (params.print_level > 0) {
        std::cout << "[HACApK] H-matrix built: "
                  << "DOF=" << m_stats.n_dof
                  << ", leaves=" << m_stats.n_leaves
                  << ", lowrank=" << m_stats.n_lowrank
                  << ", dense=" << m_stats.n_dense
                  << ", maxrank=" << m_stats.max_rank
                  << ", time=" << m_stats.build_time << "s"
                  << std::endl;
    }

    // Cache diagonal elements N_ii for Jacobi preconditioner (reused every
    // BiCGSTAB iteration). Uses virtual GetInteractionMatrixElement.
    m_diag_N.resize(m_ndof);
    ngcore::ParallelFor(ngcore::IntRange(m_ndof), [&](size_t i) {
        m_diag_N[(int)i] = GetInteractionMatrixElement((int)i, (int)i);
    });
    m_diag_cached = true;

    m_valid = true;
    return true;
}

//=========================================================================
// RadHACApKMSCManager::OnBeforeBuild
// MSC-specific precomputation switch (runs AFTER ExtractCoordinates has
// populated m_nffc, m_is_6dof/5dof/mixed_dof, and BEFORE HACApK callbacks).
//=========================================================================

void RadHACApKMSCManager::OnBeforeBuild() {
    if (!m_interaction) return;

    // Validate MSC DOF configuration
    if (m_nffc != 0 && m_nffc != 3 && m_nffc != 5 && m_nffc != 6) {
        std::cerr << "[HACApK] Warning: unexpected MSC nffc=" << m_nffc << std::endl;
    }

    // Register with callback (informational; ComputeEntry uses the manager directly)
    RadHACApKCallback::SetInteraction(m_interaction, m_n_elem, m_nffc);

    // ELF-style pre-computation for 6DOF hexahedra — shared with LU/BiCGSTAB.
    if (m_is_6dof) {
        m_interaction->PrecomputeHexaGeometry();
    }
    // ELF-style pre-computation for 3DOF tetrahedra (2025-12-26) — extracts
    // face vertices/normals for direct field computation without O(N^2)
    // SetupInteractMatrix().
    if (!m_is_6dof && !m_is_5dof && !m_is_mixed_dof) {
        PrecomputeGeometry3DOF();
    }
    // Pure wedge: precompute via radTInteraction
    if (m_is_5dof) {
        m_interaction->PrecomputeWedgeGeometry();
    }
    // FIX (2025-12-26): 3DOF tetrahedra fallback — if InteractMatrix was
    // already computed, use flat storage for O(1) element access.
    if (!m_is_6dof && !m_is_5dof && !m_is_mixed_dof && !m_geometry_3dof_ready && m_interaction->InteractMatrix != nullptr) {
        PrecomputeFlatInteractMatrix();
    }
    // Mixed DOF: precompute ALL element geometries (ComputeMixedBlockFast
    // needs each type present in the mesh).
    if (m_is_mixed_dof) {
        m_interaction->PrecomputeTetraGeometry();
        m_interaction->PrecomputeWedgeGeometry();
        m_interaction->PrecomputeHexaGeometry();
    }
}

//=========================================================================
// RadHACApKMSCManager::InitializeInvChi
// Populate m_inv_chi from material state. Matches ELF's
// initialize_chi_from_bh() for nonlinear isotropic materials (chi from the
// 2nd BH curve point) and falls back to DefineInstantKsiTensor(H=0) for
// linear materials.
//=========================================================================

void RadHACApKMSCManager::InitializeInvChi() {
    if (!m_interaction) return;
    m_inv_chi.resize(m_ndof);

    for (int i = 0; i < m_n_elem; i++) {
        radTg3dRelax* g3dRelaxPtr = m_interaction->g3dRelaxPtrVect[i];
        radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

        double chi = 1.0;
        radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
        if (NonlinMater != nullptr) {
            chi = NonlinMater->GetInitialChi_ELF_Style();
            if (chi <= 0) chi = 1.0;
        } else {
            TMatrix3d KsiTensor;
            TVector3d MrVect;
            TVector3d H_zero(0., 0., 0.);
            MaterPtr->DefineInstantKsiTensor(H_zero, KsiTensor, MrVect);
            chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
        }
        if (chi < 1.0e-6) chi = 1.0e-6;
        double inv_chi_val = 1.0 / chi;

        // All DOF on this element share the same 1/chi.
        int offset = m_dof_offset[i];
        int elem_dof = m_dof_offset[i + 1] - offset;
        for (int k = 0; k < elem_dof; k++) {
            m_inv_chi[offset + k] = inv_chi_val;
        }
    }
}

void RadHACApKBase::SetDeflation(const std::vector<int>& plaq_offsets,
                                 const std::vector<int>& dofs,
                                 const std::vector<double>& signs,
                                 double alpha) {
    m_defl_offsets = plaq_offsets;
    m_defl_dofs = dofs;
    m_defl_signs = signs;
    m_defl_alpha = alpha;
    m_defl_nplaq = plaq_offsets.empty() ? 0 : (int)plaq_offsets.size() - 1;
}

void RadHACApKBase::MatVec(const std::vector<double>& x, std::vector<double>& y) {
    if (!m_valid || !m_leafmtxp || !m_control) {
        std::fill(y.begin(), y.end(), 0.0);
        return;
    }

    int nd = HACApK_leafmtxp_get_nd(m_leafmtxp);
    HACApK_matvec_wrapper(m_leafmtxp, m_control, x.data(), y.data(), nd);

    // Matrix-free loop-mode deflation: y += alpha * L (L^T x).  L L^T acts
    // only on span(L) = the null (loop) subspace, lifting its eigenvalues
    // away from zero so the converged solution is loop-free.  O(nnz) = O(N).
    if (m_defl_nplaq > 0 && m_defl_alpha != 0.0) {
        for (int p = 0; p < m_defl_nplaq; p++) {
            double c = 0.0;
            for (int k = m_defl_offsets[p]; k < m_defl_offsets[p + 1]; k++)
                c += m_defl_signs[k] * x[m_defl_dofs[k]];
            c *= m_defl_alpha;
            for (int k = m_defl_offsets[p]; k < m_defl_offsets[p + 1]; k++)
                y[m_defl_dofs[k]] += m_defl_signs[k] * c;
        }
    }
}

void RadHACApKBase::UpdateDiagonal(const std::vector<double>& inv_chi) {
    if (!m_valid || !m_leafmtxp || !m_control) return;

    // Update stored inverse susceptibility
    m_inv_chi = inv_chi;
    RadHACApKCallback::SetInvChi(inv_chi);

    // OPTIMIZATION: Use fast diagonal update (2025-12-30)
    // Only updates true diagonal entries (i==j) using pre-computed N_ii values
    // This is O(ndof) instead of O(block_size^2 * n_diag_blocks)
    // Performance: ELF 0.39s vs old Radia 4.6s (12x slower due to slow method)
    if (m_diag_cached && m_diag_N.size() == (size_t)m_ndof) {
        HACApK_update_diagonal_fast_wrapper(m_leafmtxp, m_control,
                                             m_diag_N.data(), inv_chi.data(), m_ndof);
    } else {
        // Fallback to slow method (recomputes all entries in diagonal blocks)
        HACApK_update_diagonal_wrapper(m_leafmtxp, m_control, cHACApK_entry_ij);
    }
}

//=========================================================================
// GetInteractionMatrixElement: Access N matrix element
//
// This function returns the N(dof_i, dof_j) element of the interaction matrix.
// On-demand computation based on element DOF type:
// - 3DOF tetrahedra: Use B_comp() with PreRelax mode
// - 6DOF hexahedra: Use Yano-Sugahara MSC method (face-to-face interaction)
//
// On-demand computation is essential for H-matrix because:
// - HACApK uses ACA+ which only needs a subset of matrix elements
// - Pre-computing the full dense matrix would defeat the O(N log N) purpose
//=========================================================================

//=========================================================================
// Compute6x6Block: Calculate full 6x6 interaction block (ELF-style)
// K(face_i, face_j) = normal_i dot H_field(eval_pt_i, src_face_j)
//=========================================================================

void RadHACApKMSCManager::Compute6x6Block(int elem_i, int elem_j, double* K_mat) const {
    // Bounds check for element indices
    if (elem_i < 0 || elem_i >= m_n_elem || elem_j < 0 || elem_j >= m_n_elem) {
        std::cerr << "[HACApK] Error: Invalid element indices in Compute6x6Block: "
                  << elem_i << ", " << elem_j << " (n_elem=" << m_n_elem << ")" << std::endl;
        std::memset(K_mat, 0, 36 * sizeof(double));
        return;
    }

    radTg3dRelax* elem_row = m_interaction->g3dRelaxPtrVect[elem_i];
    radTg3dRelax* elem_col = m_interaction->g3dRelaxPtrVect[elem_j];

    if (!elem_row || !elem_col) {
        std::cerr << "[HACApK] Error: Null element pointer in Compute6x6Block: "
                  << (elem_row ? "col" : "row") << " elem" << std::endl;
        std::memset(K_mat, 0, 36 * sizeof(double));
        return;
    }

    radTPolyhedron* poly_row = dynamic_cast<radTPolyhedron*>(elem_row);
    radTPolyhedron* poly_col = dynamic_cast<radTPolyhedron*>(elem_col);

    if (!poly_row || !poly_col) {
        std::memset(K_mat, 0, 36 * sizeof(double));
        return;
    }

    // Pre-compute evaluation points and normals for row element (ELF optimization)
    TVector3d eval_pts[6];
    for (int fi = 0; fi < 6; fi++) {
        eval_pts[fi].x = 0.5 * (poly_row->FaceCenter[fi].x + poly_row->CentrPoint.x);
        eval_pts[fi].y = 0.5 * (poly_row->FaceCenter[fi].y + poly_row->CentrPoint.y);
        eval_pts[fi].z = 0.5 * (poly_row->FaceCenter[fi].z + poly_row->CentrPoint.z);
    }

    // Compute all 36 elements (source face outer loop for cache locality)
    // IMPORTANT: K_mat indexing must match GetCached6x6Element access pattern
    // K_mat[fi * 6 + fj] stores K(face_i, face_j) = normal_i dot H(eval_pt_i, src_face_j)
    for (int fj = 0; fj < 6; fj++) {
        double unit_point_charge = -1.0 * poly_col->FaceArea[fj];

        for (int fi = 0; fi < 6; fi++) {
            // Field from unit sigma on face_j
            TVector3d H_face = poly_col->FieldFromQuadFace(eval_pts[fi], fj, 1.0);
            TVector3d H_point = poly_col->FieldFromPointCharge(eval_pts[fi], unit_point_charge);

            TVector3d H_total;
            H_total.x = H_face.x + H_point.x;
            H_total.y = H_face.y + H_point.y;
            H_total.z = H_face.z + H_point.z;

            // K_ij = H_total dot normal_i
            double K_ij = H_total.x * poly_row->FaceNormal[fi].x +
                          H_total.y * poly_row->FaceNormal[fi].y +
                          H_total.z * poly_row->FaceNormal[fi].z;

            // Store K_ij / (4*pi) in row-major order: K_mat[row * 6 + col]
            // The solver will negate when building system matrix
            K_mat[fi * 6 + fj] = K_ij * RadConst::INV_FOUR_PI;
        }
    }
}

//=========================================================================
// GetCached6x6Element: ELF-style hash-based thread-local cache (NO locking!)
//
// ELF pattern from m_ppohBEM_user_func_unified.f90:
// - Single-entry cache checked first (most common case)
// - Hash-based cache (O(1) lookup, no LRU search)
// - NO global cache, NO critical sections
//
// This is THE key optimization: hash lookup + no locking = fast scaling
//
// FIX (2025-02-04): Added interaction pointer tracking to invalidate cache
// when IMA settings change between solves. Previously, stale cache values
// from non-IMA solves were incorrectly used for IMA solves.
//=========================================================================

// Hash-based cache size (must be power of 2 for fast modulo)
static constexpr int TL_HASH_SIZE = 1024;
static constexpr int TL_HASH_MASK = TL_HASH_SIZE - 1;

double RadHACApKMSCManager::GetCached6x6Element(int elem_i, int elem_j, int face_i, int face_j) const {
    // FIX (2025-02-04): Use generation counter for cache invalidation.
    // The generation counter is incremented each time a new HACApK manager builds its H-matrix.
    // This properly handles memory reuse where a new interaction object might have the
    // same pointer value as a previous (deleted) object.
    static thread_local uint64_t tl_cached_generation = 0;

    // Thread-local single-entry cache (most common case: same block accessed multiple times)
    static thread_local int tl_single_elem_i = -1;
    static thread_local int tl_single_elem_j = -1;
    static thread_local double tl_single_K_mat[36];

    // Thread-local hash-based cache (O(1) lookup, no locking!)
    static thread_local int tl_cache_elem_i[TL_HASH_SIZE];
    static thread_local int tl_cache_elem_j[TL_HASH_SIZE];
    static thread_local double tl_cache_K_mat[TL_HASH_SIZE][36];
    static thread_local bool tl_initialized = false;

    // Check if generation changed (new solve with different settings)
    // This handles IMA vs non-IMA transitions and geometry changes
    uint64_t current_gen = RadHACApKCallback::GetGeneration();
    if (tl_cached_generation != current_gen) {
        // Invalidate all caches
        tl_single_elem_i = -1;
        tl_single_elem_j = -1;
        for (int i = 0; i < TL_HASH_SIZE; i++) {
            tl_cache_elem_i[i] = -1;
            tl_cache_elem_j[i] = -1;
        }
        tl_cached_generation = current_gen;
        tl_initialized = true;  // Skip redundant initialization below
    }

    // Initialize thread-local cache on first access
    if (!tl_initialized) {
        for (int i = 0; i < TL_HASH_SIZE; i++) {
            tl_cache_elem_i[i] = -1;
            tl_cache_elem_j[i] = -1;
        }
        tl_initialized = true;
    }

    // Check single-entry cache first (fastest path)
    if (tl_single_elem_i == elem_i && tl_single_elem_j == elem_j) {
        return tl_single_K_mat[face_i * 6 + face_j];
    }

    // Compute hash index (ELF-style hash function)
    int hash_idx = ((elem_i * 73856093) ^ (elem_j * 19349663)) & TL_HASH_MASK;

    // Check hash cache (O(1) lookup!)
    if (tl_cache_elem_i[hash_idx] == elem_i && tl_cache_elem_j[hash_idx] == elem_j) {
        // Cache hit - copy to single-entry cache for repeated access
        std::memcpy(tl_single_K_mat, tl_cache_K_mat[hash_idx], 36 * sizeof(double));
        tl_single_elem_i = elem_i;
        tl_single_elem_j = elem_j;
        return tl_single_K_mat[face_i * 6 + face_j];
    }

    // Cache miss - compute the block using radTInteraction's precomputed geometry
    if (m_interaction && m_interaction->IsHexaTriDataReady()) {
        Compute6x6BlockFast(elem_i, elem_j, tl_single_K_mat);
    } else {
        Compute6x6Block(elem_i, elem_j, tl_single_K_mat);
    }
    tl_single_elem_i = elem_i;
    tl_single_elem_j = elem_j;

    // Insert into hash cache (overwrites any existing entry at this slot)
    tl_cache_elem_i[hash_idx] = elem_i;
    tl_cache_elem_j[hash_idx] = elem_j;
    std::memcpy(tl_cache_K_mat[hash_idx], tl_single_K_mat, 36 * sizeof(double));

    return tl_single_K_mat[face_i * 6 + face_j];
}

//=========================================================================
// GetInteractionMatrixElement: Optimized with O(1) lookup and LRU cache
// Supports 3DOF tetrahedra, 6DOF hexahedra, and mixed meshes
//=========================================================================

double RadHACApKMSCManager::GetInteractionMatrixElement(int dof_i, int dof_j) const {
    if (!m_interaction || dof_i < 0 || dof_i >= m_ndof || dof_j < 0 || dof_j >= m_ndof) {
        return 0.0;
    }

    // Safety check for lookup tables
    if (m_dof_to_elem.empty() || m_dof_to_local.empty()) {
        std::cerr << "[HACApK] Error: DOF lookup tables not initialized" << std::endl;
        return 0.0;
    }

    // O(1) DOF-to-element lookup (ELF-style optimization)
    int elem_i = m_dof_to_elem[dof_i];
    int local_i = m_dof_to_local[dof_i];
    int elem_j = m_dof_to_elem[dof_j];
    int local_j = m_dof_to_local[dof_j];

    // Safety check for element indices
    if (elem_i < 0 || elem_i >= m_n_elem || elem_j < 0 || elem_j >= m_n_elem) {
        return 0.0;
    }

    // Get element DOF counts
    int dof_elem_i = m_dof_offset[elem_i + 1] - m_dof_offset[elem_i];
    int dof_elem_j = m_dof_offset[elem_j + 1] - m_dof_offset[elem_j];

    // Dispatch based on DOF type of each element
    if (dof_elem_i == 6 && dof_elem_j == 6) {
        // 6DOF-6DOF: hex-hex interaction (Compute6x6BlockFast has IMA support)
        return GetCached6x6Element(elem_i, elem_j, local_i, local_j);
    } else if (dof_elem_i == 3 && dof_elem_j == 3) {
        // 3DOF-3DOF: tetra-tetra interaction
        // For mixed+IMA: Compute3x3Block_OnDemand handles IMA via B_comp() + IMA field context
        return GetCached3x3Element(elem_i, elem_j, local_i, local_j);
    } else if (dof_elem_i == 5 && dof_elem_j == 5) {
        // 5DOF-5DOF: wedge-wedge interaction (Compute5x5BlockFast has IMA support)
        return GetCached5x5Element(elem_i, elem_j, local_i, local_j);
    }

    // All cross-DOF pairs (3x5, 3x6, 5x3, 5x6, 6x3, 6x5):
    // Use unified ComputeMixedBlockFast kernel (IMA-aware, +N sign convention)
    return GetCachedMixedElement(elem_i, elem_j, dof_elem_i, dof_elem_j, local_i, local_j);
}

//=========================================================================
// GetCached3x3Element: On-demand 3x3 block computation with thread-local hash cache
// Similar to GetCached6x6Element for 6DOF hexahedra
// Uses O(1) hash lookup instead of O(n) LRU search
//
// FIX (2025-02-04): Added interaction pointer tracking to invalidate cache
// when IMA settings change between solves.
//=========================================================================

// Hash-based cache size for 3DOF (must be power of 2 for fast modulo)
static constexpr int TL_HASH_SIZE_3DOF = 1024;
static constexpr int TL_HASH_MASK_3DOF = TL_HASH_SIZE_3DOF - 1;

double RadHACApKMSCManager::GetCached3x3Element(int elem_i, int elem_j, int comp_i, int comp_j) const {
    // If flat storage is ready (pre-computed), use O(1) direct access
    if (m_flat_N_ready) {
        int64_t base_idx = ((int64_t)elem_i * m_n_elem + elem_j) * 9;
        double val = m_flat_N_data[base_idx + comp_i * 3 + comp_j];
        return val;
    }

    // On-demand computation with thread-local hash cache
    // This path is used when SetupInteractMatrix() is NOT called (HACApK optimization)

    // FIX (2025-02-04): Use generation counter for cache invalidation
    static thread_local uint64_t tl_cached_generation = 0;

    // Thread-local single-entry cache (most common case: same block accessed multiple times)
    static thread_local int tl_single_elem_i = -1;
    static thread_local int tl_single_elem_j = -1;
    static thread_local double tl_single_N_mat[9];

    // Thread-local hash cache (O(1) lookup, no locking needed!)
    static thread_local int tl_cache_elem_i[TL_HASH_SIZE_3DOF];
    static thread_local int tl_cache_elem_j[TL_HASH_SIZE_3DOF];
    static thread_local double tl_cache_N_mat[TL_HASH_SIZE_3DOF][9];
    static thread_local bool tl_initialized = false;

    // Check if generation changed (new solve with different settings)
    uint64_t current_gen = RadHACApKCallback::GetGeneration();
    if (tl_cached_generation != current_gen) {
        // Invalidate all caches
        tl_single_elem_i = -1;
        tl_single_elem_j = -1;
        for (int i = 0; i < TL_HASH_SIZE_3DOF; i++) {
            tl_cache_elem_i[i] = -1;
            tl_cache_elem_j[i] = -1;
        }
        tl_cached_generation = current_gen;
        tl_initialized = true;
    }

    // Initialize thread-local cache on first access
    if (!tl_initialized) {
        for (int i = 0; i < TL_HASH_SIZE_3DOF; i++) {
            tl_cache_elem_i[i] = -1;
            tl_cache_elem_j[i] = -1;
        }
        tl_initialized = true;
    }

    // Check single-entry cache first (fastest path)
    if (tl_single_elem_i == elem_i && tl_single_elem_j == elem_j) {
        return tl_single_N_mat[comp_i * 3 + comp_j];
    }

    // Compute hash index (same hash function as 6DOF cache)
    int hash_idx = ((elem_i * 73856093) ^ (elem_j * 19349663)) & TL_HASH_MASK_3DOF;

    // Check hash cache (O(1) lookup!)
    if (tl_cache_elem_i[hash_idx] == elem_i && tl_cache_elem_j[hash_idx] == elem_j) {
        // Cache hit - copy to single-entry cache for repeated access
        std::memcpy(tl_single_N_mat, tl_cache_N_mat[hash_idx], 9 * sizeof(double));
        tl_single_elem_i = elem_i;
        tl_single_elem_j = elem_j;

        return tl_single_N_mat[comp_i * 3 + comp_j];
    }

    // Cache miss - compute the 3x3 block on-demand
    // Use fast method with pre-computed geometry if available
    double N_mat[9];
    if (m_geometry_3dof_ready) {
        Compute3x3BlockFast(elem_i, elem_j, N_mat);
    } else {
        Compute3x3Block_OnDemand(elem_i, elem_j, N_mat);
    }
    // Both paths now return +N (physical quantity).
    // ComputeEntry() handles the sign flip to -N for the system matrix.

    // Update single-entry cache
    std::memcpy(tl_single_N_mat, N_mat, 9 * sizeof(double));
    tl_single_elem_i = elem_i;
    tl_single_elem_j = elem_j;

    // Store in hash cache (overwrites any existing entry at this hash index)
    tl_cache_elem_i[hash_idx] = elem_i;
    tl_cache_elem_j[hash_idx] = elem_j;
    std::memcpy(tl_cache_N_mat[hash_idx], N_mat, 9 * sizeof(double));

    return N_mat[comp_i * 3 + comp_j];
}

//=========================================================================
// GetCached5x5Element: On-demand 5x5 block for wedges with hash cache
// Same pattern as GetCached6x6Element / GetCached3x3Element
// Delegates to radTInteraction::Compute5x5BlockFast (IMA-aware)
//=========================================================================

static constexpr int TL_HASH_SIZE_5DOF = 512;
static constexpr int TL_HASH_MASK_5DOF = TL_HASH_SIZE_5DOF - 1;

double RadHACApKMSCManager::GetCached5x5Element(int elem_i, int elem_j, int face_i, int face_j) const {
    static thread_local uint64_t tl_cached_generation = 0;
    static thread_local int tl_single_elem_i = -1;
    static thread_local int tl_single_elem_j = -1;
    static thread_local double tl_single_K_mat[25];
    static thread_local int tl_cache_elem_i[TL_HASH_SIZE_5DOF];
    static thread_local int tl_cache_elem_j[TL_HASH_SIZE_5DOF];
    static thread_local double tl_cache_K_mat[TL_HASH_SIZE_5DOF][25];
    static thread_local bool tl_initialized = false;

    uint64_t current_gen = RadHACApKCallback::GetGeneration();
    if (tl_cached_generation != current_gen) {
        tl_single_elem_i = -1;
        tl_single_elem_j = -1;
        for (int i = 0; i < TL_HASH_SIZE_5DOF; i++) {
            tl_cache_elem_i[i] = -1;
            tl_cache_elem_j[i] = -1;
        }
        tl_cached_generation = current_gen;
        tl_initialized = true;
    }
    if (!tl_initialized) {
        for (int i = 0; i < TL_HASH_SIZE_5DOF; i++) {
            tl_cache_elem_i[i] = -1;
            tl_cache_elem_j[i] = -1;
        }
        tl_initialized = true;
    }

    if (tl_single_elem_i == elem_i && tl_single_elem_j == elem_j) {
        return tl_single_K_mat[face_i * 5 + face_j];
    }

    int hash_idx = ((elem_i * 73856093) ^ (elem_j * 19349663)) & TL_HASH_MASK_5DOF;
    if (tl_cache_elem_i[hash_idx] == elem_i && tl_cache_elem_j[hash_idx] == elem_j) {
        std::memcpy(tl_single_K_mat, tl_cache_K_mat[hash_idx], 25 * sizeof(double));
        tl_single_elem_i = elem_i;
        tl_single_elem_j = elem_j;
        return tl_single_K_mat[face_i * 5 + face_j];
    }

    // Cache miss: compute via radTInteraction::Compute5x5BlockFast (IMA-aware)
    double K_mat[25];
    if (m_interaction) {
        m_interaction->Compute5x5BlockFast(elem_i, elem_j, K_mat);
    } else {
        std::memset(K_mat, 0, 25 * sizeof(double));
    }

    std::memcpy(tl_single_K_mat, K_mat, 25 * sizeof(double));
    tl_single_elem_i = elem_i;
    tl_single_elem_j = elem_j;

    tl_cache_elem_i[hash_idx] = elem_i;
    tl_cache_elem_j[hash_idx] = elem_j;
    std::memcpy(tl_cache_K_mat[hash_idx], K_mat, 25 * sizeof(double));

    return K_mat[face_i * 5 + face_j];
}

//=========================================================================
// GetCachedMixedElement: On-demand mixed-DOF block with hash cache
// Delegates to radTInteraction::ComputeMixedBlockFast (IMA-aware)
// Handles all cross-DOF pairs: 3x5, 3x6, 5x3, 5x6, 6x3, 6x5
//=========================================================================

static constexpr int TL_HASH_SIZE_MIXED = 256;
static constexpr int TL_HASH_MASK_MIXED = TL_HASH_SIZE_MIXED - 1;
static constexpr int MAX_BLOCK_SIZE = 36;  // max DOF product: 6x6

double RadHACApKMSCManager::GetCachedMixedElement(int elem_i, int elem_j,
	int dof_i, int dof_j, int local_i, int local_j) const
{
	static thread_local uint64_t tl_cached_generation = 0;
	static thread_local int tl_single_elem_i = -1;
	static thread_local int tl_single_elem_j = -1;
	static thread_local int tl_single_dof_i = 0;
	static thread_local int tl_single_dof_j = 0;
	static thread_local double tl_single_block[MAX_BLOCK_SIZE];
	static thread_local int tl_cache_elem_i[TL_HASH_SIZE_MIXED];
	static thread_local int tl_cache_elem_j[TL_HASH_SIZE_MIXED];
	static thread_local int tl_cache_dof_i[TL_HASH_SIZE_MIXED];
	static thread_local int tl_cache_dof_j[TL_HASH_SIZE_MIXED];
	static thread_local double tl_cache_block[TL_HASH_SIZE_MIXED][MAX_BLOCK_SIZE];
	static thread_local bool tl_initialized = false;

	uint64_t current_gen = RadHACApKCallback::GetGeneration();
	if (tl_cached_generation != current_gen) {
		tl_single_elem_i = -1; tl_single_elem_j = -1;
		for (int i = 0; i < TL_HASH_SIZE_MIXED; i++) { tl_cache_elem_i[i] = -1; tl_cache_elem_j[i] = -1; }
		tl_cached_generation = current_gen; tl_initialized = true;
	}
	if (!tl_initialized) {
		for (int i = 0; i < TL_HASH_SIZE_MIXED; i++) { tl_cache_elem_i[i] = -1; tl_cache_elem_j[i] = -1; }
		tl_initialized = true;
	}

	if (tl_single_elem_i == elem_i && tl_single_elem_j == elem_j) {
		return tl_single_block[local_i * tl_single_dof_j + local_j];
	}

	int hash_idx = ((elem_i * 73856093) ^ (elem_j * 19349663)) & TL_HASH_MASK_MIXED;
	if (tl_cache_elem_i[hash_idx] == elem_i && tl_cache_elem_j[hash_idx] == elem_j) {
		int bs = tl_cache_dof_i[hash_idx] * tl_cache_dof_j[hash_idx];
		std::memcpy(tl_single_block, tl_cache_block[hash_idx], bs * sizeof(double));
		tl_single_elem_i = elem_i; tl_single_elem_j = elem_j;
		tl_single_dof_i = tl_cache_dof_i[hash_idx]; tl_single_dof_j = tl_cache_dof_j[hash_idx];
		return tl_single_block[local_i * tl_single_dof_j + local_j];
	}

	// Cache miss: compute via unified kernel
	double block[MAX_BLOCK_SIZE];
	if (m_interaction) {
		m_interaction->ComputeMixedBlockFast(elem_i, dof_i, elem_j, dof_j, block);
	} else {
		std::memset(block, 0, MAX_BLOCK_SIZE * sizeof(double));
	}

	int bs = dof_i * dof_j;
	std::memcpy(tl_single_block, block, bs * sizeof(double));
	tl_single_elem_i = elem_i; tl_single_elem_j = elem_j;
	tl_single_dof_i = dof_i; tl_single_dof_j = dof_j;

	tl_cache_elem_i[hash_idx] = elem_i; tl_cache_elem_j[hash_idx] = elem_j;
	tl_cache_dof_i[hash_idx] = dof_i; tl_cache_dof_j[hash_idx] = dof_j;
	std::memcpy(tl_cache_block[hash_idx], block, bs * sizeof(double));

	// ComputeMixedBlockFast returns +N (physical quantity) for all DOF types.
	// ComputeEntry() handles the sign flip to -N for the system matrix.
	return block[local_i * dof_j + local_j];
}

//=========================================================================
// Compute3x3Block: Compute 3x3 interaction block for tetrahedra
// N_ij = interaction matrix element between magnetization components
// Uses existing radTInteraction::InteractMatrix
//=========================================================================

void RadHACApKMSCManager::Compute3x3Block(int elem_i, int elem_j, double* N_mat) const {
    // InteractMatrix[elem_i][elem_j] returns TMatrix3df (3x3 float matrix)
    //
    // MATRIX LAYOUT FIX (2025-12-24):
    // InteractMatrix[i][j] stores TMatrix3df where:
    //   Str0 = dH/dMx = (dHx/dMx, dHy/dMx, dHz/dMx) <- COLUMN vector (response to Mx)
    //   Str1 = dH/dMy = (dHx/dMy, dHy/dMy, dHz/dMy) <- COLUMN vector (response to My)
    //   Str2 = dH/dMz = (dHx/dMz, dHy/dMz, dHz/dMz) <- COLUMN vector (response to Mz)
    //
    // For row k of the 3x3 block, we need N[i][j] element (k, l) = dH_k/dM_l
    // This requires transposed access: row k gets Str0[k], Str1[k], Str2[k]
    //
    // IMPORTANT: Radia's InteractMatrix stores positive N, but the system matrix is:
    //   A = -N - diag(1/chi) (ELF-compatible)
    // The LU solver uses: BaseMatrix = -Nij
    // For HACApK, GetInteractionMatrixElement should return -N to match.

    if (!m_interaction || !m_interaction->InteractMatrix) {
        std::memset(N_mat, 0, 9 * sizeof(double));
        return;
    }

    // Get the 3x3 block for element pair (i, j)
    // InteractMatrix is indexed by element indices, not DOF indices
    const TMatrix3df& M = m_interaction->InteractMatrix[elem_i][elem_j];

    // Extract to row-major double array with sign flip (-N) and transpose
    // Row 0 (Hx response): -dHx/dMx, -dHx/dMy, -dHx/dMz
    N_mat[0] = -static_cast<double>(M.Str0.x);  // -dHx/dMx
    N_mat[1] = -static_cast<double>(M.Str1.x);  // -dHx/dMy
    N_mat[2] = -static_cast<double>(M.Str2.x);  // -dHx/dMz
    // Row 1 (Hy response): -dHy/dMx, -dHy/dMy, -dHy/dMz
    N_mat[3] = -static_cast<double>(M.Str0.y);  // -dHy/dMx
    N_mat[4] = -static_cast<double>(M.Str1.y);  // -dHy/dMy
    N_mat[5] = -static_cast<double>(M.Str2.y);  // -dHy/dMz
    // Row 2 (Hz response): -dHz/dMx, -dHz/dMy, -dHz/dMz
    N_mat[6] = -static_cast<double>(M.Str0.z);  // -dHz/dMx
    N_mat[7] = -static_cast<double>(M.Str1.z);  // -dHz/dMy
    N_mat[8] = -static_cast<double>(M.Str2.z);  // -dHz/dMz
}

//=========================================================================
// Compute3x3Block_OnDemand: On-demand 3x3 interaction block computation
// Computes interaction directly using B_comp() without pre-computed matrix
// This is used by HACApK to avoid O(N^2) matrix pre-computation
//=========================================================================

void RadHACApKMSCManager::Compute3x3Block_OnDemand(int elem_i, int elem_j, double* N_mat) const {
    // Compute interaction from element j to observation at element i center
    // using B_comp() directly (same approach as SetupInteractMatrix)
    //
    // Returns +N (physical demagnetization tensor).
    // ComputeEntry() handles the sign flip to -N for the system matrix.

    std::memset(N_mat, 0, 9 * sizeof(double));

    if (!m_interaction || elem_i < 0 || elem_i >= m_n_elem ||
        elem_j < 0 || elem_j >= m_n_elem) {
        return;
    }

    radTg3dRelax* elem_row = m_interaction->g3dRelaxPtrVect[elem_i];
    radTg3dRelax* elem_col = m_interaction->g3dRelaxPtrVect[elem_j];
    if (!elem_row || !elem_col) return;

    TVector3d ObsPoiVect = elem_row->ReturnCentrPoint();

    radTFieldKey FieldKeyInteract;
    FieldKeyInteract.B_ = FieldKeyInteract.H_ = FieldKeyInteract.PreRelax_ = 1;

    TVector3d ZeroVect(0., 0., 0.);
    radTField Field(FieldKeyInteract, m_interaction->CompCriterium, ObsPoiVect,
                   ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
    Field.AmOfIntrctElemWithSym = m_interaction->CountRelaxElemsWithSym();

    elem_col->B_comp(&Field);

    // PreRelax mode: Field.B = dH/dMx, Field.H = dH/dMy, Field.A = dH/dMz
    // Return +N (physical): N_mat[row][col] = dH_row / dM_col
    N_mat[0] = Field.B.x;  N_mat[1] = Field.H.x;  N_mat[2] = Field.A.x;
    N_mat[3] = Field.B.y;  N_mat[4] = Field.H.y;  N_mat[5] = Field.A.y;
    N_mat[6] = Field.B.z;  N_mat[7] = Field.H.z;  N_mat[8] = Field.A.z;

    // IMA mirror contributions: B_comp above computes ONLY the direct interaction.
    // For IMA, add mirror contributions using the full IMA-aware kernel, then subtract
    // the direct part (which B_comp already computed correctly).
    if (m_interaction->IsIMAEnabled() && m_interaction->IsTetraGeomReady()) {
        double N_full[9];
        m_interaction->Compute3x3BlockFast(elem_i, elem_j, N_full);
        // N_full = direct + IMA mirrors (from radTInteraction kernel)
        // N_mat = direct (from B_comp, numerically exact for LU compatibility)
        // We want: N_mat_final = N_mat(direct) + (N_full - N_direct_kernel)
        // But N_full already contains direct+mirror, and N_mat contains direct.
        // Since both compute the same physical quantity (direct N), their difference
        // is numerical noise. So N_full - N_mat ≈ pure mirror contribution.
        // Add the mirror contribution: N_final = N_mat + (N_full - N_mat) = N_full
        // This is equivalent to just using N_full, but let's use the IMA mirror
        // extraction approach for clarity and to preserve B_comp's direct accuracy.
        //
        // Actually, the simplest correct approach: use N_full directly when IMA is on.
        // Both Compute3x3BlockFast and B_comp compute the same physics; the kernel
        // version additionally includes IMA mirrors. Any numerical difference in the
        // direct part is within discretization error.
        std::memcpy(N_mat, N_full, 9 * sizeof(double));
    }
}

//=========================================================================
// Compute3x3BlockFast: ELF-style fast 3x3 block computation for tetrahedra
// Uses pre-computed geometry arrays (no object access, no B_comp overhead)
//
// OPTIMIZATION (2025-12-31): Rewritten to compute face basis ONCE per face
// and reuse it for all 3 magnetization directions. This reduces coordinate
// transformation overhead from 12x to 4x per element pair.
//
// The 3x3 interaction matrix N[i][j] represents how magnetization M_j creates
// demagnetizing field H at element i's center:
//   H(r_i) = sum_faces [ sigma_f * integral_triangle H_field dA ] / (4*pi)
// where sigma_f = M_j dot n_f (surface charge density)
//
// For PreRelax mode, we compute dH/dM for each unit M direction.
//=========================================================================

void RadHACApKMSCManager::Compute3x3BlockFast(int elem_i, int elem_j, double* N_mat) const {
    std::memset(N_mat, 0, 9 * sizeof(double));

    if (!m_geometry_3dof_ready || elem_i < 0 || elem_i >= m_n_elem ||
        elem_j < 0 || elem_j >= m_n_elem) {
        return;
    }

    // Observation point: center of element i
    const double* obs_ptr = &m_tetra_centers[elem_i * 3];
    TVector3d obsPoint(obs_ptr[0], obs_ptr[1], obs_ptr[2]);

    // Column element centroid (for normal orientation check)
    const double* col_center = &m_tetra_centers[elem_j * 3];
    TVector3d elemCentroid(col_center[0], col_center[1], col_center[2]);

    // Unit magnetization vectors
    const TVector3d unit_M[3] = {
        TVector3d(1.0, 0.0, 0.0),
        TVector3d(0.0, 1.0, 0.0),
        TVector3d(0.0, 0.0, 1.0)
    };

    // Process each of the 4 triangular faces
    // ELF-style optimization: compute face basis ONCE per face, reuse for all 3 M directions
    for (int f = 0; f < 4; f++) {
        // Get face vertices from pre-computed geometry
        int fvIdx = (elem_j * 4 + f) * 3 * 3;
        const double* V0_ptr = &m_tetra_face_vertices[fvIdx + 0];
        const double* V1_ptr = &m_tetra_face_vertices[fvIdx + 3];
        const double* V2_ptr = &m_tetra_face_vertices[fvIdx + 6];

        TVector3d V0(V0_ptr[0], V0_ptr[1], V0_ptr[2]);
        TVector3d V1(V1_ptr[0], V1_ptr[1], V1_ptr[2]);
        TVector3d V2(V2_ptr[0], V2_ptr[1], V2_ptr[2]);

        // Compute face basis ONCE (ELF optimization - saves 2 basis computations per face)
        RadTriangleFaceBasis basis;
        RadComputeTriangleFaceBasis(V0, V1, V2, elemCentroid, basis);

        if (!basis.valid) continue;

        // Get face normal for surface charge computation
        // Use pre-computed values from m_tetra_face_normals
        int fnIdx = (elem_j * 4 + f) * 3;
        TVector3d faceNormal(
            m_tetra_face_normals[fnIdx + 0],
            m_tetra_face_normals[fnIdx + 1],
            m_tetra_face_normals[fnIdx + 2]
        );

        // Compute H field for each unit magnetization direction
        // Reuse the same face basis for all 3 directions
        for (int j = 0; j < 3; j++) {
            // Surface charge: sigma = M dot n * sign_correction
            double sigma = (unit_M[j].x * faceNormal.x +
                           unit_M[j].y * faceNormal.y +
                           unit_M[j].z * faceNormal.z) * basis.charge_sign;

            if (std::abs(sigma) < 1.0e-15) continue;

            // Compute field using pre-computed basis (fast path)
            TVector3d H = RadFieldFromTriangleFaceWithBasis(basis, sigma, obsPoint);

            // Accumulate into tensor: N_mat[i*3+j] = dH_i / dM_j (+N convention)
            // ComputeEntry() handles the sign flip to -N for the system matrix.
            N_mat[0*3 + j] += H.x;
            N_mat[1*3 + j] += H.y;
            N_mat[2*3 + j] += H.z;
        }
    }
}

double RadHACApKMSCManager::GetGenericElement(int elem_i, int elem_j, int local_i, int local_j) const {
    // Generic path: access pre-computed flat interaction matrix (+N convention)
    if (!m_interaction || m_interaction->m_flatInteractMatrix.empty()) return 0.0;

    int offset_i = m_dof_offset[elem_i];
    int offset_j = m_dof_offset[elem_j];
    int total_dof = m_interaction->m_totalDOF;

    return m_interaction->m_flatInteractMatrix[(offset_i + local_i) * total_dof + (offset_j + local_j)];
}

//=========================================================================
// End of rad_hacapk.cpp
//=========================================================================
