/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk_field.cpp
*
* Project:        RADIA
*
* Description:    Embed-in-square H-matrix field evaluation.  See
*                 rad_hacapk_field.h for the formulation.
*
* First release:  2026-07
*
-------------------------------------------------------------------------*/

#include "rad_hacapk_field.h"

#include "rad_geometry_3d.h"
#include "rad_group.h"
#include "rad_polyhedron.h"
#include "rad_application.h"
#include "rad_type_cast.h"
#include "rad_constants.h"     // RadConst::MU_0 (B = mu0 * leaf B_comp; see Compute3x3)

#include <functional>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <cmath>

// Global Radia application (handle -> object resolution), same as
// RadFieldUnified / the rest of the core.
extern radTApplication rad;

//-------------------------------------------------------------------------

RadHACApKFieldEval::RadHACApKFieldEval(int container_handle, std::vector<double> obs_points)
    : m_handle(container_handle), m_nObs(0), m_nSrc(0), m_obs(std::move(obs_points)), m_scale(1.0)
{
    if (m_obs.size() % 3 != 0)
        throw std::runtime_error("RadHACApKFieldEval: obs_points length must be a multiple of 3");
    m_nObs = (int)(m_obs.size() / 3);
}

//-------------------------------------------------------------------------
// ExtractCoordinates: resolve the container to its leaf magnetic elements,
// then build the COMBINED [obs; src] cluster-tree point set (uniform 3-DOF).
//-------------------------------------------------------------------------

void RadHACApKFieldEval::ExtractCoordinates()
{
    m_src.clear();
    m_srcM.clear();

    radThg hg;
    if (!rad.ValidateElemKey(m_handle, hg))
        throw std::runtime_error("RadHACApKFieldEval: invalid container handle "
                                 + std::to_string(m_handle));
    radTg3d* g3dPtr = radTCast::g3dCast(hg.rep);
    if (!g3dPtr)
        throw std::runtime_error("RadHACApKFieldEval: handle is not a 3D object");

    // Recursively collect leaf relaxable (magnetic) elements -- the same walk
    // as RadFieldUnified::BuildElementData (rad_field_unified.cpp:193).  NOTE:
    // ancestor-group transforms are NOT composed here, so this supports a FLAT
    // global-coordinate container (an ObjCnt of magnets); a transformed /
    // space-symmetry container would need the transform stack folded in (the
    // per-leaf B_genComp below would otherwise miss it).
    std::function<void(radTg3d*)> collect = [&](radTg3d* elem) {
        if (!elem) return;
        radTGroup* group = radTCast::GroupCast(elem);
        if (group) {
            for (auto& child : group->GroupMapOfHandlers) {
                radTg3d* childElem = radTCast::g3dCast(child.second.rep);
                if (childElem) collect(childElem);
            }
            return;
        }
        radTg3dRelax* relax = radTCast::g3dRelaxCast(elem);
        if (relax) m_src.push_back(relax);
    };
    collect(g3dPtr);

    m_nSrc = (int)m_src.size();
    if (m_nSrc == 0)
        throw std::runtime_error("RadHACApKFieldEval: container has no magnetic source elements");
    if (m_nObs == 0)
        throw std::runtime_error("RadHACApKFieldEval: no observation points");

    // Actual magnetization per source (src-slot half of the matvec vector).
    m_srcM.resize((size_t)3 * m_nSrc);
    for (int s = 0; s < m_nSrc; ++s) {
        const TVector3d& M = m_src[s]->Magn;
        m_srcM[3 * s + 0] = M.x;
        m_srcM[3 * s + 1] = M.y;
        m_srcM[3 * s + 2] = M.z;
    }

    // Combined [obs; src] point set for the square cluster tree.  Elements
    // 0..N_obs-1 are the observation points (field slots); N_obs..end are the
    // source centroids (magnetization slots).  Uniform 3-DOF.
    m_n_elem = m_nObs + m_nSrc;
    m_coordinates.resize((size_t)m_n_elem * 3);
    m_dof_offset.resize((size_t)m_n_elem + 1);

    for (int o = 0; o < m_nObs; ++o) {
        m_coordinates[3 * o + 0] = m_obs[3 * o + 0];
        m_coordinates[3 * o + 1] = m_obs[3 * o + 1];
        m_coordinates[3 * o + 2] = m_obs[3 * o + 2];
    }
    for (int s = 0; s < m_nSrc; ++s) {
        const TVector3d& c = m_src[s]->CentrPoint;
        int e = m_nObs + s;
        m_coordinates[3 * e + 0] = c.x;
        m_coordinates[3 * e + 1] = c.y;
        m_coordinates[3 * e + 2] = c.z;
    }
    for (int e = 0; e <= m_n_elem; ++e) m_dof_offset[e] = 3 * e;  // uniform 3-DOF
    m_ndof = 3 * m_n_elem;
}

//-------------------------------------------------------------------------
// OnBeforeBuild: compute the O(1) normalisation scale = max |raw entry| over a bounded sample of
// (obs, src) pairs.  The raw B-response is O(mu0*H) ~ 1e-13 for far obs; HACApK's ACA stopping is
// effectively ABSOLUTE, so without this the tiny entries never accrue rank (measured: field err stuck
// at ~0.27 until eps ~ 1e-14).  Dividing stored entries by m_scale puts them at O(1) so eps ~ 1e-8 is
// meaningful; the pybind entry()/matvec() multiply back by Scale() to return Tesla.
//-------------------------------------------------------------------------

void RadHACApKFieldEval::OnBeforeBuild()
{
    const int nOs = std::min(4, m_nObs);
    const int nSs = std::min(16, m_nSrc);
    double mx = 0.0;
    double G[9];
    for (int oi = 0; oi < nOs; ++oi) {
        int o = (nOs > 1) ? (int)((long long)oi * (m_nObs - 1) / (nOs - 1)) : 0;   // spread over obs
        for (int si = 0; si < nSs; ++si) {
            int s = (nSs > 1) ? (int)((long long)si * (m_nSrc - 1) / (nSs - 1)) : 0;
            Compute3x3(o, s, G);
            for (int k = 0; k < 9; ++k) { double a = std::fabs(G[k]); if (a > mx) mx = a; }
        }
    }
    m_scale = (mx > 0.0) ? mx : 1.0;
}

//-------------------------------------------------------------------------
// Compute3x3: B-field response of source s at observation point o, per unit
// magnetization (the exact G_ab, NOT a dipole).  Bit-consistent with rad.Fld
// (same radTFieldKey.B_ + B_genComp as RadFieldUnified::ComputeFieldSingle).
//-------------------------------------------------------------------------

void RadHACApKFieldEval::Compute3x3(int o, int s, double* G) const
{
    radTg3dRelax* src = m_src[s];
    TVector3d obs(m_obs[3 * o + 0], m_obs[3 * o + 1], m_obs[3 * o + 2]);
    TVector3d ZeroVect(0.0, 0.0, 0.0);

    // B_genComp reads elem->Magn; set it to a unit vector, evaluate, restore.
    // The parallel HACApK fill must not race on this, hence the mutex (serial
    // kernel -- correctness first; a mutation-free per-source unit-M clone /
    // analytic face kernel is the follow-up optimization).
    std::lock_guard<std::mutex> lk(m_fieldMutex);
    TVector3d saveM = src->Magn;
    try {
        for (int b = 0; b < 3; ++b) {
            TVector3d unitM(0.0, 0.0, 0.0);
            if (b == 0) unitM.x = 1.0; else if (b == 1) unitM.y = 1.0; else unitM.z = 1.0;
            src->Magn = unitM;

            radTFieldKey FieldKey;
            FieldKey.B_ = true;
            radTField Field(FieldKey, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect);
            Field.P = obs;
            src->B_genComp(&Field);

            // Leaf B_comp stores the PRE-mu0 quantity in Field.B (H_total in A/m outside; H_total+M
            // inside) -- radTg3d::B_genComp does NOT scale it; rad.Fld applies mu0 at the top level.
            // Multiply here so the H-matrix field is Tesla, bit-consistent with rad.Fld('b').
            G[0 * 3 + b] = RadConst::MU_0 * Field.B.x;   // a = 0 (Bx)
            G[1 * 3 + b] = RadConst::MU_0 * Field.B.y;   // a = 1 (By)
            G[2 * 3 + b] = RadConst::MU_0 * Field.B.z;   // a = 2 (Bz)
        }
    } catch (...) {
        src->Magn = saveM;
        throw;
    }
    src->Magn = saveM;
}

//-------------------------------------------------------------------------
// GetInteractionMatrixElement: A = [[0, K],[0, 0]] on the combined DOF space.
// Only the obs-row x src-col block is nonzero (the field response); the
// embed-zero blocks return 0 immediately (no B_genComp).
//-------------------------------------------------------------------------

double RadHACApKFieldEval::GetInteractionMatrixElement(int dof_i, int dof_j) const
{
    if (dof_i < 0 || dof_i >= m_ndof || dof_j < 0 || dof_j >= m_ndof)
        throw std::out_of_range("RadHACApKFieldEval entry index out of range: i="
                                + std::to_string(dof_i) + " j=" + std::to_string(dof_j)
                                + " ndof=" + std::to_string(m_ndof));

    // SYMMETRIC embed A = [[0, K],[K^T, 0]] (NOT [[0,K],[0,0]]): the field response fills BOTH the
    // obs-row x src-col block (K) and the src-row x obs-col block (K^T).  y = A [0; M] = [K M; 0], so the
    // obs read-out is UNCHANGED, but the matrix is now symmetric -- HACApK's cluster-tree build / fill is
    // the well-tested symmetric path (ChargeGram / MMM), whereas the one-sided [[0,K],[0,0]] mis-fills its
    // sub-blocks under a fine (multi-cluster) tree.  o=obs elem, a=obs field comp; s=src elem, b=src M comp.
    const int obsDof = 3 * m_nObs;
    int o, a, s, b;
    if (dof_i < obsDof && dof_j >= obsDof) {            // upper block: (obs row, src col) = K
        o = dof_i / 3;              a = dof_i % 3;
        s = (dof_j - obsDof) / 3;   b = (dof_j - obsDof) % 3;
    } else if (dof_i >= obsDof && dof_j < obsDof) {     // lower block: (src row, obs col) = K^T
        s = (dof_i - obsDof) / 3;   b = (dof_i - obsDof) % 3;
        o = dof_j / 3;              a = dof_j % 3;
    } else {
        return 0.0;                                    // obs-obs / src-src blocks are zero
    }

    // Per-thread single-block memo: the nine (a,b) scalars of one (o,s) block
    // share three B_genComp calls (HACApK fetches a block's entries together).
    // Keyed on the build generation too: a persistent thread_local would
    // otherwise serve a stale (o,s) block from a PRIOR instance that happened
    // to be reallocated at the same address (ABA) -- GetGeneration() bumps on
    // every BuildHMatrix, so the memo self-invalidates across builds.
    struct Memo { const void* owner = nullptr; uint64_t gen = 0; int o = -1, s = -1; double G[9]; };
    static thread_local Memo memo;
    const uint64_t gen = RadHACApKCallback::GetGeneration();
    if (memo.owner != this || memo.gen != gen || memo.o != o || memo.s != s) {
        Compute3x3(o, s, memo.G);
        memo.owner = this;
        memo.gen = gen;
        memo.o = o;
        memo.s = s;
    }
    // Stored NORMALISED (raw / m_scale) so HACApK's ACA sees O(1) entries; pybind entry()/matvec()
    // multiply by Scale() to return Tesla.
    return memo.G[a * 3 + b] / m_scale;
}
