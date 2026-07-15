/* rad_hacapk_hdiv.h -- HACApK H-matrix for the symmetric RT1 HDiv-VIM charge Gram.
 *
 * Production stores the Coulomb Gram G and applies the demagnetizing operator as
 * N = B^T G B.  Symmetric-fill HACApK plus MatVecSym makes N exactly symmetric at
 * the discrete apply level, so the physical material system
 * ((1/chi) M_mass + N) is SPD and is solved by mass-Riesz-preconditioned CG.
 *
 * Conventions (CLAUDE.md): row-major [target][source]; +N physical sign; HACApK ACA+ only.
 */
#ifndef __RAD_HACAPK_HDIV_H
#define __RAD_HACAPK_HDIV_H

#include "rad_hacapk.h"     // RadHACApKBase
#include "rad_hdiv_vim.h"   // analytic charge-potential kernels
#include "rad_hdiv_field_evaluator.h"
#include "rad_planar_charges.h"
#include <unordered_map>
#include <memory>
#include <string>
#include <utility>
#include <atomic>

// Persistent PARDISO factor of the HDiv mass for the MASS RIESZ preconditioner, cached on
// RadHACApKChargeGram across solve calls (defined in rad_hacapk_hdiv.cpp next to MassRieszPardiso).
struct RadMassRieszCache;

//-------------------------------------------------------------------------
// RadHACApKChargeGram: the charge-charge Coulomb Gram G as a HACApK H-matrix.
//-------------------------------------------------------------------------

/* The general-mesh production path.  Charges = volume and boundary charge-basis functions
 * extracted from the supported NGSolve RT1 HDiv space.
 * This manager builds the n_charge x n_charge
 * Coulomb Gram G as a HACApK H-matrix (a clean 1/r kernel over the charge centroids):
 *   G[a][b] = meas_a meas_b / (4pi |c_a - c_b|)   (a != b, centroid monopole)
 *   G[a][a] = self_energy[a]                       (the accurate sub-divided self, computed by the
 *                                                   caller per element shape -- tet/tri/hex/quad).
 * The demag operator N = B^T G B is then applied matrix-free as B^T (G-Hmatvec (B m)) with B the
 * sparse charge map; the H-matrix gives the O(N log N) Gram matvec that makes the solve scalable
 * on real geometry.  Stores +G (ComputeSystemEntry = default); MatVec/stats from RadHACApKBase. */
class RadHACApKChargeGram : public RadHACApKBase {
public:
    // MONOPOLE mode: centroids [n*3] charge centroids; measures [n] cell volumes / face areas;
    // self_energy [n] the diagonal G[a][a] (caller-computed accurate self-energy, element-shape-aware).
    RadHACApKChargeGram(std::vector<double> centroids,
                        std::vector<double> measures,
                        std::vector<double> self_energy);

    // ANALYTIC mode (M2b): the EXACT charge Gram from per-charge GEOMETRY -- matches the independent
    // analytic reference entry-by-entry.  cell_verts [n_el*12] (tets, 4 verts) then
    // face_verts [n_bf*9] (triangles, 3 verts); the n_charge charges are the n_el volume cells
    // (rho = -div M) followed by the n_bf boundary faces (sigma = M.n).  Entry
    //   G[a][b] = (1/4pi) INT_a Phi_b   (Phi_b = PhiTet/TriPotential of source b, exact analytic),
    // the outer INT_a by tet barycentric sub-points (cells) / Dunavant-5 (faces), symmetrized; the
    // diagonal is the analytic self (the Wilton/phi_tet potential is exact through the 1/r singularity).
    //
    // near_factor: the NEAR/FAR split that makes the BUILD fast.  A pair (a,b) is NEAR when
    // |c_a - c_b| <= near_factor*(size_a + size_b) and uses the expensive analytic entry; FAR pairs use
    // the cheap centroid-monopole meas_a*meas_b/(4pi r).  Physically correct: the analytic entry only
    // matters for the NEAR (non-uniform-M, div M != 0) interaction; the far field is monopole to
    // O((size/r)^2) (the validated monopole+near-correction split).  near_factor defaults to a HUGE
    // value => ALL pairs near => all-analytic (matches the independent analytic reference golden);
    // pass near_factor ~= 2 for the fast split.
    // The CELL outer-quad uses a built-in 4-pt Gauss-Duffy tet rule (64 nodes) -- the lowest-order charge is
    // constant so the inner is the EXACT analytic PhiTet and INT_T PhiTet dx is smooth, integrated to
    // ~machine precision.  (The old hardcoded equal-weight _bary_tet(3) under-integrated the volume self-
    // energy by ~6.5% -- golden-invisible because every uniform-M demag golden has div M = 0.)
    // far_quad: the FAR-pair evaluation when near_factor < inf.  0 (default) = centroid-MONOPOLE
    // (meas_a*meas_b/4pi r) -- cheap but O((size/r)^2), it slightly breaks geometric symmetry (the uniform
    // sphere transverse leak).  >0 = a low-order DOUBLE-QUADRATURE of the true 1/r kernel over a degree-2
    // rule (4-pt tet / 3-pt tri) -- ~16 cheap evals/pair (vs the full analytic's ~1e3 transcendentals), but
    // accurate to O((size/r)^4) so it reproduces the all-analytic Gram (verified: sphere transverse 7.26e-4
    // == exact 7.25e-4, vs monopole 1.19e-3).  This is the precision-preserving build speedup: near=analytic,
    // far=low-quad.  Tet/tri and polytope analytic modes both support this option.
    RadHACApKChargeGram(std::vector<double> cell_verts,
                        std::vector<double> face_verts,
                        int n_el, double near_factor = 1e30,
                        std::vector<int> image_masks = {}, std::vector<double> image_signs = {},
                        int far_quad = 0);

    // POLYTOPE analytic mode (HEX / WEDGE cells + QUAD faces): the same EXACT analytic entry as the tet/
    // triangle analytic mode above, generalized to ANY flat-faced convex cell + quad face with NO new
    // singular quadrature -- the cell Newtonian potential is the divergence-theorem sum over the cell's
    // (convex-hull) TRIANGLES of the SAME exact Wilton triangle potential (rad_hdiv::TriPotential), and a
    // quad face is two flat triangles.  Matches the independent analytic polytope reference
    // polytope path entry-by-entry.  The triangulation is supplied FROM PYTHON (which computes the convex
    // hulls / sub-triangles) so the C++ needs no hull algorithm and the two share the exact decomposition:
    //   cell_tris  : flat triangle soup (9 doubles/tri = P0,P1,P2 each xyz) of ALL cells' hull triangles;
    //   cell_troff : [n_el+1] CSR offsets into cell_tris (in TRIANGLES) -> cell c's hull tris;
    //   cell_cent  : [n_el*3] cell vertex-mean centroid (the centroid-fan apex AND the outward-normal ref);
    //   cell_meas  : [n_el]   cell volume;
    //   face_tris  : flat triangle soup of ALL boundary faces' sub-triangles (triangle->1, quad->2);
    //   face_troff : [n_bf+1] CSR offsets into face_tris; face_cent [n_bf*3]; face_meas [n_bf] area.
    // Outer quadrature: cells use the built-in 64-node Gauss-Duffy tet rule on each centroid-fan sub-tet
    // (apex = cell_cent); faces use built-in Dunavant-5 per sub-triangle -- the SAME rules the tet/tri mode
    // uses, so an all-tet/all-triangle mesh routed through here would agree with the tet mode to quadrature
    // precision (it is NOT routed here: the tet mode is kept bit-identical via cell_verts/face_verts).
    // near_factor: identical NEAR/FAR build split as the analytic mode (default 1e30 = all-analytic).
    // far_quad: same precision-preserving far option as the tet/tri ctor, here on the polytope (hex/wedge)
    // charges -- the low-order FAR rule is a degree-2 quadrature on the SAME centroid-fan sub-tets (cells,
    // 4-pt) / sub-triangles (faces, 3-pt) used for the outer quad, so it reproduces the all-analytic Gram
    // for hex/wedge far pairs at ~monopole cost.  0 (default) = centroid-monopole far.
    RadHACApKChargeGram(std::vector<double> cell_tris, std::vector<int> cell_troff,
                        std::vector<double> cell_cent, std::vector<double> cell_meas,
                        std::vector<double> face_tris, std::vector<int> face_troff,
                        std::vector<double> face_cent, std::vector<double> face_meas,
                        int n_el, double near_factor = 1e30,
                        std::vector<int> image_masks = {}, std::vector<double> image_signs = {},
                        int far_quad = 0);

    // HIGH-ORDER (order-p) mode: POLYNOMIAL charges (a monomial basis on each host element), the order-p
    // extension checked against the independent high-order analytic reference.  charge_host[c] = host element
    // index (a cell when charge_kind[c]==0, a boundary face when ==1); charge_expo[3*c+{0,1,2}] = the
    // monomial exponents in the host's REFERENCE barycentric coords (tet cell: lam1^i lam2^j lam3^k;
    // face: lam1^i lam2^j, the 3rd exponent ignored).  The SAME reference convention is used by the
    // Python charge-density map B (so B and G share the basis; N = B^T G B is then basis-invariant and
    // matches the NGSolve-L2-basis dense reference's demag).  Entry
    //   G[a][b] = (1/4pi) INT_ha INT_hb m_a(x) m_b(y)/|x-y|
    // = the monomial-WEIGHTED outer quad (m_a folded into the outer weights) x the polynomial-charge
    // inner potential PhiAt(b,.) by singularity SUBTRACTION reusing the exact PhiTet/TriPotential through
    // the 1/r singularity.  ref_tet_pts[nqt*3]/ref_tet_w[nqt] (weights sum to 1/6) + ref_tri_pts[nqr*2]/
    // ref_tri_w[nqr] (sum to 1/2) are the REFERENCE-element Gauss-Duffy rules (Python-supplied), mapped per
    // host and used for BOTH the monomial-weighted outer quad and the FIXED inner-potential subtraction
    // table.
    //
    // NEAR/FAR adaptive quadrature (build speedup, ACCURACY-PRESERVING): the per-entry cost is the nested
    // outer x inner quadrature (~O(quad^6) for vol-vol).  The expensive HIGH-quad subtraction is only needed
    // for NEAR/self pairs (through the 1/r singularity); for a well-separated (FAR) pair the kernel 1/|x-y|
    // is SMOOTH, so a CHEAP LOW-quad PLAIN double-Gauss (no PhiTet, no subtraction) is already accurate.  If
    // the optional ref_*_lo LOW-quad rules are supplied (and ho_far_factor < inf), a pair (a,b) with
    // |c_a-c_b| > ho_far_factor*(size_a+size_b) uses QuadDotFar (the low-quad plain double sum); NEAR/self
    // pairs keep the full high-quad subtraction.  This is NOT FMM/multipole (the zero-mean high-order modes
    // have zero monopole, so a monopole far is WRONG); it is just adaptive quadrature order, and the HACApK
    // ACA still compresses the well-separated low-rank blocks (now from cheap entries).  ho_far_factor
    // defaults to 1e30 (=> all pairs NEAR => the original all-high-quad behavior, golden-equivalent).
    // The production path averages both FAR directions.  A one-sided block keeps matrix symmetry when the
    // upper triangle is mirrored, but is not invariant under an explicit reflected mesh at finite quadrature
    // order.  RADIA_HDIV_HO_FAR_ONESIDED=1 is therefore diagnostic/benchmark-only.
    RadHACApKChargeGram(std::vector<double> cell_verts, std::vector<double> face_verts, int n_el,
                        std::vector<int> charge_host, std::vector<int> charge_kind,
                        std::vector<int> charge_expo,
                        std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
                        std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
                        std::vector<double> ref_tet_pts_lo = {}, std::vector<double> ref_tet_w_lo = {},
                        std::vector<double> ref_tri_pts_lo = {}, std::vector<double> ref_tri_w_lo = {},
                        double ho_far_factor = 1e30,
                        std::vector<double> ref_tet_pts_in = {}, std::vector<double> ref_tet_w_in = {},
                        std::vector<double> ref_tri_pts_in = {}, std::vector<double> ref_tri_w_in = {},
                        std::vector<int> image_masks = {}, std::vector<double> image_signs = {});

    // CURVED HIGH-ORDER mode (isoparametric P2, curve_order=2): the same monomial-charge Gram as the flat HO
    // mode above, but on a CURVED (mesh.Curve(2)) geometry -- the boundary surface charge sigma=M.n and the
    // volume charge live on the true curved element.  cell_nodes [n_cell*30] = 10 P2 nodes/tet (corners
    // 0..3, mid-edges 4=(0-1),5=(1-2),6=(2-0),7=(0-3),8=(1-3),9=(2-3)); face_nodes [n_bf*18] = 6 P2 nodes/tri
    // (corners 0..2, mid-edges 3=(0-1),4=(1-2),5=(2-0)).  charge_host/kind/expo: the monomial is in the
    // NGSolve REFERENCE frame (the Python charge map B uses the SAME reference-frame change-of-basis, so B and
    // G share the basis).  The OUTER quadrature maps the reference Gauss-Duffy points (ref_tet_pts/ref_tri_pts)
    // through the curved P2 map X(xi) + curved measure (CurvedTet/TriMapMeasure), folding xi^expo at the
    // REFERENCE point; the INNER potential is the curved Duffy (rad_hdiv::CurvedTet/TriPotential, gl/gw = an
    // nq-point Gauss-Legendre rule on [0,1]) for near/self pairs.  Optional ref_*_lo rules enable the same
    // smooth-kernel far double-Gauss path as the flat high-order operator; all low points are mapped through
    // the curved P2 geometry in C++.  No analytic curved moments and no inner-subtraction table.  ROLE
    // (de-risked 2026-06-28): curved helps NEAR-SURFACE FIELD / FLUX accuracy
    // (sigma on the true curved surface), NOT the volume-averaged demag FACTOR (which is curving-insensitive on
    // a sphere, ~3e-5).  Accuracy is the Duffy ~1e-3..1e-5 (the order>=3 conditioning caveat applies on curved).
    RadHACApKChargeGram(std::vector<double> cell_nodes, std::vector<double> face_nodes,
                        std::vector<int> cell_vertices, std::vector<int> face_vertices,
                        int n_el, int curve_order,
                        std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
                        std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
                        std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
                        std::vector<double> curve_gl, std::vector<double> curve_gw,
                        std::vector<double> ref_tet_pts_lo = {}, std::vector<double> ref_tet_w_lo = {},
                        std::vector<double> ref_tri_pts_lo = {}, std::vector<double> ref_tri_w_lo = {},
                        double ho_far_factor = 1e30,
                        std::vector<int> image_masks = {}, std::vector<double> image_signs = {});

    // CURVED POLYTOPE mode (curved hex/wedge): FULLY curved -- both the CELL volume charge AND the boundary
    // FACE surface charge live on the true mesh.Curve(2) geometry (the cell volume charge is DOMINANT for the
    // demag, NOT zero -- the lowest-order curved charge representation cannot represent uniform M exactly,
    // so div M != 0).  Each curved hex CELL
    // is split (Python, Kuhn 6-tet) into curved P2 sub-TETS, each fed through CurvedTetPotential; each curved
    // quad FACE into curved P2 sub-TRIS through CurvedTriPotential -- reusing the golden curved-tet/tri kernels.
    // cell_curved_nodes [n_cell_subtet*30] = 10 P2 nodes/sub-tet; cell_subtet_off [n_cell+1] CSR; ditto
    // face_curved_nodes [n_bf_subtri*18] + face_subtri_off [n_bf+1].  The CELL outer quad maps ref_tet_pts
    // through CurvedTetMapMeasure; the FACE outer quad maps ref_tri_pts through CurvedTriMapMeasure.
    // PhiAt(cell,.) = sum CurvedTetPotential, PhiAt(face,.) = sum CurvedTriPotential (constant charge ->
    // monomial exponent 0; curve_gl/gw the inner Duffy rule).  14 vectors + int n_el (distinct overload).
    RadHACApKChargeGram(std::vector<double> cell_curved_nodes, std::vector<int> cell_subtet_off,
                        std::vector<double> cell_cent, std::vector<double> cell_meas,
                        std::vector<double> face_curved_nodes, std::vector<int> face_subtri_off,
                        std::vector<double> face_cent, std::vector<double> face_meas,
                        std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
                        std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
                        std::vector<double> curve_gl, std::vector<double> curve_gw, int n_el);

    // HEX RT1/RT2 mode: the HDiv(hexmesh, order=p) charge Gram uses tensor Qp monomial charges
    // (8/27 per volume and 4/9 per quad face at p=1/2) over the DIRECT Q2 isoparametric geometry:
    // hex_cell_nodes [n_el*81] = 27-node triquadratic lattice (n = ix + 3*iy + 9*iz, ref (ix/2,iy/2,iz/2)),
    // quad_face_nodes [n_bf*27] = 9-node biquadratic lattice (n = iu + 3*iv), both extracted via GetTrafo at
    // the reference lattice -> ONE code path for FLAT (trilinear subset of Q2, machine-exact incl. distorted)
    // and CURVED (mesh.Curve(2)) hexes.  Quadrature follows the numpy-validated eig(M^-1 N)<=1 scheme
    // (independent prototype eig 0.998): per (target sub-simplex, source sub-simplex) pair
    // of the ref-hex 6-sub-tet / ref-quad 2-sub-tri decomposition,
    //   OUTER: near/self sub pair (|cA-cB| <= near_grade*(sA+sB), hosts near) -> a Duffy-graded product rule
    //          (gl_out/gw_out, graded toward the source sub centroid; BOTH domains graded is the key eig<=1
    //          lesson); far sub pair -> the regular symmetric rule (sym_tet_* = Keast-15 / sym_tri_* =
    //          Dunavant-7, bary pts + weights summing to the ref simplex measure 1/6 / 1/2).
    //   INNER: field point far from the source sub (> far_inner_factor*size) -> the cheap far rule
    //          (far_tet_*/far_tri_*); else the FINE corner-graded Duffy (gl_in/gw_in, graded toward the
    //          source-sub vertex nearest the field point) -- effectively exact.
    // All rules are consumed as tables (Python owns the constants).  Mixed and pyramid meshes remain
    // outside this pure-topology production constructor and fail before dispatch.
    RadHACApKChargeGram(std::vector<double> hex_cell_nodes, std::vector<double> quad_face_nodes,
                        int n_el, int n_bf,
                        std::vector<int> charge_host, std::vector<int> charge_kind,
                        std::vector<int> charge_expo,
                        std::vector<double> sym_tet_pts, std::vector<double> sym_tet_w,
                        std::vector<double> sym_tri_pts, std::vector<double> sym_tri_w,
                        std::vector<double> gl_out, std::vector<double> gw_out,
                        std::vector<double> gl_in, std::vector<double> gw_in,
                        std::vector<double> far_tet_pts, std::vector<double> far_tet_w,
                        std::vector<double> far_tri_pts, std::vector<double> far_tri_w,
                        double near_grade, double far_inner_factor,
                        std::vector<int> image_masks = {}, std::vector<double> image_signs = {});

    // 2D PLANAR mode (2026-07-03, the motor cross-section layer; memory hdiv-vim-tri-quad-motor):
    // charges rho = -div M on 2D cells (RT1: TRI P0 / QUAD Q1; RT2: TRI P1 / QUAD Q2)
    // + sigma = M.n on boundary EDGES (P1/P2), all in the NGSolve REF
    // frame with the Piola-exact extraction (the dimension-independent J-cancellation identity).  Kernel
    // = the 2D Laplace Green's function -ln(r)/(2*pi); the ln-scale shift is killed because every HDiv
    // dof's total charge is zero (divergence theorem), so N = B^T G B is scale-invariant.  Geometry =
    // Q1..Q3 polynomial maps fitted from GetTrafo (cell_type 0=tri/1=quad) -> ONE path flat/curved.
    // RT1 admits geometry through Q2; RT2 admits geometry through Q3.  Quadrature
    // (numpy-validated, C:\temp\vim2d_proto.py): REGULAR symmetric outer everywhere (the log kernel's
    // single-layer potentials are continuous, so NO graded outer is needed -- simpler than 3D); inner =
    // signed radial cones from the nearest anchor SITE for near field points, cheap far cloud otherwise.
    // Gates: eig(M^-1 N) in [0,1]; DISK demag == 1/2 exact; ellipse a:b -> N=b/(a+b); 2D Clausius-
    // Mossotti M = chi H/(1+chi/2).  The 2D-mode flag is the FIRST argument (dim2 tag) to keep the
    // overload set unambiguous.
    RadHACApKChargeGram(int dim2_tag, int geometry_order,
                        std::vector<double> cell_map, std::vector<int> cell_type,
                        std::vector<double> edge_map,
                        int n_el, int n_be,
                        std::vector<int> charge_host, std::vector<int> charge_kind,
                        std::vector<int> charge_expo,
                        std::vector<double> sym_tri_pts, std::vector<double> sym_tri_w,
                        std::vector<double> gl_quad, std::vector<double> gw_quad,
                        std::vector<double> gl_edge, std::vector<double> gw_edge,
                        std::vector<double> gl_in, std::vector<double> gw_in,
                        std::vector<double> far_tri_pts, std::vector<double> far_tri_w,
                        double near_grade, double far_inner_factor,
                        std::vector<int> image_masks = {}, std::vector<double> image_signs = {});

    // WEDGE (PRISM) RT1/RT2 mode: the HDiv(prismmesh, order=p) charge Gram uses tri-Pp (x) z-Pp
    // volume charges (6/18 per prism) over the 18-node tri-P2 (x) z-P2 lattice (node n = t + 6*iz,
    // t = the Tri6 node 0..5, iz = the z level 0..2).  Boundary faces are MIXED (2 tri caps + 3 quad
    // sides per prism): a tri face -> SurfaceL2 P1/P2 (3/6 monomials, 1 sub-tri), a quad
    // face -> SurfaceL2 Q1/Q2 (4/9 monomials, 2 sub-tris).  Faces
    // are stored in 9-node slots (a tri fills the first 6) + a per-face type array (0=tri/1=quad),
    // mirroring the 2D mode's mixed-cell storage.  The cell 3-sub-tet decomposition (WEDGEREF_TETS) + the
    // per-face-type 1/2 sub-tri decomposition drive the SAME both-domains-graded Duffy singular quadrature
    // as the hex mode (numpy de-risk eig(M_mass^-1 N) in [0,1]: 0.989 @ n=2, 0.997 @ n=3, demag_z ~ 1/3).
    // Served through the hex block memo / GetHexBlock dispatch (m_wedgemode -> QuadBlockWedge); the leaf
    // helpers (HexMonoEval, HexDuffyBary, the HexSiteRad tables, HexGetCloud / HexQuadCloud) are shared.
    // Arg 3 (vector<int> face_type) disambiguates the
    // overload from the hex ctor (whose arg 3 is int n_el).
    RadHACApKChargeGram(std::vector<double> wedge_cell_nodes, std::vector<double> face_nodes,
                        std::vector<int> face_type, int n_el, int n_bf,
                        std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
                        std::vector<double> sym_tet_pts, std::vector<double> sym_tet_w,
                        std::vector<double> sym_tri_pts, std::vector<double> sym_tri_w,
                        std::vector<double> gl_out, std::vector<double> gw_out,
                        std::vector<double> gl_in, std::vector<double> gw_in,
                        std::vector<double> far_tet_pts, std::vector<double> far_tet_w,
                        std::vector<double> far_tri_pts, std::vector<double> far_tri_w,
                        double near_grade, double far_inner_factor,
                        std::vector<int> image_masks = {}, std::vector<double> image_signs = {});
    // Q2 lattice geometry maps (PUBLIC static utilities: the file-local cloud builder uses them too).
    static void HexQ2Map(const double* nd27, const double xi[3], double X[3], double J[3][3]);
    static void QuadQ2Map(const double* nd9, const double uv[2], double X[3], double T[3][2]);
    // Values-only variants (no Jacobian): the radial inner integrands are Piola (REF measure, no |det J|),
    // so the self radial loop needs only X -- ~4x less shape work per point than the full map.
    static void HexQ2MapX(const double* nd27, const double xi[3], double X[3]);
    static void QuadQ2MapX(const double* nd9, const double uv[2], double X[3]);
    // WEDGE geometry maps: prism 18-node tri-P2 (x) z-P2 (WedgeQ2MapX), surface tri 6-node P2 (TriSurfMap,
    // 3D X since a boundary face lives in space).  Values-only (the Piola charge model never uses |det J|).
    static void WedgeQ2MapX(const double* nd18, const double xi[3], double X[3]);
    static void TriSurfMap(const double* nd18, const double uv[2], double X[3]);   // nd18 = 6 nodes x 3D
    ~RadHACApKChargeGram() override {}

    double GetInteractionMatrixElement(int a, int b) const override;

    // SYMMETRIC-FILL build (2026-07-03): the charge Gram is symmetric BY CONSTRUCTION (every entry is the
    // 0.5*(AB+BA)-symmetrized kernel), and every apply of it goes through the symmetric H-matvec, so the
    // strictly-lower H-matrix leaves are never read -- matvec_sym mirrors the upper triangle.  This shadow
    // of the (non-virtual) base build turns on cHACApK_set_sym_fill around the fill: the lower leaves stay
    // EMPTY, saving ~half the build time (ACA sampling + dense-block entry fill + the per-thread block-memo
    // rebuild of the mirror leaves) and ~half the leaf memory.  The UPPER leaves fill identically, so
    // MatVecSym is bit-identical to a full build.  Plain MatVec / MatVecTranspose are ROUTED to MatVecSym
    // (for the symmetric operator they are the same map; the base implementations would silently read the
    // empty lower leaves -- the routing makes that failure mode unrepresentable).
    bool BuildHMatrix(const RadHACApKParams& params = RadHACApKParams());
    // Initialize charge count/coordinates without allocating HACApK blocks.
    // Used by prescribed-magnetization field sources, which need the exact
    // charge geometry and B map but never apply the charge Gram matrix.
    void PrepareGeometryOnly() { ExtractCoordinates(); }
    void MatVec(const std::vector<double>& x, std::vector<double>& y) { MatVecSym(x, y); }
    void MatVecTranspose(const std::vector<double>& x, std::vector<double>& y) { MatVecSym(x, y); }

    // M3 (the iterative-solve hot kernel in C++): solve the SPD HDiv-VIM linear material system
    //   ((1/chi) M_mass + B^T G B) m = rhs
    // by Jacobi-preconditioned conjugate gradients, with G applied as THIS charge-Gram H-matvec
    // (O(N log N)) -- no dense N, no Python per-iteration glue.  This is the linear soft-iron demag
    // solve AND the symmetric Picard warmstart of the nonlinear Newton.  Sparse inputs are caller-
    // provided: B as CSR over charges (B_indptr [n_charge+1], B_indices/B_data = face columns, so
    // (B x)[charge] = sum data*x[face]); M_mass as COO (mI,mJ,mV) on the n_face DOFs; prec = the
    // Jacobi diagonal of the system (length n_face).  Returns m (length n_face); iters_out = CG iters.
    // mass_riesz=false: diagonal-Jacobi PCG (z = r/prec).  mass_riesz=true (the DEFAULT 'auto' path):
    // PCG preconditioned by a PARDISO SPD factor of the HDiv mass M_mass (z = M_mass^{-1} r, the MASS
    // RIESZ map) built from the COO (mI,mJ,mV) -- ~3-5x fewer iters, nearly mu_r-flat; `prec` is
    // then ignored.  Moves the whole linear demag solve (H-matvec + mass solve + Krylov) into C++.
    // The factor is PERSISTENT on the object (m_massRieszCache, exact-COO key): constant-mass chains --
    // the Hantila hysteresis loop (W(nu0) fixed by construction) and the C++ scalar Picard (geometry
    // M_mass; scalar inv_chi is outside the preconditioner) -- factor once.  Callers that pass a
    // per-iteration TANGENT mass (the Python nu-secant / Newton W_tan) compare-miss and refactor as before.
    // symmetric=true (DEFAULT): G is applied via the EXACTLY-symmetric H-matvec (MatVecSym, upper-tri
    // leaves define both triangles), so the +N CG operator is machine-symmetric and removes the reported
    // independently-ACA'd off-diagonal asymmetry failure mode.  symmetric=false
    // uses the general (asymmetric) MatVec (legacy / cross-check only).
    std::vector<double> SolveLinearMaterial(
        const std::vector<int>& B_indptr, const std::vector<int>& B_indices,
        const std::vector<double>& B_data, int n_face,
        const std::vector<int>& mI, const std::vector<int>& mJ, const std::vector<double>& mV,
        double inv_chi, const std::vector<double>& prec, const std::vector<double>& rhs,
        double tol, int maxit, int& iters_out, bool mass_riesz = false, bool symmetric = true,
        const std::vector<double>* x0 = nullptr);
    // Matrix-free postprocessing primitives used by the production planar and
    // three-dimensional HDiv paths.  ApplyDemagOperator computes B^T G B x
    // without materializing N.  ApplyMassRiesz computes M_mass^{-1} rhs with a
    // dedicated persistent PARDISO factor, so nonlinear material factors used
    // by SolveLinearMaterial are not evicted by local-field postprocessing.
    std::vector<double> ApplyDemagOperator(
        const std::vector<int>& B_indptr, const std::vector<int>& B_indices,
        const std::vector<double>& B_data, int n_face,
        const std::vector<double>& x, bool symmetric = true);
    std::vector<double> ApplyMassRiesz(
        const std::vector<int>& mI, const std::vector<int>& mJ,
        const std::vector<double>& mV, int n_face,
        const std::vector<double>& rhs);

    // NGSolve-style persistent operator configuration.  The charge map B is geometry/FESpace data and is
    // registered once with the Gram object; the material mass is replaced only when the constitutive tangent
    // changes.  Production matvec/solve calls then cross Python with vectors only, instead of copying CSR/COO
    // topology on every application.
    void ConfigureChargeMap(
        std::vector<int> B_indptr, std::vector<int> B_indices,
        std::vector<double> B_data, int n_face);
    void ConfigureMassMatrix(
        std::vector<int> mI, std::vector<int> mJ,
        std::vector<double> mV, int n_face);
    void ConfigureGeometryMassMatrix(
        std::vector<int> mI, std::vector<int> mJ,
        std::vector<double> mV, int n_face);
    bool HasConfiguredChargeMap() const { return m_operatorChargeConfigured; }
    bool HasConfiguredMassMatrix() const { return m_operatorMassConfigured; }
    bool HasConfiguredGeometryMassMatrix() const { return m_operatorGeometryMassConfigured; }
    int ConfiguredNFace() const { return m_operatorNFace; }
    int ConfiguredConstraintCount() const;
    std::vector<double> ApplyConfiguredDemag(
        const std::vector<double>& x, bool symmetric = true);
    void ApplyConfiguredDemag(
        const double* x, double* y, bool symmetric = true);
    void ApplyConfiguredDemagAdd(
        double scale, const double* x, double* y, bool symmetric = true);
    std::vector<double> ApplyConfiguredGeometryMass(const std::vector<double>& x);
    void ApplyConfiguredGeometryMass(const double* x, double* y);
    std::vector<double> ApplyConfiguredMassRiesz(const std::vector<double>& rhs);
    std::vector<double> SolveConfiguredLinearMaterial(
        double inv_chi, const std::vector<double>& rhs, double tol, int maxit,
        int& iters_out, bool mass_riesz = true, bool symmetric = true,
        const std::vector<double>* x0 = nullptr);
    std::vector<double> SolveConfiguredLinearMaterialAutoPrec(
        double inv_chi, const std::vector<double>& rhs, double tol, int maxit,
        int& iters_out, double& prec_min, double& prec_max,
        const std::vector<double>* x0 = nullptr);
    std::shared_ptr<rad_hdiv::HDivFieldEvaluator> CreateConfiguredFieldEvaluator(
        const std::vector<double>& magnetization,
        const rad_hdiv::FieldEvaluatorOptions& options = {}) const;
    std::shared_ptr<rad_planar_charges::PlanarFieldEvaluator> CreateConfiguredPlanarFieldEvaluator(
        const std::vector<double>& magnetization) const;
    std::vector<std::pair<std::string, double>> LastSolveTimings() const;

    // M3 (the NONLINEAR solve in C++): scalar-chi Picard for the isotropic nonlinear demag.
    // Each Picard step is a mass-Riesz SolveLinearMaterial solve of
    // ((1/chi) M_mass + B^T G B) m = H0*(M_mass mu),
    // then chi <- 0.5 chi + 0.5*chi_sec(|H|) with the closed-form saturating curve
    //   M(H) = chi0 H / (1 + chi0 |H|/Msat)   ->   chi_sec(|H|) = chi0/(1 + chi0|H|/Msat),
    // and the scalar self-consistent field H = H0 - Dscal*M_avg, Dscal = mu.(B^T G B mu)/denom,
    // M_avg = mu.(M_mass m)/denom.  Converges to the scalar fixed point M_avg = M(H0 - Dscal*M_avg)
    // -- the full nonlinear physics for an isotropic body, with NO NGSolve per iteration (the
    // per-element tensor-tangent refinement for non-uniform M stays NGSolve).  All sparse inputs as in
    // SolveLinearMaterial; Mmass_diag + N_diag are retained for API compatibility / diagnostics, but the
    // active preconditioner is the PARDISO mass-Riesz map, matching the production tet Picard path.
    struct PicardResult { std::vector<double> m; double Mavg; double chi; double Dscal; int iters; };
    PicardResult SolveNonlinearPicard(
        const std::vector<int>& B_indptr, const std::vector<int>& B_indices,
        const std::vector<double>& B_data, int n_face,
        const std::vector<int>& mI, const std::vector<int>& mJ, const std::vector<double>& mV,
        const std::vector<double>& Mmass_diag, const std::vector<double>& N_diag,
        const std::vector<double>& mu, double denom,
        double chi0, double Msat, double H0,
        int picard_iters, double cg_tol, int cg_maxit);

protected:
    void ExtractCoordinates() override;
    void OnBeforeBuild() override;
    void InitializeInvChi() override { m_inv_chi.assign(m_ndof, 0.0); }
    bool IsVariableDOF() const override { return false; }
    int  GetUniformNFFC() const override { return 1; }

private:
    void ApplyConfiguredDemagImpl(
        const double* x, double* y, double scale, bool add, bool symmetric);
    double PhiAt(int src, const double p[3]) const;   // exact analytic potential of source charge src at p
    double QuadDot(int tgt, int src) const;            // (1/4pi) sum_p w_p PhiAt(src, p) over tgt's outer quad
    // IMA image term: (1/4pi) sum_p w_p PhiAt(src, R_mask(p)) -- tgt's outer points reflected on the mask
    // axes.  Uses Phi_{R(b)}(x) = Phi_b(R(x)) (reflection isometry), so only the eval point is mirrored.
    double QuadDotRefl(int tgt, int src, int mask) const;
    // FAR low-order double-quadrature (analytic mode, far_quad>0): (1/4pi) sum_i sum_j qwf[a][i] qwf[b][j] /
    // |qpf[a][i]-qpf[b][j]| over the degree-2 far rule -- symmetric in (a,b) by construction (1/r symmetric).
    double QuadDotFarLow(int a, int b) const;

    std::vector<double> m_cent, m_meas, m_self;        // monopole mode (m_cent also = the cluster-tree points)
    int  m_n = 0;
    // analytic mode (M2b)
    bool m_analytic = false;
    int  m_n_el = 0;
    std::vector<double> m_cellV, m_faceV;              // [n_el*12], [n_bf*9]
    std::vector<std::vector<rad_hdiv::Vec3>> m_qp;     // [n] outer-quad points per charge
    std::vector<std::vector<double>>          m_qw;    // [n] outer-quad weights per charge
    std::vector<double> m_size;                        // [n] characteristic size: vol^(1/3) / area^(1/2)
    double m_near_factor = 1e30;                       // near/far split: NEAR if |c_a-c_b| <= nf*(size_a+size_b)
    int    m_far_quad = 0;                             // 0=monopole far; >0=low-order double-quad far (tet/tri)
    std::vector<std::vector<rad_hdiv::Vec3>> m_qpf;    // [n] FAR low-order quad points (built iff m_far_quad>0)
    std::vector<std::vector<double>>          m_qwf;   // [n] FAR low-order quad weights (sum = measure)
    // IMA mirror symmetry (image method): G_IMA(a,b) = G(a,b) + sum_i sign_i*0.5*(QuadDotRefl(a,b,mask_i)
    // + QuadDotRefl(b,a,mask_i)).  The 2^P-1 non-empty subsets of the P mirror planes (image_group); each
    // reflects the eval point on its axes.  Always the full analytic image (the self-on-plane image can be
    // singular -> needs the exact PhiTet/TriPotential, not a monopole far).  Empty = no IMA.
    std::vector<int>    m_image_masks;                 // [n_img] 3-bit axis mask (bit0=x,1=y,2=z) of the subset
    std::vector<double> m_image_signs;                 // [n_img] product-sign of the subset

    // POLYTOPE analytic mode (hex/wedge): per-charge source triangulation (cell hull tris / face sub-tris).
    // PhiAt(src,.) is the divergence-theorem polytope potential (cell) / sum-of-sub-triangle (face) over
    // these; the outer quadrature (m_qp/m_qw) is built in the ctor (centroid-fan / Dunavant per sub-tri).
    bool m_polytope = false;
    std::vector<std::vector<std::array<rad_hdiv::Vec3, 3>>> m_srcTris;  // [n] source triangle soup per charge

    // CURVED HIGH-ORDER (isoparametric P2) mode: m_curved sets m_highorder=true too (the QuadDot/PhiInner path
    // is shared); PhiInner -> PhiAtHO_Curved (always the curved Duffy).  m_cellNodes [n_cell*30] / m_faceNodes
    // [n_bf*18] hold the P2 high-order nodes; m_gl/m_gw the curved Duffy Gauss rule.  No analytic moments, no
    // inner-subtraction table, no near/far split (m_ho_far_factor stays 1e30).
    bool m_curved = false;
    int  m_curve_order = 0;
    std::vector<double> m_cellNodes, m_faceNodes;      // [n_cell*30] (P2 tet), [n_bf*18] (P2 tri)
    std::vector<int> m_cellVertices, m_faceVertices;   // [n_cell*4], [n_bf*3] global mesh vertex ids
    std::vector<double> m_gl, m_gw;                    // curved Duffy Gauss-Legendre rule on [0,1]
    // Symmetric reference-triangle rule for the base of each radial Duffy
    // sub-tet.  This removes reference-vertex-order dependence and reduces the
    // inner curved kernel from nq^3 to nq*ntri evaluations.
    std::vector<int> m_curvedTouchBlockIndex;          // canonical host-pair -> block slot, -1 if non-touching
    std::vector<std::vector<double>> m_curvedTouchBlocks; // precomputed symmetric touching blocks
    double m_curvedTouchBuildTime = 0.0;
    void PrecomputeCurvedTouchBlocks();
    bool CurvedTouchBlockValue(int kindA, int hostA, int localA,
                               int kindB, int hostB, int localB, double& value) const;

    // CURVED POLYTOPE mode (curved hex/wedge): both cell volume charges and boundary face charges follow
    // the curved P2 geometry.  m_srcCurvedTets[c] holds the 10-node P2 sub-tets for a cell charge;
    // m_srcCurvedTris[a] holds the 6-node P2 sub-tris for a face charge.  When m_curved_face is true,
    // PhiAt dispatches to CurvedTetPotential / CurvedTriPotential and m_qp/m_qw are mapped by the matching
    // curved outer quadrature rules.
    bool m_curved_face = false;
    std::vector<std::vector<std::array<rad_hdiv::Vec3, 6>>> m_srcCurvedTris;  // [n] curved-face P2 sub-tri nodes
    std::vector<std::vector<std::array<rad_hdiv::Vec3, 10>>> m_srcCurvedTets; // [n] curved-cell P2 sub-tet nodes

    // HEX RT1/RT2 mode (direct Q2 isoparametric geometry; see the hex ctor doc).  The charge monomial lives in
    // the HEX/QUAD REFERENCE frame (evaluated directly at ref coords -- no physical->ref inverse), geometry +
    // measure come from the Q2 lattice maps.  Entries are served block-wise (GetHexBlock); the inner is the
    // ref-frame radial decomposition (PhiInnerHexRadialVec) for near/self, cached far clouds otherwise.
    bool m_hexmode = false;
    bool m_hexCacheStatsEnabled = false;                  // opt-in block-cache stats; hot hits avoid atomics by default
    int  m_hex_n_bf = 0;
    std::vector<double> m_hexNodes, m_quadNodes;        // [n_el*81] 27-node Q2 hex, [n_bf*27] 9-node Q2 quad
    std::vector<double> m_symTetP, m_symTetW;           // regular outer tet rule (bary lam1..3; W sums 1/6)
    std::vector<double> m_symTriP, m_symTriW;           // regular outer tri rule (bary lam1..2; W sums 1/2)
    std::vector<double> m_glOut, m_gwOut;               // 1D [0,1] Gauss -> graded OUTER Duffy (near/self subs)
    std::vector<double> m_glIn, m_gwIn;                 // 1D [0,1] Gauss -> the RADIAL inner rule (PhiInnerHexRadialVec)
    std::vector<double> m_farTetP, m_farTetW;           // cheap FAR inner tet rule (bary; W sums 1/6)
    std::vector<double> m_farTriP, m_farTriW;           // cheap FAR inner tri rule (bary; W sums 1/2)
    double m_near_grade = 1.5, m_far_inner_factor = 4.0;
    std::vector<double> m_cellSubC, m_cellSubS, m_cellSubV;  // [n_el*6*3] sub-tet centroids, [n_el*6] sizes, [n_el*6*4*3] phys corners
    std::vector<double> m_faceSubC, m_faceSubS, m_faceSubV;  // [n_bf*2*3], [n_bf*2], [n_bf*2*3*3] (sub-tri)
    std::vector<unsigned char> m_hexAffineCell;              // [n_el] true when the Q2 lattice is affine
    int m_hexAffineOrder = 1, m_hexAffineMonoCount = 8, m_hexAffinePolyCount = 20;
    std::vector<double> m_hexAffineCoeff;                    // [n_el*mono*moment], RT1 8x20; RT2 27x84
    std::vector<unsigned char> m_quadAffineFace;             // [n_bf] true when the Q2 face lattice is affine
    int m_quadAffineMonoCount = 4, m_quadAffinePolyCount = 10;
    std::vector<double> m_quadAffineCoeff;                   // [n_bf*mono*moment], RT1 4x10; RT2 9x35
    bool m_hexUniformAffineCells = false;                    // same affine cell map for every cell -> translation block cache
    std::vector<int> m_hexCellLattice;                       // [n_el*3] integer lattice coordinate for uniform affine cells
    bool m_hexUniformTransHosts = false;                      // cell/face hosts are translated template copies
    std::vector<int> m_hexHostTemplate;                       // [n_el+n_bf] template id per host (cell ids and face ids are separate by kind)
    std::vector<int> m_hexHostLattice2;                       // [3*(n_el+n_bf)] half-cell lattice coordinates of host centers
    double HexMonoEval(int charge, const double xi[3]) const;   // ref-frame Q1/Q2 monomial
    // BLOCK-MEMO (the 64x co-location win): the near/far/grading decisions depend ONLY on host+sub geometry
    // (all co-located charges of a (kind,host) share m_cent/m_size), so the WHOLE directed host-pair block is
    // computed in ONE pass -- the 1/r sqrt is shared across all nT*nS monomial combos (the numpy-proto
    // block += wq*outer(Phia,inn) structure).  Bit-identical to per-entry QuadDotHex; ~64x fewer sqrt on near
    // hex-hex blocks.  m_hexLocalOf / m_cellCharges / m_faceCharges are the (kind,host)->local reverse maps.
    std::vector<int> m_hexLocalOf;                       // [n] local index of charge within its (kind,host) group
    std::vector<std::vector<int>> m_cellCharges;         // [n_el] global charge indices per cell (local order)
    std::vector<std::vector<int>> m_faceCharges;         // [n_bf] global charge indices per boundary face
    struct SolveTiming {
        double total_s = 0.0, factor_s = 0.0, prec_s = 0.0, bx_s = 0.0, gmatvec_s = 0.0;
        double btx_s = 0.0, mass_s = 0.0, dot_s = 0.0, ax_total_s = 0.0, ax_other_s = 0.0;
        double pcg_update_s = 0.0;
        double hmatvec_total_s = 0.0, hmatvec_zero_s = 0.0, hmatvec_permute_s = 0.0;
        double hmatvec_leaf_s = 0.0, hmatvec_reduce_s = 0.0, hmatvec_meta_s = 0.0;
        double hmatvec_lowrank_flop_est = 0.0, hmatvec_dense_flop_est = 0.0;
        double hmatvec_calls = 0.0, hmatvec_lowrank_leaves = 0.0, hmatvec_dense_leaves = 0.0;
        double hmatvec_mirrored_upper_leaves = 0.0, hmatvec_diagonal_leaves = 0.0;
        double hmatvec_skipped_lower_leaves = 0.0, hmatvec_last_nd = 0.0, hmatvec_last_nthr = 0.0;
        int apply_count = 0, prec_count = 0, dot_count = 0;
    };
    SolveTiming m_lastSolveTiming;
    // Persistent mass-Riesz PARDISO factor, keyed on the EXACT (n_face, mI, mJ, mV) COO arrays: reused by
    // SolveLinearMaterial reuses this when the mass is unchanged (any difference refactors), so
    // constant-mass iteration chains pay the analyze+factor ONCE.  Identical input -> identical factor ->
    // bit-identical preconditioner: timing-only.  shared_ptr keeps the .h to a forward declaration.
    std::shared_ptr<RadMassRieszCache> m_massRieszCache;
    std::shared_ptr<RadMassRieszCache> m_geometryMassRieszCache;
    std::vector<int> m_operatorBIndptr, m_operatorBIndices;
    std::vector<double> m_operatorBData;
    // Transposed charge map, built once at ConfigureChargeMap.  The production
    // operator applies B^T by row gather, avoiding an atomic scatter for every
    // H-matrix matvec and matching NGSolve's persistent BaseMatrix model.
    std::vector<int> m_operatorBTIndptr, m_operatorBTIndices;
    std::vector<double> m_operatorBTData;
    std::vector<int> m_operatorMassI, m_operatorMassJ;
    std::vector<double> m_operatorMassV;
    // Material mass in CSR, built once when the COO tangent is registered.
    // Krylov applications then use it directly instead of rebuilding sparse
    // topology at the start of every C++ solve.
    std::vector<int> m_operatorMassIndptr, m_operatorMassIndices;
    std::vector<double> m_operatorMassData;
    std::vector<int> m_operatorGeometryMassI, m_operatorGeometryMassJ;
    std::vector<double> m_operatorGeometryMassV;
    std::vector<int> m_operatorGeometryMassIndptr, m_operatorGeometryMassIndices;
    std::vector<double> m_operatorGeometryMassData;
    std::vector<unsigned char> m_operatorConstrained;
    int m_operatorNFace = 0;
    bool m_operatorChargeConfigured = false;
    bool m_operatorMassConfigured = false;
    bool m_operatorGeometryMassConfigured = false;
    // Get-or-build the persistent factor (the single shared implementation for both solve methods,
    // defined in the .cpp under HAVE_LAPACK).  Returns a PINNED shared_ptr the caller must hold for the
    // duration of its Krylov loop -- pinning makes a concurrent/nested replacement of the slot unable to
    // free a factor in use.  On a key MISS the old entry is released BEFORE the new factor is built, so
    // at most ONE factor is resident at any instant (the pre-cache peak-memory contract).  When
    // factor_s_accum is non-null it accumulates the fresh-factor wall time (a hit adds 0).
    std::shared_ptr<RadMassRieszCache> EnsureMassRieszFactor(
        const std::vector<int>& mI, const std::vector<int>& mJ, const std::vector<double>& mV,
        int n_face, const char* caller, double* factor_s_accum, bool geometry_cache = false);
    void PhiInnerHexSubVec(int kindS, int hS, int subB, const double p[3],
                           const std::vector<int>& srcG, double* inn) const;  // inner over ALL source locals (shares sqrt)
    void PhiInnerHexAffineCellSubVec(int hS, int subB, const double p[3],
                                     const std::vector<int>& srcG, double* inn) const;
    void PhiInnerHexAffineFaceSubVec(int hS, int subB, const double p[3],
                                     const std::vector<int>& srcG, double* inn) const;
    void PhiInnerHexAffineCellVec(int hS, const double p[3],
                                  const std::vector<int>& srcG, double* inn) const;
    void PhiInnerHexAffineFaceVec(int hS, const double p[3],
                                  const std::vector<int>& srcG, double* inn) const;
    std::vector<double> QuadBlockHexAffineFarProduct(
        int kindT, int hT, int kindS, int hS, int mask) const;
    // SELF inner by the tet path's PhiAtHO_Duffy RADIAL signed decomposition, ported to the REF frame:
    // anchor x0 = xiT, the outer point's OWN ref coords (the pulled-back kernel 1/|p-X(xi)| peaks there --
    // exact, no inverse), CLAMPED into the ref sub-simplex, then 4 signed radial sub-tets (3 signed
    // sub-tris on faces) from x0 with the Duffy apex AT x0: the u^2 (u) volume element kills the 1/r peak
    // exactly, and the map's warp enters only as a SMOOTH factor per quadrature point -- robust on strongly
    // distorted (|J| ratio ~0.4) and curved hexes alike, where the corner-graded-cloud / linearized-
    // subtraction schemes oscillated +-3% (eig 1.02-1.11 > 1 on the real Cubit cylinder mesh).  SELF-ONLY
    // since 2026-07-03: non-self near calls use the static-SITE radial below (the anchor-Newton branch was
    // removed with them).  m_glIn/m_gwIn is the radial 1D Gauss rule (n=5 -> 4*125 pts per cell call);
    // not cacheable here (x0 = xiT varies per outer point).
    void PhiInnerHexRadialVec(int kindS, int hS, int subB, const double p[3], const double* xiT,
                              const std::vector<int>& srcG, double* inn) const;
    // NON-SELF near inner (touch/shell, 2026-07-03): the SAME radial signed decomposition but anchored at
    // the nearest of a FIXED set of ref-space SITES (tet: 4 corners + 6 edge mids + 4 face centers +
    // centroid = 15; tri: 3+3+1 = 7).  The radial cone tiling is EXACT from ANY anchor -- the site only
    // aligns the grading with the kernel peak (p is OUTSIDE the sub-simplex here, r_min > 0), so nearest-
    // site anchoring costs only quadrature alignment, measured at the entry-drift/eig gates.  Because the
    // sites are fixed in REF space, EVERYTHING host-independent is precomputed per (sub, site) at ctor
    // time: the Q2 shape-value matrix S [nq x 27|9] (so the mapped nodes are X = S @ nodes -- no per-point
    // HexLag3), the Q1 monomial matrix M [nq x 8|4], and the signed Piola weights w [nq] (GW^3 u^2 v D,
    // orientation folded; cones whose base face contains the site are degenerate and skipped, so corner
    // sites carry 1 cone, edge 2, face-center 3, centroid 4 -- nq varies).  Per call this leaves ONE
    // nq x 27 "GEMV" + nq kernel evals: ~3-6x cheaper than the removed per-point Newton radial, with ZERO
    // per-host cloud memory (the old per-(host,sub,corner) cloud cache idea was rejected for its O(hosts)
    // thread_local footprint).
    struct HexSiteRad { int nq = 0; std::vector<double> S, M, w; };
    std::vector<HexSiteRad> m_cellSiteRad;               // [6*15] per (sub, site), ref-space tables
    std::vector<HexSiteRad> m_faceSiteRad;               // [2*7]
    std::vector<double> m_cellSiteX;                     // [n_el*6*15*3] mapped site positions (site pick)
    std::vector<double> m_faceSiteX;                     // [n_bf*2*7*3]
    void BuildHexSiteTables();                           // ctor helper: fills the four members above
    void PhiInnerHexSiteVec(int kindS, int hS, int subB, const double p[3],
                            const std::vector<int>& srcG, double* inn) const;
    std::vector<double> QuadBlockHexAffineProduct(int kindT, int hT, int kindS, int hS, int mask) const;
public:
    // Heap-stomp canary (2026-07-03 flake hunt): checksum over every hex-mode member array a block
    // computation reads, stored at ctor end.  A later mismatch proves the instance data was OVERWRITTEN
    // (the 0xc0000374-class heap corruption) rather than computed wrong.
    double HexStateChecksum() const;
    double HexStateCtorChecksum() const { return m_hex_state_sum; }
    std::vector<std::pair<std::string, double>> HexStateBreakdown() const;   // per-array (forensics)
    std::vector<std::pair<std::string, double>> HexCacheStats() const;       // block-cache hit/miss stats
    const std::vector<double>& HexStoredCellNodes() const { return m_hexNodes; }
    const std::vector<double>& HexStoredFaceNodes() const { return m_quadNodes; }
private:
    double m_hex_state_sum = 0.0;
    mutable std::atomic<long long> m_hexBlockLookups{0};
    mutable std::atomic<long long> m_hexBlockHits{0};
    mutable std::atomic<long long> m_hexBlockMisses{0};
    mutable std::atomic<long long> m_hexBlockClears{0};
    mutable std::atomic<long long> m_hexTransBlockLookups{0};
    mutable std::atomic<long long> m_hexTransBlockHits{0};
    mutable std::atomic<long long> m_hexTransBlockMisses{0};
    mutable std::atomic<long long> m_hexTransBlockClears{0};
    mutable std::atomic<long long> m_hexSymBlockLookups{0};
    mutable std::atomic<long long> m_hexSymBlockHits{0};
    mutable std::atomic<long long> m_hexSymBlockMisses{0};
    mutable std::atomic<long long> m_hexSymBlockClears{0};
    mutable std::atomic<long long> m_hexSymTransBlockLookups{0};
    mutable std::atomic<long long> m_hexSymTransBlockHits{0};
    mutable std::atomic<long long> m_hexSymTransBlockMisses{0};
    mutable std::atomic<long long> m_hexSymTransBlockClears{0};
    mutable std::atomic<long long> m_hoSymBlockLookups{0};
    mutable std::atomic<long long> m_hoSymBlockHits{0};
    mutable std::atomic<long long> m_hoSymBlockMisses{0};
    mutable std::atomic<long long> m_hoSymBlockClears{0};
    void ResetHexCacheStats();
    // mask (IMA): 0 = direct block; >0 = the mirror-image block (target host x the source host REFLECTED on
    // the 3-bit axis mask), for the reduced-symmetry (1/2,1/4,1/8) image method.  Default 0 keeps the direct
    // hex/wedge Gram byte-identical.
    std::vector<double> QuadBlockHex(int kindT, int hT, int kindS, int hS, int mask = 0) const;  // directed [nT*nS] block, INV4PI folded
    const std::vector<double>& GetHexBlock(int kindT, int hT, int kindS, int hS, int mask = 0) const;  // thread_local block cache
    const std::vector<double>& GetHexSymBlock(int kindA, int hA, int kindB, int hB, int mask = 0) const;  // cached 0.5*(AB+BA^T)

    // ---- 2D PLANAR mode (see the dim2 ctor doc) ----
    bool m_d2 = false;
    int  m_d2_n_be = 0;
    int  m_d2GeometryOrder = 1;
    int  m_d2CellMapStride = 0;
    int  m_d2EdgeMapStride = 0;
    std::vector<double> m_d2CellMap;     // polynomial map coefficients; Qq slot x 2D (tri uses Pq prefix)
    std::vector<int>    m_d2CellType;    // [n_el] 0=tri, 1=quad
    std::vector<double> m_d2EdgeMap;     // polynomial map coefficients; (q+1) x 2D
    std::vector<double> m_d2SymTriP, m_d2SymTriW;   // OUTER tri rule (bary lam1..2; W sums 1/2)
    std::vector<double> m_d2GlQ, m_d2GwQ;           // 1D [0,1] tensor outer rule for quads
    std::vector<double> m_d2GlE, m_d2GwE;           // 1D [0,1] edge outer rule
    std::vector<double> m_d2FarTriP, m_d2FarTriW;   // cheap FAR inner tri rule (bary; W sums 1/2)
    // per-sub geometry (cells: up to 4 sub-tris; edges: 1) -- centroid/size for the near test + anchor
    // SITES (tri sub: 3 corners + 3 edge mids + centroid = 7; edge: 2 ends + mid = 3), mapped positions.
    std::vector<double> m_d2CellSubC, m_d2CellSubS;   // [n_el*4*2], [n_el*4] (quad D4 centre fan)
    std::vector<double> m_d2EdgeC, m_d2EdgeS;         // [n_be*2],   [n_be]
    std::vector<double> m_d2CellSiteX;                // [n_el*4*7*2]
    std::vector<double> m_d2EdgeSiteX;                // [n_be*3*2]
    void D2CellMap(int cell_type, const double* coeff, const double xi[2], double X[2]) const;
    void D2EdgeMap(const double* coeff, double t, double X[2]) const;
    void D2EdgeTangent(const double* coeff, double t, double T[2]) const;
    // inner INT over sub subB of source (kindS,hS) of m_b(eta)*(-ln|p-X(eta)|) d(ref eta): radial cones
    // from the anchor (xiT = the outer point's own ref coords on the SELF host; else the nearest SITE)
    // for near field points, the cached-rule far cloud otherwise.
    void PhiInner2DVec(int kindS, int hS, int subB, const double p[2], const double* xiT,
                       const std::vector<int>& srcG, double* inn) const;
    std::vector<double> QuadBlock2D(int kindT, int hT, int kindS, int hS,
                                    int mask = 0) const;  // directed block, 1/(2pi) folded

    // ---- WEDGE (PRISM) RT1/RT2 mode (see the wedge ctor doc) ----  Reuses the hex-mode quadrature tables
    // (m_symTetP/W, m_symTriP/W, m_glOut/gwOut, m_glIn/gwIn, m_farTetP/W, m_farTriP/W, m_near_grade,
    // m_far_inner_factor) and the block-serving infra (m_hexLocalOf, m_cellCharges, m_faceCharges,
    // m_cent, m_size) verbatim; only the geometry (prism cells + mixed tri/quad faces) is wedge-specific.
    bool m_wedgemode = false;
    int  m_wedge_n_bf = 0;
    std::vector<double> m_wCellNodes;    // [n_el*54]  18-node tri-P2 (x) z-P2 prism lattice (n = t + 6*iz)
    std::vector<double> m_wFaceNodes;    // [n_bf*27]  9-node slots x 3D (a tri face fills the first 6 = 18)
    std::vector<int>    m_wFaceType;     // [n_bf] 0=tri (1 sub-tri, 6-node), 1=quad (2 sub-tris, 9-node)
    std::vector<double> m_wCellSubC, m_wCellSubS, m_wCellSubV;  // [n_el*3*3], [n_el*3], [n_el*3*4*3]  (3 sub-tets)
    std::vector<double> m_wFaceSubC, m_wFaceSubS, m_wFaceSubV;  // [n_bf*2*3], [n_bf*2], [n_bf*2*3*3]  (tri uses sub 0)
    std::vector<HexSiteRad> m_wCellSiteRad;     // [3*15]  per (sub-tet, site) ref-space radial tables
    std::vector<HexSiteRad> m_wFaceSiteRadTri;  // [1*7]   tri face: 1 sub-tri x 7 sites
    std::vector<HexSiteRad> m_wFaceSiteRadQuad; // [2*7]   quad face: 2 sub-tris x 7 sites
    std::vector<double> m_wCellSiteX;           // [n_el*3*15*3]  mapped site positions (nearest-site pick)
    std::vector<double> m_wFaceSiteX;           // [n_bf*2*7*3]   (tri face uses sub 0 only)
    void BuildWedgeSiteTables();                // ctor helper (fills the six members above)
    void PhiInnerWedgeSiteVec(int kindS, int hS, int subB, const double p[3],
                              const std::vector<int>& srcG, double* inn) const;   // non-self near: static-site radial
    void PhiInnerWedgeSubVec(int kindS, int hS, int subB, const double p[3],
                             const std::vector<int>& srcG, double* inn) const;    // far cloud / -> site radial
    void PhiInnerWedgeRadialVec(int kindS, int hS, int subB, const double p[3], const double* xiT,
                                const std::vector<int>& srcG, double* inn) const; // SELF exact-anchor radial
    std::vector<double> QuadBlockWedge(int kindT, int hT, int kindS, int hS, int mask = 0) const;

    // HIGH-ORDER (polynomial-charge) mode
    bool m_highorder = false;
    std::vector<int> m_host, m_kind, m_expo;           // [n] host elem, [n] 0=cell/1=face, [n*3] monomial exponents
    std::vector<int> m_nmono;                          // [n] # co-located charges per (kind,host) group -- QuadDot memo gating (skip groups of 1)
    std::vector<int> m_hoLocalOf;                      // [n] local charge index within the flat high-order host
    std::vector<std::vector<int>> m_hoCellCharges;     // [n_el] flat high-order charge indices per tetrahedron
    std::vector<std::vector<int>> m_hoFaceCharges;     // [n_bf] flat high-order charge indices per boundary triangle
    long long m_build_id = 0;                          // monotonic per-build id -> the QuadDot thread_local memo owner key (pointer-reuse-safe)
    std::vector<double> m_cellInv;                     // [n_el*9] physical->ref affine inverse per cell (row-major)
    std::vector<double> m_faceGinv;                    // [n_bf*4] 2x2 (a.a) Gram inverse per face (for 2D ref coords)
    std::vector<int> m_hoPolyDegree;                   // [n] physical-polynomial degree (flat order<=2)
    std::vector<double> m_hoPolyA, m_hoPolyB, m_hoPolyC; // [n], [n*3], [n*9]: A + B.y + y^T C y
    bool m_hoAnalyticBlock = false;                    // flat RT2: all charges covered and quadratic face modes present
    std::vector<std::vector<rad_hdiv::Vec3>> m_inP;    // [n] FIXED inner-potential Gauss points per HOST (cell/face)
    std::vector<std::vector<double>>          m_inW;   // [n] inner-potential Gauss weights (sum = host measure)
    std::vector<std::vector<double>>          m_srcval; // [n] PRECOMPUTED m_src(y_q) at the FIXED m_inP points -- bit-exact hoist of EvalMono out of the PhiAtHO inner loop (value depends only on (src,q))
    // near/far adaptive quadrature: LOW-quad tables for the cheap FAR plain double-Gauss (empty => disabled)
    double m_ho_far_factor = 1e30;                     // FAR if |c_a-c_b| > m_ho_far_factor*(size_a+size_b)
    std::vector<std::vector<rad_hdiv::Vec3>> m_qp_lo;  // [n] LOW-quad outer points (m_a folded into m_qw_lo)
    std::vector<std::vector<double>>          m_qw_lo; // [n] LOW-quad outer weights (monomial-folded)
    std::vector<std::vector<rad_hdiv::Vec3>> m_inP_lo; // [n] LOW-quad inner points (plain)
    std::vector<std::vector<double>>          m_inW_lo;// [n] LOW-quad inner weights (plain, NOT monomial-folded)
    std::vector<std::vector<double>>          m_srcval_lo; // [n] PRECOMPUTED m_src(y_q) at the FIXED m_inP_lo points (for QuadDotFar)
    double EvalMono(int charge, const double p[3]) const;   // charge's monomial at physical p (host ref-coord map)
    void InitHOPolynomialCoefficients();                    // flat order<=2 reference monomials -> physical A/B/C
    void PhiInnerHOHostVec(int kind, int host, const double p[3],
                           const std::vector<int>& charges, double* values) const;
    void PhiInnerHOCurvedHostVec(int kind, int host, const double p[3],
                                 const std::vector<int>& charges, double* values) const;
    bool CurvedHostsTouch(int kindA, int hostA, int kindB, int hostB) const;
    std::vector<double> QuadBlockHOCurvedDirect(int kindT, int hostT, int kindS, int hostS) const;
    std::vector<double> QuadBlockHOTet(int kindT, int hostT, int kindS, int hostS) const;
    const std::vector<double>& GetHOTetSymBlock(int kindA, int hostA, int kindB, int hostB) const;
    double PhiAtHO(int src, const double p[3]) const;       // polynomial-charge inner potential (subtraction, NEAR) -- superseded by PhiAtHO_Analytic for order<=2
    double PhiAtHO_Analytic(int src, const double p[3]) const; // EXACT analytic poly-charge potential (moment kernels, flat order<=2; machine precision, all pair types)
    double PhiAtHO_Duffy(int src, const double p[3]) const;    // Duffy singular-quadrature poly-charge potential (order>=3 / curved; ~1e-4)
    double PhiAtHO_Curved(int src, const double p[3]) const;   // CURVED isoparametric Duffy (rad_hdiv::CurvedTet/TriPotential at the host's P2 nodes)
    double PhiInner(int src, const double p[3]) const;        // dispatch: curved Duffy (m_curved) else analytic moments (charge deg<=2) else flat Duffy
    double QuadDotFar(int tgt, int src) const;              // cheap LOW-quad plain double-Gauss (FAR, no subtraction)
};

#endif // __RAD_HACAPK_HDIV_H
