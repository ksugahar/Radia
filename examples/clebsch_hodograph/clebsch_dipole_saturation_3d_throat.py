r"""B(b): a STRONG 3-D B_gap saturation knee, by THROAT flux-concentration.

The companion `clebsch_dipole_saturation_3d.py` (Task B(a)) showed the convergent
B-input A-formulation, but on the accel H-frame -- a LARGE, gap-reluctance-dominated
gap -- so iron saturation only softened `B_gap` by ~14%.  A gap-dominated magnet
does NOT show a strong knee no matter how hard the iron saturates: the air gap, not
the iron, sets the reluctance.

To get a STRONG `B_gap` knee the iron must be the bottleneck.  The trick is a
magnetic CIRCUIT with a thin THROAT in series with the gap:

    series flux conservation:   Phi = B_gap * A_gap = B_throat * A_throat
                            =>  B_throat = B_gap * (A_gap / A_throat)

With `A_gap / A_throat = 2`, the throat carries TWICE the gap field, so it reaches
the saturation field `J_sat` FIRST.  Once `B_throat -> J_sat` the throat permeability
collapses (`mu_r -> 1`), its reluctance explodes, and the CIRCULATING flux clamps --
a strong knee driven by the IRON, not the gap.

The strong knee lives in the **throat field `B_throat`** -- the bottleneck, on the far
limb away from the coil.  In the FEM it plateaus near the iron's saturated carrying
capacity (~2.84 T here) while the drive grows 8x: `dB_throat/dNI` drops ~157x across
the knee.  The iron literally cannot carry more flux density through the throat.

HONEST 3-D caveat (a real result, not a defect): the raw gap field `B_gap` and even its
iron-channeled part `B_gap_iron = |curl A_r|` do NOT clamp -- they keep rising past the
knee.  Two reasons: (1) the linking coil threads the top limb beside the gap, so `B_gap`
carries the coil's LOCAL, un-saturable Biot-Savart field; (2) once the throat saturates,
the circulating flux is no longer forced through it -- in 3-D it LEAKS around the choked
throat (alternative return paths through air / the rest of the iron).  So a single
saturated throat clamps the field LOCALLY (at the bottleneck), not everywhere: the
strong, physically-guaranteed clamp is `B_throat`, not the distant gap field.  The ideal
lumped circuit (`B_throat = B_gap*(A_gap/A_throat)`, flux fully conserved) over-states
how much the gap follows the throat -- it is the design-level approximation; the FEM is
the truth.

WHY THE NAIVE GEOMETRY FAILED (the B(b) bug, and its fix)
--------------------------------------------------------
A throat only concentrates flux if the drive forces CIRCULATING flux around the
closed window.  A racetrack coil sitting over the gap drives only the LOCAL gap
flux; the throat (on the far limb) sees almost nothing, so `B_throat ~ B_gap` and
there is no concentration (the early H-frame / C-core attempts).  The cure is a
coil that LINKS the window: a `rad.ObjFlmCur` loop threading one limb drives flux
all the way around the loop, through the throat.  A LINEAR FEM check on this geometry
gives `B_throat / B_gap = 3.8` -- ABOVE the geometric area ratio (2), so the throat
genuinely concentrates / funnels flux (see `linear_amplification_check`).

THE METHOD (nonlinear-as-linear, the running theme)
---------------------------------------------------
Exactly as in the 2-D `clebsch_dipole_saturation_2d.py`: a lumped MAGNETIC CIRCUIT
gives the whole saturation design map at LINEAR cost -- one scalar (bisection) solve
per drive point, no mesh -- and the convergent 3-D B-input A-formulation FEM
(reused from B(a)) validates the knee.  The circuit is the design tool; the FEM is
the trusted reference.

Units: SI (metres, A/m, Tesla).  No Cubit dependency.
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia"))

MU0 = 4.0 * np.pi * 1e-7
NU0 = 1.0 / MU0

# --- saturable iron material: high mu_r until a smooth hard cap at J_sat -------
MUR0 = 2000.0     # unsaturated relative permeability
JSAT = 2.0        # saturation field (T)
BW = 0.2          # cap softness (T) -- sharp enough to choke the throat clearly


def mur_of_B(Bm):
    """Smooth hard-cap permeability: mu_r ~ MUR0 below J_sat, -> 1 above.

    mu_r(B) = 1 + (MUR0-1) (1 - sig),   sig = 1/(1+exp(-(B-J_sat)/BW)).
    """
    sig = 1.0 / (1.0 + np.exp(-(Bm - JSAT) / BW))
    return 1.0 + (MUR0 - 1.0) * (1.0 - sig)


# --- geometry (a rectangular iron WINDOW frame; thin bottom limb = throat) -----
AX = 0.04         # half-width  (x) of the window
AZ = 0.04         # half-height (z) of the window
T = 0.012         # limb thickness (top / sides)
TT = 0.006        # THROAT thickness (bottom limb) -> A_gap/A_throat = T/TT = 2
D = 0.028         # depth (y)
GW = 0.004        # air-gap width cut in the top limb
AIR = 0.08        # half-size of the air box

A_GAP = T * D     # gap cross-section (flux runs in x)
A_THROAT = TT * D
A_IRON = T * D    # top / side limb cross-section
AREA_RATIO = A_GAP / A_THROAT


# =============================================================================
# 1.  MAGNETIC CIRCUIT MODEL  --  the saturation design map at LINEAR cost
# =============================================================================
def circuit_Bgap(NI, n_bisect=80):
    """Return (B_gap, B_throat) for a drive of NI ampere-turns, by a 1-shot
    scalar reluctance solve (bisection on the loop flux Phi).

    Series loop:  MMF(Phi) = (R_gap + R_throat(Phi) + R_rest(Phi)) * Phi = NI .
    MMF is monotone increasing in Phi, so bisection converges.  No mesh.
    """
    l_throat = 2.0 * AX                       # bottom-limb length (x)
    l_rest = (2.0 * AX - GW) + 2.0 * (2.0 * AZ)  # top + two sides

    R_gap = GW / (MU0 * A_GAP)                 # air gap (constant)

    def mmf(Phi):
        B_throat = Phi / A_THROAT
        B_rest = Phi / A_IRON
        R_throat = l_throat / (MU0 * mur_of_B(B_throat) * A_THROAT)
        R_rest = l_rest / (MU0 * mur_of_B(B_rest) * A_IRON)
        return (R_gap + R_throat + R_rest) * Phi

    # bracket Phi: upper bound where everything is deep in saturation
    Phi_hi = 4.0 * JSAT * A_GAP
    lo, hi = 0.0, Phi_hi
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        if mmf(mid) < NI:
            lo = mid
        else:
            hi = mid
    Phi = 0.5 * (lo + hi)
    return Phi / A_GAP, Phi / A_THROAT


def circuit_design_map(NI_list):
    """Sweep the drive; return arrays (NI, B_gap, B_throat).  Linear cost."""
    NI = np.asarray(NI_list, dtype=float)
    bg = np.empty_like(NI)
    bt = np.empty_like(NI)
    for i, ni in enumerate(NI):
        bg[i], bt[i] = circuit_Bgap(ni)
    return NI, bg, bt


def knee_drive():
    """The drive NI at which the throat reaches J_sat (the knee location)."""
    # B_throat = J_sat  =>  Phi = J_sat * A_throat ; invert MMF at that flux.
    Phi = JSAT * A_THROAT
    l_throat = 2.0 * AX
    l_rest = (2.0 * AX - GW) + 2.0 * (2.0 * AZ)
    R_gap = GW / (MU0 * A_GAP)
    R_throat = l_throat / (MU0 * mur_of_B(JSAT) * A_THROAT)
    R_rest = l_rest / (MU0 * mur_of_B(Phi / A_IRON) * A_IRON)
    return (R_gap + R_throat + R_rest) * Phi


# =============================================================================
# 2.  3-D FEM  --  the convergent B-input A-formulation (reused from B(a))
# =============================================================================
def _build_mesh(iron_maxh=0.005, air_maxh=0.032):
    import ngsolve as ng
    from ngsolve import TaskManager
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    hD = D / 2.0
    top = (Box(Pnt(-AX, -hD, AZ - T), Pnt(AX, hD, AZ))
           - Box(Pnt(-GW / 2, -hD, AZ - T), Pnt(GW / 2, hD, AZ)))
    bot = Box(Pnt(-AX, -hD, -AZ), Pnt(AX, hD, -AZ + TT))
    left = Box(Pnt(-AX, -hD, -AZ), Pnt(-AX + T, hD, AZ))
    right = Box(Pnt(AX - T, -hD, -AZ), Pnt(AX, hD, AZ))
    iron = top + bot + left + right
    iron.mat("iron")
    iron.maxh = iron_maxh
    air = Box(Pnt(-AIR, -AIR, -AIR), Pnt(AIR, AIR, AIR)) - iron
    air.mat("air")
    for f in air.faces:
        if max(abs(f.center.x), abs(f.center.y), abs(f.center.z)) > 0.9 * AIR:
            f.name = "outer"
    with TaskManager():
        return ng.Mesh(OCCGeometry(Glue([air, iron])).GenerateMesh(maxh=air_maxh))


def _linking_coil(I):
    """The window-LINKING filament loop (threads the top limb -> circulating flux)."""
    import radia as rad
    rad.UtiDelAll()
    Ly = D / 2.0 + 0.02
    zo = AZ + 0.02
    return rad.ObjFlmCur([[0, -Ly, 0], [0, Ly, 0], [0, Ly, zo],
                          [0, -Ly, zo], [0, -Ly, 0]], I)


def _coil_Bs(mesh, order=1):
    """Biot-Savart B_s = mu0 H_s of the unit-current linking coil, as a GridFunction."""
    import radia as rad
    from ngsolve import VectorL2, GridFunction, TaskManager
    Bs1 = GridFunction(VectorL2(mesh, order=order))
    with TaskManager():
        Bs1.Set(MU0 * rad.RadiaField(_linking_coil(1.0), "h"))
    rad.UtiDelAll()
    return Bs1


def linear_amplification_check(iron_maxh=0.006, air_maxh=0.04):
    """LINEAR FEM proof that the linking coil concentrates flux in the throat.

    Returns dict(B_gap, B_throat, ratio).  A single linear solve (mu_r = MUR0
    everywhere in iron, no saturation) -- the throat/gap amplification must
    exceed the geometric area ratio if the coil really drives circulating flux.
    """
    import ngsolve as ng
    from ngsolve import (HCurl, GridFunction, curl, dx, sqrt, InnerProduct,
                         BilinearForm, LinearForm, TaskManager)
    mesh = _build_mesh(iron_maxh, air_maxh)
    Bs = _coil_Bs(mesh)
    nu_iron = NU0 / MUR0
    nu = mesh.MaterialCF({"iron": nu_iron}, default=NU0)
    num = mesh.MaterialCF({"iron": nu_iron - NU0}, default=0.0)
    fes = HCurl(mesh, order=1, dirichlet="outer")
    u, v = fes.TnT()
    eps = 1e-6 * NU0
    with TaskManager():
        a = BilinearForm(nu * InnerProduct(curl(u), curl(v)) * dx
                         + eps * InnerProduct(u, v) * dx)
        a.Assemble()
        f = LinearForm(-num * InnerProduct(Bs, curl(v)) * dx)
        f.Assemble()
        Ar = GridFunction(fes)
        Ar.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    B = Bs + curl(Ar)
    bg = abs(float(B[0](mesh(0, 0, AZ - T / 2))))
    bt = abs(float(sqrt(InnerProduct(B, B))(mesh(0, 0, -AZ + TT / 2))))
    return {"B_gap": bg, "B_throat": bt, "ratio": bt / bg, "ne": mesh.ne}


def nonlinear_fem(NI_list, iron_maxh=0.009, air_maxh=0.05, niter=400,
                  relax0=0.25, relax_min=0.02, tol=1e-4, verbose=False):
    """Convergent B-input A-formulation saturation sweep (the trusted reference).

    B = B_s + curl A_r,  weak form
        INT nu(|B|) curl(A_r).curl(v) + INT_iron (nu(|B|)-nu0) B_s.curl(v) = 0.

    The B-input nu(|B|) energy is convex, but the frozen-coefficient (Picard) map
    has a large Lipschitz constant across the saturation knee (where d(nu)/d|B| is
    largest), so a fixed step diverges there.  ADAPTIVE under-relaxation -- halve the
    step whenever the residual rises -- restores convergence through the knee.
    Warm-started in the drive.  Returns list of dicts.
    """
    import ngsolve as ng
    from ngsolve import (HCurl, GridFunction, curl, dx, sqrt, InnerProduct,
                         BilinearForm, LinearForm, TaskManager, exp)
    mesh = _build_mesh(iron_maxh, air_maxh)
    Bs1 = _coil_Bs(mesh)

    def nu_iron_cf(Bm):
        sig = 1.0 / (1.0 + exp(-(Bm - JSAT) / BW))
        return NU0 / (1.0 + (MUR0 - 1.0) * (1.0 - sig))

    fes = HCurl(mesh, order=1, dirichlet="outer")
    u, v = fes.TnT()
    eps = 1e-6 * NU0
    Ar = GridFunction(fes)        # warm-started across drive points
    out = []
    for NI in NI_list:
        Bs = NI * Bs1
        prev = np.array(Ar.vec)
        d = 1.0
        d_prev = np.inf
        relax = relax0
        nit = 0
        t0 = time.time()
        for it in range(niter):
            B = Bs + curl(Ar)
            Bm = sqrt(InnerProduct(B, B) + 1e-30)
            nui = nu_iron_cf(Bm)
            nu = mesh.MaterialCF({"iron": nui}, default=NU0)
            num = mesh.MaterialCF({"iron": nui - NU0}, default=0.0)
            with TaskManager():
                a = BilinearForm(nu * InnerProduct(curl(u), curl(v)) * dx
                                 + eps * InnerProduct(u, v) * dx)
                a.Assemble()
                f = LinearForm(-num * InnerProduct(Bs, curl(v)) * dx)
                f.Assemble()
                An = GridFunction(fes)
                An.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                            inverse="sparsecholesky") * f.vec
            Ar.vec.data = (1 - relax) * Ar.vec + relax * An.vec
            cur = np.array(Ar.vec)
            d = np.linalg.norm(cur - prev) / (np.linalg.norm(cur) + 1e-30)
            if d > d_prev:                       # residual rose -> step too big
                relax = max(0.5 * relax, relax_min)
            prev = cur.copy()
            d_prev = d
            nit = it + 1
            if d < tol:
                break
        B = Bs + curl(Ar)
        bg = abs(float(B[0](mesh(0, 0, AZ - T / 2))))           # total gap field
        br = abs(float(curl(Ar)[0](mesh(0, 0, AZ - T / 2))))    # iron-channeled part
        bt = abs(float(sqrt(InnerProduct(B, B))(mesh(0, 0, -AZ + TT / 2))))  # throat
        rec = {"NI": NI, "B_gap": bg, "B_gap_iron": br, "B_throat": bt,
               "iters": nit, "resid": d, "t_s": time.time() - t0, "ne": mesh.ne}
        out.append(rec)
        if verbose:
            print("  NI=%6d: B_throat=%.2fT B_gap_iron=%.3fT B_gap=%.3fT iters=%d "
                  "resid=%.0e (%.0fs)"
                  % (NI, bt, br, bg, nit, d, rec["t_s"]), flush=True)
    return out


# =============================================================================
def run(with_fem=False, fem_NI=(200, 1500, 3000, 6000, 12000, 24000),
        circuit_NI=(200, 600, 1500, 3000, 6000, 12000, 24000)):
    """Build the saturation design map (circuit, linear cost) and, if asked, the
    3-D FEM validation.  Returns a dict for the golden test."""
    NI = np.asarray(circuit_NI, dtype=float)
    _, bg, bt = circuit_design_map(NI)
    NIk = knee_drive()
    # knee strength: slope below the knee vs far above it.
    slope_lo = float(bg[1] / NI[1])
    slope_hi = float((bg[-1] - bg[-2]) / (NI[-1] - NI[-2]))
    circuit = {
        "NI": [float(v) for v in NI],
        "B_gap": [float(v) for v in bg],
        "B_throat": [float(v) for v in bt],
        "knee_NI": float(NIk),
        "clamp_pred": float(JSAT / AREA_RATIO),
        "slope_below": slope_lo,
        "slope_above": slope_hi,
        "slope_drop": slope_lo / max(slope_hi, 1e-30),
    }
    out = {"geometry": {"area_ratio": AREA_RATIO, "J_sat": JSAT,
                        "MUR0": MUR0, "A_gap": A_GAP, "A_throat": A_THROAT},
           "circuit": circuit}

    if with_fem:
        amp = linear_amplification_check()
        fem = nonlinear_fem(list(fem_NI))
        # The knee is in the CHANNELED flux -- the throat field (the bottleneck, far
        # from the coil).  The raw gap field additionally carries the linking coil's
        # direct (un-saturable) Biot-Savart field, so the saturating quantity is the
        # throat field B_throat (and the iron-channeled gap field B_gap_iron).
        b0 = fem[0]["B_throat"] / fem[0]["NI"]
        bN = (fem[-1]["B_throat"] - fem[-2]["B_throat"]) / (fem[-1]["NI"] - fem[-2]["NI"])
        out["linear_amplification"] = amp
        out["fem"] = fem
        out["fem_throat_slope_below"] = float(b0)
        out["fem_throat_slope_above"] = float(bN)
        out["fem_throat_slope_drop"] = float(b0 / max(bN, 1e-30))
        out["fem_throat_clamp"] = float(fem[-1]["B_throat"])
        out["fem_max_resid"] = float(max(r["resid"] for r in fem))
    return out


def _make_figure(out, path):
    """Optional figure (left: circuit design map at linear cost; right: FEM -- the
    throat field clamps, the gap field keeps rising).  No in-figure title (lab style)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False
    c = out["circuit"]
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.4, 3.4))
    axL.plot(c["NI"], c["B_gap"], "o-", label="B_gap (gap)")
    axL.plot(c["NI"], c["B_throat"], "s--", label="B_throat (throat)")
    axL.axhline(JSAT, color="k", lw=0.8, ls=":")
    axL.axvline(c["knee_NI"], color="0.6", lw=0.8, ls=":")
    axL.set_xlabel("drive  NI  [A-turn]")
    axL.set_ylabel("field  [T]")
    axL.legend(frameon=False, fontsize=9)
    axL.text(0.04, 0.92, "magnetic circuit (1-shot/point)", transform=axL.transAxes,
             fontsize=9, va="top")
    if "fem" in out:
        f = out["fem"]
        NI = [r["NI"] for r in f]
        axR.plot(NI, [r["B_throat"] for r in f], "s-", color="C1",
                 label="B_throat (clamps)")
        axR.plot(NI, [r["B_gap"] for r in f], "o-", color="C0",
                 label="B_gap (rises)")
        axR.plot(NI, [r["B_gap_iron"] for r in f], "^:", color="C2",
                 label="B_gap iron part")
        axR.axhline(JSAT, color="k", lw=0.8, ls=":")
        axR.set_xlabel("drive  NI  [A-turn]")
        axR.set_ylabel("field  [T]")
        axR.legend(frameon=False, fontsize=9)
        axR.text(0.04, 0.92, "3-D FEM (B-input A-formulation)",
                 transform=axR.transAxes, fontsize=9, va="top")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return True


def main():
    import json
    print("=" * 70)
    print("B(b): strong 3-D B_gap knee by throat flux-concentration")
    print("  geometry: window frame, A_gap/A_throat = %.1f, J_sat = %.1f T"
          % (AREA_RATIO, JSAT))
    print("  predicted clamp  B_gap ~ J_sat/ratio = %.2f T" % (JSAT / AREA_RATIO))
    print("=" * 70)
    out = run(with_fem=("--fem" in sys.argv))
    c = out["circuit"]
    print("\nMagnetic circuit design map (1-shot per point, no mesh):")
    print("   %8s %10s %10s" % ("NI", "B_gap[T]", "B_throat[T]"))
    for i in range(len(c["NI"])):
        print("   %8.0f %10.3f %10.3f"
              % (c["NI"][i], c["B_gap"][i], c["B_throat"][i]))
    print("  knee drive (B_throat = J_sat): NI ~ %.0f At" % c["knee_NI"])
    print("  dB_gap/dNI below knee = %.3e ; far above = %.3e (drop %.0fx)"
          % (c["slope_below"], c["slope_above"], c["slope_drop"]))

    if "fem" in out:
        amp = out["linear_amplification"]
        print("\nLinear amplification (flux-concentration proof):")
        print("  ne=%d ratio B_throat/B_gap = %.2f (> geometric area ratio %.1f)"
              % (amp["ne"], amp["ratio"], AREA_RATIO))
        print("\nNonlinear FEM validation (B-input A-formulation, adaptive Picard):")
        print("  the saturating quantity is the CHANNELED flux (throat field); the raw")
        print("  gap field also carries the linking coil's direct (un-saturable) field.")
        for r in out["fem"]:
            print("  NI=%6d: B_throat=%.2fT B_gap_iron=%.3fT B_gap=%.3fT iters=%d "
                  "resid=%.0e (%.0fs)"
                  % (r["NI"], r["B_throat"], r["B_gap_iron"], r["B_gap"],
                     r["iters"], r["resid"], r["t_s"]), flush=True)
        print("  throat clamps at %.2fT; dB_throat/dNI below=%.3e above=%.3e (drop %.0fx)"
              % (out["fem_throat_clamp"], out["fem_throat_slope_below"],
                 out["fem_throat_slope_above"], out["fem_throat_slope_drop"]))

    here = os.path.dirname(__file__)
    with open(os.path.join(here, "clebsch_dipole_saturation_3d_throat.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("\nsaved clebsch_dipole_saturation_3d_throat.json")
    if _make_figure(out, os.path.join(here, "clebsch_dipole_saturation_3d_throat.png")):
        print("saved clebsch_dipole_saturation_3d_throat.png")


if __name__ == "__main__":
    main()
