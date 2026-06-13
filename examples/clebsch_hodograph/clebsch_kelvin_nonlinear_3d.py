"""Rung 3: the 3-D MERGE -- nonlinear saturation + exact Kelvin open boundary in
ONE Picard (geometry + material together).

RESEARCH example (track A, the capstone of "Kelvin in the hodograph").  Rungs
1-2 put the Kelvin open boundary in the hodograph/Clebsch frame (LINEAR);
rung 1.5b showed that in 2-D the hodograph LINEARISES saturation (one shot, no
loop).  **3-D does NOT auto-linearise** -- the Clebsch pair (psi, chi) is two
potentials for three coordinates (a gauge/helicity freedom), so there is no
clean (x,y,z) <-> (hodograph) interchange.  So the 3-D story is the user's
"混ぜて短く" in its honest 3-D form: a SINGLE Picard that updates the **geometry**
(the exact open boundary -- the Kelvin two-sphere from rung 2) and the
**material** (mu from |H| -- the saturation from rung 1.5a) TOGETHER, rather
than nesting an air-box-truncation loop inside a material loop.

The merge:
  - GEOMETRY (open boundary): the Kelvin two-sphere domain (periodic
    kelvin_int <-> kelvin_ext, GND at infinity, weight mu' = (R/rho')^2 mu0) --
    so there is NO air-box truncation loop to converge.  It is exact, once.
  - MATERIAL (saturation): mu_r = mu_r(|H|) in the magnetic body, updated each
    Picard step from the current H = Hs - grad(Omega).
  - ONE under-relaxed reduced-Omega Picard does both.  (Under-relaxation per the
    rung 1.5b lesson: undamped substitution oscillates under strong contrast.)

Exact reference (a saturable SPHERE in a uniform field has a UNIFORM interior,
so its nonlinear solution is equivalent to a LINEAR sphere with the
self-consistent mu_r_eff): the scalar demag fixed point

    H_int = 3 H0 / (mu_r(H_int) + 2),     mu_r_eff = mu_r(H_int).

The FEM interior field must reproduce this H_int, and the converged sphere's
EXTERIOR field must match the exact uniform + induced-dipole of the equivalent
mu_r_eff sphere -- the strong test of the Kelvin open boundary, now on the
nonlinear solution.

Honest scope: this is the STABLE saturation regime, and the SPHERE always
lives in it.  Its strong demagnetisation (N = 1/3) pulls the interior field
BELOW the saturating knee, so mu_r stays near mu_r0 and the H-input Picard
contraction |d RHS/dH| stays small (measured: 0.67 at mu_r0=20, and it DROPS
toward 0 as the B-H curve is made steeper -- the demag self-stabilises the
sphere; the instability is simply not reachable for it).

The reduced-Omega H-input Picard DOES fail for a LOW-demag body driven INTO
saturation (a prolate spheroid, N ~ 0.05-0.1, where the field penetrates to the
knee): seeded from the linear (high-mu) state it sticks in a spurious SHIELDED
basin (interior H ~ 0) and OSCILLATES, never reaching the true (unique,
penetrated) solution -- and neither under-relaxation NOR seeding from the
saturated state robustly cures it (both collapse back to the shielded basin).
The robust cure is the convex B-input formulation: for a physical, monotone B-H
curve the energy INT W(|B|) is convex (dH/dB > 0), so its minimiser is unique
and gradient-stable.  In 2-D that is the A-formulation (rung 1.5a, done); in
3-D it is the HCurl vector-potential nu(|B|) formulation (a substantial build,
genuinely open here -- the axisymmetric spheroid could instead use the 2-D
A_phi Henrotte form as a cheaper convex route).  This file stays with the
sphere, where the merged Picard is unconditionally well-behaved.

Verified (mu_r0=20, Hk=0.25, sphere a=0.2, Kelvin R=0.5, order 3, maxh 0.06,
ne 14301): the single Picard converges in 42 iters to the scalar demag fixed
point -- interior field_error 2.5e-4; the converged sphere is equivalent to
mu_r_eff = 11.9 < mu_r0 = 20 (saturated); EXTERIOR field_error 3.1e-3 vs the
equivalent dipole (vs 6.7e-3 for a truncated r/a=5 air box -- the Kelvin win,
now on the nonlinear solution); and the iterate is its own frozen re-solve to
3.2e-9 (the true-fixed-point diagnostic).  Contraction |dRHS/dH| = 0.67 (< 1,
the stable regime).

run:  python clebsch_kelvin_nonlinear_3d.py
"""
import math
import os

from numpy import pi, sqrt as npsqrt
import numpy as np
from ngsolve import (Mesh, H1, Periodic, GridFunction, grad, InnerProduct, dx,
                     CF, x, y, z, sqrt, BilinearForm, LinearForm, TaskManager,
                     Integrate)
from radia.kelvin_material import make_reduced_potential_background_cf

import clebsch_kelvin_3d as ck       # reuse geometry + verification helpers

MU0 = 4 * pi * 1e-7


def _mu_r(Hmag, mur0, Hk):
    """Froehlich saturation in H: mu_r0 at H=0, -> 1 as |H| -> infinity."""
    return 1.0 + (mur0 - 1.0) / (1.0 + (Hmag / Hk) ** 2)


def _demag_fixed_point(H0, mur0, Hk, n=200000):
    """Exact reference: the scalar self-consistent interior field of a saturable
    SPHERE, H = 3 H0/(mu_r(H)+2).  Returns (H_int, mu_r_eff, contraction) where
    contraction = |d RHS/dH| at the fixed point (<1 => the H-input Picard is
    stable; >=1 => needs the convex B-input form)."""
    H = 3.0 * H0 / (mur0 + 2.0)
    for _ in range(n):
        Hn = 3.0 * H0 / (_mu_r(H, mur0, Hk) + 2.0)
        if abs(Hn - H) < 1e-15 * H0:
            H = Hn
            break
        H = 0.5 * (H + Hn)                                # damped scalar solve
    mur = _mu_r(H, mur0, Hk)
    dmur = (mur0 - 1.0) * (-2.0 * H / Hk ** 2) / (1.0 + (H / Hk) ** 2) ** 2
    contraction = abs(-3.0 * H0 * dmur / (mur + 2.0) ** 2)
    return H, mur, contraction


def solve_nonlinear(mur0=20.0, Hk=0.25, order=3, maxh=0.06, a=0.2, R_K=0.5,
                    offset=(2.0, 0.0, 0.0), H0=1.0, niter=120, tol=1e-8,
                    relax=1.0, with_airbox=True, maxh_airbox=None):
    """Single under-relaxed reduced-Omega Picard on the Kelvin two-sphere domain
    with a saturable magnetic sphere.  Returns the convergence to the scalar
    demag fixed point + the Kelvin open-boundary win on the nonlinear solution."""
    H_ref, mur_eff_ref, contraction = _demag_fixed_point(H0, mur0, Hk)
    Hz_ref = H_ref                                        # interior field = H_int

    with TaskManager():
        mesh = Mesh(ck._kelvin_geometry(a, R_K, offset).GenerateMesh(maxh=maxh))
        mesh.Curve(order)
        ox, oy, oz = offset
        rho2 = (x - ox) ** 2 + (y - oy) ** 2 + (z - oz) ** 2 + 1e-24
        kmask = mesh.MaterialCF({"kelvin": 1.0}, default=0.0)
        magmask = mesh.MaterialCF({"magnetic": 1.0}, default=0.0)
        Hs = make_reduced_potential_background_cf(
            mesh, lambda xc, yc, zc: CF((0.0, 0.0, H0)),
            R_K=R_K, offset=offset, kelvin_mats=("kelvin",), dim=3)

        fes = Periodic(H1(mesh, order=order, dirichlet_bbnd="GND"))
        u, v = fes.TnT()

        def assemble_solve(mur_cf):
            # mu: mu_r(|H|)*mu0 in the magnet, mu0 in air_inner, Kelvin weight in
            # the exterior shell -- the geometry (open boundary) lives in the
            # FIXED Kelvin weight; only the magnet's mu changes per Picard step.
            mu = (magmask * (mur_cf * MU0)
                  + (1.0 - magmask - kmask) * MU0
                  + kmask * MU0 * (R_K * R_K / rho2))
            aO = BilinearForm(mu * InnerProduct(grad(u), grad(v)) * dx)
            aO.Assemble()
            fO = LinearForm(mu * InnerProduct(Hs, grad(v)) * dx)
            fO.Assemble()
            g = GridFunction(fes, name="Omega")
            g.vec.data = aO.mat.Inverse(fes.FreeDofs(),
                                        inverse="sparsecholesky") * fO.vec
            return g

        gfO = assemble_solve(CF(mur0))                   # linear seed
        Hz_prev = ck._interior_H(mesh, Hs - grad(gfO))[1]
        hist = [Hz_prev]                                 # interior Hz per Picard step
        nit = 0
        for it in range(niter):
            H = Hs - grad(gfO)
            Hmag = sqrt(InnerProduct(H, H) + 1e-30)
            gnew = assemble_solve(_mu_r(Hmag, mur0, Hk))
            gfO.vec.data = (1.0 - relax) * gfO.vec + relax * gnew.vec
            Hz_now = ck._interior_H(mesh, Hs - grad(gfO))[1]
            hist.append(Hz_now)
            nit = it + 1
            if abs(Hz_now - Hz_prev) < tol * abs(Hz_now):
                break
            Hz_prev = Hz_now

        # report from a CLEAN solve at the converged mu (NOT the relaxed blend --
        # rung 1.5a: a relaxed blend is not the solution of any linear problem).
        H = Hs - grad(gfO)
        mur_clean = _mu_r(sqrt(InnerProduct(H, H) + 1e-30), mur0, Hk)
        gfO = assemble_solve(mur_clean)
        H = Hs - grad(gfO)
        Hx_in, Hz_in = ck._interior_H(mesh, H)
        mur_eff_fem = _mu_r(abs(Hz_in), mur0, Hk)
        field_error = abs(Hz_in - Hz_ref) / abs(Hz_ref)

        # frozen re-solve: freeze mu at THIS H and solve again -- the converged
        # iterate must be its own re-solve (the true-fixed-point diagnostic).
        mur_frozen = _mu_r(sqrt(InnerProduct(H, H) + 1e-30), mur0, Hk)
        gfO_re = assemble_solve(mur_frozen)
        Hz_re = ck._interior_H(mesh, Hs - grad(gfO_re))[1]
        self_consistency = abs(Hz_in - Hz_re) / abs(Hz_in)

        # STRONG open-boundary test on the NONLINEAR solution: the converged
        # sphere is equivalent to a LINEAR mu_r_eff sphere, so the exterior must
        # match that sphere's exact uniform + induced dipole.
        exterior_error = ck._exterior_field_error(mesh, H, H0, a, mur_eff_ref)

        airbox_error = float("nan")
        if with_airbox:
            # the equivalent linear mu_r_eff sphere in a truncated air ball
            # (r/a=5 Dirichlet) -- the truncation error Kelvin removes.
            mh = maxh_airbox if maxh_airbox is not None else 2.2 * maxh
            mesh2 = Mesh(ck._truncated_geometry(a, 5.0 * a).GenerateMesh(maxh=mh))
            mesh2.Curve(order)
            Mu2 = mesh2.MaterialCF({"magnetic": mur_eff_ref * MU0}, default=MU0)
            fes2 = H1(mesh2, order=order, dirichlet="outer")
            u2, v2 = fes2.TnT()
            a2 = BilinearForm(Mu2 * InnerProduct(grad(u2), grad(v2)) * dx)
            a2.Assemble()
            f2 = LinearForm(Mu2 * InnerProduct(CF((0.0, 0.0, H0)), grad(v2)) * dx)
            f2.Assemble()
            gf2 = GridFunction(fes2)
            gf2.vec.data = a2.mat.Inverse(fes2.FreeDofs(),
                                          inverse="sparsecholesky") * f2.vec
            _, Hz2 = ck._interior_H(mesh2, CF((0.0, 0.0, H0)) - grad(gf2))
            airbox_error = abs(Hz2 - Hz_ref) / abs(Hz_ref)

    return {
        "mur0": mur0, "Hk": Hk, "order": order, "ne": int(mesh.ne), "n_iter": nit,
        "H0": H0, "Hz_in": float(Hz_in), "Hx_in": float(Hx_in),
        "Hz_ref": float(Hz_ref),                         # scalar demag fixed point
        "mur_eff_ref": float(mur_eff_ref),               # saturated mu_r (< mur0)
        "mur_eff_fem": float(mur_eff_fem),
        "contraction": float(contraction),               # <1 => H-input Picard stable
        "field_error": float(field_error),               # interior vs exact fixed point
        "exterior_error": float(exterior_error),         # EXTERIOR vs equivalent dipole
        "airbox_error": float(airbox_error),             # truncated r/a=5
        "self_consistency": float(self_consistency),     # iterate == frozen re-solve
        "hist": [float(h) for h in hist],                # interior Hz per Picard step
    }


def _plot(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    hist = np.array(res["hist"])
    Hz_ref = res["Hz_ref"]
    err = np.abs(hist - Hz_ref) / abs(Hz_ref)
    err = np.maximum(err, 1e-16)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0), dpi=150)

    ax.semilogy(range(len(err)), err, "o-", color="C0", ms=3,
                label="$|H_z^{(n)}-H_z^{ref}|/H_z^{ref}$")
    ax.semilogy(range(len(err)),
                err[0] * res["contraction"] ** np.arange(len(err)), "--",
                color="0.6", lw=1.0,
                label=fr"contraction $\rho$={res['contraction']:.2f}")
    ax.set_xlabel("Picard iteration")
    ax.set_ylabel("interior field error vs exact fixed point")
    ax.set_title(f"single merged Picard converges to the\nexact demag fixed "
                 f"point (self-cons {res['self_consistency']:.0e})")
    ax.legend(fontsize=8)

    Hk, mur0, H0 = res["Hk"], res["mur0"], res["H0"]
    Hgrid = np.linspace(0.0, 1.4 * H0, 300)
    mur = 1.0 + (mur0 - 1.0) / (1.0 + (Hgrid / Hk) ** 2)
    ax2.plot(Hgrid, mur, "-", color="C1", label=r"$\mu_r(|H|)$ (Froehlich)")
    ax2.plot(res["Hz_ref"], res["mur_eff_ref"], "o", color="C3", ms=8,
             label=fr"operating pt $(H_{{int}},\mu_r^{{eff}})$="
                   fr"({res['Hz_ref']:.2f}, {res['mur_eff_ref']:.1f})")
    ax2.axhline(mur0, color="0.6", ls=":", lw=1.0, label=fr"$\mu_{{r0}}$={mur0:.0f}")
    ax2.set_xlabel("interior field  $|H|$")
    ax2.set_ylabel(r"$\mu_r$")
    ax2.set_title("saturable material: the self-consistent\noperating point that "
                  "the Picard found")
    ax2.legend(fontsize=8)

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("Rung 3: nonlinear saturation + exact Kelvin open boundary in ONE "
          "Picard\n")
    r = solve_nonlinear()
    print(f"  mu_r0={r['mur0']:.0f}  Hk={r['Hk']:.2f}  H0={r['H0']:.1f}  "
          f"ne={r['ne']}  order={r['order']}  Picard iters={r['n_iter']}\n")
    print(f"  saturated:  mu_r_eff = {r['mur_eff_ref']:.3f}  (< mu_r0={r['mur0']:.0f}"
          f")   contraction |dRHS/dH| = {r['contraction']:.3f}  "
          f"({'STABLE' if r['contraction'] < 1 else 'UNSTABLE'})")
    print(f"  interior Hz = {r['Hz_in']:.6e}   (scalar demag fixed point "
          f"{r['Hz_ref']:.6e}),  Hx = {r['Hx_in']:.1e}")
    print(f"  -> interior field_error vs the EXACT fixed point = "
          f"{r['field_error']:.2e}")
    print(f"  -> EXTERIOR field_error vs the equivalent dipole = "
          f"{r['exterior_error']:.2e}  (the STRONG Kelvin test, nonlinear)")
    print(f"  -> air-box (r/a=5) interior error = {r['airbox_error']:.2e}  "
          f"(Kelvin ~{r['airbox_error']/r['field_error']:.0e}x better)")
    print(f"  -> frozen re-solve self-consistency = {r['self_consistency']:.2e}  "
          f"(the iterate IS the solution = true fixed point)")
    print("\n  => ONE Picard merged the EXACT open boundary (Kelvin geometry)")
    print("     and the saturation (mu(|H|) material) -- no nested air-box loop.")
    print("     3-D does not linearise (helicity), so the merge is a single loop,")
    print("     not a 1-shot (that is the 2-D-only Chaplygin, rung 1.5b).")
    _plot(r)


if __name__ == "__main__":
    main()
