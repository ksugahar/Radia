"""ellipsoid_alpha_omega_axisym.py -- full-frequency AXIAL eddy-current
polarizability alpha_c(omega) of a conducting spheroid, by an axisymmetric
FEM solve, validated against the analytic sphere and anchored at high
frequency by the analytic demag tensor.

This is the full-frequency (eddy) leg of the shape-anisotropy story.  The
analytic anchors live in ellipsoid_alpha_tensor.py (DC alpha=0, HF perfect-
conductor alpha_i=-V/(1-N_i)); here we compute the WHOLE alpha_c(omega)
curve between them with a FEM eddy-current solve, and show it is
shape-dependent.

Method (reuses the lab-verified TEAM 28 mixed phi-B axisymmetric machinery,
team28_axisym_fem.py): a conducting body of revolution sits in a UNIFORM
axial AC field imposed as the Dirichlet condition phi = B0 r^2/2 on the far
boundary (phi = r A_phi, and A_phi = B0 r/2 is the uniform-field vector
potential).  Solve (K + s N) X with that inhomogeneous BC; the induced
axial magnetic moment is

    m_z = (1/2) integral (r x J)_z dV = pi integral r^2 J_phi dr dz
        = -pi s sigma integral r phi dr dz   (J_phi = -s sigma A_phi = -s sigma phi/r)

and the axial polarizability is alpha_c(omega) = m_z mu0 / B0 (so that for a
sphere alpha_c = 4 pi a^3 G(x), the Landau response in maglev_sphere_force).

VALIDATION: for a sphere the FEM alpha_c(omega) must reproduce the analytic
4 pi a^3 G(x) -- this checks the uniform-field BC and the moment extractor.
Then prolate/oblate spheroids show alpha_c(omega) is shape-dependent and ->
the analytic HF limit -V/(1-N_c).

Run:  python ellipsoid_alpha_omega_axisym.py   (~1-2 min; a few axisym solves)
"""
import json
import os
import sys

import numpy as np
from netgen.occ import WorkPlane, MoveTo, Glue, OCCGeometry, Axis, Z
from ngsolve import (
    Mesh, H1, VectorL2, GridFunction, BilinearForm, Integrate, CF, x, dx,
    grad, TaskManager, ngsglobals,
)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "sphere"))  # maglev_sphere_force
from maglev_sphere_force import G_exact, delta  # noqa: E402
import ellipsoid_alpha_tensor as ET  # noqa: E402
from _validation_output import validation_output  # noqa: E402

mu0 = 4 * np.pi * 1e-7
SIGMA = 5.8e7          # copper
A_SPH = 5e-3           # equatorial semi-axis used everywhere (r-direction)
B0 = 1.0               # imposed uniform axial field amplitude (alpha is linear in B0)


def build_mesh(a_eq, c, box=0.10, maxh=0.01, maxh_cond=None):
    """Axisymmetric (r,z) mesh: a half-spheroid conductor (semi-axes a_eq in
    r, c in z) centred at the origin, inside a large air box [0,box]x[-box,box].
    Far boundary named 'outer' (Dirichlet phi=B0 r^2/2 imposed there)."""
    if maxh_cond is None:
        maxh_cond = min(a_eq, c) / 6.0
    boxf = MoveTo(0, -box).Rectangle(box, 2 * box).Face()
    boxf.edges.name = "outer"
    # OCC requires major >= minor.  Oblate/sphere (a_eq >= c): major = a_eq
    # along r.  Prolate (c > a_eq): build major = c along the default x, then
    # rotate 90 deg so the long axis lands on z.
    if c > a_eq:
        ell = WorkPlane().MoveTo(0, 0).Ellipse(c, a_eq).Face()
        ell = ell.Rotate(Axis((0, 0, 0), Z), 90)
    else:
        ell = WorkPlane().MoveTo(0, 0).Ellipse(a_eq, c).Face()
    leftrect = MoveTo(-a_eq, -c).Rectangle(a_eq, 2 * c).Face()   # r<0 half
    cond = ell - leftrect
    cond.faces.name = "cond"
    cond.maxh = maxh_cond
    air = boxf - cond
    air.faces.name = "air"
    shape = Glue([air, cond])
    return Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh))


def alpha_c(mesh, omega, p=2):
    """Axial polarizability alpha_c(omega) [m^3] of the 'cond' body in a unit
    uniform axial AC field, via the mixed phi-B axisym solve."""
    fesPhi = H1(mesh, order=p, dirichlet="outer", complex=True)
    fesB = VectorL2(mesh, order=p - 1, complex=True)
    fes = fesPhi * fesB
    (u, B), (v, dB) = fes.TnT()

    def Bop(phi):
        g = grad(phi)
        return CF((-1 / x * g[1], 1 / x * g[0]))

    nu = 1.0 / mu0                                    # non-magnetic everywhere
    sig = mesh.MaterialCF({"cond": SIGMA}, default=0.0)
    s = 1j * omega
    a = BilinearForm(nu * x * (B * dB - Bop(u) * dB + Bop(v) * B) * dx)
    a += (v * (s * sig * (u / x))) * dx
    with TaskManager():
        a.Assemble()
        gf = GridFunction(fes)
        # inhomogeneous Dirichlet: phi = B0 r^2/2 on 'outer' (axis r=0 -> 0)
        gf.components[0].Set(CF(B0 * x * x / 2), definedon=mesh.Boundaries("outer"))
        res = gf.vec.CreateVector()
        res.data = -a.mat * gf.vec
        gf.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * res
        # induced axial moment m_z = -pi s sigma int r phi dr dz over 'cond'
        mom = Integrate(x * gf.components[0] * dx(mesh.Materials("cond")), mesh)
    m_z = -np.pi * s * SIGMA * mom
    alpha = complex(m_z) * mu0 / B0
    # The FEM solver uses s = +j*omega (e^{+j omega t} engineering convention);
    # the analytic G_exact (maglev_sphere_force) uses the e^{-j omega t}
    # physics convention, so its imaginary part has the opposite sign.  The
    # physical content (|alpha|, Re[alpha], loss magnitude) is identical; we
    # conjugate the FEM result to report alpha_c in the SAME convention as the
    # rest of this folder (so it equals 4 pi a^3 G(x) for a sphere, not its
    # conjugate).  Re[alpha] -- the only part the levitation force uses -- is
    # unaffected by the conjugation.
    return np.conj(alpha)


def sphere_analytic(omega, a=A_SPH):
    """Analytic axial polarizability of a conducting sphere: 4 pi a^3 G(x)."""
    return 4 * np.pi * a**3 * G_exact(omega)


def main():
    ngsglobals.msg_level = 0
    print(f"Axisymmetric FEM eddy-current alpha_c(omega); Cu sigma={SIGMA:.2e}, "
          f"a_eq={A_SPH*1e3:.1f} mm")

    # ---- check 1: VALIDATE the pipeline on a sphere vs analytic 4 pi a^3 G(x) ----
    print("\n[check 1] sphere: FEM alpha_c(omega) vs analytic 4 pi a^3 G(x)")
    mesh_s = build_mesh(A_SPH, A_SPH, box=0.10, maxh=0.012)
    print(f"   sphere mesh: {mesh_s.ne} elements, materials={mesh_s.GetMaterials()}")
    print("   f[Hz]   a/delta   FEM Re,Im [mm^3]        analytic Re,Im [mm^3]   relerr")
    worst = 0.0
    for f in (2e3, 1e4, 5e4, 2e5):
        w = 2 * np.pi * f
        afem = alpha_c(mesh_s, w)
        aan = sphere_analytic(w)
        rel = abs(afem - aan) / abs(aan)
        worst = max(worst, rel)
        print(f"  {f:7.0f}  {A_SPH/delta(w):6.2f}  "
              f"({afem.real*1e9:+8.2f},{afem.imag*1e9:+8.2f})  "
              f"({aan.real*1e9:+8.2f},{aan.imag*1e9:+8.2f})  {rel*100:5.2f}%")
    print(f"   worst relerr vs analytic sphere = {worst*100:.2f}%")
    assert worst < 0.05, f"FEM sphere alpha_c off by {worst*100:.1f}% from analytic G(x)"

    # ---- check 2: HF limit -> analytic demag-tensor value -V/(1-N_c) ----
    print("\n[check 2] HF limit alpha_c(omega->hi) -> -V/(1-N_c) (analytic demag)")
    f_hi = 5e5
    w_hi = 2 * np.pi * f_hi
    shapes = {
        "sphere   1:1": (A_SPH, A_SPH),
        "prolate  1:2": (A_SPH, 2 * A_SPH),
        "oblate   2:1": (A_SPH, 0.5 * A_SPH),
    }
    print("   shape         N_c     -V/(1-N_c)[mm^3]   FEM Re[alpha_c]   relerr  |a_c|/V")
    hf = {}
    for name, (ae, c) in shapes.items():
        N = ET.demag_tensor((ae, ae, c))
        V = 4 / 3 * np.pi * ae * ae * c
        kappa_hf = -V / (1.0 - N[2])
        m = build_mesh(ae, c, box=0.10, maxh=0.012)
        afem = alpha_c(m, w_hi)
        rel = abs(afem.real - kappa_hf) / abs(kappa_hf)
        hf[name] = (kappa_hf, afem.real, N[2], V)
        print(f"   {name}   {N[2]:.4f}   {kappa_hf*1e9:+9.2f}      "
              f"{afem.real*1e9:+9.2f}    {rel*100:5.2f}%   {abs(afem.real)/V:.3f}")
        assert rel < 0.08, f"{name}: HF alpha_c off {rel*100:.1f}% from -V/(1-N_c)"
    # the pure SHAPE signature is the per-volume response |alpha_c|/V = 1/(1-N_c):
    # flatter-along-the-field (oblate, large N_c) excludes more flux per volume
    # than elongated-along-the-field (prolate, small N_c).  (Absolute |alpha_c|
    # also tracks volume; per-volume isolates the shape.)
    pv = {k: abs(hf[k][1]) / hf[k][3] for k in shapes}
    assert pv["prolate  1:2"] < pv["sphere   1:1"] < pv["oblate   2:1"], \
        "per-volume axial response must increase prolate < sphere < oblate"

    # ---- check 3: full-frequency curve is shape-dependent ----
    print("\n[check 3] alpha_c(omega) full curve is shape-dependent "
          "(DC->0, HF->-V/(1-N_c))")
    freqs = np.array([5e2, 2e3, 1e4, 5e4, 2e5, 5e5])
    curves = {}
    for name, (ae, c) in shapes.items():
        m = build_mesh(ae, c, box=0.10, maxh=0.012)
        vals = [alpha_c(m, 2 * np.pi * f) for f in freqs]
        curves[name] = vals
        re = " ".join(f"{v.real*1e9:7.1f}" for v in vals)
        print(f"   {name}  Re[alpha_c] (mm^3): {re}")
    # DC end: |alpha_c| smallest at lowest freq for every shape (eddy onset)
    for name in shapes:
        amp = [abs(v) for v in curves[name]]
        assert amp[0] < amp[-1], f"{name}: |alpha_c| must grow from DC to HF"

    # ---- persist ----
    data = {
        "sigma": SIGMA, "a_eq_m": A_SPH, "B0": B0,
        "sphere_validation_worst_relerr": worst,
        "hf_anchor": {k: {"N_c": hf[k][2], "kappa_analytic_m3": hf[k][0],
                          "alpha_fem_re_m3": hf[k][1]} for k in shapes},
        "freqs_Hz": freqs.tolist(),
        "alpha_c_re_m3": {k: [v.real for v in curves[k]] for k in shapes},
        "alpha_c_im_m3": {k: [v.imag for v in curves[k]] for k in shapes},
    }
    output = validation_output("ellipsoid_alpha_omega_axisym_results.json", HERE)
    with open(output, "w") as fp:
        json.dump(data, fp, indent=2, allow_nan=False)
    print(f"\n wrote {output}")
    plot(freqs, curves, hf, shapes)


def plot(freqs, curves, hf, shapes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "font.family": "serif", "axes.linewidth": 0.8})
    fig, ax = plt.subplots(figsize=(8 / 2.54, 6 / 2.54))
    styles = {"sphere   1:1": ("C0", "-"), "prolate  1:2": ("C1", "--"),
              "oblate   2:1": ("C3", "-.")}
    for name in shapes:
        col, ls = styles[name]
        re = np.array([-v.real * 1e9 for v in curves[name]])   # plot -Re (lift)
        ax.semilogx(freqs, re, col + ls[0] if False else col, linestyle=ls,
                    lw=1.2, label=name.replace("  ", " "))
        ax.axhline(-hf[name][0] * 1e9, color=col, ls=":", lw=0.6)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel(r"$-\mathrm{Re}[\alpha_c]$ (mm$^3$)")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout(pad=0.3)
    out = os.path.join(HERE, "ellipsoid_alpha_omega_axisym.png")
    fig.savefig(out, dpi=300)
    print(f" wrote {out}")


if __name__ == "__main__":
    main()
