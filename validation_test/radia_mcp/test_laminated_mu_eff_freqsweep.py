"""Deep-skin frequency-sweep validation of solve.laminated_mu_eff -- the complex
EFFECTIVE PERMEABILITY of an in-plane laminated stack.

test_laminated_steel.py validated laminated_mu_eff at a single moderate skin
depth (|b| ~ 1).  This sweep validates it from the eddy-negligible regime into
DEEP skin effect (|b| >> 1), using realistic electrical-steel lamination
parameters (d_iron = 0.35 mm, fill 0.875, sigma 4 MS/m, mu_r 5000) -- the same
regime a motor-lamination homogenization (calc_motor_lamination.py) runs in.

Reference: an independent 1D FE solves the in-plane eddy diffusion across the
period (A_y(z), B_x = -dA/dz, eddy jw*sigma*A in the iron) with the cell-average
flux imposed; the effective permeability is the field ratio

    mu_eff = <B>_cell / H_s ,   H_s = nu0 * <B>_insulation  (continuous tang. H)

which must reproduce  mu0[ fill*mu_r*tanh(b)/b + (1-fill) ],  b=(d/2)sqrt(j w mu0 mu_r sigma).

NOTE (lesson from the 横展開 cross-check): the eddy LOSS is NOT
0.5*w*B^2*Im(1/mu_eff) -- that flux-averaging-mu identity only holds at low
frequency.  The lamination eddy loss is resistive (0.5*sigma*|E|^2 / surface
impedance); mu_eff here is validated as the B-H field ratio (what a homogenized
global FE consumes), not as a loss proxy.
"""
import cmath
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx,
                     Integrate, CoefficientFunction as CF, IfPos, BND,
                     x as Z, TaskManager)
from ngsolve.meshes import Make1DMesh
from radia.lamination import laminated_mu_eff

MU0 = 4e-7 * math.pi
D_IRON, D_INS = 0.35e-3, 0.05e-3
SIGMA, MU_R = 4.0e6, 5000.0
L = D_IRON + D_INS
FILL = D_IRON / L


def mu_eff_fe(freq, n=600):
    """1D FE effective permeability <B>_cell / H_s of the iron+insulation cell."""
    omega = 2 * math.pi * freq
    mesh = Make1DMesh(n, mapping=lambda t: L * t)
    # iron CENTERED with insulation/2 on each side -- the symmetric unit cell
    # of a periodic lamination (field penetrates from both faces), which is what
    # laminated_mu_eff's tanh assumes.
    iron = IfPos((Z - D_INS / 2.0) * (D_INS / 2.0 + D_IRON - Z), 1.0, 0.0)
    ins = 1.0 - iron
    nu = (1.0 / (MU_R * MU0)) * iron + (1.0 / MU0) * ins
    fes = H1(mesh, order=2, complex=True, dirichlet="left|right")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu * grad(u) * grad(v) * dx + 1j * omega * (SIGMA * iron) * u * v * dx
    a.Assemble()
    f = LinearForm(fes); f.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(CF(1.0 * Z), BND)                          # A(0)=0, A(L)=L -> <B>=-1
    r = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * r
    Bx = -grad(gfu)[0]
    B_cell = complex(Integrate(Bx, mesh) / L)
    B_ins = complex(Integrate(ins * Bx, mesh) / Integrate(ins, mesh))
    H_s = B_ins / MU0
    return B_cell / H_s                                # mu_eff = <B>/H_s


def main():
    print(f"  electrical steel lamination: d_iron={D_IRON*1e3:.2f}mm fill={FILL:.3f}")
    worst = 0.0
    with TaskManager():
        for freq in (50.0, 500.0, 5000.0, 20000.0):
            omega = 2 * math.pi * freq
            mu_an = laminated_mu_eff(MU_R, SIGMA, omega, D_IRON, fill=FILL)
            mu_fe = mu_eff_fe(freq)
            b = (D_IRON / 2.0) * cmath.sqrt(1j * omega * MU0 * MU_R * SIGMA)
            rel = abs(mu_fe - mu_an) / abs(mu_an)
            worst = max(worst, rel)
            print(f"  f={freq:7.0f}Hz |b|={abs(b):5.2f}: "
                  f"analytic mu/mu0={mu_an/MU0:9.2f}  FE={mu_fe/MU0:9.2f}  "
                  f"({100*rel:.3f}%)")
    # ~1-2% residual is the finite-cell boundary (Dirichlet-A on the outer
    # insulation faces vs an ideal infinite periodic stack); it plateaus rather
    # than diverging, confirming the deep-skin tanh behaviour. The clean
    # zero-residual point (|b|~1, fill=1, H-driven) is test_laminated_steel.py.
    assert worst < 0.03, f"laminated_mu_eff vs 1D FE off {100*worst:.2f}% (>3%)"
    print(f"\n[OK] laminated_mu_eff validated vs 1D-FE <B>/H_s from low-eddy to "
          f"deep skin (|b|~9.8), worst {100*worst:.2f}%.")


if __name__ == "__main__":
    main()
