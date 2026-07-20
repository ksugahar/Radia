"""Validates the FEMM laminated-steel material model (radia_ngsolve.solve):

  AC  : laminated_mu_eff(mu_r, sigma, omega, d, fill) -- the complex effective
        permeability of in-plane laminations.  Checked against a direct 1D
        finite-element solve of the eddy-current diffusion ACROSS one lamination
        (H'' = j w mu sigma H, H = Hs at both faces): the FE volume-average
        <B>/(mu0 H) must reproduce the closed form  mu_r * tanh(b)/b.

  DC  : laminated_reluctivity_tensor(...) -- the anisotropic (fill-factor)
        reluctivity.  Checked by a uniform field NORMAL to the stack: a
        geometry-RESOLVED steel/air layer stack and the single-block HOMOGENIZED
        tensor must give the same magnetic energy (= 1/2 B^2 V / mu_perp).
        Done for both stacking directions to pin both tensor entries (and catch
        an x/y orientation swap).
"""
import cmath
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

MU0 = 4e-7 * math.pi


def check_ac_complex_mu():
    from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx,
                         Integrate)
    from ngsolve.meshes import Make1DMesh
    from radia_mcp.radia_ngsolve.solve import laminated_mu_eff

    mu_r, sigma, freq, d = 1000.0, 2.0e6, 1000.0, 0.5e-3
    omega = 2 * math.pi * freq
    mu = MU0 * mu_r
    Hs = 1.0

    # Reference element [0,1] mapped to the lamination thickness [-d/2, d/2].
    # Weak form for the deviation h = H - Hs (Dirichlet h=0 at both faces):
    #   int h' v'  +  B2 int h v  =  -B2 Hs int v ,    B2 = j w mu sigma d^2
    B2 = 1j * omega * mu * sigma * d * d
    mesh = Make1DMesh(200)
    fes = H1(mesh, order=3, complex=True, dirichlet="left|right")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx + B2 * u * v * dx
    f = LinearForm(fes)
    f += -B2 * Hs * v * dx
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec

    mean_h = Integrate(gfu, mesh)             # int over [0,1] = average
    mu_eff_fe = mu * (Hs + mean_h) / Hs       # mu * <H>/Hs
    mu_eff_formula = laminated_mu_eff(mu_r, sigma, omega, d, fill=1.0)

    rel = abs(mu_eff_fe - mu_eff_formula) / abs(mu_eff_formula)
    fe, fo = mu_eff_fe / MU0, mu_eff_formula / MU0
    print(f"  [AC]  mu_eff/mu0  FE      = {fe.real:.4f} {fe.imag:+.4f}j")
    print(f"        mu_eff/mu0  formula = {fo.real:.4f} {fo.imag:+.4f}j    "
          f"rel diff = {100*rel:.3f}%")
    assert rel < 1e-2, f"laminated_mu_eff off by {100*rel:.2f}% vs 1D FE (>1%)"


def _energy_uniform_field(nu, lam_resolved, lam_normal, mu_r, fill,
                          N=10, W=1.0, order=2):
    """Magnetic energy of a uniform field NORMAL to the stack, computed on a
    geometry-resolved steel/air stack (lam_resolved=True) or on the single-block
    homogenized anisotropic tensor (lam_resolved=False).  Returns (energy, mesh).
    """
    from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction, grad,
                         dx, InnerProduct, BND, x as xc, y as yc)
    from netgen.occ import OCCGeometry, MoveTo, Glue
    from radia_mcp.radia_ngsolve.solve import (reluctivity,
                                               laminated_reluctivity_tensor)

    # Uniform B normal to the sheets: stack along x -> B=(B0,0) (A = B0*y);
    # stack along y -> B=(0,B0) (A = -B0*x).
    B0 = 1.0
    A_bnd = B0 * yc if lam_normal == "x" else -B0 * xc

    if lam_resolved:
        period = W / N
        wsteel = fill * period
        faces = []
        for i in range(N):
            lo = i * period
            if lam_normal == "x":      # strips stacked along x (vertical strips)
                s = MoveTo(lo, 0).Rectangle(wsteel, W).Face()
                g = MoveTo(lo + wsteel, 0).Rectangle(period - wsteel, W).Face()
            else:                       # strips stacked along y (horizontal)
                s = MoveTo(0, lo).Rectangle(W, wsteel).Face()
                g = MoveTo(0, lo + wsteel).Rectangle(W, period - wsteel).Face()
            s.faces.name = "steel"
            g.faces.name = "air"
            faces += [s, g]
        shape = Glue(faces)
        shape.edges.name = "outer"
        mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=period / 3))
        nu_use = reluctivity(mesh, {"steel": mu_r})
    else:
        shape = MoveTo(0, 0).Rectangle(W, W).Face()
        shape.faces.name = "lam"
        shape.edges.name = "outer"
        mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=W / 12))
        nu_use = laminated_reluctivity_tensor(
            mesh, {"lam": (mu_r, fill)}, lam_normal=lam_normal)

    fes = H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += nu_use * grad(u) * grad(v) * dx
    a.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(A_bnd, BND)
    res = gfu.vec.CreateVector()
    res.data = -a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res
    tmp = gfu.vec.CreateVector()
    tmp.data = a.mat * gfu.vec
    energy = 0.5 * InnerProduct(gfu.vec, tmp)
    return energy, W * W


def check_static_anisotropy():
    mu_r, fill = 1000.0, 0.5
    mu_perp = MU0 / (fill / mu_r + (1.0 - fill))
    e_exact = 0.5 * (1.0 / mu_perp) * 1.0   # 1/2 * B0^2 * V / mu_perp (B0=V=1)
    for lam_normal in ("x", "y"):
        e_res, _ = _energy_uniform_field(None, True, lam_normal, mu_r, fill)
        e_hom, _ = _energy_uniform_field(None, False, lam_normal, mu_r, fill)
        rel = abs(e_hom - e_res) / abs(e_res)
        rel_exact = abs(e_hom - e_exact) / e_exact
        print(f"  [DC {lam_normal}] energy resolved={e_res:.6e} "
              f"homogenized={e_hom:.6e} exact={e_exact:.6e}  "
              f"hom-vs-res={100*rel:.3f}%  hom-vs-exact={100*rel_exact:.3f}%")
        assert rel < 1e-2, f"homogenized vs resolved off {100*rel:.2f}% (>1%)"
        assert rel_exact < 1e-2, f"homogenized vs exact off {100*rel_exact:.2f}%"


def main():
    from ngsolve import TaskManager
    with TaskManager():
        check_ac_complex_mu()
        check_static_anisotropy()
    print("\n[OK] laminated-steel model validated: AC complex mu_eff vs 1D eddy "
          "FE (<1%), DC anisotropic tensor vs resolved stack (<1%).")


if __name__ == "__main__":
    main()
