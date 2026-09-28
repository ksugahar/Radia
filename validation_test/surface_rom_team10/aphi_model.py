"""A-phi model (HCurl nograds x H1(steel)) shared by the FOM and the ROM.

Energy per backward-Euler step (phi is scaled by dt):
  E(A,phi) = int_steel sigma/(2dt) |A - A_n + grad phi|^2 + w(|curl A|)
           + int_rest nu_K/2 |curl A|^2 + int eps nu0/2 |A|^2 + int_steel reg/2 phi^2 - i(t) <f, A>
DOF classes of the compound vector [A | phi]:
  gamma : A DOFs on the steel surface (the only steel DOFs coupled to the air)
  int   : all other steel DOFs, i.e. interior A DOFs and every phi DOF
  air   : A DOFs outside the steel
"""
import os, time
import numpy as np
from ngsolve import *
from team13_model import *


class ModelAphi:
    def __init__(self, maxh=.006, order=3, dt=0.002, gauge=None):
        gauge = float(os.environ.get('TEAM_GAUGE', '1e-6')) if gauge is None else gauge
        self.gauge = gauge; self.dt = dt; self.order = order; self.maxh = maxh
        self.mesh = m = build_mesh_half(maxh_steel=maxh)
        sym_y = 'sym_y' in m.GetBoundaries()
        self.Vh = Periodic(HCurl(m, order=order, nograds=True, dirichlet='sym_y' if sym_y else '', dirichlet_bbnd='GND'))
        self.Qh = H1(m, order=order + 1, definedon='steel', dirichlet='sym_y' if sym_y else '')
        self.fes = fes = self.Vh * self.Qh
        (u, p), (v, q) = fes.TnT(); w, _ = energy_density()
        nu = make_kelvin_nu_cf(m, R_K, OFFSET); steel = m.Materials('steel'); rest = ~steel
        s = SIGMA_PLACEHOLDER; reg = 1e-12 * s / dt
        self.gn = GridFunction(fes); An = self.gn.components[0]
        self.a = BilinearForm(fes, symmetric=True)
        self.a += Variation(0.5 * s / dt * (u - An + grad(p)) * (u - An + grad(p)) * dx(steel))
        self.a += Variation(w(sqrt(1e-12 + curl(u) * curl(u))) * dx(steel, bonus_intorder=2))
        self.a += Variation(0.5 * nu * curl(u) * curl(u) * dx(rest, bonus_intorder=2))
        self.a += Variation(0.5 * gauge * NU_0 * u * u * dx)
        self.a += Variation(0.5 * reg * p * p * dx(steel))
        _, Jp = source_linear_form(self.Vh, m, 1.0, order)
        self.f = LinearForm(Jp * v * dx('coil', bonus_intorder=2)).Assemble()
        self.fs = [self.f]                      # one unit-current linear form per independently driven coil
        if 'coil2' in m.GetMaterials():
            _, Jp2 = coil2_source(self.Vh, m, 1.0, order)
            self.fs.append(LinearForm(Jp2 * v * dx('coil2', bonus_intorder=2)).Assemble())
        # exact linear pieces used by the ROM
        self.k_rest = BilinearForm(nu * curl(u) * curl(v) * dx(rest, bonus_intorder=2) + gauge * NU_0 * u * v * dx(rest), symmetric=True).Assemble()
        self.k_steel_lin = BilinearForm(gauge * NU_0 * u * v * dx(steel) + reg * p * q * dx(steel), symmetric=True).Assemble()
        self.m_sig = BilinearForm(s * (u + grad(p)) * (v + grad(q)) * dx(steel), symmetric=True).Assemble()
        self.c_steel = BilinearForm(curl(u) * curl(v) * dx(steel), symmetric=True).Assemble()
        self.gfu = GridFunction(fes)
        n = fes.ndof; nA = self.Vh.ndof; free = fes.FreeDofs(); st = fes.GetDofs(steel); gam = fes.GetDofs(m.Boundaries('steel_surf'))
        isA = np.arange(n) < nA
        fr = np.array([free[i] for i in range(n)]); stv = np.array([st[i] for i in range(n)]); gv = np.array([gam[i] for i in range(n)])
        self.idx_gamma = np.where(fr & stv & gv & isA)[0]
        self.idx_int = np.where(fr & stv & ~(gv & isA))[0]
        self.idx_air = np.where(fr & ~stv)[0]
        self.idx_steelA = np.where(stv & isA)[0]
        self.isA = isA

    def bfield(self, g): return curl(g.components[0])

    def currents(self, i):
        c = np.atleast_1d(np.asarray(i, dtype=float))
        if len(c) != len(self.fs): raise ValueError(('currents', len(c), 'sources', len(self.fs)))
        return c

    def energy(self, vec, i):
        return self.a.Energy(vec) - sum(ci * InnerProduct(fs.vec, vec) for ci, fs in zip(self.currents(i), self.fs))

    def residual(self, vec, i, out):
        self.a.Apply(vec, out)
        for ci, fs in zip(self.currents(i), self.fs): out.data -= ci * fs.vec

    def steel_B_norm2(self, d):
        g = self.gfu.vec.CreateVector(); g.FV().NumPy()[:] = d; h = g.CreateVector()
        self.c_steel.mat.Mult(g, h); return float(np.dot(d, h.FV().NumPy()))


def run_fom(m, fn, steps, max_newton=40, tol=1e-10, log=print):
    g = m.gfu; g.vec[:] = 0; m.gn.vec[:] = 0; dt = m.dt
    r = g.vec.CreateVector(); du = g.vec.CreateVector(); tr = g.vec.CreateVector()
    snaps = [g.vec.FV().NumPy().copy()]; hist = []; solves = 0; t0 = time.time()
    for k in range(1, steps + 1):
        i = fn(k * dt); m.gn.vec.data = g.vec
        for it in range(40):
            if it == max_newton: break
            m.a.AssembleLinearization(g.vec); m.residual(g.vec, i, r)
            du.data = m.a.mat.Inverse(m.fes.FreeDofs(), inverse='sparsecholesky') * r; solves += 1
            dec = abs(InnerProduct(du, r))
            if it == 0: d0 = max(dec, 1e-300)
            if dec <= tol * d0 or dec < 1e-24: break
            E0 = m.energy(g.vec, i); tau = 1.
            while True:
                tr.data = g.vec - tau * du
                if m.energy(tr, i) <= E0 or tau < 1e-4: break
                tau *= .5
            g.vec.data = tr
        else: raise RuntimeError(('Newton', k))
        snaps.append(g.vec.FV().NumPy().copy())
        pos = team10_sections(m.bfield(g), m.mesh)
        hist.append(dict(t=k * dt, current_AT=[float(c) for c in m.currents(i)], newton=it, B_pos=pos))
        if k % 10 == 0: log(f'aphi {k}/{steps} newton {it} S1 {pos["S1"]:.3f} {time.time()-t0:.0f}s')
    return np.array(snaps), hist, dict(seconds=time.time() - t0, linear_solves=solves)
