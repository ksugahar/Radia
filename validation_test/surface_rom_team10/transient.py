"""Nonlinear transient eddy-current FOM on the TEAM-13 geometry (Kelvin open boundary).

Backward Euler, A-formulation (no scalar potential; the steel plates are the only
conductors), energy-minimising Newton per step:
  E(A) = int_steel sigma/(2dt)|A-A_n|^2 + w(|curl A|)
       + int_rest nu_K/2 |curl A|^2 + gauge nu0/2 |A|^2  -  i(t) <f_unit, A>
Snapshots of the full DOF vector are saved for the ROM.

Placeholders (TEAM 10 values pending): SIGMA_PLACEHOLDER, waveform parameters.
Usage: python transient.py <maxh_steel> <order> <wave> [steps]
waves: rise (i=NI(1-exp(-t/tau))), rise_fast, pulse (rise then decay), reverse
"""
import sys, json, time, platform, hashlib
from pathlib import Path
import numpy as np
import ngsolve
from ngsolve import *
from team13_model import *

ROOT = Path(__file__).resolve().parent
import os as _os
T_END = float(_os.environ.get('TEAM_T_END', '0.12'))   # override only for truncated FOM studies
NI0 = 1000.0


def waveform(name):
    if name == 'rise': return lambda t: NI0 * (1 - np.exp(-t / .05))
    if name == 'rise_fast': return lambda t: NI0 * (1 - np.exp(-t / .02))
    if name == 'pulse': return lambda t: 1.2 * NI0 * (1 - np.exp(-t / .015)) * np.exp(-t / .08)
    if name == 'reverse': return lambda t: NI0 * (1 - np.exp(-t / .03)) * np.cos(2 * np.pi * t / .12)
    if name == 'half': return lambda t: .5 * NI0 * (1 - np.exp(-t / .05))
    if name == 'team10': return lambda t: TEAM10_TURNS * TEAM10_IM * (1 - np.exp(-t / TEAM10_TAU))
    # saturation-state tests (the permeability distribution differs from the training states)
    if name == 'team10_x0.3': return lambda t: 0.3 * TEAM10_TURNS * TEAM10_IM * (1 - np.exp(-t / TEAM10_TAU))
    if name == 'team10_x2': return lambda t: 2.0 * TEAM10_TURNS * TEAM10_IM * (1 - np.exp(-t / TEAM10_TAU))
    if name == 'team10_off': return lambda t: TEAM10_TURNS * TEAM10_IM * (1 - np.exp(-t / TEAM10_TAU)) * (t <= 0.06 + 1e-12)
    # two independently driven coils (TEAM_COIL2=1): currents (main coil, web coil) in AT
    NI1, NI2 = TEAM10_TURNS * TEAM10_IM, 400.0
    rise = lambda t, a: a * (1 - np.exp(-t / TEAM10_TAU))
    pulse = lambda t, a: 1.2 * a * (1 - np.exp(-t / .015)) * np.exp(-t / .08)
    two = {'c1_rise': lambda t: np.array([rise(t, NI1), 0.0]),
           'c1_pulse': lambda t: np.array([pulse(t, NI1), 0.0]),
           'c2_rise': lambda t: np.array([0.0, rise(t, NI2)]),
           'c2_pulse': lambda t: np.array([0.0, pulse(t, NI2)]),
           'mix_add': lambda t: np.array([rise(t, NI1), rise(t, NI2)]),
           'mix_oppose': lambda t: np.array([rise(t, NI1), -rise(t, NI2)]),
           # candidate pool for the automatic pattern selection (run_pattern_select.py); none equals a test
           'mixp_add': lambda t: np.array([pulse(t, NI1), pulse(t, NI2)]),
           'mixp_oppose': lambda t: np.array([pulse(t, NI1), -pulse(t, NI2)]),
           'c2_pulse_x2': lambda t: np.array([0.0, pulse(t, 2 * NI2)]),
           'c1_rise_half': lambda t: np.array([rise(t, .5 * NI1), 0.0]),
           'c1_rise_neg': lambda t: np.array([-rise(t, NI1), 0.0])}
    if name in two: return two[name]
    raise ValueError(name)


class Model:
    """Shared FOM operators (also used by the ROM)."""
    def __init__(self, maxh=.006, order=1, dt=T_END / 60, gauge=None, half=False):
        if gauge is None: gauge = float(_os.environ.get('TEAM_GAUGE', '1e-6'))   # air mass regularisation / nu0
        self.gauge = gauge
        self.half = half
        self.mesh = (build_mesh_half if half else build_mesh)(maxh_steel=maxh); m = self.mesh
        self.fes = spaces(m, order); self.dt = dt; self.order = order; self.maxh = maxh
        u, v = self.fes.TnT(); w, _ = energy_density()
        nu = make_kelvin_nu_cf(m, R_K, OFFSET)
        self.An = GridFunction(self.fes)
        steel = m.Materials('steel'); rest = ~steel
        self.a = BilinearForm(self.fes, symmetric=True)
        self.a += Variation(0.5 * SIGMA_PLACEHOLDER / dt * (u - self.An) * (u - self.An) * dx(steel))
        self.a += Variation(w(sqrt(1e-12 + curl(u) * curl(u))) * dx(steel, bonus_intorder=2))
        self.a += Variation(0.5 * nu * curl(u) * curl(u) * dx(rest, bonus_intorder=2))
        self.a += Variation(0.5 * gauge * NU_0 * u * u * dx(rest))
        self.f, _ = source_linear_form(self.fes, m, 1.0, order)
        # linear "rest" operator (air extension) and steel curl mass (error metric)
        self.k_rest = BilinearForm(self.fes, symmetric=True)
        self.k_rest += nu * curl(u) * curl(v) * dx(rest, bonus_intorder=2) + gauge * NU_0 * u * v * dx(rest)
        self.k_rest.Assemble()
        self.c_steel = BilinearForm(curl(u) * curl(v) * dx(steel), symmetric=True).Assemble()
        self.gfu = GridFunction(self.fes)
        free = self.fes.FreeDofs()
        st = self.fes.GetDofs(steel); gam = self.fes.GetDofs(m.Boundaries('steel_surf'))
        self.idx_gamma = np.array([i for i in range(self.fes.ndof) if free[i] and st[i] and gam[i]])
        self.idx_int = np.array([i for i in range(self.fes.ndof) if free[i] and st[i] and not gam[i]])
        self.idx_air = np.array([i for i in range(self.fes.ndof) if free[i] and not st[i]])

    def energy(self, vec, i):
        return self.a.Energy(vec) - i * InnerProduct(self.f.vec, vec)

    def residual(self, vec, i, out):
        self.a.Apply(vec, out); out.data -= i * self.f.vec

    def steel_B_norm2(self, d):
        """d^T C_steel d for a numpy DOF vector (|curl|^2 over steel)."""
        g = self.gfu.vec.CreateVector(); g.FV().NumPy()[:] = d; h = g.CreateVector()
        self.c_steel.mat.Mult(g, h); return float(np.dot(d, h.FV().NumPy()))


def run_fom(model, wave, steps=60, tol=1e-10, max_newton=40):
    """max_newton < 40 gives cheap, inexact training snapshots: the step is accepted
    after max_newton Newton updates (max_newton=1: linearly implicit Euler)."""
    fn = waveform(wave); m = model; dt = T_END / steps; assert abs(dt - m.dt) < 1e-15
    gfu = m.gfu; gfu.vec[:] = 0; m.An.vec[:] = 0
    r = gfu.vec.CreateVector(); du = gfu.vec.CreateVector(); trial = gfu.vec.CreateVector()
    snaps = [gfu.vec.FV().NumPy().copy()]; hist = []; t0 = time.time()
    for k in range(1, steps + 1):
        i = fn(k * dt); m.An.vec.data = gfu.vec
        for it in range(40):
            if it == max_newton: break                      # inexact training snapshot
            m.a.AssembleLinearization(gfu.vec); m.residual(gfu.vec, i, r)
            du.data = m.a.mat.Inverse(m.fes.FreeDofs(), inverse='pardiso') * r
            dec = abs(InnerProduct(du, r))
            if it == 0: dec0 = max(dec, 1e-300)
            if dec <= tol * dec0 or dec < 1e-24: break
            E0 = m.energy(gfu.vec, i); tau = 1.
            while True:
                trial.data = gfu.vec - tau * du
                if m.energy(trial, i) <= E0 or tau < 1e-4: break
                tau *= .5
            gfu.vec.data = trial
        else: raise RuntimeError(('Newton', k, dec / dec0))
        snaps.append(gfu.vec.FV().NumPy().copy())
        B = curl(gfu)
        pos = team10_sections(B, m.mesh) if CASE == 'team10' else average_B_positions(B, m.mesh, n=(3, 6))
        Jp = -SIGMA_PLACEHOLDER * (gfu - m.An) / dt
        jc = Jp(m.mesh(0.0, 0.0, .03))   # centre plate probe
        hist.append(dict(t=k * dt, current_AT=float(i), newton=it, B_pos=pos, J_probe=[float(c) for c in jc]))
        if k % 10 == 0:
            b1 = pos['S1'] if isinstance(pos, dict) else pos[0]
            print(wave, k, steps, 'newton', it, 'B1', round(b1, 3), round(time.time() - t0, 1), flush=True)
    return np.array(snaps), hist, time.time() - t0


if __name__ == '__main__':
    maxh, order, wave = float(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    steps = int(sys.argv[4]) if len(sys.argv) > 4 else 60
    half = len(sys.argv) > 5 and sys.argv[5] == 'half'
    max_newton = int(sys.argv[6]) if len(sys.argv) > 6 else 40
    import os
    if os.environ.get('NGS_THREADS'): SetNumThreads(int(os.environ['NGS_THREADS']))
    with TaskManager():
        model = Model(maxh, order, T_END / steps, half=half)
        snaps, hist, sec = run_fom(model, wave, steps, max_newton=max_newton)
    out = ROOT / 'results' / (('team10_' if CASE == 'team10' else '') + f'fom_h{maxh*1000:g}_p{order}_n{steps}' + ('_half' if half else '') + ('_quarter' if QUARTER else '')
                              + (f'_nw{max_newton}' if max_newton < 40 else '') + (f'_T{T_END*1000:g}ms' if abs(T_END - 0.12) > 1e-12 else '') + (f'_g{model.gauge:g}' if model.gauge != 1e-6 else '')); out.mkdir(parents=True, exist_ok=True)
    np.save(out / f'{wave}.npy', snaps)
    res = dict(case=f'{CASE}_transient_{wave}', team_case=CASE, sigma=SIGMA_PLACEHOLDER, half=half, quarter=QUARTER, gauge=model.gauge, max_newton=max_newton,
               linear_solves=int(sum(min(h['newton'] + 1, max_newton) for h in hist)), maxh_steel=maxh, order=order, steps=steps, dt=T_END / steps,
               sigma_placeholder=SIGMA_PLACEHOLDER, ndof=model.fes.ndof, n_gamma=len(model.idx_gamma),
               n_interior=len(model.idx_int), n_air=len(model.idx_air), solve_s=sec, history=hist,
               host=platform.node(), ngsolve=ngsolve.__version__,
               sources={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in ['team13_model.py', 'transient.py']})
    (out / f'{wave}.json').write_text(json.dumps(res, indent=2))
    print('done', wave, 'ndof', model.fes.ndof, 'gamma', len(model.idx_gamma), 'int', len(model.idx_int), 'sec', round(sec, 1))
