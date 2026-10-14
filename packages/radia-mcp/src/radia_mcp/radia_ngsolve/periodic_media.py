"""1-D periodic (Bragg / multilayer) media at normal incidence -- the photonic-bandgap /
distributed-Bragg-reflector member of the high-frequency family.

A bi-layer unit cell (n1 over length L1, n2 over L2) has the Abeles characteristic matrix per
layer  M_i = [[cos d_i, (i/eta_i) sin d_i], [i eta_i sin d_i, cos d_i]],  d_i = n_i omega L_i / c,
eta_i ~ 1/n_i (non-magnetic). For the infinite periodic stack the Bloch wavenumber K obeys

    cos(K Lambda) = 1/2 tr(M1 M2)
                  = cos d1 cos d2 - 1/2 (z + 1/z) sin d1 sin d2 ,   z = n1/n2 ,

(Lambda = L1+L2). Where |1/2 tr| <= 1 the cell passes a propagating Bloch mode (allowed band);
where |1/2 tr| > 1, K is complex -> the wave is evanescent (a STOP BAND / photonic bandgap).

A QUARTER-WAVE stack (n1 L1 = n2 L2 = lambda0/4) is the canonical distributed Bragg reflector:
at the design frequency f0 both layers are a quarter wave (d1 = d2 = pi/2), so 1/2 tr(f0) =
-1/2 (z + 1/z), which is < -1 for any index contrast -> f0 sits at the CENTRE of the widest
stop band. No contrast (n1 = n2) gives |1/2 tr| <= 1 everywhere (no gap, as it must).
"""
import numpy as np

C_LIGHT = 299792458.0


def bragg_half_trace(n1, L1, n2, L2, freq):
    """Half-trace 1/2 tr(M1 M2) = cos(K Lambda) of the bi-layer unit cell at frequency ``freq``
    [Hz] (scalar or array). |value| <= 1 -> allowed (propagating) band; > 1 -> stop band."""
    w = 2.0 * np.pi * np.asarray(freq, dtype=float)
    d1 = n1 * w * L1 / C_LIGHT
    d2 = n2 * w * L2 / C_LIGHT
    z = n1 / n2
    return np.cos(d1) * np.cos(d2) - 0.5 * (z + 1.0 / z) * np.sin(d1) * np.sin(d2)


def bragg_bandgaps(n1, L1, n2, L2, fmin, fmax, npts=200000):
    """Stop bands (photonic bandgaps) of the infinite bi-layer stack in [fmin, fmax], as a list of
    (f_lo, f_hi) intervals where |1/2 tr| > 1 (evanescent Bloch mode). Sampled on ``npts`` points."""
    f = np.linspace(fmin, fmax, npts)
    forbidden = np.abs(bragg_half_trace(n1, L1, n2, L2, f)) > 1.0
    gaps = []
    in_gap = False
    start = 0.0
    for i in range(len(f)):
        if forbidden[i] and not in_gap:
            in_gap, start = True, f[i]
        elif not forbidden[i] and in_gap:
            in_gap = False
            gaps.append((float(start), float(f[i - 1])))
    if in_gap:
        gaps.append((float(start), float(f[-1])))
    return gaps


def quarter_wave_stack_layers(f0, n1, n2):
    """Layer thicknesses (L1, L2) of the quarter-wave (Bragg) stack designed at ``f0`` [Hz]:
    L_i = lambda0 / (4 n_i) with lambda0 = c / f0. Returns ``(L1, L2)`` [m]."""
    lam0 = C_LIGHT / f0
    return lam0 / (4.0 * n1), lam0 / (4.0 * n2)
