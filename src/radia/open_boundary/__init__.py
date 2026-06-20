# -*- coding: utf-8 -*-
"""radia.open_boundary -- exact exterior DtN open boundary as a Cauer Ladder
Network (the "Zs-DtN-CLN" open boundary).

A radia CORE method (ships in the radia wheel; pure Python on numpy/scipy).  It
realises the EXACT exterior Dirichlet-to-Neumann symbol of a SEPARABLE truncation
(sphere / circle) as a finite Cauer ladder -- exact at n+1 stages, well-conditioned,
passive -- for magneto-quasistatic (eddy / diffusion) open boundaries where you want
an exact, DC-well-conditioned boundary AND a compact passive circuit / ROM.

Read `dtn_cln` (the module docstring) for the SCOPE (where this beats a CFS-PML and
where it does NOT) and the PRIOR ART (Grote-Keller / Hagstrom-Warburton / Warburg-Cauer
/ Kameari CLN -- this is a verified reusable operator, not a novelty claim).  See
docs/open_boundary/OPEN_BOUNDARY_MAP.md for the open-boundary selector.

Quick start::

    import radia.open_boundary as ob
    # exact eddy/diffusion DtN of multipole n=2 at a sphere R0, s = i*omega
    g = ob.eddy_dtn(2, 1j * 50.0, R0=0.1, mu_sigma=4e-7 * 3.14159 * 1e6)
    # the Cauer ladder is EXACT at n+1 stages and reproduces it:
    stages = ob.cauer_ladder(2)
    g_cln = ob.eval_ladder(stages, 1j * 50.0, R0=0.1, mu_sigma=...)
    # time-domain Robin realisation: one auxiliary ODE per companion pole (Re<0):
    poles = ob.companion_poles(2)
"""
from .dtn_cln import (  # noqa: F401
    reverse_bessel_theta,
    reverse_bessel_roots,
    eddy_dtn,
    eddy_dtn_rational_q,
    wave_dtn,
    cauer_ladder,
    eval_ladder,
    companion_poles,
    sqrt_s_passive_ladder,
    eval_sqrt_ladder,
)

__all__ = [
    "reverse_bessel_theta",
    "reverse_bessel_roots",
    "eddy_dtn",
    "eddy_dtn_rational_q",
    "wave_dtn",
    "cauer_ladder",
    "eval_ladder",
    "companion_poles",
    "sqrt_s_passive_ladder",
    "eval_sqrt_ladder",
]
