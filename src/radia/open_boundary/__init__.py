# -*- coding: utf-8 -*-
"""radia.open_boundary -- exact exterior DtN open boundaries.

A radia CORE method (ships in the radia wheel; pure Python on numpy/scipy).  For a
SEPARABLE truncation (sphere / circle) `dtn_exact` evaluates the EXACT exterior
Dirichlet-to-Neumann symbol per multipole -- a reverse-Bessel rational function in
s (wave) or q = sqrt(s) (magneto-quasistatic eddy / diffusion) -- together with its
passive, stable companion poles for a transient Robin boundary and a passive
Foster-form realisation of the sqrt(s) diffusion memory.  `kelvin_dtn` BUILDS the
DtN of a non-separable or material exterior by Kelvin-transformed FEM and reduces
it over a band with a low-degree rational fit.

Read the `dtn_exact` module docstring for the SCOPE (where this beats a CFS-PML and
where it does NOT) and the prior art (Grote-Keller / Hagstrom-Warburton).  See
docs/open_boundary/OPEN_BOUNDARY_MAP.md for the open-boundary selector.

Quick start::

    import radia.open_boundary as ob
    # exact eddy/diffusion DtN of multipole n=2 at a sphere R0, s = i*omega
    g = ob.eddy_dtn(2, 1j * 50.0, R0=0.1, mu_sigma=4e-7 * 3.14159 * 1e6)
    # time-domain Robin realisation: one auxiliary ODE per companion pole (Re<0):
    poles = ob.companion_poles(2)

`cauer_ladder` / `eval_ladder` (module `dtn_cln`) are still exported and are
scheduled for removal with the repository's Cauer-ladder code.
"""
from .dtn_exact import (  # noqa: F401
    reverse_bessel_theta,
    reverse_bessel_roots,
    eddy_dtn,
    eddy_dtn_rational_q,
    wave_dtn,
    companion_poles,
    sqrt_s_passive_ladder,
    eval_sqrt_ladder,
)
from .dtn_cln import (  # noqa: F401
    cauer_ladder,
    eval_ladder,
)
from .kelvin_dtn import (  # noqa: F401
    kelvin_fem_radial_dtn,
    kelvin_dtn_matrix,
    steklov_spectrum,
    band_rational_fit,
)

__all__ = [
    # dtn_exact -- exact closed-form separable DtN
    "reverse_bessel_theta",
    "reverse_bessel_roots",
    "eddy_dtn",
    "eddy_dtn_rational_q",
    "wave_dtn",
    "companion_poles",
    "sqrt_s_passive_ladder",
    "eval_sqrt_ladder",
    # dtn_cln -- Cauer realisation, scheduled for removal
    "cauer_ladder",
    "eval_ladder",
    # kelvin_dtn -- Kelvin-built material-aware / non-separable DtN
    "kelvin_fem_radial_dtn",
    "kelvin_dtn_matrix",
    "steklov_spectrum",
    "band_rational_fit",
]
