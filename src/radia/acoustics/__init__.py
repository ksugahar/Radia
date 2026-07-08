"""radia.acoustics -- analytic acoustic scattering references (application lane).

Acoustic Helmholtz scattering is an APPLICATION lane built on ngsolve.bem (the
Helmholtz kernel), NOT part of the Laplace-only radia core (see the "Green's
Function: Laplace Kernel Only" policy).  This subpackage holds the analytic
partial-wave sphere-scattering references ported from the readable MATLAB
acoustic FEM/BEM teaching lane (matlab-acoustic-fembem), used as the gold
standard for validating numerical acoustic BEM (ngsolve.bem) solves -- the
"Complement NGSolve" strategy: NGSolve/ngsolve.bem provide the numerical
Helmholtz FEM/BEM, radia.acoustics provides the analytic truth to check them.

Validated 3-way (validation_test/acoustics/): Python analytic == MATLAB analytic
to ~1e-14, and analytic == ngsolve.bem numerical (soft sphere) to ~2e-5.

Public API (analytic scattering references, pure numpy/scipy):
  soft_sphere_scattering     -- sound-soft (p = 0)
  rigid_sphere_scattering    -- sound-hard (dp/dn = 0)
  fluid_sphere_scattering    -- penetrable fluid sphere (Anderson 1950)
  elastic_sphere_scattering  -- solid elastic sphere (Faran 1951)

Submodule (needs NGSolve, import explicitly):
  radia.acoustics.fsi        -- elastic fluid-structure interaction coupled solve
                                (NGSolve VectorH1(order=p) interior + spherical DtN
                                exterior); fsi.fsi_dtn_solve, fsi.sphere_mesh.
"""

from .scattering import (
    soft_sphere_scattering,
    rigid_sphere_scattering,
    fluid_sphere_scattering,
    elastic_sphere_scattering,
)

__all__ = [
    "soft_sphere_scattering",
    "rigid_sphere_scattering",
    "fluid_sphere_scattering",
    "elastic_sphere_scattering",
]
