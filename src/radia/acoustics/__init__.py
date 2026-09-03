"""Radia's NGSolve-based acoustic application and validation lane.

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

Public API (readable SciPy analytical references, intentionally not native):
  soft_sphere_scattering     -- sound-soft (p = 0)
  rigid_sphere_scattering    -- sound-hard (dp/dn = 0)
  fluid_sphere_scattering    -- penetrable fluid sphere (Anderson 1950)
  elastic_sphere_scattering  -- solid elastic sphere (Faran 1951)

Submodules (need NGSolve / ngsolve.bem, import explicitly):
  radia.acoustics.fsi        -- elastic fluid-structure interaction coupled solve
                                (NGSolve VectorH1(order=p) interior + spherical DtN
                                exterior); fsi.fsi_dtn_solve, fsi.sphere_mesh.
  radia.acoustics.cq         -- Lubich convolution-quadrature time-domain sound-soft
                                BEM (ngsolve.bem Helmholtz single layer at complex
                                kappa); cq.cq_soft_sphere_scattering,
                                cq.soft_sphere_scattering_complex_k.

The analytical references remain outside Radia's C++/pybind11/MEX kernels so
they provide an implementation-independent check of the numerical solvers.
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
