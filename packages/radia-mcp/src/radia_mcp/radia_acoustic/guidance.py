"""NGSolve/ngsolve.bem-first guidance for ``radia.acoustics``."""

from __future__ import annotations


_TOPICS = {
    "overview": """# Radia acoustics

`radia.acoustics` is the production acoustic application lane. NGSolve owns
finite-element spaces, assembly, geometry, orientation, quadrature, and mapped
evaluation; `ngsolve.bem` owns Helmholtz boundary operators. Radia supplies
readable analytical scattering references, FEM-BEM/FSI coupling conventions,
spherical DtN closure, and Lubich convolution quadrature. The references remain
outside C++/pybind11/MEX so they can judge the numerical path independently.
This server does not execute or
manage MATLAB and does not duplicate the educational MATLAB solver.
""",
    "ngsolve_bem": """# Helmholtz BEM

Use `ngsolve.bem.HelmholtzSL` with an NGSolve surface space and caller-owned
`TaskManager`. Validate signs, normals, boundary labels, and frequency-domain
fields against `radia.acoustics.soft_sphere_scattering` or another analytic
reference before interpreting a plot. Complex wavenumbers are supported for
CQ Laplace nodes. Human visualization belongs in `netgen.webgui`; durable
headless evidence should be JSON plus mesh/result artifacts.
""",
    "fsi": """# Acoustic FSI

Use `radia.acoustics.fsi.fsi_dtn_solve`. The solid is an NGSolve `VectorH1`
elasticity FEM; normal displacement couples to exterior pressure. The current
exact exterior closure is a spherical Helmholtz DtN operator, not PML and not a
Kelvin transform. Before solving, verify the boundary name, spherical radius
deviation, material wave speeds (`cL > sqrt(2)*cT` for positive lambda), density,
frequency, FE order, and DtN harmonic capacity. Cross-check sphere cases with
`radia.acoustics.elastic_sphere_scattering` and the stiff limit with the rigid
sphere reference.
""",
    "cq": """# Lubich convolution quadrature

Use `radia.acoustics.cq.cq_soft_sphere_scattering`. Pin the method (BDF1/BDF2),
time step, number of samples, damping radius, convention `s=delta(zeta)/dt`, and
complex wavenumber `kappa=i*s/c`. Each CQ node is an independent
`ngsolve.bem.HelmholtzSL` solve. Validate frequency nodes against
`soft_sphere_scattering_complex_k`, then check conjugate symmetry, small
imaginary leakage after inverse FFT, causality, and time-step refinement.
""",
    "validation": """# Validation order

1. Inspect materials, boundaries, and finite-element spaces before solving.
2. Check analytic sphere solutions.
3. Check independent NGSolve/ngsolve.bem formulations. The pure double-layer
   validation solves `(1/2 I + K) mu = -u_inc` and compares `D mu` with the
   sound-soft sphere in `validation_test/acoustics/test_double_layer_bem.py`.
4. For FSI, refine mesh/order and test the stiff rigid-sphere limit.
5. For CQ, test every complex-frequency solve before the inverse FFT, then test
   reality, causality, and refinement.
Heavy convergence/timing studies belong in `validation_test/acoustics` on an
idle compute host; fast API contracts belong in `tests/acoustics`.
""",
}


def acoustic_usage(topic: str = "overview") -> str:
    """Return production acoustic guidance by topic."""
    key = topic.strip().lower().replace("-", "_")
    if key not in _TOPICS:
        raise ValueError(f"unknown topic {topic!r}; choose from {sorted(_TOPICS)}")
    return _TOPICS[key]


def acoustic_capabilities() -> dict:
    """Describe ownership and stable public entry points."""
    return {
        "schema": "radia-mcp.radia-acoustic-capabilities/v1",
        "owner": "radia.acoustics",
        "numerical_backends": ["NGSolve", "ngsolve.bem"],
        "matlab_runtime": False,
        "topics": sorted(_TOPICS),
        "apis": {
            "analytic": ["soft_sphere_scattering", "rigid_sphere_scattering", "fluid_sphere_scattering", "elastic_sphere_scattering"],
            "fsi": ["radia.acoustics.fsi.sphere_mesh", "radia.acoustics.fsi.fsi_dtn_solve"],
            "cq": ["radia.acoustics.cq.bdf_delta", "radia.acoustics.cq.cq_soft_sphere_scattering", "radia.acoustics.cq.soft_sphere_scattering_complex_k"],
        },
        "education_solver": "ksugahar/matlab-acoustic-fembem (separate)",
    }
