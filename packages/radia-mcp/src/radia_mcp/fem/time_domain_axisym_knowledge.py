"""Time-domain FEM, axisymmetric (Henrotte), harmonic balance, high-frequency."""

TIME_DOMAIN_FEM = r"""
# Time-domain FEM (TD-FEM)

For transient EM analysis: solve A(t) directly with time-stepping.

## Time stepping schemes

| Scheme | Order | Implicit? | Stability |
|--------|-------|-----------|-----------|
| Backward Euler | 1 | yes | unconditionally stable |
| Crank-Nicolson | 2 | yes | unconditionally stable, may oscillate |
| BDF2 (Backward Differentiation) | 2 | yes | A-stable |
| RK4 | 4 | no (explicit) | CFL-limited |
| Newmark-beta | varied | yes/no | tunable damping |

## In MQS context

For the discrete system
    M_sigma dA/dt + K(A) A = f(t)
with M_sigma = mass matrix in conductor, K(A) = nonlinear stiffness.

Implicit Euler:
    (M_sigma + Δt * K(A_new)) A_new = M_sigma A_old + Δt * f
    → solve at each time step (Picard or Newton for nonlinear)

## Lab use

- `radia.calc_motor_transient` (transient motor solver)
- Hollaus MSFEM with TD coupling
- See `radia_mcp.motor.calc_motor_transient`

## References (this folder)

[LOCAL] `10_FEM_定式化/22_時間領域FEM/` (2 files)
"""


HENROTTE_AXISYM = r"""
# Henrotte basis for axisymmetric magnetic (★ lab CORE)

Reference: Henrotte 1993 IEEE TMAG axisymmetric — see
`radia_mcp.radia_ngsolve.axifemm_documentation` and CLAUDE.md
"Axisymmetric FE: Henrotte for Magnetic, Standard H1 for Scalar".

## The problem

Axisymmetric magnetic A_phi (curl-curl) has integrand:

    B_z = (1/r) ∂(r A_phi)/∂r

This **1/r** kernel cannot be integrated accurately by standard Gauss
quadrature near the axis r=0. Result: O(h) instead of O(h^2)
convergence in standard FE, and outright wrong B near axis.

## Henrotte's solution

Use basis {1, r^2, z} (instead of standard {1, r, z}) and substitute
s = r^2. Then the 1/r kernel **disappears analytically**, integration
becomes exact.

## API in radia

```python
import radia.axifem as ax

mesh = Mesh(...)                       # axis-aligned (r, z) mesh
fes = ax.H1Henrotte(mesh, order=p)     # p = 1 (Q1) or p = 2 (Q2)

a_mag = BilinearForm(fes, symmetric=True)
a_mag += ax.AxiHenrotteStiffnessBFI(mu_cf)
a_mag += ax.AxiHenrotteSigmaMassBFI(sigma_cf)   # eddy-current term
```

## Why NOT use Henrotte for axisymmetric scalars

Per CLAUDE.md POLICY (2026-05-10): scalar Laplacians use **standard H1
with 2*pi*r weighting** (FEMM canonical). The Jacobian is smooth, no
1/r kernel, so standard FE is fine and matches FEMM 4.2 production
convention.

## References

- CLAUDE.md "Axisymmetric FE" section (POLICY 2026-05-10)
- `docs/axifemm/FORMULATION.md` (lab repo) — full derivation
- `memory/reference_femm_source_axisym_conventions.md`
- [LOCAL] `10_FEM_定式化/21_軸対称/` (1 file)

## Cross-reference

- `radia_mcp.radia_ngsolve.axifemm_documentation` — complete API + examples
- `radia_mcp.mathematica.recipes.HENROTTE` — symbolic derivation
"""


HARMONIC_BALANCE = r"""
# Harmonic Balance (HB) method — periodic steady-state without time-stepping

For periodic steady-state with nonlinear material (e.g. iron under AC
excitation), HB avoids time-stepping convergence by solving for
**Fourier coefficients directly**.

## Formulation

Expand A(t) = Σ_k A_k e^{i k ω t}, substitute into nonlinear equation:

    Σ_k M_sigma · (ikω A_k) e^{ikωt} + K_nonlin(A) = Σ_k f_k e^{ikωt}

Apply Galerkin in time (or balance harmonics) → block matrix system for
{A_0, A_1, ..., A_N}.

## Advantages

- ✓ Direct steady-state (no transient ramp-up)
- ✓ Capture multiple harmonics (lab: 1st + 3rd for iron saturation)
- ✓ Continuation in amplitude possible (gradual approach to saturation)

## Disadvantages

- ✗ Block matrix is (N+1) × (N+1) sub-blocks → grows fast with N
- ✗ Convergence harder than transient (continuation often needed)
- ✗ Implementation more complex

## Lab status

- Reference layer; not yet a production panel
- Could integrate with `radia.calc_motor_transient` for steady-state
- See research papers in `10_FEM_定式化/23_ハーモニックバランス法/`

## References (this folder)

[LOCAL] `10_FEM_定式化/23_ハーモニックバランス法/` (2 files, 70 MB)
"""


HIGH_FREQUENCY = r"""
# High-frequency FEM (Helmholtz / full Maxwell)

Note: per CLAUDE.md POLICY, lab is **Laplace kernel only** (MQS/Darwin).
High-frequency Helmholtz / full Maxwell is OUT OF SCOPE for radia.

The 31 files / 3.7 GB in `10_FEM_定式化/24_高周波/` are reference layer
only, for users who need full-wave analysis (use ngsolve.bem or
specialized HF tools).

## Why HF is out of scope

1. **Computational cost**: HF mesh requires h < λ/10 = c/(10 f)
   At 1 GHz: λ = 30 cm → mesh ~ 3 cm → millions of DOFs for any object
2. **Kernel choice**: Helmholtz e^{-jkr}/r vs Laplace 1/(4πr)
   Radia chose Laplace for efficiency in DC-MHz range
3. **Tooling**: ngsolve.bem already supports Helmholtz; no reason to
   duplicate

## When users need HF

→ Use `ngsolve.bem` with `kernel='helmholtz'`, k = ω/c
→ Or specialized HF tools (HFSS, FEKO, CST)
→ See `radia_mcp.bem.surface_ie.cfie` for HF surface IE theory

## References (this folder, reference only)

[LOCAL] `10_FEM_定式化/24_高周波/` (31 files, 3.7 GB)
        — PML, ABC, EFIE/MFIE/CFIE, scattering, antennas
"""


CIRCUIT_COUPLING = r"""
# FE-circuit coupling

For motor / transformer simulation, the FE field problem must be
coupled to the external circuit (winding voltage/current, switches,
load impedance).

## Approaches

| Approach | When |
|----------|------|
| **Strong coupling** (monolithic) | Field unknowns + circuit currents in same system |
| **Weak coupling** (staggered) | Iterate between FEM solve and circuit solve |
| **Equivalent-circuit FEM** | Reduce field to L(I, theta) and Lookup table |
| **Schur complement** | Eliminate field DOFs at coil interface → impedance |

## Hanser 2025 Schur voltage drive (lab core)

For Schöbinger-Hanser voltage-driven coil modeling:
- Compute effective coil impedance Z_eff(ω) via FEM Schur complement
- Couple Z_eff to circuit time-stepping
- See `radia.calc_motor_lamination v3+` and `radia_mcp.motor.hollaus_eddy`

## References (this folder)

[LOCAL] `10_FEM_定式化/25_回路連携/` (2 files)

## Cross-reference

- `radia_mcp.motor` — Hanser 2025 Schur, Hollaus DEIM
- `radia_mcp.mor` — model order reduction for circuit coupling
- `radia_mcp.peec` — circuit-side topology for PEEC+FEM hybrid
"""


def get_time_domain_axisym_knowledge(topic: str = "henrotte_axisym") -> str:
    """Dispatch time-domain + axisymmetric + harmonic balance + HF topics.

    Topics:
        time_domain_fem    - TD-FEM time-stepping schemes
        henrotte_axisym    - ★ Henrotte basis for axisym magnetic (lab CORE)
        harmonic_balance   - HB method for periodic steady-state
        high_frequency     - HF Helmholtz (out of lab scope, reference)
        circuit_coupling   - FE-circuit coupling (Hanser Schur)
        all                - Everything
    """
    topic = topic.lower().strip()
    if topic in ("time_domain_fem", "td_fem", "tdfem", "time_domain"):
        return TIME_DOMAIN_FEM
    if topic in ("henrotte_axisym", "henrotte", "axisym", "axisymmetric"):
        return HENROTTE_AXISYM
    if topic in ("harmonic_balance", "hb"):
        return HARMONIC_BALANCE
    if topic in ("high_frequency", "hf", "helmholtz"):
        return HIGH_FREQUENCY
    if topic in ("circuit_coupling", "circuit", "coupling", "schur"):
        return CIRCUIT_COUPLING
    if topic == "all":
        return "\n\n".join([TIME_DOMAIN_FEM, HENROTTE_AXISYM, HARMONIC_BALANCE,
                            HIGH_FREQUENCY, CIRCUIT_COUPLING])
    return (f"Unknown topic '{topic}'. Available: time_domain_fem, "
            "henrotte_axisym, harmonic_balance, high_frequency, "
            "circuit_coupling, all.")
