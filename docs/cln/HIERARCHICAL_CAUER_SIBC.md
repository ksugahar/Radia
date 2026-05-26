# Hierarchical Cauer: Warburg–Schur Termination for SIBC

Academic reference for the single-conductor hierarchical Cauer
construction that closes a bulk Cauer Ladder Network (CLN) at the
non-rational $\sqrt{s}$ Surface-Impedance-Boundary-Condition (SIBC)
asymptote by terminating it with a Warburg element via Schur
complement elimination.

Companion to the foundational [CAUER_LADDER_NETWORK.md](CAUER_LADDER_NETWORK.md)
(bulk CLN, Tanimoto–Kameari method) and to the multi-conductor
extension [BEM_CLN.md](BEM_CLN.md).

This document tracks the Paper 1 manuscript:
**"Hierarchical Cauer: Warburg–Schur Termination of CLN for SIBC"**
(Sugahara, Nagamine, Hane, 2026), and the IGTE 2026 1-page digest.

## 1. The CLN termination gap

A finite-stage Cauer ladder must terminate at its $N$-th rung either
inductively ($|Y|\!\sim\!f^{-1}$, L-term) or resistively
($|Y|\!\sim\!f^{0}$, R-term) [3].  Neither matches the SIBC tail
$|Y|\!\sim\!K_{\rm SIBC} f^{-1/2}$ which arises from surface
boundary-layer diffusion ── a process distinct from the bulk
diffusion the CLN ladder discretises.

The mismatch is **fundamental, not numerical**:
$\sqrt{s}\!\notin\!\mathbb{Q}(s)$ is non-rational, so no rational
ROM (Padé, Foster sum, Cauer ladder, AAA) of finite order can reach
the SIBC asymptote.  Increasing the rank $N$ improves accuracy near
DC but the post-$f_N$ behaviour stays anchored to one of the
$f^{-1}$ / $f^{0}$ types.

## 2. The Warburg–Schur construction (Theorem 1)

### 2.1 Statement

Let $K_r(s) = Q^\top K(s) Q$ be the bulk CLN reduced matrix on
Krylov basis $Q$ of size $N$, $b_r = Q^\top b$ the projected
source, and $Y_{\rm CLN}(s) = \sigma\,b_r^\top K_r^{-1} b_r$ the
bulk CLN admittance.

Append a single Warburg degree of freedom with impedance
$z(s) = (s+d)/(K_{\rm SIBC}\sqrt{s})$ and form the augmented system
$K_{\rm aug} = \mathrm{diag}(K_r, z)$, $b_{\rm aug} = (b_r, 1)^\top$.
Schur-eliminating the Warburg DOF onto the bulk system yields

$$
Y_R(s) = \sigma\, b_{\rm aug}^\top K_{\rm aug}^{-1} b_{\rm aug}
       = Y_{\rm CLN}(s) + \frac{K_{\rm SIBC}\sqrt{s}}{s+d}
$$

with $K_{\rm SIBC} = S\sqrt{\sigma/\mu}$, where $S$ is the perimeter
(2D) or surface area (3D), and $d > 0$ is the single bulk–surface
crossover constant.

### 2.2 Properties

| Limit | Value | Interpretation |
|---|---|---|
| $Y_R(0)$ | $Y_{\rm CLN}(0)$ exactly | bulk DC preserved (Warburg vanishes via $(s+d)$) |
| $Y_R(s\to\infty)$ | $K_{\rm SIBC}/\sqrt{s}$ | SIBC asymptote recovered |
| Continued-fraction structure | two-level Cauer | bulk ladder ($N$ rungs of $L$–$R$) terminated by Warburg = input impedance of semi-infinite uniform $R$–$C$ ladder |

The clean **DC/AC split** is the structural insight: bulk CLN
exclusively owns the DC limit; the Warburg block carries the
AC $\sqrt{s}$ component of SIBC.  This is not an empirical fit ── the
$(s+d)$ factor ensures $Y_W(0)\!=\!0$ as a matter of pole structure.

### 2.3 Minimality

Theorem 1 of the Paper 1 manuscript states that this is the
**minimal** ROM (in DOF count) that simultaneously preserves $Y(0)$
and reaches the $\sqrt{s}$ SIBC asymptote.  Replacing the single
Warburg DOF with multiple $\sqrt{s}$-blocks adds redundant DOFs
without improving the wall-band residual (verified empirically:
canonical CLN3+Schur basis 4 gives 17% wall-band peak on Cu cylinder;
Multi-K $M_0\!=\!3+4$pt+Schur basis 8 gives 17.5%; Multi-K $M_0\!=\!3+7$pt+Schur basis 11 gives 17.0% — same floor).

The wall-band peak residual is intrinsic to **any single-DOF
$\sqrt{s}$ closure on this geometry**, not a limitation of the
specific construction.

## 3. Time-domain analysis

The non-rational $\sqrt{s}$ would normally preclude finite-dimensional
state-space realisation.  Two complementary routes resolve this.

### 3.1 Closed-form step response

The Warburg block's step response is

$$
H(t) = \mathcal{L}^{-1}\!\left\{\frac{\sqrt{s}}{s(s+d)}\right\}
     = \frac{2}{\sqrt{\pi d}}\,\mathrm{dawsn}(\sqrt{dt})
$$

where $\mathrm{dawsn}(x) = e^{-x^2}\int_0^x e^{u^2}\,du$ is the Dawson
function.  Limits: $H(0) = 0$, peak at $t\approx 1/d$,
$H(\infty) = 0$.

### 3.2 Diffusive Foster quantisation

For numerical state-space realisation, expand

$$
\frac{\sqrt{s}}{s+d} \;\approx\; \sum_{k=1}^{N_\xi} \frac{c_k}{s + \xi_k}
$$

where $\xi_k$ are log-spaced around $d$ over 8 decades and $c_k$
are least-squares residues fit on the $j\omega$ axis.  $N_\xi = 50$
gives frequency-domain max rel-err $1.1\times 10^{-4}$, time-domain
step-response max rel-err $1.8\times 10^{-3}$ vs the Dawson
closed-form.

This is the operational route in `examples/hierarchical_cauer_sibc/`:
`diffusive_quadrature.py` constructs the expansion, downstream scripts
treat the Warburg block as a 50-state LTI subsystem indistinguishable
from a rational ROM at the engineering level.

## 4. Verification

| Geometry | Reference | Result | Script |
|---|---|---|---|
| 2D Cu cylinder ($a\!=\!5$ mm) | exact Bessel | basis 4: 17% wall band, sub-0.01% asymptote | `frequency_domain/circle_cpe_schur_plot.py` |
| 3D Cu sphere | Bessel-like closed form | $K_{\rm SIBC} = 4\pi R^2\sqrt{\sigma/\mu}$ recovered | `frequency_domain/3d_sphere.py` |
| 3D Cu cuboid | Mellin asymptote ($c_0$, $c_1$, $c_2$) | $K_{\rm SIBC} = 2(ab\!+\!bc\!+\!ca)\sqrt{\sigma/\mu}$ + edge corrections | `frequency_domain/3d_cuboid.py` |
| 3D cuboid via NGSolve FEM | itself (production) | Hierarchical Cauer extracted from full FEM Y(s) sweep | `frequency_domain/ngsolve_cuboid_Y.py`, `engineering/production_ngsolve_cube_kelvin.py` |
| PWM transient (cylinder) | independent 1D FDM PDE | rel-err $7\times 10^{-8}$ port; $7\times 10^{-4}$ internal $H_z(r,t)$ | `engineering/pwm_transient_field.py` vs `engineering/pde_pwm_reference.py` |

## 5. Field reconstruction

The same state variables used to compute port admittance also yield
the full spatial field via linear combination of mode functions:

$$
H_z(\mathbf{r}, t) = H_{\rm ext}(t) - \sum_{k=1}^{N} m_k(t)\,\phi_k(\mathbf{r}) - \sum_{j=1}^{N_\xi} n_j(t)\,\psi_j(\mathbf{r})
$$

where $\phi_k$ are the **bulk volume modes** (Bessel
$J_0(j_{0,k}r/a)$ for cylinder; FEM Krylov basis columns $q_k$ for
general 3D) and $\psi_j$ are the **Warburg surface modes**
$\psi_j(r) = e^{-(a-r)/\delta_j}$ with $\delta_j = 1/\sqrt{\xi_j\mu\sigma}$.

This distinguishes the hierarchical Cauer construction from standard
ROMs that return port quantities only.  See
[`examples/hierarchical_cauer_sibc/field_reconstruction/`](../../examples/hierarchical_cauer_sibc/field_reconstruction/)
for cylinder, sphere, and 3D FEM-Krylov demos.

## 6. Curvature and edge corrections (3D)

For 3D geometries the leading $K_{\rm SIBC}\!=\!S\sqrt{\sigma/\mu}$ is
universal, but second-order accuracy on curved or sharp-cornered
bodies requires HOIBC-style corrections:

- **Curvature**: Senior–Mitzner tower
  $\gamma_1\!=\!-H$, $\gamma_2 = (K_{\rm Gauss}\!-\!H^2)/2$,
  $\gamma_3 = -H^3 + H K_{\rm Gauss} + \lambda_3 \Delta_S H$
  (Rytov principal-direction form).  Implemented in
  `src/radia/netgen_mesh_curvature.py`.
- **Sharp edges**: Mellin-transform asymptote for cuboid right-angle
  dihedral gives $c_2 = +48/(\pi \mu^{3/2}\sqrt{\sigma})$, edge-length
  independent.  Implemented in supplement Mathematica notebooks.

These corrections enter as multiplicative factors / additive terms
on $K_{\rm SIBC}$, preserving the overall Warburg–Schur structure.

## 7. Scope and limitations

**Port-driven vs volume-source**: Theorem 1's $r(f)\!\to\!1$ claim
implicitly assumes that $Y_{\rm CLN}(s\to\infty)\to 0$, which holds
for port-driven admittance problems (Paper 1 §V cuboid benchmark)
but **not** for volume-source problems with uniform $J_0$ forcing
(e.g., Hiruma 2023 §III.A).  Volume-source problems leave an FE
residual offset $Y_{\rm CLN}(\infty)\!=\!\sigma(\text{area} - b_r^\top M_r^{-1} b_r)$
that decays only algebraically with mesh refinement.  For such
problems, FE-level treatment (XFEM, Hiruma 2023) is recommended
over the augmented-CLN ROM.

**High-stage rank**: Foster pole orthogonalisation breaks down past
$N\!\approx\!4$ in float64 (Hankel matrix condition number exceeds
$10^{15}$, see [CAUER_LADDER_NETWORK.md §6](CAUER_LADDER_NETWORK.md)).
The hierarchical Cauer construction is most useful at low-$N$ +
Schur block: basis 4 (CLN3 + 1 Warburg) is the sweet spot.

**Single-DOF $\sqrt{s}$ floor**: The 17% wall-band residual on the
Cu cylinder cannot be reduced below this floor by adding Krylov
expansion points; it is intrinsic to the geometric form and is the
price paid for the asymptote-correct rational+irrational
factorisation.  Multi-K is a verified equivalent (see
`frequency_domain/circle_multik_schur_plot.py` for the empirical
floor demonstration).

## 8. References

[1] D. Givoli, J. B. Keller, "A Finite Element Method for Large Domains," Comput. Methods Appl. Mech. Engrg. 76(1), 41–66, 1989.
[2] A. Kameari et al., "Cauer Ladder Network Representation of Eddy-Current Fields for Model Order Reduction Using FEM," IEEE Trans. Magn. 54(11), 7202804, 2018.
[3] T. Matsuo, "Representation of fractional-power-law frequency dependence in magnetic properties using a Cauer circuit," J. Magn. Magn. Mater. Art. 174229, 2026.
[4] J. E. B. Randles, "Kinetics of Rapid Electrode Reactions," Discuss. Faraday Soc. 1, 11–19, 1947.
[5] H. Köster, A. König, O. Bíró, "Cauer Ladder Network as Proper Generalised Decomposition Realisation," IEEE Trans. Magn. 57(6), 7401204, 2021.
[6] R. Kuriyama et al., "Multi-expansion-point Krylov for Eddy-Current ROMs," IEEJ Trans. PE 139(7), 555–562, 2019.
[7] H. Sugahara, H. Nagamine, Y. Hane, "Hierarchical Cauer: Warburg–Schur Termination of CLN for SIBC," IGTE 2026 (1-page digest); Paper 1 manuscript in preparation.

## 9. See also

- [CAUER_LADDER_NETWORK.md](CAUER_LADDER_NETWORK.md) — bulk CLN foundations (Tanimoto–Kameari)
- [BEM_CLN.md](BEM_CLN.md) — multi-conductor BEM-style assembly (Paper 2)
- [CLN_3D_CUBOID.md](CLN_3D_CUBOID.md) — 3D cuboid Cu benchmark
- [`examples/hierarchical_cauer_sibc/`](../../examples/hierarchical_cauer_sibc/) — canonical reference implementations
- `radia-mcp` topic [`cln_sibc_orthogonal`](../../packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/cln_sibc_orthogonal.py) — machine-readable knowledge

---
**Status**: Cylinder verified (digest), 3D sphere / cuboid verified
(supplement scripts), Paper 1 manuscript in preparation.
**Last updated**: 2026-05-26.
