# Hierarchical Cauer + Warburg-Schur SIBC Termination

Reference implementation for the **Warburg-Schur termination of CLN
for SIBC** ── the canonical hierarchical Cauer ladder that takes the
bulk Cauer Ladder Network (CLN) of eddy-current diffusion and closes
it with a single non-rational $\sqrt{s}$ Warburg block via Schur
complement elimination.

## Theory

Given a bulk CLN reduced matrix $K_r(s)$ on Krylov basis $Q$ of size
$N$, appending one DOF carrying the Warburg impedance
$z(s) = (s+d)/(K_{\rm SIBC}\sqrt{s})$ and Schur-eliminating it yields

$$
Y_R(s) = Y_{\rm CLN}(s) + \frac{K_{\rm SIBC}\sqrt{s}}{s+d}
$$

with $K_{\rm SIBC} = S \sqrt{\sigma/\mu}$ ($S$ = surface area), and $d
> 0$ the single bulk-surface crossover constant. The Warburg block
carries the **AC $\sqrt{s}$ component** of SIBC (vanishing at DC via
the $(s+d)$ factor, leaving DC to bulk CLN) — clean DC/AC split with
exact preservation of $Y_R(0) = Y_{\rm CLN}(0)$ and asymptote
$Y_R(s) \to K_{\rm SIBC}/\sqrt{s}$ as $s \to \infty$.

Time-domain analysis becomes tractable by expanding the Warburg block
into a **diffusive Foster quantisation**
$\sum_{k=1}^{N_\xi} c_k/(s+\xi_k)$ ($N_\xi=50$) ── finite-dimensional
LTI realisation matching the analytic step response
$(2/\sqrt{\pi d})\,\mathrm{dawsn}(\sqrt{dt})$ within $1.8\times 10^{-3}$.

## Directory layout

| Subfolder | Contents | What it demonstrates |
|---|---|---|
| `frequency_domain/` | cylinder / sphere / cuboid 3D Y(s) verification | Hierarchical Cauer matches exact Bessel / FEM across 8 decades |
| `time_domain/` | Dawson closed form + diffusive Foster | Warburg expansion → time-domain LTI |
| `field_reconstruction/` | $H_z(\mathbf{r},t)$, $J(\mathbf{r},t)$ from same state variables | Same ROM gives port + spatial field |
| `multi_conductor/` | $N$-conductor mutual coupling via H-matrix | Paper 2 application |
| `engineering/` | PWM transient + NGSolve production + IH workpiece | Engineering applications |
| `nonlinear_esim/` | ESIM-coupled Karl iteration | Paper 2: nonlinear steel workpiece |

## Canonical reference: cylinder (`frequency_domain/circle_cpe_schur_plot.py`)

The minimal demo. Reproduces digest Fig. 1:
- Bessel exact $Y_{\rm cyl}(s) = \pi a^2 \sigma \cdot 2I_1(\gamma a)/[\gamma a I_0(\gamma a)]$
- Schur-augmented CLN3 (basis $4$ = 3 Krylov + 1 Warburg DOF, $d/2\pi = 1.58\times 10^5$ Hz)
- Wall-band peak rel-err $17\%$, sub-$0.01\%$ at DC and along the SIBC asymptote

```bash
cd frequency_domain && python circle_cpe_schur_plot.py
# → circle_cpe_schur.pdf
```

## Related publications

- **Sugahara et al., "Hierarchical Cauer: Warburg-Schur Termination of CLN for SIBC"**, IGTE 2026 digest (cylinder verification, 1-page)
- **Paper 1**: Theorem 1 (Schur uniqueness) + cuboid 3D verification + curvature/edge corrections
- **Paper 2**: Multi-conductor BEM-CLN coupling + IH applications + ESIM nonlinear

## References

[1] D. Givoli, J. B. Keller, "A Finite Element Method for Large Domains," CMAME 76(1), 41-66, 1989.
[2] A. Kameari et al., "Cauer Ladder Network Representation of Eddy-Current Fields for MOR Using FEM," IEEE TMag 54(11), 7202804, 2018.
[3] T. Matsuo, "Fractional-Power-Law Frequency Dependence in Magnetic Properties via Cauer Circuit," JMMM Art. 174229, 2026.
[4] J. E. B. Randles, Discuss. Faraday Soc. 1, 11 (1947).
