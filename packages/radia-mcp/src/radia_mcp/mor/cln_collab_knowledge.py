"""CLN collaborator work + form-conversion + application examples.

Covers W:\\30_CauerLadderNetwork\\{2021_01_22_CauerI_to_CauerII,
2022_10_09_遠藤@法政, 2023_08_08_松本@法政, 2026_04_01_長方形CLN}
plus the 2017_05_17_インバータ回路.docx application.

External collaborators are HOSEI UNIVERSITY -- Endo (2022) and
Matsumoto (2023). 2021 self-work was the CauerI <-> CauerII form
conversion (Lanczos tridiagonalization route).

Companion to `cln_knowledge.py` (canonical CLN method) and the
upcoming `cln_systematic_knowledge.py` (folders 01-16).
This file is Group E in the Sugahara-lab CLN corpus absorption.
"""

OVERVIEW = """
# Group E -- CLN form conversion, external collaborators, applications

This module collects four threads of CLN work that do not fit the
"core method" (cln_knowledge.py) and do not fit the systematic
problem catalog (01_単一電源 ... 16_FP_CLN, cln_systematic).

| Year | Folder                            | Topic                                            | Role        |
|------|-----------------------------------|--------------------------------------------------|-------------|
| 2017 | 2017_05_17_インバータ回路.docx     | CLN -> SPICE inverter coupling (application)     | self        |
| 2021 | 2021_01_22_CauerI_to_CauerII      | CauerI <-> CauerII form conversion via Lanczos   | self        |
| 2022 | 2022_10_09_遠藤@法政              | 2D vector-potential CLN, COMSOL coupling          | Endo @ Hosei|
| 2023 | 2023_08_08_松本@法政              | (slim) follow-up handoff slides                  | Matsumoto @ Hosei|
| 2026 | 2026_04_01_長方形CLN              | Rectangular-cross-section CLN (planned)          | TBD         |

## Why Group E matters

Most of the Sugahara-lab CLN corpus is solver method (cln_knowledge)
plus a catalog of problem types (cln_systematic).  Group E is the
*infrastructure layer* underneath both:

1. **Form conversion** (2021).  All CLN derivations naturally produce
   a "CauerI" ladder (series R, parallel L, E-mode first).  Many
   downstream uses -- SPICE export, dual representations, magnetic
   energy bookkeeping -- need the equivalent "CauerII" form (series L,
   parallel R, H-mode first).  The conversion is a *Lanczos*
   tridiagonalization in a non-standard inner product.  Without this
   technology, the rest of the lab's CLN work cannot interoperate
   with the SPICE / circuit-simulator world.

2. **External collaborators** (2022 Endo, 2023 Matsumoto).  Both are
   from Hosei University.  Endo brought a 2D multi-conductor +
   surrounding cylinder problem, which became the lab's first
   *COMSOL-coupled* CLN demo (MATLAB drives COMSOL via the
   LiveLink-MATLAB API).  This is also the only CLN folder in which
   the "source reconstruction" mode -- recovering A0 from the CLN
   coefficients -- is exercised in 2D vector-potential form.

3. **Application coupling** (2017 inverter).  The CLN's selling point
   is "milliseconds in the time domain for any input waveform".  The
   inverter docx pairs the CLN as a passive subcircuit fed by an IGBT
   bridge.  This is the prerequisite for the 07_電気力学 stream in
   the systematic catalog.

4. **TBD: rectangular CLN** (2026).  Folder is empty.  Listed here
   so the next agent will not miss it.
"""

CAUER_I_VS_II = r"""
# CauerI vs CauerII -- WHY two ladder forms exist

## Continued-fraction expansion of the impedance

For a 1D skin-effect problem (plate of thickness d, conductivity sigma,
permeability mu) the per-unit-width impedance is

    Z(s) = j omega mu * (2 / (k d)) * tan(k d / 2) * d
    k    = sqrt(-j omega mu sigma) = sqrt(s mu sigma)         (s = j omega)

There are *two* canonical continued-fraction expansions of any rational
(or asymptotically rational) Z(s):

* **CauerI** -- expansion around `s = 0` (DC).
  The first element is whatever is finite at s=0 (a resistance R for a
  resistive skin problem).  The structure is

      Z(s) = R0 + 1 / (s L1 + 1 / (R2 + 1 / (s L3 + 1 / (R4 + ...))))

  Realized as a *series R - parallel L* ladder:

        R0      R2      R4
      o--+---+--+---+--+---+--- ...
             |      |      |
            L1     L3     L5
             |      |      |
      o-----+------+------+--- ...

  Even-indexed elements are resistors, odd-indexed are inductors.
  CauerI is what falls out *directly* from the lab's recursive FEM
  formulation when starting with the E-mode (electric basis E0 first):
  R is the "first thing you see" at DC.

* **CauerII** -- expansion around `s = infty` (high frequency).
  The first element is whatever dominates as s -> infty (an inductance
  L for the same skin problem).  The structure is

      Z(s) = s L0 + 1 / (1/R1 + 1 / (s L2 + 1 / (1/R3 + 1 / (s L4 + ...))))

  Realized as a *series L - parallel R* ladder:

        L0      L2      L4
      o--+---+--+---+--+---+--- ...
             |      |      |
            R1     R3     R5
             |      |      |
      o-----+------+------+--- ...

  Even-indexed are inductors, odd-indexed are resistors.
  CauerII is what naturally appears starting from the H-mode (magnetic
  basis H0 first): L is the "first thing you see" at high frequency.

## When to use which

| Scenario                                       | Preferred form |
|------------------------------------------------|----------------|
| Magnetic energy bookkeeping                    | CauerII        |
| SPICE export with an inductor at the input     | CauerII        |
| Resistive DC bookkeeping / Joule loss          | CauerI         |
| SPICE export with a resistor at the input      | CauerI         |
| Best low-frequency approximation @ small N     | CauerI         |
| Best high-frequency approximation @ small N    | CauerII        |
| Match a pre-existing measured Y-parameter file | depends on form|

The "best at small N" rows are the substantive point: truncating
CauerI to N stages preserves the moments around s=0 (so the DC and
low-frequency asymptotes are exact), and truncating CauerII preserves
the moments around s=infty (so the high-frequency asymptote is exact).
Whichever range you care about, pick the corresponding form.

## Why they are equivalent at fixed N (untruncated)

Both ladder forms are realizations of the same rational function Z(s)
of degree N.  Two rational functions of the same numerator/denominator
degree that match on infinitely many points must be identical.  The
recursion in either form is therefore a Euclidean-style division of
two polynomials, and the two forms are linked by a change of basis
on the underlying state-space matrices (K, N) -- which is exactly
what the Lanczos algorithm computes.
"""

CAUER_FORM_CONVERSION_CODE = r"""
# CauerI -> CauerII conversion -- Lanczos in non-standard inner product

## The matrix problem

Start from the CauerI matrix system at frequency s:

    (RR + s LL) * I = V,         V = [1, 0, 0, ...]^T
    Z(s) = 1 / I(1)

where RR is *diagonal* (CauerI: R0, R2, R4, ...) and LL is *tridiagonal*
(CauerI: after the +/-1 difference operator U is applied to diag(L1,
L3, L5, ...), then U^T diag(L) U is tridiagonal symmetric).

The CauerII representation is the same impedance function with the
roles swapped -- *RR'* tridiagonal, *LL'* diagonal:

    (R_tri + s L_diag) * I' = V',   Z(s) = 1 / (Q I')(1)

The change-of-basis matrix Q is *not* orthogonal in the usual sense;
it is orthogonal *in the K inner product*:

    Q^T K Q = R_diag    (diagonal)
    Q^T N Q = L_tridiag (tridiagonal)
    Q^T Q   != I        (NOT identity)

with K = LL (CauerI inductance matrix) and N = RR (CauerI resistance
matrix).  This is exactly the "two-matrix simultaneous reduction"
that the Lanczos three-term recurrence solves.

## The recurrence (MATLAB, from W:\\30_...\\MATLAB\\lanczos.m)

```matlab
% K is to be diagonalized; N is to be tridiagonalized; v0 is starting vector
w_o(:,1)   = K \ v0;
alpha(1)   = 1 / (w_o(:,1)' * K * w_o(:,1));
w_e(:,1)   = -alpha(1) * w_o(:,1);

for n = 1:N_iter
    beta(n)    = w_e(:,n)' * N * w_e(:,n);
    w_o(:,n+1) = w_o(:,n) + (K \ (N * w_e(:,n))) / beta(n);
    alpha(n+1) = 1 / (w_o(:,n+1)' * K * w_o(:,n+1));
    w_e(:,n+1) = w_e(:,n) - alpha(n+1) * w_o(:,n+1);
end

% Reading off:
%   alpha == CauerII series inductances (L0, L2, L4, ...)
%   beta  == CauerII parallel resistances (R1, R3, R5, ...)
%   Q(:,n) = alpha(n) * w_o(:,n)  -> change-of-basis matrix
```

## Mapping back to CauerII

After the recurrence,

    L_CauerII(k) = alpha(k)        for k = 1, 2, ..., N
    R_CauerII(k) = beta(k)         for k = 1, 2, ..., N

and the impedance can be computed either from the transformed matrices

    zz = L_tridiag + s * R_diag      % NB: roles swapped from CauerI
    zz * II = vv                     % vv = [1; 0; ...]
    Z(s) = 1 / (Q * II)(1)

or directly by walking the CauerII ladder structure with L_CauerII /
R_CauerII as element values.

## Equivalence verification (untruncated)

`verify_cln_i_ii_equivalence.py` in
W:\\30_CauerLadderNetwork\\2021_01_22_CauerI_to_CauerII\\ verifies
the round-trip CauerI -> Lanczos -> CauerII -> impedance against the
direct CauerI impedance:

| N stages | Max relative error in Z(s) over frequency sweep |
|----------|-------------------------------------------------|
| 3        | 1.4e-11   (machine precision)                   |
| 5        | 3.4e-11   (machine precision)                   |
| 7        | 6.1e-11   (machine precision)                   |
| 10       | 6.4e-06   (drifting -- Lanczos loss-of-orthog)  |
| 15       | divergent (catastrophic loss of orthogonality)  |

This is the *classical* finite-precision Lanczos pathology -- the K
inner product loses orthogonality as N grows.  Two practical fixes
are well known:

* re-orthogonalize at every step against all previously computed
  w_o / w_e (full reorthogonalization, O(N^2 * dim) work);
* use the Q-route: compute Q explicitly and report
  `R_diag = Q^T K Q`, `L_tridiag = Q^T N Q` -- the projection cancels
  most of the accumulated drift even when alpha / beta themselves
  drift.  The verification script confirms the Q-route stays at
  ~1e-11 even at N = 15 where reading alpha / beta directly diverges.

## Practical limit

For lab work the conversion is reliable up to N ~ 7.  For larger N
prefer one of:

(a) Run CauerI directly (no conversion).
(b) Build CauerII directly from the FEM recursion starting in the
    H-mode rather than via Lanczos conversion.
(c) Reorthogonalize at every step (slow but stable).

The script set in 2021_01_22_... documents (a) and (b) -- e.g.
`Analytical_CLN_for_H_mode_for_iron_sheet.m` derives the CauerII
parameters symbolically from the H-mode recursion, bypassing the
numerical Lanczos route entirely.

## CauerI parameter formulas for the lab's two canonical conductors

Lifted verbatim from `CauerI_to_CauerII.m`:

**Plate (thickness d, E-mode CauerI):**
```
R(1) = 0        % or 1e-16 to avoid division
R(n) = (4n - 5) * 4 / (sigma * d)     for n >= 2
L(n) = d * mu / (4n - 3)              for n >= 1
Z_exact(s) = s mu * (2 / (k d)) * tan(k d / 2) * d,  k = sqrt(s mu sigma)
```

**Round wire (radius a, E-mode CauerI):**
```
R(n) = (2n - 1) / (sigma * pi * a^2)
L(n) = mu / (8 n pi)
Z_exact(s) = k / (2 pi a) * J0(k a) / J1(k a) / sigma
```

The plate formula is what you would derive by hand from the H-mode
diffusion equation in the slab; the wire formula uses Bessel J0/J1
of complex argument.  Both are *closed-form CauerI ladders* -- no FEM
needed -- and are the gold reference for verifying any new CauerI
solver in the lab.

## Thermal-diffusion CauerI (subfolder 熱伝導/)

`熱伝導/E0_mode_analytical_Hiruma.m` shows that the *same* Lanczos /
ladder machinery applies to the heat-diffusion equation on a disk of
radius a with Dirichlet T(a)=0 and symmetry dT(0)=0:

    -alpha * (1/r) d/dr [r dT/dr] = source

The R/L pair is replaced by a *thermal* R = 1/lambda(2n-1) and
*thermal* L = lambda(2n).  This is the Hiruma-method connection
between CLN (electromagnetic) and Cauer-ladder reduced-order
thermal models -- relevant to the 12_比留間法 folder of the
systematic catalog.
"""

ENDO_HOSEI = r"""
# 2022_10_09 -- ENDO @ Hosei University collaboration

## What problem

A 2D quasi-static eddy-current configuration in vacuum:

* Four square aluminum conductors of side 30 mm placed at
  +/- 30 mm in x and +/- 35 mm in y (square pattern with a vertical
  asymmetry to allow a moving-source offset study).
* A central circular aluminum conductor of radius `a` (a = 100 mm or
  120 mm in different variants).  Surrounded by an air shell of
  radius `a_out = 1.2 * a`.
* Material: sigma = 37e6 S/m (Endo problem) or 57e6 S/m
  (`CLN_VectorPotential_solver.m`).  Both are "aluminum-like".

## What was computed

Five passes per geometry, all driven from MATLAB:

1. **`CLN_VectorPotential_comsol.m`** -- builds the geometry,
   physics (Helmholtz / induction-current), and mesh in COMSOL
   via the LiveLink-MATLAB API.  Saves a COMSOL `.mph` file.

2. **`VectorPotential.m`** -- analytical vector potential A_z(x, y)
   from a pair of uniform-current rectangular conductors (closed-form
   in terms of logs and arctangents).  Used as the *incident*
   source field A_0 that the CLN must reproduce / reconstruct.

3. **`CLN_VectorPotential_solver.m`** -- main driver:
   * loads precomputed CLN basis fields A_z^{(n)} and current
     densities J_z^{(n)} from `CLN/CLN_VectorPotential.mat`
     (these come from the FEM-based CLN recursion on the same mesh
     as the COMSOL model);
   * for each frequency f in {1, 50, 150, 300, 1000} Hz and each
     vertical offset of the moving source in {-20, -15, ..., +20} mm,
     projects the source A_z onto each basis J_z^{(n)} to get the
     CauerI "voltage source" coefficients V_n;
   * runs the CLN in the frequency domain (`cln = f_domain(cln,f,V)`)
     for nStage in {2, 3, 4, 5};
   * reconstructs A_z(x, y, f) = A_z^{source} - j omega A_z^{CLN}
     and saves the field for plotting.

4. **`COMSOL_f_domain.m`** -- runs the same configuration directly
   as a frequency-domain COMSOL solve, varying offset and frequency,
   to provide the *reference* against which the CLN reconstruction
   is compared.

5. **`plot_contour_A_*.m`** -- 2D contour overlays of:
   * A0 source reconstruction (analytical input vs CLN basis fit)
   * A_z from CLN (any frequency, any nStage, any offset)
   * A_z from COMSOL frequency-domain
   At identical color scales so visual error is the side-by-side
   comparison.

The PDF set `回路図(1).pdf` ... `回路図(6).pdf` documents the
SPICE-style circuit diagram of the CLN at each nStage, and
`4段_CLN_I.pdf` / `4段_CLN_II.pdf` show the four-stage CauerI and
CauerII variants for the same problem.

The pptx set (`2022_11_03_CLN.pptx`, `2022_11_16_CLN.pptx`,
`CLN法の計算手順_0728.pptx`) is the working draft of the
calculation-procedure handout that was eventually given to Endo.

## What was the research question

The lab side wanted to test:

* whether the CLN -- which had been validated only on the 1-source
  cases up to that point -- could handle a *moving* multi-source
  configuration without re-deriving the basis;
* whether the CauerII form (preferred by Endo's circuit-simulator
  workflow) gave the same answer as CauerI for the multi-source 2D
  case;
* whether ~5 stages were sufficient to recover the 50 Hz - 1 kHz
  response of the central conductor at 5% accuracy.

The MATLAB / COMSOL trip-up loop above is the answer:
overlay the contour plots, eyeball at each offset / frequency, count
where CLN diverges from COMSOL.

## Key files (Endo)

| File                                    | Purpose                                |
|-----------------------------------------|----------------------------------------|
| `CLN_VectorPotential_solver.m`          | Main driver; projects source onto CLN  |
| `CLN_VectorPotential_comsol.m`          | COMSOL geometry / physics / mesh build |
| `VectorPotential.m`                     | Closed-form A_z of moving source pair  |
| `COMSOL_f_domain.m`                     | Reference frequency-domain COMSOL run  |
| `plot_contour_A_CLN.m`                  | A_z reconstruction from CLN coefficients|
| `plot_contour_A_COMSOL_f_domain.m`      | A_z from COMSOL (reference)            |
| `plot_contour_A0_source_reconstruction.m`| Visual check of A_0 basis projection  |
| `2022_11_16_CLN.pptx`                   | Handout for Endo (latest version)      |
| `回路図(1..6).pdf`                       | CLN circuit diagrams 1-6 stages        |
| `4段_CLN_I.pdf` / `4段_CLN_II.pdf`       | Side-by-side I / II form at 4 stages   |
| `CLN/offset_*/`                         | Saved A_z(.mat) per offset, nStage, f  |
| `COMSOL/offset_*/`                      | Saved A_z(.mat) from COMSOL reference  |

## Lessons inherited from this collaboration

1. **External-collaborator interface convention.**  The lab provides
   the closed-form A_0 source (here `VectorPotential.m`) AND the
   precomputed CLN basis as a `.mat` -- the collaborator does not
   re-derive either.  They only run the projection and the f-domain
   solve.

2. **COMSOL LiveLink-MATLAB is the chosen coupling**, not Java API or
   `mphsave`/`mphload` round-trip.  All geometry / physics / mesh /
   solve commands are issued from MATLAB so the parametric sweep
   (5 frequencies x 9 offsets x 4 nStage = 180 cases) is one
   `for`-nest, not 180 GUI clicks.

3. **The CLN basis is mesh-locked.**  `.mat` carries (Jz, Az, Rn, Ln)
   tied to one specific COMSOL mesh.  If Endo re-meshes for any
   reason, the projection step has to be redone (the basis must live
   on the same nodes as the new source field).  This is the source of
   most "wrong answer" complaints from external collaborators -- the
   .mat looks correct but the mesh underneath has changed.

4. **The "source reconstruction" mode** (just project A_0 onto the
   CLN basis without solving) is a cheap sanity check: if A_0 cannot
   be represented by the basis at the chosen nStage, the f-domain CLN
   answer is meaningless.  Run it first, look at
   `plot_contour_A0_source_reconstruction.m` output, then proceed.
"""

MATSUMOTO_HOSEI = r"""
# 2023_08_08 -- MATSUMOTO @ Hosei University handoff

## What is here

The folder is intentionally thin:

* `補足スライド.pdf` -- supplementary slide deck (handoff document).
* `Mesh_Node243717_Element485552.txt` -- mesh dump, 243,717 nodes /
  485,552 elements.
* `Mesh_Node370882_Element_739882.txt` -- larger mesh dump, 370,882
  nodes / 739,882 elements.

No MATLAB or COMSOL source.  No new method.

## What it represents

Matsumoto picked up the Endo project line in 2023.  The role of this
folder is the *baseline mesh transfer* from the lab to Hosei plus a
slide deck explaining the calculation procedure (a follow-on of
`2022_10_09_遠藤@法政\CLN法の計算手順_0728.pptx`).

The two mesh sizes (~244k nodes / ~371k nodes) are the medium-
and high-density variants of the multi-conductor 2D problem from
Endo's setup, scaled up so the CLN basis fields can be tested at
finer resolution.  No CLN .mat or solver outputs were preserved
in this folder -- those live on Matsumoto's side.

See W:\\30_CauerLadderNetwork\\2023_08_08_松本@法政\\補足スライド.pdf
for the actual content (this module deliberately does not duplicate
slide bullet points -- the PDF is short enough to read directly).

## Why it is a separate topic

Even though there is no new code, this folder defines the *contract*
between the lab and Hosei for the second academic year of the
collaboration: which meshes are canonical, what node count is the
upper limit Hosei could run, and what calculation procedure they
were expected to follow.  Future "Hosei CLN" inquiries should be
checked against this folder first.

## Key files (Matsumoto)

| File                                     | Purpose                              |
|------------------------------------------|--------------------------------------|
| `補足スライド.pdf`                       | Calculation procedure handout        |
| `Mesh_Node243717_Element485552.txt`      | Medium-density baseline mesh         |
| `Mesh_Node370882_Element_739882.txt`     | High-density baseline mesh           |
"""

RECTANGULAR_CLN_TBD = r"""
# 2026_04_01 -- Rectangular-cross-section CLN (TBD, placeholder)

## Status

**Empty folder.**  Created 2026-04-01 as a placeholder; no files yet.

## What we expect

The historical CLN method derivations (both Kameari 2018 base paper
and the Sugahara-lab follow-ups) start from one of two canonical
1D cross-sections:

* **Slab / plate** (thickness d, infinite width) -- gives the
  `(4n-5)*4 / sigma d` / `d mu / (4n-3)` CauerI formula above.
* **Round wire / circular cylinder** (radius a) -- gives the
  Bessel J0/J1 closed form.

Real induction-heating and motor-winding cross-sections are
*rectangular* (litz strands, busbars, rectangular slot conductors).
The CLN for a rectangular cross-section does not separate into a
single transverse coordinate -- it is genuinely 2D in (x, y) within
the conductor and requires either:

(a) A 2D FEM-based CLN recursion (so the basis fields E^{(n)}(x, y),
    H^{(n)}(x, y) live on a 2D mesh, similar to the Endo problem
    but with the conductor as the *primary* CLN domain rather than
    as a passive observer).

(b) A double Fourier-series expansion (closed-form analytic) that
    generalizes the slab formula to rectangular geometry -- analogous
    to how the Wakao-Igarashi-Fujiwara-Kameari Part-1 series treats
    the 2D rectangular bar (radia_mcp.radia_ngsolve.analytical_formulas
    topic `rectangular_bar_2d`).

(c) A *separable approximation* using the slab formula in each
    transverse direction, treating the cross-section as a product
    of two slabs -- valid only at the "thin both ways" limit.

The 2026 plan (per lab roadmap, not yet committed in this folder) is
to attempt option (a) first -- because the FEM-based CLN
infrastructure already exists in cln_systematic and only the
conductor 2D domain needs to be plugged in.

## Why this matters

If option (a) lands successfully, the lab gains *the* CLN that works
for litz strand bundles and motor slot conductors directly, without
the strand-by-strand FEM brute force.  This would close out the
"CLN cannot handle rectangular conductors" criticism that has come
up repeatedly from industrial collaborators since 2018.

## What to do when this folder gains contents

1. Update RECTANGULAR_CLN section title from "TBD" to the actual
   scope.
2. Add the cross-section formula (if closed-form) or the FEM-recursion
   pseudocode (if option (a) won).
3. Add a citation to the JIEEJ / IEEE Trans Mag paper that documents
   the result -- by 2026 this should be in submission.
4. Cross-reference into `cln_knowledge.CLN_BASIC_RECURSION` so the
   reader of the core method knows that the rectangular case follows
   the same recursion with a different basis.
"""

INVERTER_APPLICATION = r"""
# 2017_05_17_インバータ回路.docx -- CLN as inverter subcircuit

## What the docx covers

A first-cut design study from 2017 that pairs the CLN with an
inverter circuit:

* CLN models the eddy-current load (a magnetic actuator or the
  primary side of an induction heater) as a small CauerII ladder.
* The inverter is a standard H-bridge (4 IGBTs with antiparallel
  diodes) driven by a PWM waveform.
* The CLN ladder is exported as a SPICE subcircuit (each L_k and
  R_k becomes a named SPICE element) and dropped into the H-bridge
  netlist.
* The combined simulation runs in LTspice (or any SPICE-compatible
  simulator) at the system level -- one transient run gives
  switch losses + conduction losses + eddy-current losses in the
  load, in milliseconds.

Without the CLN this is intractable: the eddy-current load alone is
a 1-7 hour transient FEM solve per switching cycle, and the H-bridge
needs to be resolved at the IGBT switching frequency (10-100 kHz).

## Why the docx is short

It is a *design memo*, not a paper -- it lays out:

1. The intended netlist topology (which CauerII node connects to
   which inverter midpoint).
2. The recipe for emitting the CLN as SPICE: one `L` line per
   inductor, one `R` line per resistor, named consistently with the
   stage index.
3. A worked example with N = 5 stages.

The numbers / waveforms shown are illustrative, not validated -- the
validated work is in the 07_電気力学 folder of the systematic catalog
(`cln_systematic.electrodynamic_coupling`).

## Pairing with `07_電気力学`

The 2017 docx is the *prototype* for the entire 07_電気力学 stream:

* 2017 docx: H-bridge + CLN, design memo, no measurements.
* 07_電気力学 systematic folder: H-bridge + CLN + bench measurements
  on real induction-heater hardware, with the CLN expanded to handle
  multi-stage non-uniform input.

Anyone trying to understand "why is CLN in `cln_systematic` paired
with an inverter at all" should read the 2017 docx first to see the
motivation, then jump to the 07_電気力学 docs for the full story.

## Lessons from the inverter coupling

1. **SPICE name collisions.**  If you emit `L1`, `R1`, ..., `L5`,
   `R5`, you will clash with the inverter sub-circuit's own L1 / R1
   (gate resistor + parasitic inductor).  Always prefix the CLN
   names: e.g. `LCLN1`, `RCLN1`, ... so the netlist is unambiguous.

2. **Numbering follows the form.**  CauerI: R is even (R0, R2, ...)
   and L is odd (L1, L3, ...).  CauerII: L is even (L0, L2, ...) and
   R is odd (R1, R3, ...).  Pick one and stick with it through the
   netlist or the inverter waveform plots will be impossible to
   sanity-check by hand.

3. **Initial conditions matter.**  At t=0 the inverter is energizing;
   if the CLN ladder is initialized to zero current on every inductor,
   the first half-cycle wastes itself charging the ladder.  For
   periodic-steady-state runs, *precharge* the CauerII inductors with
   the DC-bias current first (`uic` flag in SPICE), then start the
   PWM waveform.

4. **Stage count needs only to cover the PWM bandwidth.**  Switching
   frequency f_sw = 20 kHz means the CLN only needs to represent Z(s)
   accurately up to ~200 kHz (10th harmonic of f_sw).  N = 5
   typically suffices for the slab / wire formulas; N = 7 for safety.
   Do not over-stage -- more stages slow down SPICE for no accuracy
   gain.
"""

LAB_LESSONS = r"""
# Group-E lessons (form conversion, collaborators, applications)

## L1.  Pick the form by where you need accuracy

CauerI is "tight at DC", CauerII is "tight at infinity".  Truncating
to N stages, the *moments* of Z(s) match exactly around the expansion
point (s=0 for CauerI, s=infty for CauerII) and lose accuracy as you
move away.  Choosing the wrong form at small N gives a 30 percent
error in the wrong frequency band that looks like a bug.

## L2.  Lanczos works only at small N -- below ~7 stages

The non-orthogonality (in the standard inner product) of the Lanczos
basis lets numerical noise grow stage by stage.  The lab measured
machine precision at N <= 7, ~1e-6 drift at N = 10, divergence at
N = 15.  Above N ~ 7 either reorthogonalize, or derive the target
form *directly* from its own FEM recursion (e.g.
`Analytical_CLN_for_H_mode_for_iron_sheet.m` for CauerII from H-mode).

## L3.  Q-route is always more stable than the alpha / beta route

`Z_recovered(s) = 1 / [Q (alpha + s beta)^{-1} V](1)` is more stable
than `Z_recovered(s) = 1 / I_first_node(s)` computed from the bare
alpha / beta values.  The projection Q^T (...) Q absorbs most of the
accumulated Lanczos drift.  When in doubt, store Q and recover Z via
the projection, not via the reduced ladder.

## L4.  External collaborators get *one* deliverable: precomputed .mat

In the Endo collaboration the lab shipped a single `.mat` file
(`CLN_VectorPotential.mat`) containing (Jz, Az, Rn, Ln) for the
agreed mesh.  Endo's job was projection + f-domain solve, *not*
re-running the CLN recursion.  This makes the boundary contract
clean: "if the .mat doesn't reproduce your A_0, tell us; we re-export
.mat" rather than "Endo accidentally re-meshed and broke the basis".

## L5.  Mesh-lock is the silent killer in shared CLN work

The CLN basis fields are defined on a specific mesh.  Any re-mesh on
the collaborator side -- adaptive refinement, edge of the air shell
moved, even a different node numbering after a fresh COMSOL import --
invalidates the basis.  Verify the mesh first (node count, element
count, file hash if possible) before debugging the impedance answer.

## L6.  Source reconstruction first, f-domain solve second

`plot_contour_A0_source_reconstruction.m` shows whether the chosen
nStage can even *represent* the source A_0.  If the reconstructed A_0
visibly differs from the analytical A_0, no f-domain answer at any
frequency will be correct.  Always run this sanity check first.

## L7.  SPICE coupling: name everything `LCLN<k>`, `RCLN<k>`

Prefix all CLN ladder elements with `CLN` in any SPICE export.  The
inverter, the load resistor, the snubbers, and stray-parasitic models
in the surrounding netlist all want to be `L1`, `R1`, etc., and you
will spend an afternoon debugging a name collision otherwise.

## L8.  Precharge CauerII before periodic-steady-state runs

The CauerII inductor at the input carries the DC bias current.
Starting from zero, the inverter wastes one half-cycle charging it
up before steady state begins.  Set `.IC I(LCLN0)=I_dc` and use
the `uic` simulator flag.  This collapses transient-to-steady-state
time from 50+ cycles to ~5.
"""

KEY_FILES = r"""
# Key files referenced in this module

## Cauer I <-> Cauer II (2021_01_22)

| Path                                              | Role                                |
|---------------------------------------------------|-------------------------------------|
| `MATLAB/lanczos.m`                                | Two-matrix Lanczos (Q-route)        |
| `MATLAB/CauerI_to_CauerII.m`                      | Main driver, plate + wire examples  |
| `MATLAB/CauerI_to_CauerII_analytic.m`             | Symbolic round-trip verification    |
| `MATLAB/calc_cln_i_params.m`                      | Plate / wire CauerI parameter helper|
| `MATLAB/Analytical_CLN_for_H_mode_for_iron_sheet.m`| Direct CauerII from H-mode (no Lanczos) |
| `MATLAB/check_CLN_with_analytical.m`              | Z(s) vs `s mu tan(kd/2)/(kd/2) d`   |
| `verify_cln_i_ii_equivalence.py`                  | Python equivalence test, N=3..15    |
| `verify_lanczos_unified.py`                       | U=V=Q verification                  |
| `draw_cln_circuits.py`                            | SchemDraw -> cln_type_{i,ii}.pdf    |
| `cln_type_i.pdf`, `cln_type_ii.pdf`               | Reference circuit diagrams          |
| `docs/2021_01_22_CauerIからCauerIIの変形.docx`     | Original derivation memo            |
| `docs/2021_04_23_CauerIからCauerIIの変形.docx`     | Updated derivation (post analytic)  |
| `熱伝導/E0_mode_analytical_Hiruma.m`               | Thermal CauerI on disk (Hiruma)     |

## Endo @ Hosei (2022_10_09)

| Path                                       | Role                              |
|--------------------------------------------|-----------------------------------|
| `CLN_VectorPotential_solver.m`             | Main MATLAB driver                |
| `CLN_VectorPotential_comsol.m`             | COMSOL geometry + physics build   |
| `VectorPotential.m`                        | Closed-form A_z source            |
| `COMSOL_f_domain.m`                        | COMSOL reference solve            |
| `plot_contour_A_CLN.m`                     | CLN A_z contour                   |
| `plot_contour_A_COMSOL_f_domain.m`         | COMSOL A_z contour                |
| `plot_contour_A0_source_reconstruction.m`  | Source-projection sanity check    |
| `plot_contour_A0_VectorPotential_check.m`  | A_0 vs analytical sanity check    |
| `4段_CLN_I.pdf`, `4段_CLN_II.pdf`           | 4-stage I/II form circuit diagrams|
| `回路図(1..6).pdf`                          | Stage-1..6 circuit diagrams       |
| `CLN法の計算手順_0728.pptx`                 | Procedure handout (early)         |
| `2022_11_03_CLN.pptx`, `2022_11_16_CLN.pptx`| Procedure handout (later versions)|
| `CLN/CLN_VectorPotential.mat`              | Precomputed CLN basis (the .mat)  |
| `CLN/offset_*/`                            | CLN reconstruction results        |
| `COMSOL/offset_*/`                         | COMSOL reference results          |
| `V_offset.png`                             | Headline result figure            |

## Matsumoto @ Hosei (2023_08_08)

| Path                                     | Role                          |
|------------------------------------------|-------------------------------|
| `補足スライド.pdf`                       | Procedure handout (Matsumoto) |
| `Mesh_Node243717_Element485552.txt`      | Medium baseline mesh          |
| `Mesh_Node370882_Element_739882.txt`     | High-density baseline mesh    |

## Inverter application (2017)

| Path                                                       | Role                |
|------------------------------------------------------------|---------------------|
| `W:\\30_CauerLadderNetwork\\2017_05_17_インバータ回路.docx` | Design memo (CLN + H-bridge)|

## Rectangular CLN (2026)

| Path                                              | Role                       |
|---------------------------------------------------|----------------------------|
| `W:\\30_CauerLadderNetwork\\2026_04_01_長方形CLN\\`| Empty placeholder, see TBD |
"""

CITATIONS = r"""
# Citations

## CauerI and CauerII canonical forms

* W. Cauer, "Synthesis of linear communication networks", McGraw-Hill,
  1958 (English translation of Cauer's 1930s German work).  The
  continued-fraction expansions around s=0 (CauerI) and s=infty
  (CauerII) of a rational impedance Z(s) are the historical
  definitions used here.

* The two ladder topologies are also covered in any standard text on
  network synthesis (Van Valkenburg, "Introduction to Modern Network
  Synthesis", 1960; Weinberg, "Network Analysis and Synthesis", 1962).

## Lanczos tridiagonalization route to CauerII

* C. Lanczos, "An iteration method for the solution of the eigenvalue
  problem of linear differential and integral operators", J. Res.
  Natl. Bur. Stand. 45 (1950) 255-282.

* In the lab's two-matrix non-standard inner-product form, the
  procedure follows Saad's treatment (Y. Saad, "Iterative methods
  for sparse linear systems", 2nd ed., SIAM 2003, Ch. 6.6).  The
  finite-precision pathology of Lanczos (loss of orthogonality
  beyond N ~ 7-10 stages without reorthogonalization) is a textbook
  result; see Parlett, "The Symmetric Eigenvalue Problem", SIAM 1998,
  Ch. 13.

## CLN base paper

* A. Kameari, H. Ebrahimi, K. Sugahara, Y. Shindo, T. Matsuo,
  "Cauer Ladder Network Representation of Eddy-Current Fields for
  Model Order Reduction Using Finite-Element Method", IEEE Trans.
  Magn. 54(3): 7201804, 2018.  This is the foundational CLN paper.

* K. Sugahara, A. Kameari, H. Ebrahimi, Y. Shindo, T. Matsuo,
  "Finite-Element Analysis of Unbounded Eddy-Current Problems Using
  Cauer Ladder Network Method", IEEE Trans. Magn. 54(3): 7200704,
  2018.  Companion paper handling open-boundary problems.

* T. Matsuo, A. Kameari, K. Sugahara, Y. Shindo, "Matrix Formulation
  of the Cauer Ladder Network Method for Efficient Eddy-Current
  Analysis", IEEE Trans. Magn. 54(11): 7205805, 2018.  The matrix
  form that the Lanczos route in this module operates on.

## Closed-form impedance references (CauerI parameter derivations)

* For the plate / slab `(4n-5)*4 / sigma d`, `d mu / (4n-3)` formulas:
  this falls out of the standard skin-effect series expansion of
  `Z = (s mu / k) tan(kd/2)` around s=0, k = sqrt(s mu sigma).  See
  any classical EM textbook (e.g. Stratton, "Electromagnetic Theory",
  McGraw-Hill 1941, sec. 5.18).

* For the round-wire `(2n-1) / (sigma pi a^2)`, `mu / (8 n pi)`:
  same approach starting from `Z = (k / (2 pi a sigma)) J0(ka)/J1(ka)`.

## Endo / Matsumoto @ Hosei University

* If a journal paper was produced from these collaborations, it
  is not yet identified in the lab archive.  TBD -- to be added
  when the citation surfaces.  The folders themselves are the
  primary record.

## Inverter application

* The 2017 inverter docx is an internal design memo; it has not
  been published.  The validated bench measurements live in the
  systematic catalog folder `07_電気力学` (covered separately in
  `cln_systematic_knowledge.py`).

## Thermal CLN

* Hiruma method (lab-internal name) -- the algebraic CLN-style
  reduction of the heat-diffusion equation; documented in the
  systematic-catalog folder `12_比留間法` and partially mirrored
  in `2021_01_22_CauerI_to_CauerII\\熱伝導\\`.
"""

_TOPICS = {
    "overview": OVERVIEW,
    "cauer_i_vs_ii": CAUER_I_VS_II,
    "cauer_form_conversion_code": CAUER_FORM_CONVERSION_CODE,
    "endo_hosei": ENDO_HOSEI,
    "matsumoto_hosei": MATSUMOTO_HOSEI,
    "rectangular_cln_tbd": RECTANGULAR_CLN_TBD,
    "inverter_application": INVERTER_APPLICATION,
    "lab_lessons": LAB_LESSONS,
    "key_files": KEY_FILES,
    "citations": CITATIONS,
    "all": "",
}
_TOPICS["all"] = "\n\n".join(v for k, v in _TOPICS.items() if k != "all")


def get_cln_collab_documentation(topic: str = "all") -> str:
    """Return the CLN collaborator / form-conversion / application doc for ``topic``.

    Available topics:
        overview, cauer_i_vs_ii, cauer_form_conversion_code, endo_hosei,
        matsumoto_hosei, rectangular_cln_tbd, inverter_application,
        lab_lessons, key_files, citations, all.
    """
    if topic in _TOPICS:
        return _TOPICS[topic]
    return ("Unknown topic '" + topic + "'. Available: "
            + ", ".join(sorted(_TOPICS.keys())))
