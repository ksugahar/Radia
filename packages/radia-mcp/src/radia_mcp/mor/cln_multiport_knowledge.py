"""CLN multi-port / multi-expansion / 3D extensions.

Covers public-safe curated corpus, 04_展開点, 10_タイル,
11_3次元}. Source-side and structural extensions of the basic
single-port linear CLN: vector port impedance Z_ij, multi-expansion-
point convergence acceleration (Kuriyama 2019), periodic/tiled
domains, and 3D extension.

This module is Group B of the Sugahara Lab CLN absorption.  See
`cln_knowledge.py` for the canonical single-port background (Kameari-
Ebrahimi-Sugahara-Shindo-Matsuo 2018 base paper, basic recursion,
Cauer ladder topology, single-expansion convergence).  The four
topics here are the lab's specific extensions that turn the textbook
single-port CLN into a tool that handles real engineering systems
(transformers, multi-conductor coils, Litz wires, WPT, 3D iron).
"""

OVERVIEW = """
# CLN extensions: multi-port, multi-expansion, tiled, 3D

The vanilla CLN of Kameari et al. 2018 builds a SCALAR ladder that
matches a single-port impedance Z(jw) of an eddy-current FEM problem.
That is enough for a one-winding inductor on a textbook geometry.
Real engineering systems require four extensions that this module
covers:

  1. MULTI-PORT (public-safe curated corpus)
     -- N-conductor windings, transformers, multi-source coils
     -- Output: a Y_ij(jw) (or Z_ij(jw)) MATRIX, not a scalar
     -- Two formulations exist: "voltage-controlled-current-source"
        coupled CLNs (Sugahara 2017) and "matrix-form Kameari"
        (Matsuo 2017 / 2018c); the matrix form converges at lower
        ladder order so it is the production choice.

  2. MULTI-EXPANSION POINT (public-safe curated corpus)
     -- Add s_0 != 0 (typically s_0 = j*omega_0 of the operating
        frequency) so the same ladder works DC -> high-f with few
        stages.
     -- Replace the static operator C^T nu C by K = C^T nu C + s_0
        sigma (or K = C rho C^T + s_0 mu for T-formulation).  Same
        recursion, K^{-1} replaces (C^T nu C)^{-1}.
     -- Variants: A-formulation (.docx 有限要素法展開点付.doc),
        3D variant (有限要素法展開点付3d.doc), T-formulation
        (有限要素法展開点付t.doc), A-K formulation
        (有限要素法展開点付ak.doc).
     -- Published as Kuriyama-Kameari-Ebrahimi-Fujiwara-Sugahara-
        Shindo-Matsuo 2019 (cited in cln_knowledge.py).

  3. TILED / PERIODIC (public-safe curated corpus)
     -- A graphical representation of ladder networks where each
        circuit element becomes a rectangular tile; total impedance
        = aspect ratio of the assembled tile pattern.
     -- A pedagogical / debugging tool, NOT a separate solver --
        helps the lab visualize convergence and termination at a
        glance, especially for multi-port (multi-color-tile) cases.

  4. 3D EXTENSION (public-safe curated corpus)
     -- Vector A_phi formulation with HCurl/Nedelec edge elements
        + H1 nodal Lagrange multiplier (gauge); FreeFEM++ for mesh
        prototyping (test.edp), Gridap.jl for production
        (A_phi_Gridap.jl), 3D mesh experiments via FreeFEM msh3 +
        SurfaceHex+Sphere builders (2013_03_03 folder).
     -- The dominant blocker historically has been HCurl gauge
        treatment (kernel of curl, ungauged null space) and 3D mesh
        density vs. 2D-axisym alternatives; lab pragmatic choice is
        often 2D axisym (calc_inductance.py style) when the device
        admits it, and 3D CLN is held in reserve.

Together these four turn CLN from a 2D textbook demo into a tool
that produces SPICE-compatible N-port equivalent circuits for the
lab's induction-heating coils, multi-turn air-core inductors, Litz
wires, transformers, and WPT pads -- all from one FEM solve.
"""

MULTIPORT_THEORY = """
# Multi-port CLN: theory (folder 03_複数電源)

## From scalar Z to matrix Y_ij

The single-port CLN of Kameari 2018 produces one ladder generating
the scalar impedance Z(jw) seen at one terminal pair.  For an
N-port device (transformer, multi-coil WPT, multi-turn winding
with several externally addressable taps) the analog is a Y matrix:

    [ I_1 ]   [ Y_11 Y_12 ... Y_1N ] [ V_1 ]
    [ I_2 ] = [ Y_21 Y_22 ... Y_2N ] [ V_2 ]
    [ ... ]   [   ...              ] [ ... ]
    [ I_N ]   [ Y_N1 ...      Y_NN ] [ V_N ]

The lab has TWO derivations, both in the 03_複数電源 folder.

## Path A: coupled CLNs with voltage-controlled current sources
(Sugahara 2017, 2017_02_27_複数電源のKameari法.docx)

Idea: build ONE single-port CLN per port (Y_ii diagonal) using the
ordinary Kameari recursion -- short all other ports, apply 1 V at
port i.  The off-diagonal Y_ij comes from a side-channel reading:
during the generation of CLN_i, compute the mode current G_2n that
each mode injects into port j (the shorted port).  G_2n is one of
the unknown currents already produced by the recursion, so this
adds no extra solve.

The coupled circuit (Fig. 4 of the docx) is two single-port CLNs
laid in parallel, joined by voltage-controlled current sources
G_2n_1 (sensing R_2n_2 terminal voltage, injecting into port 1)
and G_2n_2 (sensing R_2n_1, injecting into port 2).

Y-parameter formulae (from the derivation):

    Y_11 = (admittance of CLN generated by 1 V at port 1,
            port 2 shorted)               -- diagonal, ordinary CLN
    Y_22 = (admittance of CLN generated by 1 V at port 2,
            port 1 shorted)               -- diagonal, ordinary CLN
    Y_12 = sum_n G_2n_2 / Z_2n_(stage 1)
    Y_21 = sum_n G_2n_1 / Z_2n_(stage 2)

The Joule loss inside the resistors of the two-port CLN equals the
Joule loss of the original FEM problem -- field reconstruction is
only needed if the spatial distribution of loss is required.

Once generated, the two-port CLN drops into LTSpice and any other
external circuit (load impedance, switching converter, etc.) can be
attached without regenerating the CLN.  This is the WHY of the
two-port extension: a single-port CLN must be re-built whenever the
load changes; the two-port CLN is load-independent.

## Path B: matrix-form Kameari (Matsuo 2017, 2018c)
(2017_12_30_行列版Kameari法.docx)

For low-frequency power-electronics problems the two ports often
have ISOLATED grounds; the natural network description is the
INDUCTANCE MATRIX, not Z (which forces a common ground for a T
equivalent).  Sugahara's evaluation (folder 行列版) concluded that
the matrix-form Kameari "produces accurate results with fewer
ladder stages than the controlled-source two-port CLN" -- so it is
the production choice in the lab.

The matrix form generalizes equations (11)-(13) of the A-form
Kameari recursion (single port):

    C^T nu C (a_{n+1} - a_n) = (1/lambda_{2n}) sigma e_{2n}    (11)
    C (e_{2n+2} - e_{2n})    = -(1/lambda_{2n+1}) C a_{n+1}    (12)
    lambda_{2n}   = e_{2n}^T sigma e_{2n}                       (13a)
    lambda_{2n+1} = a_{n+1}^T C^T nu C a_{n+1}                  (13b)

into a MATRIX recursion where a_n and e_n carry one column per port
(so a_n is N_dof x N_port, not N_dof x 1), and lambda becomes an
N_port x N_port symmetric positive matrix.  See section MULTIPORT_CODE
for the actual MATLAB / FreeFEM++ code.

## Practical lessons (from 03 corpus)

- For N ports the matrix lambda has N(N+1)/2 independent entries
  per stage.  The 1-source two-port (folder 2_port_1-wire) shows
  that 5-stage matrix CLN matches FEM from DC to 1 MHz.

- Termination choice (R_terminate vs L_terminate, n even vs n odd)
  matters at higher frequency; the lab convention is L-termination
  with odd n=N (=2k+1) for best mid-band fit, R-termination with
  even n=N for best DC fit.

- The grounding topology (T equivalent vs PI equivalent vs raw
  matrix) is dictated by the application -- see the matrix-form
  docx for the rule of thumb: PI equivalent works when port
  grounds are isolated (Y matrix), T equivalent forces common
  ground (Z matrix), inductance matrix works either way.
"""

MULTIPORT_CODE = """
# Multi-port CLN: code (folder 03_複数電源)

## Reference geometry

Two-port, parallel round-trip conductors (Fig. 1 of the docx):
- conductor radius a = 0.5 mm
- conductor spacing d = 0.1 mm
- conductivity sigma = 5.76e7 S/m (copper)
- region: x in [0, 15 mm], y in [-8 mm, 8 mm] for the small case,
  x in [0, 50 mm], y in [-25 mm, 25 mm] for the multi-turn case

The mesh is built in FreeFEM++ (check.edp, Multi-turnLadderSeries.edp);
see KEY_FILES below.

## FreeFEM++ recursion driver (Multi-turnLadderSeries.edp)

The axisymmetric A-formulation Kameari recursion with multi-turn
series connection.  The key recursion body:

    // First stage, n = 0:
    for(int i=0; i<nxy; ++i){
        real j0n = int2d(Th, conRegion[i])(e0*sigma/x);
        vn[i] = 1./j0n;        // per-turn unit voltage
        V += vn[i];
    }
    for(int i=0; i<nxy; ++i) vn[i] /= V;   // normalize series
    EE[0] = sum_i ((Th(x,y).region==conRegion[i])) * e0 * vn[i];
    lambda[0] = int2d(Th)(sigma * EE[0]*EE[0] / x) * 2*pi;
    R = 1/lambda[0];

    // n-th stage, n = 1..nmax:
    Vdc J = sigma * EE[n-1] / lambda[n-1];
    Vh   A, w;
    solve Prob1(rA, w) = int2d(Th)( (1/mu) * (dx(w)*dx(rA) +
                                              dy(w)*dy(rA)) / x )
                       - int2d(Th)(w * J / x)
                       + on(OUTER, rA = 0);
    A = rA;
    if(n > 1) A = AA[n-2] + A;            // accumulate
    lambda[n] = int2d(Th)((1/mu) * (dx(A)^2 + dy(A)^2) / x) * 2.*pi;
    AA[n] = A;

    ++n;                                   // odd-step electric
    EE[n] = EE[n-2];
    for(int i=0; i<nxy; ++i){
        EE[n] -= AA[n-1]/lambda[n-1] * (Th(x,y).region==conRegion[i]);
        real h = int2d(Th, conRegion[i])(sigma*EE[n]*EE[0]/x) / norm[i];
        EE[n] -= h * EE[0] * (Th(x,y).region==conRegion[i]);
    }
    lambda[n] = int2d(Th)(sigma * EE[n]*EE[n] / x) * 2*pi;

Notes:
- The `2*pi` factor and division by x (= r in cylindrical coords)
  are the axisymmetric volume integration measure dV = 2*pi*r*dr*dz.
- The Gram-Schmidt step `EE[n] -= h * EE[0] * (region indicator)`
  preserves per-conductor orthogonality of the electric modes.
- The lambda alternation (even = R-stage, odd = L-stage) is the
  same as the single-port recursion.

## Impedance reconstruction (continued fraction)

After lambda[0..N] are computed, the impedance at any frequency is
a closed-form continued fraction (no extra FEM solve):

    real omega = 2*pi*freq;
    complex omegaj = 1i*omega;
    complex Z;
    int n = nmax;
    Z = 1/lambda[n];                       // termination
    while (n > 0) {
        Z = 1/Z;  --n;
        Z = Z + 1/(omegaj * lambda[n]);    // L stage
        --n;
        Z = 1/Z;
        Z = Z + 1/lambda[n];               // R stage
    }
    // Z is the port impedance at frequency `freq`

## MATLAB post-processing: Y-parameter assembly (transformer/CLN_scalar.m)

For each pair of ports (n_port, m_port) and each frequency:

    s = j*2*pi*F(nf);
    Z = diag(R);                                    % R from FEM lambda
    for n = 1:N-1
        Z(n:n+1, n:n+1) = Z(n:n+1, n:n+1) + ...
            s * [[L(n), -L(n)]; [-L(n), L(n)]];     % L stage
    end
    Z(N, N) = Z(N, N) + s*L(N);                     % termination
    V = [1; zeros(N-1, 1)];                         % unit drive at port 1
    II = Z \\ V;                                     % LU solve
    Y(nf) = (Imn .* R) * II;                        % off-diagonal Y_mn

The cross-coupling current Imn comes from the side-channel readings
during CLN generation (see MULTIPORT_THEORY above): file
`m_port/Im.dat` written by the FreeFEM++ pass holds `J[n] = ...`.

## LTSpice integration

Each two-port CLN can be exported as a SPICE .cir netlist
(transformer/CLN_scalar.m, cases 1-7):

    R{port}_{n} n{port}_{n} n{port}_{n+2} {R(stage)}
    L{port}_{n} n{port}_{n+1} g{port} {L(stage)} Rser=0
    G{port}_{mport}_{n} n{port}_0 g{port} n{mport}_{n} n{mport}_{n+2}
                                              {coupling_gain}

with each port group sharing its own ground node g{port}.  LTSpice
.ac sweep then reproduces the FEM impedance across the band.
"""

MULTI_EXPANSION_THEORY = """
# Multi-expansion-point CLN: theory (folder 04_展開点)

## Why one expansion point hits a wall

Kameari 2018 expands around s = 0 (DC, the static magnetic
operator C^T nu C).  Each stage is a STATIC FEM solve.  The series

    e_2n+2 = e_2n - (1/lambda_{2n+1}) a_{n+1}                  (Kameari eq. 15)
    a_{n+1} = a_n + (1/lambda_{2n}) (C^T nu C)^{-1} sigma e_2n

converges geometrically with radius proportional to |s/s_0|^{-1}.
For s_0 = 0 the "distance" to a target operating frequency
omega_op is |omega_op|, which can be very large -- 85 kHz WPT,
500 kHz IH, MHz-band converters -- and the series converges slowly
or not at all in the practical N <= 60 stages.

## The trick (Kuriyama 2019)

Introduce a NON-ZERO expansion point s_0 = j*omega_0.  Replace the
static operator C^T nu C by the COMPLEX-SHIFTED operator

    K = C^T nu C + s_0 sigma                                    (eq. 17)

and keep the recursion structure intact (eq. 14, 15 of
有限要素法展開点付.doc):

    K (a_{n+1} - a_n) = (1/lambda_{2n}) sigma e_2n              (14)
    e_{2n+2} - e_{2n} = -(1/lambda_{2n+1}) a_{n+1}              (15)
    lambda_{2n}   = e_{2n}^T sigma e_{2n}                       (16a)
    lambda_{2n+1} = a_{n+1}^T K a_{n+1}                         (16b)

K is invertible at any s_0 not equal to the system's resonance, so
K^{-1} replaces (C^T nu C)^{-1} -- exactly the same one-time LU
factorization cost, with PARDISO etc.  Boundary conditions and
forcing terms are absorbed into a_0 / e_0 (the static seed) and the
recursion ignores them after stage 0.

The price: a_n and e_n become COMPLEX vectors (since K is complex),
lambda is complex, and the resulting ladder has complex L and R
elements.  The Sugahara Lab fix (展開点付き20170919.docm, section 3)
is to IGNORE the imaginary residual at stage 0 and re-introduce it
through a single "drift correction" at the end -- the bulk of the
recursion stays real and the ladder is interpretable as ordinary
L/R/G elements.

## What changes vs. single-expansion

| Aspect                  | s_0 = 0 (Kameari)   | s_0 != 0 (Kuriyama)    |
|-------------------------|---------------------|------------------------|
| LHS operator            | C^T nu C            | C^T nu C + s_0 sigma   |
| Cost per stage          | one static FEM      | one COMPLEX FEM (same  |
|                         |                     | sparsity)              |
| LU re-use across stages | yes (one factor)    | yes (one COMPLEX       |
|                         |                     | factor)                |
| Convergence radius      | |omega_op|          | |omega_op - omega_0|   |
| Stages to <1% at 85 kHz | >100                | ~6 (Kuriyama 2019      |
|                         |                     | result, WPT case)      |

## Multi-stage expansion (switch s_0 mid-way)

The most powerful variant: use s_0 = 0 for the first M stages
(captures DC limit exactly), then SWITCH to s_0 = j*omega_op for
stages M+1..N (captures operating frequency exactly).  The
orthogonality auxiliary `p_{2i-1}` terms (Kuriyama eq. 11-12)
and one extra inductor L* appear at the switch point, but no
re-factorization of the FIRST M stages is needed.

The 2017_07_02_展開点付.docx and 展開点付き20170919.docm record
the lab's experimental sweep: 0-stage vs 2-stage-from-switch vs
3-stage-from-switch.  At a target omega_0 in the kHz band the
2-stage-from-switch case typically wins.

## T-formulation variant (有限要素法展開点付t.doc)

The T-formulation (current vector potential T, scalar magnetic
potential omega) admits the same expansion-point trick with

    K_T = C rho C^T + s_0 mu                                    (eq. 17, T-form)

where rho = 1/sigma is the resistivity matrix and mu the
permeability matrix.  Boundary handling is cleaner (no Coulomb
gauge needed) but the formulation is restricted to conductor-only
domain, so air must be handled externally.

## 3D variant (有限要素法展開点付3d.doc)

The 3D case has the same K = C^T nu C + s_0 sigma form with
Whitney 1-form (edge) elements and Whitney 2-form (face) elements.
The added complication is the kernel of curl (Coulomb gauge or
tree-cotree decomposition) -- the lab uses a Lagrange multiplier
in H1 (see 11_3次元 / A_phi_Gridap.jl, THREE_D_EXTENSION below).

## A-K formulation (有限要素法展開点付ak.doc)

An "A and K" hybrid -- specifics in `public-safe curated corpus
04_展開点\\有限要素法展開点付ak.doc`.  Lab notes mark this as
research-stage; production uses A-form 3D.
"""

MULTI_EXPANSION_CODE = """
# Multi-expansion-point CLN: code (folder 04_展開点)

## Reference geometry (One Conductor_CLN_s0.edp)

A simplest possible test: ONE circular conductor of radius a inside
a circular air bubble of radius b = 1.1*a.  Skin depth at f_0 = 100 kHz
for copper (sigma = 5.76e7, mu_0) is delta ~ 0.21 mm -- the conductor
radius 0.5 mm is ~ 2.4 delta, deep into eddy-current regime.

    real f0 = 1.e5;       // expansion-point frequency (Hz)
    real s0 = 2.*pi*f0;   // expansion-point angular frequency
    int N  = 8;           // ladder stages

The expansion-point operator is built once per stage:

    Vh Az, w;
    solve Prob(Az, w) =
        int2d(Th)( (1/mu)*(dx(Az)*dx(w) + dy(Az)*dy(w))
                 + s0 * sig * Az * w )                          // K = ...
      - int2d(Th)(w * J)                                        // forcing
      + on(Outer, Az = 0);                                      // BC

This single FreeFEM++ `solve` block is the heart of the
expansion-point recursion -- it replaces the static `solve
... int2d(Th)((1/mu)*grad(Az).grad(w))` of the vanilla CLN.
For n = 1 (the bootstrap stage) the same operator is used with the
boundary condition `on(Outer, Az = -V)` (V = 1 unit voltage drive)
instead of an interior J forcing:

    if (n == 1) {
        solve Prob(Az, w) =
            int2d(Th)( (1/mu)*(dx(Az)*dx(w)+dy(Az)*dy(w)) + s0*sig*Az*w )
            + on(Outer, Az = -V);                               // V seed
        E = Az;
        lambda[n] = 1. / int2d(Th)( (1./mu)*(dy(E)^2 + dx(E)^2)
                                  + s0*sig*E*E );
        AA[n] = lambda[n] * E;
    }

## Impedance reconstruction with non-zero s_0

The continued fraction has a residual-correction term `s/(s - s_0)`
on every R stage (because the residual was non-zero at s_0):

    complex s = 2i*pi*freq;
    complex Z = 0.;
    int n = N;
    if (N % 2 == 1) { Z += s*Ls(n); --n; }   // odd termination = L
    while (n >= 2) {
        Z += s/(s - s0) * Rs(n);             // R stage with correction
        --n;
        Z = 1./Z;
        Z += 1./(s*Ls(n));                   // L stage
        Z = 1./Z;
        --n;
    }
    // Z is the port impedance at frequency `freq`

The `s/(s - s_0)` factor reduces to 1 when s_0 = 0 (recovering the
vanilla CLN reconstruction) and to 0 when s -> s_0 (the residual
vanishes at the expansion point, as designed).

## Triple-recurrence verification (有限要素法展開点付.doc, eq. 21, 30)

The lab routinely sanity-checks orthogonality by re-deriving the
three-term recurrence from the recursion equations:

    lambda_{2n+1} e_{2n+2}
    + [(1/lambda_{2n}) K^{-1} sigma - (lambda_{2n+1} + lambda_{2n-1}) I] e_{2n}
    + lambda_{2n-1} e_{2n-2}                                  = 0    (eq. 21)

    lambda_{2n+2} a_{n+2}
    + [(1/lambda_{2n+1}) K^{-1} sigma - (lambda_{2n+2} + lambda_{2n}) I] a_{n+1}
    + lambda_{2n} a_n                                          = 0    (eq. 30)

Compute and print

    e_{2k}^T sigma e_{2n+2}    for all k <= n         (should be 0)
    a_{k+1}^T K a_{n+2}        for all k <= n         (should be 0)

at the end of the recursion -- non-zero values (above floating
roundoff ~ 1e-12) flag a bug in the K^{-1} sparse-solve precision,
a Dirichlet BC leak, or a numerical instability in the
Gram-Schmidt step.
"""

PERIODIC_TILED = """
# Tiled / periodic CLN representation (folder 10_タイル)

## What the "tile" representation is

タイル張り.docx introduces a graphical convention for ladder
networks: each circuit element becomes a rectangle whose

  - horizontal extent represents current,
  - vertical extent represents voltage,
  - aspect ratio = impedance,
  - area = power.

A pure RESISTOR is a fixed-aspect-ratio tile (Ohm's law independent
of frequency).  A pure INDUCTOR is a tile whose horizontal extent
scales with s (so the aspect ratio scales with s = j*omega) -- at
high omega it becomes tall and thin (looks like a high-Z element),
at low omega it becomes wide and short (low-Z, "shorted out").
SERIES connection of two elements = stack the tiles vertically.
PARALLEL connection = lay them side by side.  Total impedance =
aspect ratio of the assembled rectangle.

For a ladder

    R_0 - L_0 - R_1 - L_1 - ... - termination (R_n or L_n)

the tile diagram is a sequence of tall thin L tiles inter-leaved
with square R tiles, with termination as a final tile (R-square or
L-rectangle depending on parity).

## Why the lab uses it

The tile picture makes three things instantly visible:

  1. WHERE the ladder converges:  far past the point where each
     L_k tile becomes much wider than its corresponding R_k tile
     (or vice versa), the remaining stages contribute negligible
     area -- this is the visual analog of the
     |Z_n - Z_{n-1}|/|Z_n| < tol termination criterion.

  2. TERMINATION choice (R vs L):  a final R-tile produces a
     finite DC limit, a final L-tile produces an open-circuit DC
     limit.  The right choice depends on whether the physical
     device is DC-conductive or DC-isolated (transformer secondary).

  3. MODE INSPECTION:  in the lab's multi-port plots
     (mcp-server-document figures), each ladder stage is drawn as
     a tile of a different color -- it is immediately obvious which
     mode dominates the response at which frequency.

## How it appears in practice

The 10_タイル folder ships ~98 .emf rendered tile diagrams at
f = 100 kHz and f = 1 MHz (49 stages each) -- these are stage-by-
stage visualizations of the multi-turn air-core coil from
03_複数電源 (マルチターン空芯コイルのＣＬＮ解析.docx, see Fig. 5
"タイル図" referenced in that document).  They are pedagogical /
debugging artifacts, not solver inputs.

## Relation to Bloch / Floquet periodic CLN

NOTE: The folder name "タイル" (tile) is the GRAPHICAL tile, NOT
a Bloch-Floquet periodic-structure CLN.  A true periodic-cell CLN
(unit cell + Bloch phase) is in scope for `Homogenization-by-Biot
.docx` in folder 03 (see HOMOGENIZATION_BIOT_SAVART below), not
in the タイル folder.  The lab knowledge base treats these as
distinct concepts.

## Homogenization by Biot-Savart CLN (03/Homogenization-by-Biot.docx)

A spiritually-similar but technically different idea: a unit cell
of a multi-conductor structure (e.g. a Litz wire bundle) is solved
ONCE with a Biot-Savart-derived boundary A condition, replacing
the full air-region mesh:

  - B_1 = 1 T uniform inside the cell (no mesh needed for L_1).
  - L_1 = (1/mu_0) * V_cell.
  - E_2 = electric mode satisfying sigma E_2 = curl H boundary,
    solved on conductors only.
  - B_3 (next magnetic mode) requires A_3 on a surface surrounding
    the conductors; A_3 from Biot-Savart over the known J_2 = sigma E_2:

        A_3(r) = (mu_0 / 4*pi) integral over conductors of
                 J_2(r') / |r - r'| dV'                          (eq. 7)

    The k-periodic repetition of the cell is handled by SUMMING
    Biot-Savart contributions from periodic images at
    d_k = (kx*dx, ky*dy, kz*dz) with k in a finite Brillouin set.

This is the lab's path to "homogenized periodic CLN" -- a cousin
to but not the same as the 10_タイル graphical tiles.
"""

THREE_D_EXTENSION = """
# 3D CLN extension (folder 11_3次元)

## State of the art in the lab

Folder 11_3次元 contains three artifacts that together sketch the
3D-CLN frontier in the Sugahara Lab:

  1. test.edp (FreeFEM++) -- the 3D mesh-build smoke test.  Builds
     an 8x8x8 layered hex mesh on the unit cube + plain 3D Stokes
     solve.  Used to confirm FreeFEM++ msh3 / buildlayers /
     medit pipeline works on the lab Windows machines.  NOT a CLN
     run; it is the "is the 3D mesh tool chain alive?" test.

  2. 2013_03_03_3次元メッシュ/Laplace3d.edp + MeshSurface.idp -- a
     3D mesh-construction library (SurfaceHex + Sphere) with
     metric-driven sphere meshing via adaptmesh.  Builds a
     boundary mesh of a cube-with-ball assembly for downstream
     tetrahedral mesh by `tetg`:

         mesh3 ThH = SurfaceHex(N, B, L, 1);
         mesh3 ThS = Sphere(0.5, h, 7, 1);
         mesh3 ThHS = ThH + ThS;
         savemesh(ThHS, "Hex-Sphere.mesh");
         // then tetg(ThHS, switch="pqaAAYYQ", ...) for volume mesh

     Still a mesh / topology pipeline test, no CLN content.

  3. 20231120_A_phi_Lukas定式/A_phi_Gridap.jl -- the LIVE 3D-CLN
     implementation in Gridap.jl (Julia FEM library), A-phi
     formulation.

The Gridap.jl script is the production-grade 3D-CLN; the FreeFEM++
files are legacy / scratch.

## A-phi formulation (Gridap.jl, A_phi_Gridap.jl)

The continuous problem:

    curl(nu * curl A) + j*omega*sigma*(A + grad phi) = J   in Omega
    div(sigma*(A + grad phi)) = 0                          in conductor
    n x A = 0                                              on S_boundary

The Lagrange multiplier phi (scalar potential in H1) imposes the
Coulomb gauge `div sigma (A + grad phi) = 0` weakly via the
saddle-point form.  Gridap encodes this as a MultiFieldFESpace:

    Vcurl = TestFESpace(model, ReferenceFE(nedelec, Float64, order);
                        conformity=:Hcurl, dirichlet_tags="S_boundary")
    VH1   = TestFESpace(model, ReferenceFE(lagrangian, Float64, order);
                        conformity=:H1, dirichlet_tags="S_boundary")
    Y = MultiFieldFESpace([Vcurl, VH1])
    X = MultiFieldFESpace([Ugcurl, UgH1])

The CLN recursion is the standard A-form one (Kameari eq. 11-13),
but assembled in 3D.  Each stage is a saddle-point solve:

    a((a, phi), (a', phi')) =
        int( (1/mu)*curl(a) . curl(a') ) dOmega
      + int( grad(phi) . a' ) dOmega
      + int( a . grad(phi') ) dOmega
    b((a', phi')) = int( a' . J ) dOmega

    op = AffineFEOperator(a, b, Y, X)
    xh = solve(BackslashSolver(), op)
    uh, ph = xh        # uh = A (HCurl), ph = phi (H1)

After stage `nStage`:

    Ao = uh
    if nStage == 1:    A = R * Ao
    else:              A = A + R * Ao
    L = sum( int(R*J*A) dOmega )
    J = J - sigma * A / L                  # next-stage electric mode

Output is per-stage VTK files for ParaView, exporting

    Ao, A, B = curl(A), phi, J

with cellfields tagged for inspection of each ladder mode.

## What works, what doesn't (lab notes)

WORKS:
- Order-1 Nedelec + Order-1 H1 multiplier on tetrahedral meshes
  (Gmsh v2 or v4).  Geometry: cylinder_Gmsh_2D5_v2.msh (a 2.5D
  cylinder, code is true 3D but mesh is extruded).
- Backslash direct solver for small N_dof (<~50k).  Memory blows
  up around 100-200k DOF on a desktop machine -- larger problems
  need an iterative solver (BiCGStab + AMS preconditioner from
  `radia.sparsesolv_ngsolve`, see CLAUDE.md "Compact HX for HCurl
  Problems").

DOES NOT WORK / known issues:
- The "S_boundary" Dirichlet tag is hard-coded in the .jl -- a
  different mesh with a different boundary name needs the script
  edited.  No production-grade tag-discovery yet.
- Coreform Cubit v2 .msh files (commented-out `mesh/input_v2.msh`,
  `mesh/input_v4.msh`) are still in flux; the lab uses Gmsh-
  produced meshes as the reference.
- Multi-expansion-point + 3D is NOT yet combined in production --
  the A-form 3D recursion uses s_0 = 0.  Promoting Kuriyama 2019's
  K = C^T nu C + s_0 sigma to the saddle-point 3D form is open
  research (rough sketch in 有限要素法展開点付3d.doc).

## Why 3D is held in reserve

For lab-typical devices (axisymmetric coils, IH workpieces, WPT
pads) the 2D-axisymmetric path (calc_inductance.py + the
multi-port 2D Kameari) is FASTER, BETTER VALIDATED, and produces
the same Z_ij(jw) to engineering accuracy.  3D is reserved for
non-axisymmetric geometries (transformers with asymmetric
windings, motor rotors, fully 3D iron yokes).  The Gridap.jl
prototype is the lab's option for when 3D becomes necessary --
not the default.
"""

LAB_LESSONS = """
# Hard-won lab lessons (multi-port / multi-expansion specific)

These are pitfalls the lab has hit and fixed -- not in any paper,
just in folder docx notes and ML discussions.

## L1. Always verify orthogonality at the end (every multi-port run)

After the recursion, compute the matrix

    M_ik = e_{2i}^T sigma e_{2k}    for all i <= k

and print it.  Off-diagonal entries above 1e-10 indicate that the
Gram-Schmidt step lost orthogonality -- usually a sign that K^{-1}
was solved with insufficient precision (PARDISO single-precision,
or BiCGStab tol > 1e-12).  See `Multi-turnLadderSeries.edp` lines
142-153 for the check.  Single-precision K^{-1} is enough for
N <= 5 stages; for N > 20 the lab uses double-precision PARDISO.

## L2. Port-current normalization matters for multi-conductor series

In a multi-turn coil where all turns are connected in series
externally, `vn[i]` is the per-turn voltage that ensures TOTAL
current = 1 A in the series circuit.  The lab convention
(Multi-turnLadderSeries.edp lines 87-93):

    for(int i=0; i<nxy; ++i){
        real j0n = int2d(Th, conRegion[i])(e0*sigma/x);
        vn[i] = 1./j0n;     // per-turn unit-current voltage
        V += vn[i];          // sum of per-turn voltages = port V
    }
    for(int i=0; i<nxy; ++i) vn[i] /= V;   // normalize series

Forgetting the second loop produces an N-conductor PARALLEL
two-port instead of an N-conductor SERIES winding -- impedance
off by ~N^2.  Bug seen and fixed in the 2017_02_21 -> 2017_02_27
revision of the docx (folder 03/2_port_1-wire/).

## L3. Matrix-form lambda is symmetric but NOT necessarily positive

In the matrix form (folder 行列版), the per-stage lambda is an
N_port x N_port symmetric matrix.  It is GUARANTEED symmetric
(by the integration-by-parts identity) but can have negative
eigenvalues if the chosen ports do not span an admissible
N-port (e.g. two ports that physically short each other through
a perfectly-coupled magnetic core).  The lab fix: check
`eig(lambda)` at every stage; if any eigenvalue < 0, the port
definition is degenerate and must be re-specified.

## L4. Termination parity is application-dependent (NOT one-size-fits-all)

The lab's rule of thumb (2017_12_30_行列版Kameari法.docx):

    odd N (terminate on L)  -> good DC limit for INDUCTORS
    even N (terminate on R) -> good DC limit for RESISTIVE devices

A transformer secondary terminated on R will show a non-zero DC
admittance even though the physical secondary is DC-open --
choose L termination.  A multi-turn resistive coil terminated on
L will show an inductive DC component that is not in the FEM
data -- choose R termination.  Get this wrong and the multi-port
Y matrix has the right shape but the wrong dc asymptote.

## L5. PI vs T equivalent has GROUNDING implications

PI equivalent (Y matrix) works for isolated grounds (transformer,
optocoupler, isolation amplifier).  T equivalent (Z matrix)
FORCES a common ground.  The matrix-form Kameari paper
(2017_12_30) specifically calls this out: if your two ports
have isolated grounds, do NOT use a T equivalent -- use the
inductance matrix directly (or convert to PI).

## L6. K = C^T nu C + s_0 sigma is COMPLEX -- watch arithmetic precision

The expansion-point operator K is a COMPLEX symmetric (NOT
Hermitian) sparse matrix.  PARDISO's complex-symmetric LU
factorization is the natural fit (mtype = 6 in MKL PARDISO).
COCR (Conjugate Orthogonal CR) is the lab default iterative
backend for K -- see `radia.sparsesolv_ngsolve.COCRSolver`.
Casual use of CG (which assumes Hermitian SPD) gives wrong
answers on K -- check the matrix type before choosing a solver.

## L7. Stage-0 residual is NOT zero at non-zero s_0 -- account for it

In single-expansion (s_0 = 0), the static seed a_0 satisfies
C^T nu C a_0 = j_0 exactly, so the stage-0 residual is exactly
zero.  At s_0 != 0, the seed satisfies K a_0 = j_0 BUT
K = C^T nu C + s_0 sigma INCLUDES the conductive term -- so the
residual at frequency s != s_0 is sigma * (s - s_0) * a_0 which
is NOT zero.  This residual is what the `s/(s - s_0)` factor in
the impedance reconstruction (MULTI_EXPANSION_CODE) is correcting
for.  Forgetting it produces a constant offset on every Y_ij at
high frequency.

## L8. Tile diagrams are pedagogy, not solver state

The 10_タイル diagrams (49 stages, 100 kHz and 1 MHz) are
SAVED OUTPUT of a post-processing run -- they are not consumed
by any solver.  Their value is for thesis writeups and reviewer
explanations.  Do not try to "import" them into the recursion;
their information is already implicit in the lambda array.

## L9. 3D + multi-expansion + multi-port is unsolved (research-grade)

The triple combination (3D HCurl + s_0 != 0 + N port matrix) is
NOT yet validated in the lab.  The Gridap.jl prototype handles
3D + 1 port + s_0 = 0; matrix-form 3D is sketched in
有限要素法展開点付3d.doc as future work.  For now, use the 2D
axisym matrix CLN for production multi-port + multi-expansion;
hold 3D in reserve for asymmetric geometries.

## L10. PEEC + CLN is a natural hybrid (not yet productized)

PEEC (radia.lanczos_reduction.PRIMA) and CLN are both
Krylov-based MOR methods.  PEEC excels at filament-coil + air
(no iron); CLN excels at FEM-mesh + iron + eddy current.  A
hybrid -- PEEC for the coil L/R, CLN for the workpiece impedance
delta-L -- is exactly the calc_inductance.py / FEM-SIBC split
documented in CLAUDE.md.  Treat them as complementary, not
competing.
"""

KEY_FILES = """
# Key files cited in this module

The 03_複数電源 folder is the biggest (17 files + 4 subfolders);
04_展開点 is concise (10 files); 10_タイル has 1 docx + ~100 EMF
renderings; 11_3次元 has 6 files in 3 subfolders.

Top files ranked by educational value:

| Path                                                          | Type | Why it matters                                                  |
|---------------------------------------------------------------|------|-----------------------------------------------------------------|
| 03_複数電源/Multi-turnLadderSeries.edp                        | .edp | Live FreeFEM++ axisym multi-turn CLN driver (179 lines)         |
| 03_複数電源/2_port_1-wire/2017_02_27_複数電源のKameari法.docx | .docx| THE multi-port derivation, 2-port voltage-controlled CS form    |
| 03_複数電源/2_port_1-wire/2017_12_30_行列版Kameari法.docx     | .docx| Matrix-form Kameari (the production choice)                     |
| 03_複数電源/2-Port-CLNReport-_derivation.docx                 | .docx| English derivation, Y_ij from short/open tests                  |
| 03_複数電源/Homogenization-by-Biot.docx                       | .docx| Biot-Savart unit-cell CLN for Litz / multi-conductor periodic   |
| 03_複数電源/transformer/CLN_scalar.m                          | .m   | MATLAB SPICE .cir netlist generator for 2-port CLN              |
| 03_複数電源/transformer/check_source_Y.m                      | .m   | MATLAB Y_ij assembly from per-port FEM lambda (4-port)          |
| 03_複数電源/transformer_hassan/Compare_CLN_FEM.m              | .m   | Production CLN-vs-FEM comparison driver (439 lines)             |
| 03_複数電源/マルチターン空芯コイルのＣＬＮ解析.docx           | .docx| Multi-turn 100-conductor coil case study (with mode plots)      |
| 03_複数電源/多導体コイルのＣＬＮ.docx                         | .docx| Multi-conductor axisym coil + mode dimensions (R_n, L_n table)  |
| 04_展開点/One Conductor_CLN_s0.edp                            | .edp | THE expansion-point recursion in FreeFEM++ (1 conductor, K op)  |
| 04_展開点/有限要素法展開点付.doc                              | .doc | Full A-form expansion-point derivation (eq. 11-30)              |
| 04_展開点/有限要素法展開点付t.doc                             | .doc | T-formulation variant (current vector potential)                |
| 04_展開点/有限要素法展開点付3d.doc                            | .doc | 3D variant of expansion-point (research-stage)                  |
| 04_展開点/展開点付き20170919.docm                             | .docm| Comprehensive notebook -- Krylov view, Neumann series, equiv-ct |
| 10_タイル/タイル張り.docx                                     | .docx| Tile-graphical representation of ladder networks                |
| 11_3次元/20231120_A_phi_Lukas定式/A_phi_Gridap.jl             | .jl  | Live 3D CLN (Gridap.jl, A-phi formulation, Nedelec + H1)        |
| 11_3次元/2013_03_03_3次元メッシュ/Laplace3d.edp               | .edp | 3D mesh assembly (cube + sphere) via FreeFEM++ msh3             |
| 11_3次元/2013_03_03_3次元メッシュ/MeshSurface.idp             | .idp | SurfaceHex / Sphere mesh-builder library                        |

Subfolder structure of 03_複数電源:
- 2_port_1-wire/{1,2,LTSpice,行列版,直列_順方向,直列_逆方向,並列_順方向,電源1,電源2}
- 2_port_1-wire_open/{IABC,Kelvin_NG} -- open-port two-port + Kelvin BC research
- transformer/{1_single,2,2_single,3_single,4_single,電源1..4_single} -- 4-port T case
- transformer_hassan/ -- Hassan Ebrahimi's transformer demo + LTSpice .asc

Subfolder structure of 11_3次元:
- 2013_03_03_3次元メッシュ/ -- legacy FreeFEM++ 3D mesh smoke tests
- 20231120_A_phi_Lukas定式/ -- production Gridap.jl A-phi CLN
"""

CITATIONS = """
# Citations

Authoritative published references for the four extensions covered
here.  See `cln_knowledge.py` CLN_OVERVIEW for the canonical
foundational papers (Kameari 2018a/b, Matsuo 2018c, Ebrahimi 2020).

## Multi-port CLN

@article{Sugahara2017MultiPort,
  author  = {Sugahara, K. and Kameari, A.},
  title   = {Kameari-Hou wo Mochita Fukusuu Dengen no Kaiseki
             (Analysis of Multiple Power Sources Using the Kameari
             Method)},
  journal = {Lab note},
  year    = {2017},
  note    = {Folder: public-safe curated corpus
             2_port_1-wire\\2017_02_27_複数電源のKameari法.docx},
}

@article{Matsuo2017MatrixCLN,
  author  = {Matsuo, T. and Kameari, A. and Sugahara, K. and Shindo, Y.},
  title   = {Matrix Formulation of the Cauer Ladder Network Method
             for Efficient Eddy-Current Analysis},
  journal = {IEEE Trans. Magn.},
  volume  = {54},
  number  = {11},
  pages   = {7205805},
  year    = {2018},
  doi     = {10.1109/TMAG.2018.2840145},
}

## Multi-expansion point CLN

@article{Kuriyama2019MultiExpansion,
  author  = {Kuriyama, S. and Kameari, A. and Ebrahimi, H. and
             Fujiwara, K. and Sugahara, K. and Shindo, Y. and
             Matsuo, T.},
  title   = {Cauer Ladder Network With Multiple Expansion Points for
             Efficient Model Order Reduction of Eddy-Current Field},
  journal = {IEEE Trans. Magn.},
  volume  = {55},
  number  = {6},
  pages   = {7203404},
  year    = {2019},
  doi     = {10.1109/TMAG.2019.2900441},
}

@article{Sugahara2017ExpansionDerivation,
  author  = {Sugahara, K.},
  title   = {Yuugen Yousoshikuukan jou de no Kameari-Hou ni Tsuite
             (Tenkaiten-Tsuki) (On the Kameari Method in the Finite
             Element Space (with Expansion Points))},
  journal = {Lab note},
  year    = {2017},
  note    = {Folder: public-safe curated corpus
             有限要素法展開点付.doc and 3d/t/ak variants},
}

## Homogenization / periodic (relevant to TILED section)

@article{HirumaIgarashi2017Macroperm,
  author  = {Hiruma, S. and Igarashi, H.},
  title   = {Macro Toujiritsu wo Fukumu Kei no Toukai Kaeshi Goushin
             (Equivalent Circuit Synthesis of Systems with Macroscopic
             Permeability)},
  journal = {IEEJ Joint Technical Meeting on Static Apparatus and
             Rotating Machinery},
  number  = {SA-17-64, RM-17-95},
  year    = {2017},
  note    = {Cited in public-safe curated corpus
             マルチターン空芯コイルのＣＬＮ解析.docx as the source
             of the homogenization framing.},
}

## 3D CLN reference (NOT yet published, in-house)

@unpublished{Sugahara2023ApjGridap,
  author  = {Sugahara, K. and Lukas, [given name TBD]},
  title   = {A-phi Formulation 3D CLN in Gridap.jl},
  year    = {2023},
  note    = {Folder: public-safe curated corpus
             20231120_A_phi_Lukas定式\\A_phi_Gridap.jl.  Production
             3D-CLN script (Julia, Gridap.jl, Nedelec+H1
             saddle-point form).},
}

## Background sibling MOR -- PRIMA (Loop-Star PEEC)

@article{Odabasioglu1998PRIMA,
  author  = {Odabasioglu, A. and Celik, M. and Pileggi, L. T.},
  title   = {PRIMA: Passive Reduced-Order Interconnect
             Macromodeling Algorithm},
  journal = {IEEE Trans. CAD},
  volume  = {17},
  number  = {8},
  pages   = {645--654},
  year    = {1998},
  doi     = {10.1109/43.712097},
  note    = {Used in radia.lanczos_reduction.PRIMASchurExtractor
             for PEEC circuit MOR.  Cousin of CLN -- both use
             Krylov/Cauer expansion but PRIMA targets circuit-level
             PEEC matrices while CLN targets eddy-current FEM.},
}
"""


_TOPICS = {
    "overview": OVERVIEW,
    "multiport_theory": MULTIPORT_THEORY,
    "multiport_code": MULTIPORT_CODE,
    "multi_expansion_theory": MULTI_EXPANSION_THEORY,
    "multi_expansion_code": MULTI_EXPANSION_CODE,
    "periodic_tiled": PERIODIC_TILED,
    "three_d": THREE_D_EXTENSION,
    "lab_lessons": LAB_LESSONS,
    "key_files": KEY_FILES,
    "citations": CITATIONS,
    "all": "",
}
_TOPICS["all"] = "\n\n".join(v for k, v in _TOPICS.items() if k != "all")


def get_cln_multiport_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "overview"               - Big picture: four extensions and
                                 why they exist
      "multiport_theory"       - Folder 03: N-port Y_ij matrix,
                                 voltage-controlled CS form vs
                                 matrix-form Kameari
      "multiport_code"         - Folder 03: FreeFEM++ recursion,
                                 MATLAB SPICE export, impedance
                                 reconstruction
      "multi_expansion_theory" - Folder 04: K = C^T nu C + s_0 sigma,
                                 single vs multi expansion, A/T/3D
                                 variants
      "multi_expansion_code"   - Folder 04: One Conductor_CLN_s0.edp,
                                 impedance reconstruction with
                                 s/(s-s_0) correction
      "periodic_tiled"         - Folder 10: tile-graphical
                                 representation + Biot-Savart unit-
                                 cell homogenization
      "three_d"                - Folder 11: A-phi Gridap.jl
                                 production 3D CLN, FreeFEM++
                                 legacy mesh tools
      "lab_lessons"            - 10 hard-won pitfalls (orthogonality,
                                 normalization, K-complex, etc.)
      "key_files"              - File-by-file citation table
      "citations"              - BibTeX-ish reference list
    """
    t = topic.lower().strip()
    if t in _TOPICS:
        return _TOPICS[t]
    return ("Unknown topic '" + topic + "'. Available: "
            + ", ".join(sorted(_TOPICS.keys())))
