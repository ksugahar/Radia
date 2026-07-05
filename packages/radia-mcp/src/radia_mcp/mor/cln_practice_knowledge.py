"""CLN practice / foundations: 4 source folders + 3 root reference docs.

Sugahara lab's Cauer Ladder Network MATLAB+COMSOL implementation
practice corpus -- absorbed from public-safe curated corpus

Scope (Group A: Foundations + Practice):
  - 01_single_source   (public-safe curated corpus)
  - 02_boundary_cond   (public-safe curated corpus)
  - 09_series_expand   (public-safe curated corpus)
  - linear-CLN-practice (public-safe curated corpus)
  - A-phi-formulation-2.pdf (root)
  - 二次元準静磁場解析の再考.docx (root)
  - 2016_10_02_ベッセルの不等式.docx (root)

Companion to cln_knowledge.py (theoretical / publication-paper level).
This module is the IMPLEMENTATION-LEVEL knowledge: the actual MATLAB
class, the COMSOL `withsol` trick, the Robin BC implementation, the
Kelvin-vs-Infinite-Element-vs-Dirichlet trade-offs, the analytical
verification machinery.
"""

OVERVIEW = """
# CLN practice and foundations (Sugahara lab MATLAB+COMSOL corpus)

The lab's public-safe curated corpus tree contains roughly a decade of
Cauer Ladder Network implementation experiments.  This module absorbs
the foundational and 2020 single-frequency-input practice layer.

## What is in this module vs cln_knowledge.py

- cln_knowledge.py = the PAPER-LEVEL CLN: big-picture, canonical
  Kameari-Ebrahimi-Sugahara-Shindo-Matsuo 2018, recursion, multiple
  expansion points, applications.
- cln_practice_knowledge.py (this file) = the CODE-LEVEL CLN: the
  actual 71-line MATLAB CLN class, the COMSOL withsol() recurrence
  pattern, the analytical Legendre-form verification formulae up to
  order 9, the three competing open-boundary treatments (Dirichlet on
  truncated disk, Infinite Element layer, Kelvin transformation),
  and the lab-internal lessons learned from the practice scripts.

## Source folders absorbed

| Folder                                              | What it contains                                 |
|-----------------------------------------------------|--------------------------------------------------|
| 01_単一電源                                          | Single-source theoretical writeups + COMSOL .m  |
| 02_境界条件                                          | BC comparison: Kelvin vs IABC/IE vs Dirichlet   |
| 09_級数展開                                          | Series expansion theory (Kameari amendments)    |
| 2020_11_04_線形のCLNの練習                            | The "linear CLN practice" -- CLN.m + drivers    |
| A-phi-formulation-2.pdf                             | A-phi formulation reference                     |
| 二次元準静磁場解析の再考.docx                          | 2D quasi-static rethink                         |
| 2016_10_02_ベッセルの不等式.docx                       | Bessel inequality foundation                    |

## How to use

Topics (each is a self-contained markdown block):

  overview              - this section
  matlab_class          - CLN.m, the 71-line class explained
  e0_mode               - E0 (electric) mode practice scripts
  h1_mode               - H1 (magnetic source) mode + 3 BC variants
  boundary_conditions   - Kelvin vs Infinite vs Dirichlet trade-offs
  robin_bc              - Check_Robin_BC.m verification pattern
  series_expansion      - 09_級数展開 + Bessel inequality
  a_phi                 - A-phi formulation summary
  quasi_static_2d       - 2D MQS rethink doc
  lab_lessons           - hard-won pitfalls from lab practice
  key_files             - top 15 .m files with role + line count
  citations             - paper + docx references
"""


CLN_MATLAB_CLASS = """
# The CLN.m MATLAB class (71 lines, canonical lab implementation)

Location: public-safe curated corpus

This is the SINGLE canonical Sugahara-lab MATLAB class that consumes
the per-stage CLN basis fields {R_n, L_n, J_n, A_n} produced by a
finite-element step (COMSOL via LiveLink) and turns them into:

- f-domain impedance Z(jw) and per-frequency JJ, AA reconstruction
- t-domain transient response via ode23 integration of the ladder ODE

## Class definition (verbatim)

```matlab
classdef CLN
    properties
        nStage = [];
        R = [];        % per-stage resistances R_0..R_{N-1}
        L = [];        % per-stage inductances L_0..L_{N-1}
        J = [];        % per-stage J_n basis (cell array of vectors)
        A = [];        % per-stage A_n basis (cell array of vectors)
        Rmat = [];     % diag(R)
        Lmat = [];     % tridiagonal L block
        f = [];        % swept frequencies
        Z = [];        % Z(jw) at each f
        II = [];       % stage currents (per f or per t)
        JJ = {};       % reconstructed J(r) at each f or t step
        AA = {};       % reconstructed A(r) at each f or t step
        t = []         % t-domain time samples
    end
    methods
        function cln = CLN(R, L, J, A, nStage)
            cln.R = R(1:nStage);
            cln.L = L(1:nStage);
            cln.J = J(1:nStage);
            cln.A = A(1:nStage);
            cln.nStage = nStage;
        end
        function cln = calc_LR_mat(cln)
            cln.Rmat = diag([cln.R]);
            cln.Lmat = zeros(size(cln.Rmat));
            for n = [1:size(cln.Lmat,1)-1]
                cln.Lmat(n:n+1,n:n+1) = cln.Lmat(n:n+1,n:n+1) ...
                    + [[cln.L(n),-cln.L(n)];[-cln.L(n),cln.L(n)]];
            end
            cln.Lmat(end,end) = cln.Lmat(end,end) + cln.L(end);
        end
        function cln = f_domain(cln,f,V)
            cln.f = f;
            cln.Z = zeros(size(f));
            cln.JJ= {};
            cln.AA = {};
            for nf = [1:length(f)]
                cln.II  = [cln.Rmat + j*2*pi*f(nf)*cln.Lmat]\\[V];
                cln.Z(nf) = cln.II(1);
                cln.JJ{nf} = zeros(size(cln.J{1}));
                for n = [1:cln.nStage]
                    cln.JJ{nf} = cln.JJ{nf} + cln.II(n)*cln.R(n)*cln.J{n};
                end
                cln.AA{nf} = zeros(size(cln.A{1}));
                for n = [1:cln.nStage-1]
                    cln.AA{nf} = cln.AA{nf} + (cln.II(n)-cln.II(n+1))*cln.A{n};
                end
            end
        end
        function cln = t_domain(cln,tspan,V,I0)
            Rmat = cln.Rmat;
            Lmat = cln.Lmat;
            N = cln.nStage;
            dIdt = @(t,II) [Lmat]\\[[V(t);zeros(N-1,1)]-Rmat*II];
            [cln.t,cln.II] = ode23(dIdt,tspan,I0);
            for nt = [1:length(cln.t)]
                cln.JJ{nt} = zeros(size(cln.J{1}));
                cln.AA{nt} = zeros(size(cln.A{1}));
                for n = [1:cln.nStage]
                    cln.JJ{nt} = cln.JJ{nt} + cln.II(nt,n)*cln.R(n)*cln.J{n};
                end
                for n = [1:cln.nStage-1]
                    cln.AA{nt} = cln.AA{nt} + (cln.II(nt,n)-cln.II(nt,n+1))*cln.A{n};
                end
            end
        end
    end
end
```

## What each method does

### Constructor `CLN(R, L, J, A, nStage)`

Truncates input vectors/cells to the first `nStage` stages.  Inputs:
- `R` = vector of stage resistances [Ohm or per-unit per the model]
- `L` = vector of stage inductances [H or per-unit]
- `J` = cell array of per-stage basis current densities J_n(x)
- `A` = cell array of per-stage basis vector potentials A_n(x)

These come from the FEM-side recursion (calc_LR loop in
CLN_E0_mode.m / CLN_H1_mode.m -- see e0_mode / h1_mode topics).

### `calc_LR_mat(cln)` -- build Rmat and Lmat

The system matrix is built by hand from the ladder topology:

  Rmat = diag(R_0, R_1, ..., R_{N-1})

  Lmat = tridiagonal:
    [ L_0,        -L_0,                       ]
    [-L_0,  L_0+L_1,  -L_1,                   ]
    [        -L_1,   L_1+L_2, -L_2,           ]
    [                                ...      ]
    [                       ...,  L_{N-2}+L_{N-1} ]

The last-row diagonal entry has `+ cln.L(end)` from the trailing loop.
This matches the Cauer ladder MNA: each inductor L_n connects nodes n
and n+1, with a final L at the bottom branch.

### `f_domain(cln, f, V)` -- frequency-domain solve

For each frequency f_k in the sweep:
  II = (Rmat + j*2*pi*f_k*Lmat) \\ V
  Z(f_k) = II(1)   <- driving-point impedance at the first node
  JJ{k} = sum_n II_n * R_n * J_n(x)
  AA{k} = sum_n (II_n - II_{n+1}) * A_n(x)

The `JJ`/`AA` formulae are the standard CLN reconstruction:
- Current density J_total = sum of stage_current * stage_R * basis_J_n
  (the R weighting comes from the basis normalisation in the FEM step)
- Vector potential A_total = sum of (II_n - II_{n+1}) * basis_A_n,
  which is the "current through inductor L_n" times A_n -- physically
  this is the magnetic-flux contribution of stage n.

V is the input excitation vector, typically `V = [1; 0; ...; 0]` for
a unit-voltage drive at the first port (driving-point impedance run).

### `t_domain(cln, tspan, V, I0)` -- time-domain ODE solve

The ladder ODE is:
  L_mat * dI/dt + R_mat * I = V(t) e_1     (e_1 = [1;0;0;...])

Reformulated as:
  dIdt = @(t,II) [Lmat] \\ [[V(t); zeros(N-1,1)] - Rmat*II];

Then integrated with `ode23` (Bogacki-Shampine RK 2/3, default tspan,
default tolerances).  After the ODE solve, JJ and AA reconstructions
are done time-step by time-step, identical pattern to f_domain.

A commented stub `dIdt = @(t,II) [Lmat]\\[ここだけ変える]-Rmat*II];`
is a hand-note from Sugahara saying "if you want a different input
profile, this is the line to edit" -- the rest of the class is
input-agnostic.

## Why this class is the lab's canonical artefact

- 71 lines total, single file, no external deps beyond MATLAB base
- consumes only {R, L, J, A, nStage} -- the FEM is decoupled
- handles both f-domain (Z(jw) sweep) and t-domain (transient)
- can be re-invoked with different excitation/BCs by just changing V
  or I0 -- no FEM resolve

The downstream Check_CLN_class.m demonstrates usage against a
closed-form Bessel-function reference (see lab_lessons / key_files).
"""


E0_MODE_PRACTICE = """
# E0 mode practice: cylindrical conductor + COMSOL HelmholtzEquation hack

File: public-safe curated corpus
File: public-safe curated corpus

## The geometry and physics

Round copper cylinder of radius a=3 mm, sigma=57e6 S/m, mu_r=1, in 2D.
"E0 mode" = the m=0 (axisymmetric) excitation, where:
- J_z is the unknown current density (axial, only z component)
- A_z is the unknown vector potential (axial)
- The recursion alternates between J_z and A_z stages.

## The HelmholtzEquation hack (why COMSOL needs this)

COMSOL has no built-in "compute one stage of a Cauer recursion" primitive.
The lab uses the **HelmholtzEquation** physics with c=0, a=1, f=<custom>
which reduces the Helmholtz equation
  -div(c grad u) + a*u = f
to simply  u = f.  This lets us assemble an arbitrary source-side
expression as `f` and have COMSOL evaluate it on the FE mesh as a new
field -- effectively a projection operator.

Snippets from CLN_E0_mode.m:
```matlab
model.component('comp1').physics.create('hzeq1', 'HelmholtzEquation', 'geom1');
model.component('comp1').physics('hzeq1').field('dimensionless').field('Jz');
model.component('comp1').physics('hzeq1').feature('heq1').setIndex('c', [0 0 0 0], 0);
model.component('comp1').physics('hzeq1').feature('heq1').setIndex('a', 1, 0);
model.component('comp1').physics('hzeq1').feature('heq1').setIndex('f', 0, 0);
model.component('comp1').physics('hzeq1').prop('ShapeProperty').set('order', 7);

model.component('comp1').physics.create('hzeq2', 'HelmholtzEquation', 'geom1');
model.component('comp1').physics('hzeq2').field('dimensionless').field('Az');
% ... same c/a setup ...
```

So `hzeq1` is the J_z storage field and `hzeq2` is the A_z storage
field.  The Poisson physics `poeq1` (At = mu0*mur*Jz with Dirichlet
BC) gives the auxiliary "tilde A" needed for the orthogonalisation
step.

Note `order=7` shape functions -- the basis fields are smooth
polynomials (Legendre-like) inside the disk, so high-order shape
functions are necessary for accuracy at later stages.

## The `withsol` solution-chaining trick

The recursion needs each new stage to access the previous stage's
J_z, A_z, and At values.  COMSOL's `withsol('sol2', expr)` evaluates
`expr` using the stored solution `sol2` (the previous step) without
re-solving.  This is the key idiom for stitching stages together:

```matlab
% Stage n+1 J_z source = J_z_n - sigma * A_z_n / L_n
model.component('comp1').physics('hzeq1').feature('heq1').setIndex('f', ...
   'withsol(''sol2'',Jz) - sig*withsol(''sol2'',Az)/withsol(''sol2'',L)', 0);

% Stage n+1 A_z source = A_z_n
model.component('comp1').physics('hzeq2').feature('heq1').setIndex('f', ...
   'withsol(''sol2'',Az)', 0);
```

After each `model.sol('sol1').runAll;` the solution is stored back
to `sol2`, ready for the next iteration.  This is a single-solution
in-place chaining pattern.

(CLN_E0_mode_ver_2.m is a refinement that gives each stage its OWN
solution name `sol0, sol1, sol2, ...` so you can post-evaluate any
intermediate stage without re-running.  See ver_2 for the multi-sol
pattern.)

## R_n / L_n probe definitions (Domain Integral probes)

```matlab
model.component('comp1').probe.create('G', 'Domain');
model.component('comp1').probe('G').set('type', 'integral');
model.component('comp1').probe('G').set('expr', 'Jz^2/sig');
model.component('comp1').probe.create('L', 'Domain');
model.component('comp1').probe('L').set('type', 'integral');
model.component('comp1').probe('L').set('expr', 'Jz*Az/G');
```

So:
  G = integral( J_z^2 / sigma ) dV         (intermediate)
  R_n = 1 / G                              (per-stage resistance)
  L_n = integral( J_z * A_z / G ) dV       (per-stage inductance)

Pulled back into MATLAB via:
```matlab
Rn(n) = 1/mphglobal(model,'withsol(''sol2'',G)');
Ln(n) = mphglobal(model,'withsol(''sol2'',L)');
```

## Analytical Legendre-form verification (built into the script)

CLN_E0_mode.m hard-codes the closed-form expressions for A_n(r) and
J_n(r) for n = 1..9 in the round cylinder, derived by symbolic
integration in MATLAB Symbolic Toolbox.  Example for n=1, n=2:

  A_1(r) = mu0*mur*(a^2 - r^2) / (4*a^2*pi)
  A_2(r) = -mu0*mur*(a^4 + 3*r^4 - 4*a^2*r^2) / (8*a^4*pi)
  A_3(r) = mu0*mur*(a^6 - 10*r^6 + 18*a^2*r^4 - 9*a^4*r^2) / (12*a^6*pi)
  ...
  A_9(r) = mu0*mur*(a^18 - 24310*r^18 + 115830*a^2*r^16
                  - 231660*a^4*r^14 + 252252*a^6*r^12
                  - 162162*a^8*r^10 + 62370*a^10*r^8
                  - 13860*a^12*r^6 + 1620*a^14*r^4
                  - 81*a^16*r^2) / (36*a^18*pi)

(See the source file for the full nine formulae; they are the lab's
verification ground truth.)

Reference values for R_n and L_n in this case (n-indexing from 1):
  R_o(n) = (2n-1) / (a^2 * sigma * pi)        -- analytical R_n
  L_o(n) = mu0*mur / (8*pi) / n               -- analytical L_n

The driver prints the % error per stage:
```matlab
disp({sprintf('R%d',2*(n-1)),
      sprintf('%20.16e',Ro(n)),
      sprintf('%20.16e',Rn(n)),
      sprintf('%.4f%%',(Rn(n)-Ro(n))/Ro(n)*100)});
```

A passing run shows < 0.01% error up to n=5 at hmax=100e-6 mesh.

## Verifying the multi-sol pattern (CLN_E0_mode_ver_2.m)

The v2 driver creates 2N+1 separate study steps (sol0, sol1, ...,
sol{2N}), each running either the J_z stage (stat1, even n) or the A_z
stage (stat2, odd n).  This makes EACH stage post-evaluable, at the
cost of more disk for the stored solutions.  Probes are also created
per stage: G0, G2, G4, ... and L1, L3, L5, ...

The R/L collection becomes:
```matlab
case 0  % even n -> R stage
    Rn = [Rn, 1/mphglobal(model, sprintf('withsol(''sol%d'',G%d)',n,n))];
    Ro = [Ro, (n+1)/(a^2*sig*pi)];
case 1  % odd n -> L stage
    Ln = [Ln, mphglobal(model, sprintf('withsol(''sol%d'',L%d)',n,n))];
    Lo = [Lo, (mu0*mur)/(8*pi)/((n+1)/2)];
```

The commented blocks at the bottom of ver_2 contain hooks for the
Galerkin orthogonalisation correction terms (AzJz_{n,m}, AtAz_{n,m},
JzJz_n, AzAz_n) which are needed when basis orthogonality drifts at
higher stages.  See lab_lessons for when these matter.
"""


H1_MODE_PRACTICE = """
# H1 mode practice: magnetic source variant + 3 BC choices

Files in public-safe curated corpus
  - CLN_H1_mode.m                            (Robin BC, single disk)
  - CLN_H1_mode_Infinite_Element.m           (COMSOL Infinite Element)
  - CLN_H1_mode_Kelvin_NG.m                  (Kelvin transformation -- FAILED)
  - CLN_rectangle_H1_mode_Infinite_Element.m (rectangular conductor variant)

## What "H1 mode" means

E0 = m=0 axisymmetric mode (J on cylinder axis, A on cylinder axis).
H1 = m=1 mode: the conductor sees a transverse uniform-field-like
excitation.  Source term has an `(x/sqrt(x^2+y^2))` factor at stage 0,
so the basis is asymmetric (sin(theta)-like) rather than purely
axisymmetric.

The lab needed H1 to handle "the cylinder in an external transverse
field" geometry that approximates the standard SIBC / induction-
heating cell problem.

## The three BC strategies (and which one works)

| Strategy            | File                                           | Status |
|---------------------|------------------------------------------------|--------|
| Robin (-a/mur*n.grad on disk edge) | CLN_H1_mode.m                  | OK     |
| COMSOL Infinite Element layer      | CLN_H1_mode_Infinite_Element.m | OK     |
| Kelvin transformation (twin-disk)  | CLN_H1_mode_Kelvin_NG.m        | FAILED (NG) |

The filename `_Kelvin_NG.m` is itself the lab's pithy verdict: NG = "no
good".  See boundary_conditions topic for WHY.

## CLN_H1_mode.m: Robin BC pattern (working)

The trick is the Robin BC on the auxiliary At Poisson problem:
```matlab
model.component('comp1').physics('poeq1').create('dir1', 'DirichletBoundary', 1);
model.component('comp1').physics('poeq1').feature('dir1').selection.all;
model.component('comp1').physics('poeq1').feature('dir1').setIndex('r', ...
   '-a/mur*(nx*Atx+ny*Aty)', 0);
```

The Dirichlet "r" coefficient is set to `-a/mur * n . grad(At)` which
turns the apparently-Dirichlet condition into a Robin-style coupling
between At and its normal derivative at the disk boundary -- this
mimics what infinity would impose on a transverse-field problem
without needing an exterior domain.

(See Check_Robin_BC.m for a focused verification of this exact BC
trick on a single Poisson solve -- robin_bc topic.)

## Stage-0 source (the m=1 trigger)

```matlab
model.component('comp1').physics('hzeq1').feature('heq1').setIndex('f', ...
   'sig*2*mur/(mur+1)*(sqrt(x^2+y^2)/a)*(x/sqrt(x^2+y^2))', 0);
```

This is sigma * 2*mur/(mur+1) * (r/a) * cos(theta), the seed
m=1-symmetric source in cylindrical coordinates.

## Stage recursion (identical pattern to E0)

```matlab
model.component('comp1').physics('hzeq1').feature('heq1').setIndex('f', ...
   'withsol(''sol2'',Jz) - withsol(''sol2'',Az)/withsol(''sol2'',L)', 0);
% (note: no sigma factor here -- it's already absorbed into L definition)
```

The probe G uses `Jz^2/sig` (same as E0) but probe L uses
`Jz*Az/G/sig` -- the extra `/sig` is the m=1 normalisation.

## H1 mode analytical reference (Legendre-form, mu_r-dependent)

For Robin BC m=1 mode, A_n(r) has both r-powers AND mu_r-mixed terms:

  A_1(r) = sig*(mu0*r*(a^2*mur - mur*r^2 + 3*a^2 - r^2)) / (4*a^3*pi)
  A_2(r) = -sig*(mu0*r*(a^4*mur + 2*mur*r^4 + 11*a^4 + 14*r^4
                       - 27*a^2*r^2 - 3*a^2*mur*r^2)) / (4*a^5*pi)
  ... (up to n=9)

These hard-coded formulae let the script auto-verify each iteration.

Reference values:
  R_o(n) = n*(mur + 2n^2 - 1)^2 / (a^2 * mur^2 * sigma * pi)
  L_o(n) = mu0*(mur + (2n^2 - 1))*(mur + (2(n+1)^2 - 1)) / ((16n+8)*mur*pi)

Notice the explicit mu_r dependence: the H1 mode formulae are NOT
mu_r-independent (unlike E0).  This is one of the lab's surprises:
the radial structure of the basis fields depends on mu_r through the
Robin BC's `-a/mur` factor.

## CLN_H1_mode_Infinite_Element.m: alternative BC

Same H1 setup, but the BC is provided by COMSOL's built-in Infinite
Element coordinate scaling:
```matlab
model.component('comp1').coordSystem.create('ie1', 'InfiniteElement');
model.component('comp1').coordSystem('ie1').selection.set([1]);   % outer ring
model.component('comp1').coordSystem('ie1').set('ScalingType', 'Cylindrical');
```

A 2-disk geometry is used:
  - inner disk (radius a)  = conductor
  - middle ring (a to a_1) = real air
  - outer ring (a_1 to a_2) = Infinite Element compressed exterior

Shape function order is bumped to 3 for the IE layer (lower-order
gives meshing artefacts in the compressed exterior).

## CLN_rectangle_H1_mode_Infinite_Element.m

Same Infinite Element pattern, but the inner conductor is rectangular
(size sqrt(pi)*0.8*a x sqrt(pi)/0.8*a, equal-area to the disk).
Used as a sanity check: H1 basis on a rectangle should still produce
sensible R_n / L_n via FEM even if no closed form is available.

## Why CLN_H1_mode_Kelvin_NG.m failed

See boundary_conditions topic.  Briefly: the lab built a twin-disk
Kelvin transformation (real disk + image disk at y = 2.1*a, connected
by GeneralExtrusion coupling so the field is matched across them).
The Magnetic Potential coupling `mp1, mp2` with `genext1(Ax)` etc.
SHOULD impose the kelvin pair identification, but in the H1 recursion
the higher-stage A_n source terms broke the coupling consistency and
the iteration diverged.  The current best path for open H1 is Robin
or Infinite Element, NOT Kelvin -- at least with this implementation.
"""


BOUNDARY_CONDITIONS = """
# Boundary conditions for open-domain CLN: 3 strategies compared

Source: public-safe curated corpus
Files:
  - 2016_12_01_2導体モデルに亀有法を適用_iabc.docx  (IABC on 2-conductor)
  - 2017_05_21_境界条件.docx                       (BC overview)
  - AJ積分.docx                                    (A.J integral identity)
  - IABC境界条件に関して.docx                       (IABC BC theory)
  - IABC層内の磁気エネルギー.docx                    (IABC layer magnetic energy)
  - ac_iabc/                                       (AC IABC numerical test)
  - busbar/                                        (busbar full comparison)

The single biggest unknown when you start CLN in 2D is: how do you
truncate the infinite exterior?  Three competitors:

| BC strategy         | Idea                                                          | Pros                          | Cons (lab observed) |
|---------------------|---------------------------------------------------------------|-------------------------------|---------------------|
| Dirichlet on big disk | A = 0 on disk of radius R >> a                              | trivial to implement, COMSOL native | exterior must be VERY large for low-order modes; mesh cost; perturbs L_n by O((a/R)^2) |
| Robin (IABC) on disk  | -a/mur * n . grad(A) = A on conductor edge                  | exact for m=0,1 transverse field; no exterior needed | mu_r-dependent; coupling derivation per mode |
| COMSOL Infinite Element | coordinate stretch on outer annulus (Cylindrical Scaling) | conceptually "exact"; COMSOL native | shape order ~3 minimum; numerical conditioning poor at later stages |
| Kelvin (twin-disk)    | mirror the disk across plane y=offset, identify boundaries  | open-boundary representation natural | identification coupling drifts at higher stages -> divergence (CLN_H1_mode_Kelvin_NG.m) |

## What the lab found per BC

### Dirichlet on truncated big disk

Simplest, cheapest to set up.  Always converges -- but R/L values
drift away from the analytical reference for the same mode count.
Acceptable if the goal is f-domain response near the source band and
you can afford a 10-20a exterior.

### Robin BC (IABC-style, the lab's main result)

The "Impedance Absorbing Boundary Condition" coupling
  -a/mur * n . grad(A) = A
is exact for the m=1 mode on the disk boundary -- it represents the
1/r decay of the m=1 multipole at infinity.  Generalises to higher m
with appropriate -a/(m*mur) scaling but the lab's working scripts
only use m=0/1 so far.

The 02_境界条件\\IABC境界条件に関して.docx documents the derivation
in detail (binary docx, see file).  Key idea: the IABC is an
operator-based BC, not a point-wise one -- it relates A to its normal
derivative through a frequency- and mu_r-dependent coefficient.

For DC / quasi-static (omega -> 0) the IABC becomes the simple
"-a/mur * n.grad" Robin form that CLN_H1_mode.m uses.

### COMSOL Infinite Element

Uses COMSOL's built-in `coordSystem('ie1', 'InfiniteElement')` with
`ScalingType='Cylindrical'`.  Effectively an exact open-boundary
representation via coordinate stretching to infinity, but requires:
- a 3-region geometry (conductor / real-air / IE-air annulus)
- shape function order >= 3 (order 1-2 give mesh artifacts)
- careful selection of the IE outer-radius ratio (1.5x-2.0x of the
  inner air worked in CLN_H1_mode_Infinite_Element.m)

### Kelvin transformation (the "NG" path)

Geometry: conductor disk + mirror "Kelvin disk" placed at y = 2.1*a
distance.  GeneralExtrusion couplings (`genext1`, `genext2`) map field
values between the two disks' boundaries, imposing the Kelvin
condition that A(r) outside is the inverted A(R^2/r) inside.

In CLN_H1_mode_Kelvin_NG.m the lab built this with full COMSOL
InductionCurrents physics (`mf`) plus the HelmholtzEquation pair --
but the H1 recursion source terms, when projected through the
Kelvin coupling, introduced inconsistencies that grew per stage and
the iteration diverged.

The verdict (`NG.m` suffix): for the H1 m=1 mode with single source,
DO NOT use Kelvin coupling.  Use Robin or Infinite Element instead.

The lab DID get Kelvin to work for separate (single-physics) tests
in Kelvin_2D_check_Je_source.m and Kelvin_2D_check_fixed_bc.m --
these verify the Kelvin coupling itself is mechanically correct.
The breakage is specific to the CLN recursion's source-term
projection, NOT the Kelvin BC itself.

## The busbar case study (02_境界条件\\busbar\\)

The lab put a busbar (rectangular conductor) through all four BC
strategies and stored the resulting per-stage R/L vs. an FEMM
reference.  Output files:
  - IABC_2_dl=0.02_8_stage.txt_{L,R}.emf   (IABC at 8 stages)
  - Kelvin_dl=0.02_16_stage.txt_{L,R}.emf  (Kelvin needs 16 stages)
  - Robin_dl=0.02_8_stage.txt_{L,R}.emf    (Robin matches IABC at 8)
  - FEMM_dl=0.01/                          (FEMM reference)

Takeaway: IABC/Robin converge in 8 stages, Kelvin needs ~16 to reach
similar accuracy in this case, and Dirichlet (not shown in the
filenames here) needs even more.

## The AJ integral identity (AJ積分.docx)

A practical identity used in the lab's derivations:
  integral_Omega (A . J) dV = integral_dOmega (A x H) . n dS
                            + flux-linkage terms
which is the magnetic-energy / co-energy identity.  Used to verify
that the per-stage L_n = integral(J_n * A_n / G_n) dV matches the
magnetic-energy interpretation (1/2 L I^2 storage).

## Quick decision guide

For 2D round-conductor CLN with frequency sweep 0 - few MHz:
- Use Robin BC if you can derive the -a/(m*mur) coupling for your mode
- Use Infinite Element if you cannot, and you can afford order-3 FEM
- Use Dirichlet only as a last resort, with R/a >= 20
- Do NOT use Kelvin coupling in the CLN recursion (lab finding 2020)
"""


ROBIN_BC = """
# Check_Robin_BC.m: focused verification of the Robin BC

File: public-safe curated corpus

A standalone driver that ONLY verifies the COMSOL Robin BC setup
against an analytical Poisson solve.  Useful when chasing whether a
mismatched H1 mode result is due to the BC implementation or the
stage recursion.

## What it does

1. Build a 2D circle of radius a in COMSOL.
2. Set up `PoissonEquation` physics with field `At`.
3. Apply `DirichletBoundary` on ALL boundaries with `r = a` (i.e.
   the prescribed Dirichlet value IS `a` literally -- not zero!).
4. Source `f = 1` everywhere.
5. Solve.

Code:
```matlab
model.component('comp1').physics.create('poeq1', 'PoissonEquation', 'geom1');
model.component('comp1').physics('poeq1').field('dimensionless').field('At');
model.component('comp1').physics('poeq1').create('dir1', 'DirichletBoundary', 1);
model.component('comp1').physics('poeq1').feature('dir1').selection.all;
model.component('comp1').physics('poeq1').feature('dir1').setIndex('r', 'a', 0);
model.component('comp1').physics('poeq1').feature('peq1').set('f', '1');
```

This is the "Robin-equivalent Dirichlet" used in CLN_H1_mode.m,
isolated and tested against the analytical solution:

## Symbolic reference

```matlab
syms AA(r) a;
DA = diff(AA);
assume(a>0);
assume(r>0 & r<a);

cond1 = DA(0) == 0;        % regular at origin
cond2 = AA(a) == a;        % prescribed Dirichlet value
ode = -diff(r*diff(AA,r),r)/r == sym('1');
A = dsolve(ode,[cond1 cond2]);
```

Then `subs(A, {r,a}, {r_, a_})` is plotted against the COMSOL
`At` field -- a 1-line visual check that the BC behaves as the
symbolic computation predicts.

## Why this matters

The Robin BC in CLN_H1_mode.m uses
  setIndex('r', '-a/mur*(nx*Atx+ny*Aty)', 0)
which interprets the Dirichlet RHS string AS a Robin term (the `r`
field becomes a function of nx, ny, and the gradient components,
making it implicit in At).  COMSOL's Dirichlet implementation handles
this correctly only if the gradient discretisation matches the field
discretisation.

Check_Robin_BC.m is the focused test: same physics, same BC class,
constant `r = a` (no gradient coupling) -- if THIS gives the wrong
answer, your COMSOL install or physics setup is broken.  If it gives
the right answer, you can trust the Robin BC in the H1 driver.

## When to run

Run Check_Robin_BC.m FIRST when:
- a fresh COMSOL install
- after any LiveLink/COMSOL version upgrade
- when CLN_H1_mode.m results stop matching the analytical Legendre
  formulae and you cannot tell whether the bug is in the BC, the
  recursion, or the basis-evaluation chain
"""


SERIES_EXPANSION_THEORY = """
# Series expansion theory + Bessel's inequality (foundation)

Sources:
  - public-safe curated corpus
  - public-safe curated corpus
  - public-safe curated corpus

(All binary docx -- the descriptions below are the lab's standing
notes on what these documents establish.)

## Bessel's inequality role in CLN

The 2016_10_02 doc establishes that the CLN basis fields {J_n} and
{A_n} satisfy Bessel's inequality with respect to the
conductivity-weighted inner product:
  sum_n |<f, J_n>_sigma|^2 <= ||f||_sigma^2
where the weighted inner product is
  <f, g>_sigma = integral_Omega sigma * f * g dV

This is a fundamental sanity check on the orthogonalisation step:
Bessel's inequality MUST hold; if a numerical CLN run violates it
(typically due to mesh anisotropy or order mismatch), the basis has
ceased to be orthogonal and the Cauer circuit will not represent
Z(jw) correctly.

The practical check (not in the docs but standard in the lab):
- after N stages, compute Parseval residual = sum_n (I_n * R_n)^2
  vs. integral J(r)^2 / sigma dV at the chosen working frequency
- the residual should monotonically decrease to zero as N grows
- if it INCREASES at some stage, orthogonality has broken -- usually
  fixed by raising the FE shape function order or refining mesh

## Series convergence of the CLN basis

The 10-04_old and 亀有追記0921 documents prove (and amend) the
convergence theorem: under the assumption that the conductor cross-
section is Lipschitz and sigma is uniformly positive, the partial
sum  sum_{n=0}^{N-1} I_n * R_n * J_n  converges to the true J(r) in
L^2_sigma(Omega) as N -> infinity, with rate determined by the
spectral gap of the demagnetisation operator.

The lab's experience (see lab_lessons):
- For smooth m=0 mode on a disk, convergence is geometric in N
  (5 stages give 1e-6 in L_inf)
- For boundary-layer-heavy problems (skin depth << conductor),
  convergence is much slower; many more stages needed and/or
  multiple expansion points (cln_knowledge.py)

## Why this is in the foundations module

The Bessel-inequality and convergence theorems are the THEORETICAL
guarantee that CLN works.  Without them, CLN is just a heuristic
truncation.  The lab's docx files contain the explicit proofs that
make Sugahara's lab THE authority on CLN's mathematical foundations.

The Kameari amendment (亀有追記0921) tightens the convergence rate
proof for the specific case of axisymmetric conductors -- the
n-dependent rate matches the (mu_r + 2n^2 - 1) coefficient seen
empirically in the H1 mode analytical R_o formula.

## TBD content

Specific theorem statements / proofs not transcribed here (binary
docx).  See:
  public-safe curated corpus
  public-safe curated corpus
  public-safe curated corpus
"""


A_PHI_FORMULATION = """
# A-phi formulation reference

File: public-safe curated corpus

(Binary PDF, see file for full content.)

## What it covers

The A-phi formulation is the standard low-frequency electromagnetic
formulation:
  curl(nu curl A) + sigma (dA/dt + grad(phi)) = J_0
  div(sigma (dA/dt + grad(phi))) = 0

with A = magnetic vector potential, phi = electric scalar potential.
This is one of the two CLN seed formulations (the other is plain A
with Coulomb gauge).

## Why the lab keeps it as a top-level reference

CLN was originally derived in the A-only formulation (Kameari 2018
TMAG canonical paper), but multi-conductor and current-vs-voltage
driven setups benefit from the A-phi separation:
- A captures the magnetic energy / Cauer L's
- phi captures the lumped-port voltage drop -- making the connection
  to circuit-equation drives cleaner

The reference document derives:
- The space-time separation E(t,x) = sum v_n(t) E_n(x), H(t,x) =
  sum i_n(t) H_n(x) in the A-phi formulation
- The R_n / L_n integral definitions in terms of phi-gauge fields
- The recursion that produces J_n / A_n / phi_n simultaneously
- The matching to the Cauer ladder

## Connection to the practice scripts

CLN_E0_mode.m and CLN_H1_mode.m as shipped use plain A formulation
(no phi).  The A-phi extension is needed when:
- The driven port has a prescribed voltage (not current)
- Multiple conductors share a common ground reference
- DC current sharing among parallel paths must be tracked

The lab has draft scripts for A-phi CLN (not in 2020_11_04 folder;
in 01_単一電源\\1-wire\\ as text data, not as MATLAB code yet).
Promotion of A-phi to MATLAB-class form is a TBD upgrade path.

## TBD

Specific equations / derivations not transcribed.  See the PDF:
public-safe curated corpus
"""


QUASI_STATIC_2D_RETHINK = """
# 2D quasi-static rethink (二次元準静磁場解析の再考.docx)

File: public-safe curated corpus

(Binary docx, see file for full content.)

## What this document does

A retrospective look (~2020-2022, exact date TBD from file metadata)
at the lab's standard 2D MQS (magneto-quasi-static) formulation,
asking the question:

"When we set up `HelmholtzEquation` with c=0 in 2D, what is the
INSIDE-OUT view of the operator?"

The conclusion (per lab-internal summary):
- The standard 2D MQS equation
    -nu * laplacian(A_z) + jw*sigma*A_z = J_0z
  is mathematically equivalent to the m=0 Helmholtz equation in
  cylindrical coords:
    (1/r) d/dr (r dA/dr) + (something with sigma, mu) = J
- The HelmholtzEquation-with-c=0 trick (E0_mode / H1_mode) is in
  effect a "projection" operator applied AFTER the FE solve, not
  before -- it builds the recursion's right-hand sides from the
  previous step's field.

## Why it matters

The doc clarifies a common misconception:
- "HelmholtzEquation with c=0" is NOT solving a Helmholtz equation
- It is using COMSOL's PDE framework as a way to express
  "u = f on the FE mesh" with the FE shape functions of choice
- The actual PHYSICS is in the source term `f`, which the lab
  constructs via `withsol('sol2', ...)` to encode the recursion.

This view enables the 2D MQS CLN code to be lifted to:
- 2D MQS in non-circular geometries (already done -- rectangles via
  CLN_rectangle_H1_mode_Infinite_Element.m)
- Axisymmetric (r,z) MQS -- needs the cylindrical-divergence form
  with the r * Jacobian (see Sugahara lab Henrotte work for the
  general axisym pattern)
- 3D vector A -- requires HCurl basis instead of scalar Helmholtz,
  not yet implemented in the practice code

## TBD

Specific equations / numerics not transcribed.  See the docx:
public-safe curated corpus
"""


LAB_LESSONS = """
# Lab lessons (pitfalls learned from 01_単一電源 + 2020_11_04 practice)

These are the specific "do this / don't do that" findings that came
out of running the scripts, not from reading the papers.  All are
direct, lab-internal observations.

## 1. Use Stationary, not Frequency Domain, for the basis stages

Each CLN stage is a STATIC FEM problem in the recursion variables --
the time/frequency dependence is in the EXTERNAL ladder ODE, not in
each stage.  Use `model.study('std1').create('stat', 'Stationary');`
NOT `Frequency Domain`.  All practice scripts (CLN_E0_mode.m,
CLN_H1_mode.m, etc.) do this correctly; the mistake to AVOID is
copying from a Frequency-Domain template.

## 2. Segregated vs Fully Coupled solver

Default in the lab scripts is **Segregated** for E0 mode
(`create('se1', 'Segregated');`) and **Fully Coupled** for H1 mode
(`create('fc1', 'FullyCoupled');`).  Reason: H1 mode mixes Jz, Az,
and the InductionCurrents physics's Ax/Ay/Az -- the multiple
strongly-coupled fields converge slower with segregated, and the
Robin BC's Dirichlet/gradient mix is unstable under segregated.

Rule of thumb: if your stage involves only HelmholtzEquation + Poisson
projection (E0 mode), Segregated is OK.  If it involves
`InductionCurrents` physics or Magnetic Potential coupling
(Kelvin / H1-rectangle), use Fully Coupled.

## 3. Shape function order MUST be high enough for late stages

E0 mode driver sets `order=7` on hzeq1/hzeq2/poeq1.  H1 mode default
is `order=1` (the scripts often have `order=3` commented out).  At
stage n the basis fields are polynomial of degree ~2n; if shape
function order p < 2n then the FE approximation is at most piecewise
quadratic-in-stage and the recursion diverges OR the analytical
Legendre formulae stop matching.

Quick fix: bump `prop('ShapeProperty').set('order', 7)` for the
HelmholtzEquation pair (cheap, no nonlinearity) and ~3 for the
InductionCurrents (expensive, but only needs 3 for typical 5-10
stages).

## 4. Probe definition order matters

Probes `G` and `L` MUST be defined BEFORE the recursion loop starts,
NOT inside the loop.  If you create probes inside the loop, COMSOL
re-evaluates them on each iteration and the `withsol('sol2', G)`
reference goes stale.  All practice scripts do this correctly.

## 5. `withsol` references must point to the COMPLETED stage

The lab idiom is:
  - solve sol1 (stage 0)
  - reference 'sol2' in stage 1's source term
  - solve sol1 (which updates sol2 due to StoreSolution)
  - reference 'sol2' in stage 2
  - ...

This depends on `model.sol('sol1').create('su1', 'StoreSolution');`
being present.  If you remove StoreSolution, the recursion
silently uses the wrong solution.

## 6. Mesh density on the conductor BOUNDARY drives accuracy

CLN_E0_mode.m uses `hmax = 1000e-6` (1 mm) on a 3 mm disk -> ~6
elements across the diameter -- this is enough for stages 1-3 but
errors blow up at stage 5+.  CLN_E0_mode_ver_2.m bumps to
`hmax = 100e-6` (0.1 mm) for stable stages 1-9.  For H1 mode the
boundary curvature matters more; use `hmax = 100e-6` on a 3mm disk
from the start.

## 7. The `_NG` filename suffix is a real artifact

`CLN_H1_mode_Kelvin_NG.m` is the lab's pithy "this approach failed"
marker.  Do NOT try to fix Kelvin-coupled H1 CLN -- multiple postdocs
have already tried.  Use Robin or Infinite Element instead.  If
someone proposes Kelvin-CLN coupling in a future paper, the lab's
standing answer is "we tried; see CLN_H1_mode_Kelvin_NG.m".

## 8. The analytical Legendre formulae are TRUSTED REFERENCES

The 9 hard-coded A_n(r) and J_n(r) formulae in CLN_E0_mode.m are the
lab's gold standard for verifying any new CLN implementation
(MATLAB-class refactor, port to Python, port to NGSolve, etc.).
A new implementation MUST reproduce A_1(r) through A_9(r) to within
1e-6 absolute on the unit disk -- otherwise it has a bug.

The formulae for H1 mode are analogous but mu_r-dependent (the
mu_r + 2n^2 - 1 family).

## 9. `CLN_practice_parameter_sweep.m` is a smoke test, not a workflow

The parameter-sweep driver is a 67-line demo of how to use COMSOL's
parametric solver to vary the source amplitude `nn` across stages.
It does NOT do a real CLN recursion (the `nn=1,2,3,4` is just a
linearity check on the Helmholtz solve).  Don't mistake it for the
canonical run; use CLN_E0_mode.m / CLN_H1_mode.m for actual CLN.

## 10. State separation between f-domain and t-domain

In CLN.m, the `f_domain` and `t_domain` methods both overwrite
`cln.II`, `cln.JJ`, `cln.AA`, `cln.t`.  If you run f_domain THEN
t_domain on the same `cln` object, the f-domain results are LOST.
Save the f-domain output explicitly:
  cln_f = f_domain(cln, f, V);   Z_save = cln_f.Z;
  cln_t = t_domain(cln, [0,t1], Vt, I0);

(MATLAB's value-class semantics save you from cross-contamination
within a single method call, but successive calls overwrite.)

## 11. Bessel inequality residual is the best convergence diagnostic

When you don't know if N stages is enough:
  R_N = ||J - J_N||_sigma^2  /  ||J||_sigma^2
where J = the true (e.g. FEM) current density and J_N = the N-stage
CLN reconstruction.  Geometric decrease in R_N vs N means you are
converging; oscillation or growth means you broke orthogonality
(see lessons 3, 6, 8).

This is from the 2016_10_02_ベッセルの不等式.docx theoretical
backing -- the inequality guarantees R_N -> 0 as N -> infinity for
correct basis fields.
"""


KEY_FILES = """
# Top 15 key MATLAB / docx files (with role + line count)

| Lines | File                                                  | Role                                                  |
|-------|-------------------------------------------------------|-------------------------------------------------------|
|    71 | 2020_11_04\\CLN.m                                       | The canonical 71-line CLN class (f_domain / t_domain) |
|   216 | 2020_11_04\\CLN_E0_mode.m                                | E0 (axisymmetric) practice with analytical reference |
|   268 | 2020_11_04\\CLN_E0_mode_ver_2.m                          | E0 with per-stage solutions (sol0..sol{2N})           |
|   179 | 2020_11_04\\CLN_H1_mode.m                                | H1 (m=1) practice with Robin BC                       |
|   225 | 2020_11_04\\CLN_H1_mode_Infinite_Element.m               | H1 with COMSOL Infinite Element BC                    |
|   275 | 2020_11_04\\CLN_H1_mode_Kelvin_NG.m                      | H1 with Kelvin (FAILED -- documented)                 |
|   226 | 2020_11_04\\CLN_rectangle_H1_mode_Infinite_Element.m     | H1 on rectangular conductor + IE BC                   |
|    72 | 2020_11_04\\Check_Robin_BC.m                             | Robin BC isolated verification                        |
|   138 | 2020_11_04\\Kelvin_2D_check_Je_source.m                  | Kelvin coupling test (Je source variant)              |
|   133 | 2020_11_04\\Kelvin_2D_check_fixed_bc.m                   | Kelvin coupling test (fixed BC variant)               |
|    72 | 2020_11_04\\Check_CLN_class.m                            | End-to-end test of CLN.m vs Bessel reference          |
|    72 | 2020_11_04\\case01.m                                     | Smoke test, single Poisson with analytic source func  |
|   256 | 2020_11_04\\case02.m                                     | Smoke test, COMSOL.mph-style scripted driver           |
|    67 | 2020_11_04\\CLN_practice_parameter_sweep.m               | Parametric COMSOL sweep demo (not a CLN run)          |
|   239 | 2020_11_04\\関数読み込みの練習.m                          | "function-loading practice" -- helper drafts           |

## Subfolders also absorbed

| Path                                                     | Notes                                              |
|----------------------------------------------------------|----------------------------------------------------|
| 2020_11_04\\解析解\\                                       | Analytical references (CLN_z.m, E0_mode_analytical_*)|
| 2020_11_04\\解析解\\symbolic_{Dirichlet,Hiruma,Robin}\\    | Per-BC symbolic-toolbox derivations               |
| 01_単一電源\\1-wire\\                                       | Single-wire dataset (.dat files, .asc, .raw)      |
| 01_単一電源\\1-wire_open\\                                  | Open-BC single-wire dataset                        |
| 01_単一電源\\2-wire_+{+,-}\\                                | Two-wire pair (parallel / opposite polarity)       |
| 01_単一電源\\2-wire_++_寄り_NG\\                            | "Close 2-wire failed" -- another _NG marker        |
| 01_単一電源\\circle\\                                       | Single circle test data                            |
| 01_単一電源\\many_wire\\                                    | Multi-wire generalisation tests                    |
| 02_境界条件\\busbar\\                                       | Busbar full BC comparison (IABC vs Kelvin vs Robin)|
| 02_境界条件\\ac_iabc\\                                      | AC IABC convergence study                           |

## docx writeups (binary, see files)

| File                                                  | What it documents                            |
|-------------------------------------------------------|----------------------------------------------|
| 01_単一電源\\160913.docx                                | Initial 2016-09-13 notes on the formulation |
| 01_単一電源\\2016_12_01_2導体モデルに亀有法を適用.docx  | 2-conductor Kameari method application      |
| 01_単一電源\\CLN(Cauer Ladder Network) method with A-_ formulation.docx | A-phi formulation overview |
| 01_単一電源\\Foster.docx                               | Foster (vs Cauer) network comparison         |
| 01_単一電源\\Kameari法における直交性について.doc        | Orthogonality of Kameari's method            |
| 01_単一電源\\A法とA-Phi法.docx                         | A vs A-phi formulation                       |
| 01_単一電源\\亀有法の理解のための簡単なモデル検討.docx   | Simple-model study for understanding         |
| 01_単一電源\\円柱渦電流Cauer梯子回路の計算例.docx        | Cylindrical eddy-current Cauer example       |
| 01_単一電源\\円柱渦電流の（円筒Ｅモード）のCauer梯子回路の表現.docx | Cylinder E-mode Cauer expression  |
| 01_単一電源\\回路方程式2.docx                          | Circuit equations writeup                     |
| 01_単一電源\\有限要素法.doc                            | FEM background notes                           |
| 01_単一電源\\有限要素空間上でのKameari法について.doc    | Kameari method on FE spaces                   |
| 01_単一電源\\疑問点.docx                               | Open questions list                            |
| 01_単一電源\\CLN回路方程式 <-> Kameari法 漸化式.docx    | CLN circuit eqn <-> Kameari recursion link    |
| 02_境界条件\\IABC境界条件に関して.docx                  | IABC BC theory                                 |
| 02_境界条件\\IABC層内の磁気エネルギー.docx              | Magnetic energy in IABC layer                  |
| 02_境界条件\\AJ積分.docx                               | A.J integral identity                          |
| 02_境界条件\\2017_05_21_境界条件.docx                  | 2017 BC overview                                |
| 09_級数展開\\10-04_old.docx                            | Series expansion theory (pre-amendment)        |
| 09_級数展開\\亀有追記0921.docx                          | Kameari's amendment to convergence proof       |

## PDF references

| File                                              | Notes                                  |
|---------------------------------------------------|----------------------------------------|
| A-phi-formulation-2.pdf  (root)                   | A-phi formulation reference            |
| 01_単一電源\\intern2.pdf                            | Intern report                          |
| 02_境界条件\\busbar\\Canova_Giaccone2009_pre_print_version.pdf | Busbar BC paper (2009)      |

## Folders explicitly NOT in this module's scope

(For Group B / future absorption -- not in public-safe curated corpus
Group A list)

  - 03_複数電源       (multi-source -- needs separate module)
  - 04_展開点         (expansion points -- belongs with multi_expansion in cln_knowledge.py)
  - 05_高周波         (high-frequency)
  - 06_非線形         (nonlinear -- belongs with cln_knowledge.py nonlinear topic)
  - 07_電気力学       (electrodynamics)
  - 08_終端問題       (termination problem)
  - 10_タイル         (tile)
  - 11_3次元          (3D)
  - 12_比留間法        (Hiruma's method)
  - 13_BEM+FEM        (BEM+FEM hybrid)
  - 14_ブロッホ方程式  (Bloch equation)
  - 15_誤差論@長嶺氏   (error theory by Nagamine)
  - 16_FP_CLN          (FP-CLN variant)
  - 2017_05_17_インバータ回路.docx
  - 2020_12_07_高周波CLNの練習
  - 2021_01_22_CauerI_to_CauerII
  - 2022_10_09_遠藤@法政
  - 2023_08_08_松本@法政
  - 2026_04_01_長方形CLN

These can be absorbed in a Group B / Group C follow-up module.
"""


CITATIONS = """
# Citations (paper + docx + script references)

## Canonical CLN papers (Sugahara lab)

@article{kameari2018cln,
  author = {Kameari, A. and Ebrahimi, H. and Sugahara, K.
            and Shindo, Y. and Matsuo, T.},
  title  = {Cauer Ladder Network Representation of Eddy-Current
            Fields for Model Order Reduction Using Finite-Element
            Method},
  journal = {IEEE Trans. Magn.},
  volume = {54},
  number = {3},
  pages  = {7201804},
  year   = {2018}
}

@article{sugahara2018unbounded,
  author = {Sugahara, K. and Kameari, A. and Ebrahimi, H.
            and Shindo, Y. and Matsuo, T.},
  title  = {Finite-Element Analysis of Unbounded Eddy-Current
            Problems Using Cauer Ladder Network Method},
  journal = {IEEE Trans. Magn.},
  volume = {54},
  number = {3},
  pages  = {7200704},
  year   = {2018}
}

@article{matsuo2018matrix,
  author = {Matsuo, T. and Kameari, A. and Sugahara, K. and Shindo, Y.},
  title  = {Matrix Formulation of the Cauer Ladder Network Method
            for Efficient Eddy-Current Analysis},
  journal = {IEEE Trans. Magn.},
  volume = {54},
  number = {11},
  pages  = {7205805},
  year   = {2018}
}

@article{kuriyama2019multiple,
  author = {Kuriyama, K. and Kameari, A. and Ebrahimi, H.
            and Fujiwara, K. and Sugahara, K. and Shindo, Y.
            and Matsuo, T.},
  title  = {Cauer Ladder Network With Multiple Expansion Points
            for Efficient Model Order Reduction of Eddy-Current
            Field},
  journal = {IEEE Trans. Magn.},
  volume = {55},
  number = {6},
  pages  = {7203404},
  year   = {2019}
}

@article{ebrahimi2020modal,
  author = {Ebrahimi, H. and Sugahara, K. and Matsuo, T.
            and Kaimori, H. and Kameari, A.},
  title  = {Modal Decomposition of 3-D Quasi-Static Maxwell
            Equations by Cauer Ladder Network Representation},
  journal = {IEEE Trans. Magn.},
  volume = {56},
  number = {3},
  pages  = {7513004},
  year   = {2020}
}

## Background / sibling references

@article{canova2009busbar,
  author = {Canova, A. and Giaccone, L.},
  title  = {(See 02_境界条件\\busbar\\Canova_Giaccone2009_pre_print_version.pdf)},
  journal = {(TBD -- pre-print version filed in busbar/)},
  year   = {2009},
  note   = {Used in the lab's busbar BC comparison study}
}

## Internal lab documents (docx, not formally published)

Cited inline in this module by full path under public-safe curated corpus
(see key_files topic for the complete table).  Authors are lab members
(Sugahara, Kameari attribution noted on '亀有追記' files).  Dates
embedded in filenames (2016, 2017, 2020 ranges) where available.

## Script-level credits

  - CLN.m (71 lines)                        - lab canonical, Sugahara
  - CLN_E0_mode.m / CLN_E0_mode_ver_2.m     - lab, 2020-11-04 dated
  - CLN_H1_mode.m and variants              - lab, 2020-11-04 dated
  - Kelvin_2D_check_*.m                     - lab, Kelvin verification
  - Check_Robin_BC.m                        - lab, BC isolation test
  - Check_CLN_class.m                       - lab, end-to-end check
  - 解析解\\symbolic_*\\*.m                  - lab, symbolic-toolbox
                                              derivations per BC

## Related radia-mcp modules

  - radia_mcp.mor.cln_knowledge             - paper-level CLN (sibling)
  - radia_mcp.mor.systematic_knowledge      - MOR systematic overview
  - radia_mcp.mor.bibliography_index_knowledge  - MOR bibliography

## TBD citations

- The 2D quasi-static rethink docx author + date
- The Bessel inequality docx 2016-10-02 author (likely Kameari)
- The A-phi formulation PDF -- pre-print of which paper?
- The intern report intern2.pdf -- whose intern, what year?
"""


_TOPICS = {
    "overview": OVERVIEW,
    "matlab_class": CLN_MATLAB_CLASS,
    "e0_mode": E0_MODE_PRACTICE,
    "h1_mode": H1_MODE_PRACTICE,
    "boundary_conditions": BOUNDARY_CONDITIONS,
    "robin_bc": ROBIN_BC,
    "series_expansion": SERIES_EXPANSION_THEORY,
    "a_phi": A_PHI_FORMULATION,
    "quasi_static_2d": QUASI_STATIC_2D_RETHINK,
    "lab_lessons": LAB_LESSONS,
    "key_files": KEY_FILES,
    "citations": CITATIONS,
    "all": "",
}
_TOPICS["all"] = "\n\n".join(v for k, v in _TOPICS.items() if k != "all")


def get_cln_practice_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      all                   - everything concatenated
      overview              - module scope and source folders
      matlab_class          - CLN.m 71-line class explained
      e0_mode               - E0 (axisymmetric) COMSOL practice
      h1_mode               - H1 (m=1) COMSOL practice (3 BC variants)
      boundary_conditions   - Kelvin vs Infinite Element vs Dirichlet
      robin_bc              - Check_Robin_BC.m verification pattern
      series_expansion      - Bessel inequality + convergence theory
      a_phi                 - A-phi formulation reference
      quasi_static_2d       - 2D MQS rethink doc
      lab_lessons           - hard-won pitfalls (11 lessons)
      key_files             - top 15 .m files with role / line count
      citations             - papers + docx + script credits
    """
    topic = topic.lower().strip()
    if topic in _TOPICS:
        return _TOPICS[topic]
    return (
        "Unknown topic '" + topic + "'. Available: "
        + ", ".join(sorted(_TOPICS.keys()))
    )
