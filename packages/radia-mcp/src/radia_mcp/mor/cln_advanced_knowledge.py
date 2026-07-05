"""CLN advanced extensions: HF / nonlinear / circuit / Bloch /
finite-period (FP).

Covers public-safe curated corpus, 06_非線形, 07_電気力学,
14_ブロッホ方程式, 16_FP_CLN, 2020_12_07_高周波CLNの練習}.

Sister module to `radia_mcp.mor.cln_knowledge` (basic linear CLN).
This module is the *advanced physics* layer: where the basic Cauer
ladder breaks (high frequency / nonlinear saturation / hysteresis /
mechanical coupling) and which lab-developed extension recovers it.
"""

OVERVIEW = """
# CLN advanced extensions -- where basic Cauer breaks

The basic CLN (cln_knowledge.py) assumes:
  (a) low frequency (quasi-static, no displacement current),
  (b) linear magnetic material (mu = const),
  (c) the EM problem is the only physics in the loop, and
  (d) iron sheets / cores can be lumped into ONE bulk stage.

Each of those assumptions has a counter-example in real engineering
problems, and each has its own Sugahara-Lab WIP folder in
public-safe curated corpus  This module catalogues the FOUR
extensions and the ONE cross-application:

| Extension | Folder | Status (as observed) | What breaks without it |
|-----------|--------|----------------------|------------------------|
| High frequency / Helmholtz | 05_高周波, 2020_12_07_高周波CLNの練習 | symbolic prototypes | wave skin effect at GHz, displacement current |
| Nonlinear / saturation     | 06_非線形 | three generations of attempts | iron at saturation, induction heating |
| Electrodynamic / circuit   | 07_電気力学 | FreeFEM-based magnetic-bearing study | force vs frequency, mechanical coupling |
| FP-CLN (Fixed-Point)       | 16_FP_CLN | CEFC 2024 paper + production MATLAB classes | hysteresis, harmonic balance, time-domain nonlinear |
| Bloch equation (NMR/MRI)   | 14_ブロッホ方程式 | one design document (2019) | application sketch, not a working CLN extension |

The single biggest message:

  FP-CLN (Fixed-Point CLN) is the lab's working solution for
  nonlinear + hysteresis problems.  It linearizes the LHS with a
  virtual constant permeability mu_FP and pushes ALL nonlinearity
  into a controlled-source term on the RHS.  Basis functions and
  Cauer R, L matrices stay CONSTANT.  CEFC 2024:
  Sugahara-Tobita-Matsuo-Takahashi, "Fixed-Point Cauer Ladder
  Network Method for Eddy-current Problems with Hysteresis."

The earlier nonlinear CLN attempts (06_非線形, 2017-2020) all
parameterized the circuit elements themselves, which is what FP-CLN
explicitly avoids.  See section NONLINEAR for that lineage.
"""

HIGH_FREQUENCY = """
# High-frequency CLN (Helmholtz regime)

Folders:
  public-safe curated corpus               (15 files, 14 .m)
  public-safe curated corpus  (9 .m + 4 .mph)

## Motivation -- where basic CLN runs out

Basic CLN is built on the magnetic-quasi-static (MQS) eddy-current
equation
   curl(nu curl A) + sigma * dA/dt = J_s
which drops displacement current.  At HF (above a few hundred MHz
for a 1 mm conductor in copper) the displacement current term
   epsilon * d^2 A / dt^2
is no longer negligible.  The PDE becomes the Helmholtz equation
   curl(nu curl A) + sigma dA/dt + epsilon d2A/dt2 = J_s
and the Cauer recursion changes: each "magnetic step" now uses the
Helmholtz Green's function (Bessel J_0, J_1 of a complex argument)
instead of the Laplace Green's function used in MQS-CLN.

This is the "Helmholtz regime" extension.  Note that the Radia
project itself has an explicit POLICY against shipping a Helmholtz
kernel in the C++ core (see CLAUDE.md "Green's Function: Laplace
Kernel Only (MQS/Darwin)").  HF-CLN is therefore a research /
benchmarking layer in MATLAB, not a feature of the Radia C++
solver.

## What the lab actually wrote

Two-pass symbolic approach in MATLAB Symbolic Math Toolbox:

  - GenKameariBasis.m  -- recursion for J_n(r), A_n(r) basis on a
    circular conductor (radius a).  Generates `kameari.m` with
    explicit J(1)..J(11), A(1)..A(11) up to order r^20.
  - GenPowerSeries.m   -- alternative recursion producing
    `power_series.m` / `power_series2.m` for two boundary
    treatments (nType=1 vs nType=2).
  - bessel_taylor.m    -- Taylor expansion of (a*k*J_0(kr)) /
    (2*J_1(ka)) -- the Bessel ground truth that R_n, L_n converge
    to.
  - CLN.m, CLN2.m      -- symbolic Cauer R_n, L_n recursion for
    small N=3 + Taylor expansion of input admittance h(s) around
    s = -k_d^2 / (sigma mu d^2).  `switch 'H'/'E'` selects the
    boundary trimming (set H or set E zero at surface), giving
    different L_n tails.
  - kameari.m, plot_kameari.m, power_series.m, power_series2.m,
    plot_power_series1.m, plot_power_series2.m, check_function*.m
                        -- generated code + plotting drivers,
    compared against the lossless Bessel limit.

## High-frequency lossless-ring driver (2020_12_07 練習)

`check_CLN_high_freq_lossless.m` is the most useful HF prototype:
SYMBOLIC `dsolve` of the Helmholtz equation on a coaxial ring
(radii a, b):
    d/dr( r dAA/dr ) / r  - m^2 AA / r^2  =  - mu_0 * mu_r * J / cos(m t)
two BC sets:
  case 4:  AA(a) = 0, AA(b) = 0      (closed coax)
  case 5:  AA(a) = 0, AA(b) free     (open at outer radius)

J_{n+1} update uses epsilon (displacement) instead of sigma (eddy):
   J_{n+1} = simplify( J_n - epsilon * A_n / L_n )
That is the LITERAL switch from MQS-CLN to Helmholtz-CLN.  Cross-
coupling matrices C_{n1,n2}, L_{n1,n2} are integrated symbolically
and saved to `high_freq_lossless_{4,5}.mat`.

`check_CLN_high_freq_infinity.m` runs the same recursion in the
b -> infinity limit (radiating mode, N=2, f=300 MHz).

## Lab findings (what worked, what did not)

- Convergence: HF Cauer recursion in symbolic form converges to the
  Bessel-ratio AC impedance only when the boundary trimming is
  consistent with the chosen formulation (set-H zero or set-E zero
  on the surface).  Mixing them gives non-orthogonal basis.

- Conditioning: symbolic expansion at order r^20 starts producing
  coefficients in the 10^30 range.  Numerical evaluation needs
  variable-precision arithmetic (`vpa`); IEEE double precision
  loses accuracy.

- No FEM driver yet: 05_高周波 stays at the analytical-circular-
  wire toy problem and never crosses over to 2D/3D FEM.  The lab
  did not promote HF-CLN into a working FEM solver -- the relevant
  application (induction heating workpieces) lives in `radia.peec`
  and `radia.esim_cell_problem` instead.

## Companion document

`互除法とPade近似.pptx` (in 05_高周波) -- discusses the relation
between the Cauer ladder continued fraction and Pade approximation
via the Euclidean algorithm.  TBD content; not opened in this scan.
"""

HIGH_FREQUENCY_CODE = """
# High-frequency CLN code excerpts

## The Helmholtz-Cauer recursion (check_CLN_high_freq_lossless.m)

The key driver loop -- note the use of `dsolve` to handle the
Helmholtz ODE symbolically at each stage, and the J update with
epsilon (displacement current) instead of the usual MQS sigma:

```matlab
syms AA(r) A mur mu0 sig t a r
DA = diff(AA);
% ... boundary conditions cond1, cond2 set per case ...

s = j*2*pi*f;
k = 2*pi*f*sqrt(mur*mu0*ep0);
rr = linspace(ra,rb,100);

for n = [1:N]
    R(n) = 1/double(int(int(J(n)^2*r,t,0,2*pi),r,ra,rb));
    ode = diff(r*diff(AA,r),r)/r - m^2*AA/r^2 ...
          == -mu0*mur*J(n)/cos(m*t);
    if n == 1
        A(n) = simplify( R(n)*dsolve(ode,[cond1 cond2])*cos(m*t) );
    else
        A(n) = simplify( A(n-1) ...
              + R(n)*dsolve(ode,[cond1 cond2])*cos(m*t) );
    end
    L(n) = double( ep0*R(n)*int(int(J(n)*A(n)*r,t,0,2*pi),r,ra,rb) );
    J(n+1) = simplify( J(n) - ep0*A(n)/L(n) );   % <-- ep0, not sig
end
```

Contrast with the standard *low-freq* CLN recursion in
`cln_knowledge.py` ("recursion" topic): the only structural change
is `ep0` replacing `sig` in the J update, and the recursion using
the Helmholtz operator (curl curl + k^2) instead of the Laplace
operator.

## Gram-Schmidt construction of the radial basis (GramShmidt.m)

The high-freq basis is built by Gram-Schmidt orthogonalization
against the conductor cross-section measure 2*pi*r dr:

```matlab
syms x r a
N = 5;
w = 2*pi*sym('r');
for n = [0:N]
    sum = sym('0');
    for k = [0:n-1]
        sum = sum + simplify( ...
            int(w*Pd(k+1)*r^(2*n),r,0,a) ...
            / int(w*Pd(k+1)^2,r,0,a) * Pd(k+1) );
    end
    Pd(n+1) = r^(2*n) - simplify(sum);
end
```

This produces the "axisymmetric Legendre-like" polynomials on
[0, a] that the Cauer-Kameari basis turns out to be (modulo
boundary trimming).
"""

NONLINEAR = """
# Nonlinear CLN (saturation + hysteresis)

Folder: public-safe curated corpus
  Subfolders (chronological generations):
   - 2017_11_19_FreeFEM_非線形の練習         (FEMM driver, GenFEMM.m)
   - 2017_12_12_非線形の案_01                 (first symbolic + numeric attempt)
   - 2019_10_16_非線形ランチョス              (Lanczos-style nonlinear CLN)
   - 2020_08_06_京大資料_非線形              (Kyoto Univ slides, PDF)
   - 2020_08_09_非線形CLNの練習              (COMSOL LiveLink + tanh BH)
   - 2023_09_26_飛田モデル_jω法              (Tobita's jω model, the
                                              pre-FP-CLN production line)

## The hardness

The basic CLN relies on the constant-mu assumption to factor the
linear operator ONCE and re-use the L_n, R_n matrices across all
time steps / frequency points.  When mu = mu(|B|), the operator is
no longer constant and the obvious choices are:

  (a) Re-build the operator at every step.  Defeats the MOR purpose.
  (b) Linearize around an operating point (local mu_eff).  Each
      operating point needs its own CLN; recombination is awkward.
  (c) Use a fixed-point split: linear left, nonlinear right.

The lab tried (a) and (b) for several years before settling on (c),
which is the FP-CLN of section FP_CLN below.

## Lineage of attempts (the 06_非線形 folder)

### Generation 1: 2017_12_12_非線形の案_01 (parameterized basis)

First attempt parameterized the basis: a new J_n(r; sigma, mu)
rebuilt at every step from local mu(B).  `non_linear_CLN.m` (200
lines) precomputes 8 R_n, L_n and basis tables H, E from
`1..8.mat`, then time-steps a Cauer ODE with either
  case 1: linear mu (BB/mu -> HH), or
  case 2: tanh saturation (HH = atanh(BB/2)/mu)
and compares J(r,t) reconstruction to the Bessel reference.
Companion: `genBH.m` (synthetic BH curve B = 2 tanh(mu H / 2)).

Companion docs: `2017_{09_26_BH空間,12_05_非線形,12_12_非線形}_Kameari
法.docx`, `円形導体の渦電流問題.docx`, `14-02.pdf` (TBD digest).

### Generation 2: 2019_10_16_非線形ランチョス (Lanczos / fixed point)

Pivots to Lanczos-style basis generation:
  ./E0.m, gen_BH.m, gen_tV.m, plotV.m, main.m, compare_rJ.m
  + subfolders CLN_transient, CLN_transient_fixed_point_{1,2},
    基底の構築 (basis construction), 静磁場 (magnetostatic).
The driver `main.m` loops over (I1, I2, I3) currents and caches
FreeFEM++ solutions in `mat/<i0>_<i2>_<i4>.mat`.  This is the
proto-FP-CLN -- the IDEA of linear-LHS + nonlinear-RHS is here,
but the basis is still rebuilt each iteration.

### Generation 3: 2020_08_09_非線形CLNの練習 (COMSOL LiveLink)

Validation against COMSOL FEM on a 2D circular conductor:
  case0.m   COMSOL pipeline: 8 analytic functions a_n_m(x,y)
            = J_n * A_m on circle radius a, FreeQuad mesh,
            integrate vs symbolic S0(n,m) -- orthogonality check.
  case1.m   only A_n exported, S0(n) = int(2 pi r A_n) check.
  case2.m   the BIG one: sweeps (i0, i2, i4), computes tanh-
            saturated H(B):
              Az = i0*A_0 + i2*A_2 + i4*A_4
              Bx, By = curl Az;  |B| = sqrt(Bx^2 + By^2)
              Hx, Hy = atanh(|B|)/|B| * (Bx, By)
              Jz = dHy/dx - dHx/dy
            integrates J_n * Jz per stage.  Saves to mat/.
  compare_freefem.m -- cross-check vs FreeFEM++ at
                       public-safe curated corpus

This is the ONLY 06_非線形 file that drives nonlinear CLN through
a full FEM (COMSOL) end-to-end.  Companion .mph: `1D_FEM.mph`,
`2020_08_11_comsol_test.mph`, `2021_04_06_1D_FEMの練習.mph`,
`comsol_test.mph`, `FEM_1D.mph`, `積分.mph`.  Pptx
`2023_03_03_FP-CLN.pptx` is the moment the FP-CLN naming lands;
work then moves to 16_FP_CLN.

### Generation 4: 2023_09_26_飛田モデル_jω法 (Tobita's jω model)

Immediate FP-CLN predecessor.  Files: `Gen_FEMM_model.m`,
`read_mesh.m`, six `compare_results_*.m` (CLN vs f-domain vs
contour vs line).  Companion: `2023_10_03_委員会資料_飛田法.pptx`.
The "jω" nomenclature is the frequency-domain harmonic-balance
treatment that 16_FP_CLN/HaromonicBalance takes over.

## Lab lesson (from the entire 06_非線形 arc)

Five+ years of attempts to "build the nonlinear basis directly"
were structurally complicated AND made the resulting CLN circuit
elements time-varying, which negates the MOR benefit.  The eventual
fix (FP-CLN, 16_FP_CLN, 2023 CEFC paper) was to KEEP THE BASIS
LINEAR with a virtual mu_FP and push the nonlinearity into a
controlled-source term -- the lesson is that the "MOR-friendly"
linearization wins even when the virtual permeability is "wrong",
because the iteration still converges and the circuit stays small.
"""

NONLINEAR_CODE = """
# Nonlinear CLN code excerpts (06_非線形)

## Non-linear CLN driver (2017_12_12, non_linear_CLN.m)

The first-generation attempt -- B-input with tanh saturation,
iterating over a Cauer ladder with time-stepping:

```matlab
% Pre-loaded R(n), L(n), basis tables H(:,n), E(:,n) from 1..8.mat
mu0 = 4*pi*1e-7;  sig = 5.8e7;  a = 0.001;
N = 7;

% Synthetic BH curve: linear vs tanh-saturated
%  B = mu*H            (linear)
%  B = 2*tanh(mu*H/2)  (saturating)

% Cauer block U(n,m) = 1 if n <= m, else 0
U = triu( ones(N,N) );

% Driver loop -- case 2 = saturating
b = zeros(N,1);
for nt = [1:length(t)]
    BB = mu*H*b;            % H is the precomputed basis matrix
    HH = atanh(BB/2)/mu;    % invert tanh saturation
    for n = [1:N]
        h(n,1) = trapz(r, 2*pi*mu*r.*H(:,n).*HH/L(n));
    end
    e(1,1) = R(1)*e0*cos(omega*t(nt));
    for n = [2:N+1]
        e(n,1) = R(n) * ( e(n-1,1)/R(n-1) - h(n-1) );
    end
    dbdt = U*e([2:N+1],1);
    b = b + diag(L) \\ dbdt * dt;
end
```

Problem: `H` (the basis matrix) was built assuming constant mu, so
when mu changes effectively due to saturation, this driver does NOT
update the basis -- it only updates the *coefficients*.  That's the
proto-form of the fixed-point split, but without the proper FP-CLN
derivation that 16_FP_CLN provides.

## Per-mu polynomial basis library (2019_10_16, E0.m)

When mu is the *only* nonlinear parameter and the geometry is the
canonical circular conductor, the basis is closed-form in (r, a,
sigma, mu_r, mu0).  E0.m provides the library up to N=9:

```matlab
function CLN = Hn(r,t,a,sig,mur,mu0)
    CLN(1).R = 1/(a^2 * sig * pi);
    CLN(2).R = 3/(a^2 * sig * pi);
    CLN(3).R = 5/(a^2 * sig * pi);
    % ... 2*n - 1 pattern up to n = 9
    CLN(1).L = (mu0 * mur)/(8 * pi);
    CLN(2).L = (mu0 * mur)/(16 * pi);
    % ... 8*n pattern up to n = 9
    CLN(1).J = sig;
    CLN(2).J = -(sig.*(a.^2 - 2.*r.^2))./a.^2;
    CLN(3).J = (sig.*(a.^4 + 6.*r.^4 - 6.*a.^2.*r.^2))./a.^4;
    % ... polynomials of degree 2(n-1) up to n = 9
end
```

Same J_n polynomials reappear in the FP-CLN `E0.m` for the cylinder
problem (16_FP_CLN/TimeDomain/円筒_*).  These closed-form bases are
the only reason the entire CLN program is benchmarkable -- 1D round-
conductor problems have analytic Bessel ground truth.

## COMSOL LiveLink basis-check (2020_08_09, case0.m)

```matlab
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('Model');
model.modelNode.create('comp1', true);

% Declare the analytic basis functions a_n_m(x,y) = J_n * A_m
for n = [1:N]
for m = [1:M]
    model.func.create(sprintf('a_%d_%d',n,m), 'Analytic');
    model.func(sprintf('a_%d_%d',n,m)).model('comp1');
    model.func(sprintf('a_%d_%d',n,m)).set('expr', ...
        char(simplify(subs( J(n)*A(m), 'r', sqrt(x^2+y^2) ))));
    model.func(sprintf('a_%d_%d',n,m)).set('args', {'x' 'y'});
end
end

% Mesh & integrate -- compare against symbolic ground truth S0(n,m)
for n = [1:N]; for m = [1:M]
    model.result.evaluationGroup.create(...)
    S(n,m) = mphtable(model, sprintf('eg_%d_%d(x,y)',n,m)).data;
    e(n,m) = (S(n,m) - S0(n,m))/S0(n,m)*100;
end; end
```

This is the most useful reference for understanding *how the lab
validates a CLN basis on a real FEM mesh* -- precompute the basis
symbolically, push it into an FE tool as analytic functions, do the
integrations there, and compare to the symbolic ground truth.
"""

ELECTRODYNAMIC_CIRCUIT = """
# Electrodynamic / circuit-coupled CLN (07_電気力学)

Folder: public-safe curated corpus
  Main artifacts:
   - shindo.edp                            (FreeFEM++ driver, 437 lines)
   - 20170509.docm                         (Word doc-with-macros, TBD)
   - Kameari法による振動解析r2.docx        (Kameari method for vibration)
   - 回路連成問題のFEM解析へのKameari法の適用.docx
                                            (Kameari method applied to
                                             circuit-coupled FEM)
   - FreeFemによる節点力.docx              (FreeFEM nodal force)
   - 簡単にz方向の力を求める方法.docx       (simple z-force formula)
   - 運動を含む電磁場.docx                  (EM fields with motion)
   - mian.m                                 (one-line Biot-Savart test)
   - 進藤問題.xlsx                          (Shindo problem data sheet)
   - model.dxf                              (geometry)
   - FEM_1.fem, FEM_1.ans                   (FEMM input + result)
   - アニメーション ジュール損分布 構成a(45°).wmv  (Joule loss animation)
   - EMSolution/input.ems, inputAC.ems, pre_geom2D.neu

## Scope -- what this folder actually is

This is the lab's "Kameari method for a circuit-coupled axisym
magnetic actuator" study.  The point is NOT to derive a new CLN
recursion; the point is to *use* the CLN as a circuit-equivalent
block that plugs into a larger SPICE-like AC/DC analysis of an
electromagnet driving a movable disk through air gap.

The model (shindo.edp) is the canonical "Shindo problem": an
axisymmetric inner+outer yoke pair, a coil wound between them, an
iron disk on the under side separated by a 1 mm air gap, all
embedded in an outer air region.  Material constants: mu_r = 500
for yoke and disk; sigma_coil = 45.6 MS/m (space-factor corrected);
sigma_iron = 8 MS/m.

The CLN recursion (the 11-stage E0/A1/E2/A3/.../E10 chain inside
shindo.edp) extracts:
  - lambda(0)..lambda(10) -- the per-stage R, L "amplitudes"
  - dF(0)..dF(4)          -- the per-stage z-force coefficients
    via the Maxwell stress over the half-air-gap surface
  - FDC                   -- the static-force constant

Then it runs a frequency sweep:

```matlab
for freq = 0.001 ... 1000  (decade per 10 steps)
    omega = 2*pi*freq;
    h(n) = sum_j Y(n,j) * CoilDCVoltage;   % per-stage current
    Force = sum_n dF(n) * h(n);             % combined force
    ForceI = Force * (Vdc / sum h);         % current-input variant
    log(freq, |ForceI|, arg ForceI)
end
```

## Why this is the "electrodynamic / circuit" CLN

Three physics flavours in one recursion:
  (1) FEM eddy-current in iron disk + yoke (standard CLN A-form),
  (2) circuit constraint: coil sees fixed DC voltage (2.450 V)
      across 1.040 ohm DC resistance -- effective AC impedance
      depends on the magnetic-reaction shift handled by CLN,
  (3) per-mode force via Maxwell stress on a thin band at the gap
      centre, dF(n) contributions.
Output: per-frequency force transfer function ForceI(jomega) for
controller design / mechanical-bandwidth plots.  This is the
canonical "use CLN as a SPICE block" pattern.

## Lab finding (from filenames; .docx are binary)

The lab established that Kameari/CLN exports as a small (NN=5
stage) RL block + per-mode dF for a 1-DOF mechanical equation --
"magnetic actuator as high-order RL circuit driving a lumped
mass".  The single .wmv ("Joule loss distribution, 45 deg phase")
animates J^2 on the disk mesh showing where AC losses concentrate.

## Status

This is a one-off study (2017, Shindo problem) rather than a
generally-applicable extension.  The technique generalises -- any
axisym CLN-derived [R, L] + per-mode force coefficients dF can be
inserted into SPICE -- but the lab did not productionize a "circuit
export" path the way it did for FP-CLN.

Related earlier paper:
  - Shindo-Kameari-Sugahara-Matsuo, "Dynamical model of an
    electromagnet by Cauer ladder network representation of eddy-
    current fields," LDIA 2017 (referenced in
    16_FP_CLN/Digests/references.bib).

The mover (translation) case is treated in the separate paper:
  - Shindo-Yamamoto-Sugahara-Matsuo-Kameari,
    "Equivalent Circuit in Cauer Form for Eddy Current Field
    Including a Translational Mover," IEEE Trans. Magn. 56(12), 2020.
And independently:
  - Sabariego-Vanbroekhoven-Gyselinck-Kuo-Peng, "RL-Ladder Circuit
    Models for Eddy-Current Problems With Translational Movement,"
    IEEE Trans. Magn. 58(9), 2022.
Neither of these has source code in public-safe curated corpus, they appear only
as citations.
"""

BLOCH_EQUATION = """
# Bloch equation / NMR-MRI cross-application (14_ブロッホ方程式)

Folder: public-safe curated corpus
  Single file: 2019_08_11_ブロッホ方程式.docx

## Honest scope statement

This folder is ONE design document (2019-08-11) and nothing else.
No MATLAB, no .mph, no driver code, no figures.  The file is binary
.docx so its contents were not opened during this scan; only the
existence and the title were observed.

## What it is intended to do (inferred from the title)

The Bloch equation
    dM/dt = gamma * M x B - (1/T1)(M_z - M_eq) e_z - (1/T2)(M_x e_x + M_y e_y)
governs nuclear magnetization in an NMR / MRI experiment.  It is a
*first-order ODE in time* with a NONLINEAR coupling term (M x B)
between the magnetization vector M and an externally controlled RF
field B.

There is a natural mathematical analogy between (a) the Bloch
relaxation operator with T1, T2 time constants and (b) a CLN ladder
where each stage has its own R_n / L_n time constant tau_n = L_n /
R_n.  A multi-tau Bloch model can be viewed as a CLN ladder driven
by an RF field.

## Status -- speculative

The document very likely sketches THE MAPPING: which Bloch-equation
quantity corresponds to which CLN stage variable.  But there is no
implementation, no validation, and no follow-up file in the same
folder.  This is a "potential future direction" placeholder rather
than a working extension.

For an MCP user asking "can I use CLN for NMR/MRI?", the honest
answer is: the lab has thought about it (2019), the document is at
public-safe curated corpus,
and there is no code yet.  See instead `radia_mcp.nmr_mri` for the
lab's actually-implemented NMR-MRI knowledge layer.
"""

FP_CLN = """
# FP-CLN -- Fixed-Point Cauer Ladder Network (16_FP_CLN)

Folder: public-safe curated corpus
  Subfolders:
   - HaromonicBalance\\       (harmonic-balance / frequency-domain
                                FP-CLN)
   - TimeDomain\\             (time-stepping FP-CLN, including the
                                cylinder and steel-sheet benchmarks)
   - Digests\\                (CEFC 2024 paper drafts -- the
                                canonical "FP-CLN" reference)
   - 疑問点\\                  (open questions / iteration concerns)

## What "FP" stands for

FP = Fixed-Point (固定点) -- as in the *fixed-point iteration*
technique for nonlinear field problems.  This is the lab's
*production* nonlinear CLN, finalised in CEFC 2024 as
   K. Sugahara, M. Tobita, T. Matsuo, Y. Takahashi
   "Fixed-Point Cauer Ladder Network Method for Eddy-current
    Problems with Hysteresis"
   CEFC 2024, JSPS KAKENHI JP20K04443 + Power Academy grant
   (file: 16_FP_CLN/Digests/CEFC_digest.tex, .pdf)

It is NOT "finite period CLN" or "fractional polynomial CLN" --
those are mis-readings.

## The single key idea

Start from the nonlinear eddy-current A-formulation
    curl( nu(B) curl A ) + sigma dA/dt = J_s.
The reluctivity nu(B) is the problem -- it makes the LHS operator
time-varying.

Fixed-Point split with a *virtual* constant reluctivity nu^FP:
    curl( nu^FP curl A^(k) ) + sigma dA^(k)/dt
      = J_s
      - curl( H(B^(k)) - nu^FP B^(k) )                  (*)

  - Left-hand side is now LINEAR with constant nu^FP, so all the
    standard CLN basis functions, R matrix, L matrix can be derived
    ONCE and re-used.
  - All nonlinearity is in the RHS, expanded as a controlled source
    term:
       zeta^(k)_{2n} = (1/R_{2n}) * integral over Omega_C of
         { j_{2n} . (1/sigma) curl[H(B^(k)) - nu^FP B^(k)] dV }

  - The fixed-point iteration loop is:
       (1) solve the linear CLN system
              (R + d/dt L) * i^(k) = v - zeta^(k)
       (2) reconstruct B^(k+1) from i^(k) via the (linear) basis
       (3) update zeta^(k+1) by integrating (H(B^(k+1)) - nu^FP B^(k+1))
           against the test functions j_{2n}
       (4) repeat until ||zeta^(k+1) - zeta^(k)|| < tol

The mu^FP value is a tuning parameter -- the paper demonstrates
convergence for mu^FP_r = 2, 5, 8 on the same steel-sheet benchmark.

## What this avoids

Compared to the parameterized-basis nonlinear CLN tried in
06_非線形 (Gen 1-4), FP-CLN:
  - keeps the basis functions CONSTANT (one factorisation)
  - keeps R and L matrices CONSTANT (Cauer circuit topology fixed)
  - works the same way in time-domain and harmonic-balance
  - handles HYSTERESIS directly via the H(B) constitutive map
    (no additional "play model" / "Preisach" wrapping)
  - is implemented as a small MATLAB OO class (CLN.m) in only
    ~90 lines

## Two production drivers (the canonical "use FP-CLN" pattern)

### Harmonic Balance: 16_FP_CLN/HaromonicBalance/HB_t_B_x.m

Three nCASE switches inside one driver:
  nCASE = 1: linear material  (mur=100, BH = @(B) B/mu)
  nCASE = 2: nonlinear sat    (BH from B = 2 tanh(H/10e3) + mu0 H)
  nCASE = 3: HYSTERESIS       (BH model with hysteresis, this is
                                the headline case in the CEFC paper)
For each nCASE the driver sweeps mu^FP_r over {2, 5, 8} and N_f over
{1, 3, 5} harmonics, demonstrating the result is INDEPENDENT of
mu^FP within the convergence tolerance.

Output: animated B(x, t) plots saved to
  HB_t_x_mur_FP_{linear, nonlinear, hysteresis}/mur_FP=*.png + .mat
(also exported as PDF figures in the Digests folder).

### Time-Domain: 16_FP_CLN/TimeDomain/nonlin_timedomain_*

Multiple problem geometries, each with its own subfolder:
  - 鋼板問題_解析解\\        (steel sheet, analytical 1D
                              -- the CEFC benchmark)
  - 鋼板問題_NGSolve\\       (steel sheet, NGSolve FEM)
  - 円筒_磁場面に平行問題_解析解\\   (cylinder, B parallel to surface)
  - 円筒_磁場面に垂直問題_解析解\\   (cylinder, B perpendicular)
  - TBI_飛田モデル\\         (Tobita-Brand-Iron model, NGSolve)
  - TBI_飛田モデル_ギャップ無\\      (gap-less variant)

Each folder has:
  - `CLN.m`               -- the OO class (same 90-line definition,
                              with t_domain() and f_domain() methods)
  - `E0.m`                -- the closed-form per-mu basis (R, L,
                              A, B, H, J polynomial library)
  - `lglnodes.m`          -- Legendre-Gauss-Lobatto quadrature
                              nodes for the 1D radial / through-
                              thickness integrations
  - `nonlin_timedomain_current_input_CLN.m`        -- main driver
  - `nonlin_timedomain_current_input_CLN_restart.m` -- restart from
                                                       saved state

## The CLN MATLAB class (canonical 90 lines)

The OO class is identical across HaromonicBalance and TimeDomain
folders.  See FP_CLN_CODE topic for the full listing.  Three methods:
  - calc_LR_mat(cln)        Build diag(R) and tridiag(L) matrices.
  - f_domain(cln,f,V)       Solve [R + j 2*pi*f L]*I = V at each f.
                            Reconstruct JJ, AA, HH, BB on the radial
                            quadrature grid.
  - t_domain(cln,tspan,V,I0)   Solve d/dt(L I) + R I = V(t) using
                                MATLAB `ode23`.

## NGSolve driver (the FEM-coupled production path)

`TimeDomain/TBI_飛田モデル/CLN/NGSolve_CLN.m` -- couples FP-CLN to
NGSolve via the .mat exchange:
  - NGSolve generates 2D mesh (rectangular iron sheet, periodic BC)
  - exports BxI, ByI, VHx, VHy, Rmat, mag, nonmag, mu_FP
    as `NGSolve_CLN_TYPE=0_mur_FP=*.mat`
  - MATLAB FP-CLN driver loads .mat, runs Newton or AD-based
    fixed-point iteration (see `nonlin_autodiff_prototype.m` for
    the autodiff variant using `myAD`)
  - Output: B(x,y) snapshots overlaid on NGSolve mesh

This is the closest thing to a "production deployed" FP-CLN
pipeline -- though it still requires MATLAB on the FEM side, and
has NOT yet been ported into a radia.* Python module.

## Open questions (16_FP_CLN/疑問点)

  - 2023_09_19_反復計算.docx     -- "iteration computation"
  - 2導体モデルに単純化.docx     -- "simplify to 2-conductor model"

(These two docx files appear to track unresolved convergence /
formulation questions; not opened in this scan.)

## Companion theory references (see KEY_FILES topic for full list)

`FixedPoint法.pdf`, `hys_lamin_equations2.pdf`, `Section6.pdf`,
`プレイモデル.pdf` (play model), and Word docs `A法のMaxwell方程式.docx`,
`積層鋼板.docx` (laminated steel sheet) in 16_FP_CLN root.  HB
theory in `HaromonicBalance/`: `BH.pdf`, `Fixed_Point法.pdf`,
`2023_11_27_HBの行列.docx`, `定式化.docx`.

## Citation (cite this for FP-CLN)

```bibtex
@inproceedings{sugahara_fpcln_cefc2024,
  author = {Sugahara, Kengo and Tobita, Miwa and Matsuo, Tetsuji
            and Takahashi, Yasuhito},
  title  = {Fixed-Point {Cauer} Ladder Network Method for
            Eddy-current Problems with Hysteresis},
  booktitle = {IEEE CEFC 2024},
  year   = {2024},
  note   = {JSPS KAKENHI JP20K04443, JP23K03826;
            Power Academy Research Grant}
}
```
"""

FP_CLN_CODE = """
# FP-CLN code excerpts (16_FP_CLN)

## The CLN MATLAB class (canonical, ~90 lines)

Production OO class shared across HaromonicBalance, TimeDomain
folders.  Found verbatim at HaromonicBalance/CLN.m,
TimeDomain/円筒_*/CLN.m, TimeDomain/鋼板問題_解析解/CLN.m, with a
variant at TimeDomain/TBI_飛田モデル/CLN/NGSolve_CLN.m.

```matlab
classdef CLN
    properties
        nStage; R; L; J; A; H; B;
        Rmat; Lmat; f; Z; II; JJ; AA; HH; BB; t;
    end
    methods
        function cln = CLN(R, L, A, B, H, J, nStage)
            cln.R = R(1:nStage);  cln.L = L(1:nStage);
            cln.J = J(1:nStage);  cln.A = A(1:nStage);
            cln.H = H(1:nStage);  cln.B = B(1:nStage);
            cln.nStage = nStage;
        end
        function cln = calc_LR_mat(cln)
            cln.Rmat = diag([cln.R]);
            cln.Lmat = zeros(size(cln.Rmat));
            for n = [1:size(cln.Lmat,1)-1]
                cln.Lmat(n:n+1, n:n+1) = cln.Lmat(n:n+1, n:n+1) ...
                  + [[cln.L(n), -cln.L(n)]; [-cln.L(n), cln.L(n)]];
            end
            cln.Lmat(end,end) = cln.Lmat(end,end) + cln.L(end);
        end
        function cln = f_domain(cln, f, V)
            % Solve [R + j 2*pi*f L] I = V at each f, reconstruct
            % JJ{nf}, AA{nf}, HH{nf}, BB{nf} on the quadrature grid
            % ...
        end
        function cln = t_domain(cln, tspan, V, I0)
            % Time-domain via ode23:
            dIdt = @(t,II) [cln.Lmat] \\ ...
                   [[V(t); zeros(cln.nStage-1,1)] - cln.Rmat*II];
            [cln.t, cln.II] = ode23(dIdt, tspan, I0);
        end
    end
end
```

## Harmonic-Balance FP-CLN driver (HB_t_B_x.m, hysteresis case)

```matlab
mu0 = 4*pi*1e-7;  mur_FPs = [2,5,8];  Nfs = [1,3,5];  nStage = 4;
for mur_FP = mur_FPs
    mu_FP = mur_FP * mu0;  sig = 1.0e7;  d = 5e-3;
    [x,w,P] = lglnodes(100-1);
    x = d/2*transpose(x);  w = d/2*transpose(w);
    f = 50;  omega = 2*pi*f;

    cln = E0(x, d, sig, mu_FP);                    % LINEAR basis
    cln = CLN(cln.R, cln.L, cln.A, cln.B, cln.H, cln.J, nStage);
    cln = calc_LR_mat(cln);
    t = transpose(linspace(0, 1/f, 100));
    V0 = zeros(Nfs(end), nStage);   V0(1,1) = 2.5;

    for niter = [0:1000]
        % (1) linear solve per harmonic (zeta = V1 from prev iter)
        for nf = Nfs
            V(nf,:)  = V0(nf,:) + V1(nf,:);
            II(nf,:) = transpose( ...
              [cln.Rmat + j*omega*nf*cln.Lmat] \\ ...
              [transpose(V(nf,:))]);
        end
        % (2) reconstruct B(x,t) from harmonics II
        % (3) H_FP(t,x) = sign(B)*BH(|B|) - B/mu_FP    <-- nonlinear
        % (4) project H_FP onto basis cln.H{n} -> hn(t,n)
        % (5) Cauer-sum hn to get jn(t,n) (zeta in time)
        % (6) FFT jn(t,n) over t -> j_re(nf,n), j_im(nf,n)
        % (7) V1 = -complex(j_re, j_im); check norm(V - Vp)/norm(V)
        if err < 1e-3; break; end
    end
end
```

Key features:
  - basis (cln.B, cln.H, cln.R, cln.L) computed ONCE per mu_FP,
    OUTSIDE the iteration;
  - iteration only updates zeta + harmonic currents II(nf,:);
  - nonlinear physics is ONLY in `H_FP = sign(B)*BH(|B|) - B/mu_FP`;
  - hysteresis handled by BH function with its own internal state
    (play model / Preisach -- external to FP-CLN).

## Time-Domain FP-CLN with autodiff Newton (nonlin_autodiff_prototype.m)

Forward-mode AD on the FP residual via `myAD`:

```matlab
f_AD = @(II) func_AD(II, BxI, ByI, mu_, nu_, dnudB_, mu0, mu_FP, ...
                     VHx, VHy, V0, Rmat, mag, nonmag);
for niter = [1:1000]
    II_AD = myAD(II);
    f_    = f_AD(II_AD);
    f     = getvalues(f_);   df = getderivs(f_);
    II    = II - df \\ f;     % Newton step
    if norm(f) < 1e-6; break; end
end
```

Much faster than finite-difference Jacobian; `func_AD` is
identical to `func` except intermediate variables are AD numbers.
"""

LAB_LESSONS = """
# Hard-won lessons across the advanced CLN regimes

1) **The MOR-friendly linearisation wins, even when "wrong".**
   Five+ years of attempts (06_非線形, Gen 1 through Gen 4) to make
   the CLN basis itself nonlinear all ran into the same trap: time-
   varying R, L matrices defeat the MOR speedup.  The FP-CLN
   solution (16_FP_CLN, CEFC 2024) keeps the LHS LINEAR with a
   *deliberately wrong* mu^FP and pushes the nonlinearity to a RHS
   controlled source.  The fixed-point iteration converges for a
   wide range of mu^FP values (validated mu^FP_r = 2, 5, 8) and the
   final solution does not depend on mu^FP.

2) **Symbolic prototypes scale to ~N = 10 stages; beyond that, VPA.**
   The HF symbolic recursion (05_高周波) produces basis coefficients
   in the 10^30 range by order r^20.  IEEE double precision loses
   accuracy; use MATLAB `vpa` (variable-precision arithmetic) for
   any HF expansion beyond N = 8 or so.  For ground-truth checks
   against Bessel functions, use `bessel_taylor.m` with `taylor(...,
   k, 2*N)` to align orders.

3) **Boundary trimming is NOT optional.**
   In the recursive basis construction (both MQS and Helmholtz), the
   choice of which surface quantity to set to zero (H_tangential or
   E_normal) at the conductor boundary changes the L_n tail.  Mixing
   them in different stages breaks orthogonality.  CLN.m (05_高周波)
   exposes this as the `switch 'H'` / `case 'E'` selector -- pick
   ONE and stay with it.

4) **The Cauer L matrix is TRIDIAGONAL, the R matrix is DIAGONAL.**
   This is consistently enforced in `calc_LR_mat` across all
   FP-CLN folders:
     R = diag(R_n)
     L = tridiag(...,[L_n -L_n; -L_n L_n]...) + diag(L_end at corner)
   When porting to a new language, do NOT collapse L into a full
   matrix -- the tridiagonal structure is the entire point of "this
   is a LADDER".

5) **Hysteresis is handled by the constitutive map, not by CLN.**
   FP-CLN treats H(B) as a black-box function with internal state.
   For play-model hysteresis the state advances each time the BH
   function is called.  CLN itself stores nothing hysteretic --
   only the linear circuit history (currents I_n at the previous
   time step).  This separation is the source of the lab's "FP-CLN
   handles hysteresis" claim.

6) **Restart your time-domain CLN integrator at every nonlin
   convergence failure.**
   The `_restart.m` files in 16_FP_CLN/TimeDomain show the pattern:
   if the FP iteration fails to converge inside a time step, drop
   back to the last converged state, halve the time step, and
   retry.  Hard-coded `if niter == 2000; error('failed'); end` is
   used as a safety valve.  Do NOT silently advance with a non-
   converged iterate.

7) **The 1D radial benchmark (round wire) is the universal sanity
   check.** Every CLN extension (HF, nonlinear, FP) ships with a
   round-wire driver and compares against the analytical Bessel
   solution at the end.  When porting to a new geometry, ALWAYS
   include a round-wire regression test first.

8) **For HF, switch sigma -> epsilon in the J update.**  This is
   the *only* structural change needed to convert the MQS recursion
   to the Helmholtz recursion.  See `check_CLN_high_freq_lossless.m`:
     J(n+1) = simplify( J(n) - ep0 * A(n) / L(n) );
   vs the MQS form:
     J(n+1) = simplify( J(n) - sig * A(n) / L(n) );

9) **Save Cauer matrices to .mat at each (mu_FP, geometry) pair.**
   The lab keeps `NGSolve_CLN_TYPE=*_mur_FP=*.mat` as the
   interchange format between MATLAB FP-CLN and NGSolve FEM
   meshing.  These files contain (BxI, ByI, VHx, VHy, Rmat,
   mag, nonmag, mu_FP) and the MATLAB driver just `load`s them.
   Same pattern is recommended when porting to Python: pickle /
   .npz the precomputed Cauer matrices, separate from the
   nonlinear iteration logic.

10) **Cauer matrix for steel sheet has a recursive pattern.**
    For an iron sheet of thickness d, the Cauer recursion (from
    16_FP_CLN/TimeDomain/鋼板問題_解析解/E0.m) produces:
      R_2n   = 4 (4 n - 1) / (sigma d)
      L_2n+1 = mu^FP d / (4 n + 1)
    This is the closed form referenced in Shindo-Miyazaki-Matsuo
    2015 (Legendre expansion for a magnetic sheet) and matches
    Eq. (anal) of the CEFC 2024 paper.

11) **A circuit-coupled CLN (07_電気力学) needs both Cauer matrices
    AND per-mode force coefficients.**  The Shindo-electromagnet
    pipeline extracts:
      [R, L]       -- standard Cauer matrices (5 stages)
      dF(0)..dF(4) -- per-stage z-force from Maxwell stress
      FDC          -- static force
    and feeds them as a small RL block + per-mode force-input
    matrix into the mechanical 1-DOF equation.  Without dF(n) the
    mechanical bandwidth cannot be computed; without [R, L] the
    drive impedance cannot.  Both are required for any controller
    design.
"""

KEY_FILES = """
# Top files referenced (advanced CLN corpus)

Base path: public-safe curated corpus

## Theory references (PDFs)

  16_FP_CLN\\Digests\\CEFC_digest.{pdf,tex} + references.bib
      -- Sugahara-Tobita-Matsuo-Takahashi, CEFC 2024.
         THE canonical FP-CLN reference.
  16_FP_CLN\\{FixedPoint法.pdf, hys_lamin_equations2.pdf,
              プレイモデル.pdf}
      -- fixed-point + hysteresis-laminate + play-model theory.
  16_FP_CLN\\HaromonicBalance\\Efficient Circuit Representation of
              Eddy-Current Field.pdf, Fixed_Point法.pdf
      -- HB theory + FP-in-HB.
  16_FP_CLN\\TimeDomain\\a130372_IEEJ-20220307X01601.pdf
      -- IEEJ 2022 conference paper (TD FP-CLN).
  06_非線形\\2020_08_06_京大資料_非線形\\2020-07-31-nolin-cln.pdf,
              2020-09-15-slide.pdf
      -- Kyoto Univ slides on early nonlinear CLN.
  06_非線形\\2017_12_12_非線形の案_01\\14-02.pdf
      -- 2017 conference digest (TBD).

## Production code

  16_FP_CLN\\HaromonicBalance\\CLN.m            -- 90-line OO class.
  16_FP_CLN\\HaromonicBalance\\HB_t_B_x.m       -- HB driver, 3
       cases (linear/nonlinear/hysteresis) x 3 mu^FP x 3 harmonics.
  16_FP_CLN\\HaromonicBalance\\E0.m             -- per-mu basis lib.
  16_FP_CLN\\TimeDomain\\nonlin_static_CLN.m    -- static FP driver.
  16_FP_CLN\\TimeDomain\\nonlin_autodiff_prototype.m
                                                  -- AD-Newton FP.
  16_FP_CLN\\TimeDomain\\鋼板問題_解析解\\非線形解析\\
       nonlin_timedomain_current_input_CLN{,_restart}.m
                                                  -- CEFC 2024
                                                     benchmark.
  16_FP_CLN\\TimeDomain\\TBI_飛田モデル\\CLN\\
       NGSolve_CLN.m, nonlin_timedomain_CLN_NGSolve.m
                                                  -- FEM-coupled
                                                     FP-CLN.

## Symbolic HF prototypes

  05_高周波\\GenKameariBasis.m, bessel_taylor.m, CLN.m, CLN2.m
       -- symbolic Cauer recursion + Bessel ground-truth.
  2020_12_07_高周波CLNの練習\\check_CLN_high_freq_lossless.m,
       check_CLN_high_freq_infinity.m, GramShmidt.m
       -- coaxial-ring Helmholtz Cauer + Gram-Schmidt basis.

## Electrodynamic-circuit + nonlinear lineage + Bloch

  07_電気力学\\shindo.edp              -- 437-line FreeFEM++ driver
       (11-stage CLN + Maxwell-stress dF(n) + frequency sweep).
  06_非線形\\2017_12_12_非線形の案_01\\non_linear_CLN.m, genBH.m
                                        -- Gen 1: tanh-saturated BH.
  06_非線形\\2020_08_09_非線形CLNの練習\\case{0,1,2}.m
                                        -- COMSOL LiveLink test.
  06_非線形\\2023_09_26_飛田モデル_jω法\\compare_results_*.m,
       Gen_FEMM_model.m, read_mesh.m   -- Tobita's jω method.
  14_ブロッホ方程式\\2019_08_11_ブロッホ方程式.docx
                                        -- single design doc, 2019.
"""

CITATIONS = """
# Citations -- advanced CLN

These are the journal / conference references cited inside
16_FP_CLN/Digests/references.bib (the CEFC 2024 paper's bib),
filtered to the advanced-physics extensions.  For the basic linear
CLN canon (Kameari 2018 group), see `cln_knowledge.py` -- those
papers are repeated below ONLY when they are also cited from the
advanced literature.

## FP-CLN (the headline paper for THIS module)

```bibtex
@inproceedings{sugahara_fpcln_cefc2024,
  author={Sugahara,K. and Tobita,M. and Matsuo,T. and Takahashi,Y.},
  title={Fixed-Point Cauer Ladder Network Method for Eddy-current
         Problems with Hysteresis},
  booktitle={IEEE CEFC 2024}, year={2024},
  note={JSPS KAKENHI JP20K04443 and Power Academy Research Grant}}
```

## Earlier nonlinear-CLN attempts + fixed-point method

```bibtex
@article{nonlin_cln_eskandari2020,
  author={Eskandari,H. and Matsuo,T.},
  title={Cauer Ladder Network Representation of a Nonlinear
         Eddy-Current Field Using a First-Order Approximation},
  journal={IEEE Trans. Magn.}, volume={56}, number={2}, year={2020},
  doi={10.1109/TMAG.2019.2946354}}

@article{nonlin_cln_tobita2022,
  author={Tobita,M. and Matsuo,T.},
  title={Nonlinear Model Order Reduction of Induction Motors Using
         Parameterized Cauer Ladder Network Method},
  journal={IEEE Trans. Magn.}, volume={58}, number={9}, year={2022},
  doi={10.1109/TMAG.2022.3171743}}

@article{dlala_fixedpoint_2007,
  author={Dlala,E. and Belahcen,A. and Arkkio,A.},
  title={Locally Convergent Fixed-Point Method for Solving
         Time-Stepping Nonlinear Field Problems},
  journal={IEEE Trans. Magn.}, volume={43}, number={11}, year={2007},
  pages={3969-3975}, doi={10.1109/TMAG.2007.904819}}

@article{dlala_fixedpoint_fast_2008,
  author={Dlala,E. and Belahcen,A. and Arkkio,A.},
  title={A Fast Fixed-Point Method for Solving Magnetic Field
         Problems in Media of Hysteresis},
  journal={IEEE Trans. Magn.}, volume={44}, number={6}, year={2008},
  pages={1214-1217}, doi={10.1109/TMAG.2007.916673}}
```

## CLN with translational mover (07_電気力学 generalisation)

```bibtex
@article{shindo_cln_mover_2020,
  author={Shindo,Y. and Yamamoto,R. and Sugahara,K. and Matsuo,T.
          and Kameari,A.},
  title={Equivalent Circuit in Cauer Form for Eddy Current Field
         Including a Translational Mover},
  journal={IEEE Trans. Magn.}, volume={56}, number={12},
  year={2020}, doi={10.1109/TMAG.2020.3023234}}

@article{sabariego_rlladder_mover_2022,
  author={Sabariego,R.V. and Vanbroekhoven,B. and Gyselinck,J.
          and Kuo-Peng,P.},
  title={RL-Ladder Circuit Models for Eddy-Current Problems With
         Translational Movement},
  journal={IEEE Trans. Magn.}, volume={58}, number={9},
  year={2022}, doi={10.1109/TMAG.2022.3170996}}

@inproceedings{shindo_dynamical_ldia2017,
  author={Shindo,Y. and Kameari,A. and Sugahara,K. and Matsuo,T.},
  title={Dynamical model of an electromagnet by Cauer ladder
         network representation of eddy-current fields},
  booktitle={LDIA 2017}, year={2017},
  doi={10.23919/LDIA.2017.8097231}}
```

## Cauer-circuit for laminated steel sheet (1D analytic benchmark)

```bibtex
@article{shindo_legendre_2015,
  author={Shindo,Y. and Miyazaki,T. and Matsuo,T.},
  title={Cauer circuit representation of the homogenized eddy-
         current field based on the Legendre expansion for a
         magnetic sheet},
  journal={IEEE Trans. Magn.}, volume={52}, number={3}, year={2015}}
```

## MOR alternatives + Cauer-PGD (comparison set, not CLN itself)

```bibtex
@article{pierquin_pod_2015, author={Pierquin,A. and Henneron,T.
  and Clenet,S. and Brisset,S.},
  title={Model-Order Reduction of Magnetoquasi-Static Problems
         Based on POD and Arnoldi-Based Krylov Methods},
  journal={IEEE Trans. Magn.}, volume={51}, number={3}, year={2015}}

@article{henneron_pod_dei_2015, author={Henneron,T. and Clenet,S.},
  title={Model-Order Reduction of Multiple-Input Non-Linear
         Systems Based on POD and DEI Methods},
  journal={IEEE Trans. Magn.}, volume={51}, number={3}, year={2015}}

@article{sato_pade_2016, author={Sato,Y. and Igarashi,H.},
  title={Generation of Equivalent Circuit From Finite-Element
         Model Using Model Order Reduction},
  journal={IEEE Trans. Magn.}, volume={52}, number={3}, year={2016},
  doi={10.1109/TMAG.2015.2476370}}

@article{hollaus_msfem_2018,
  author={Hollaus,K. and Schoberl,J.},
  title={Some 2-D Multiscale Finite-Element Formulations for the
         Eddy Current Problem in Iron Laminates},
  journal={IEEE Trans. Magn.}, volume={54}, number={4}, year={2018},
  doi={10.1109/TMAG.2017.2777395}}

@article{koster_pgd_cln_2021,
  author={Koster,N. and Konig,O. and Biro,O.},
  title={Proper Generalized Decomposition With Cauer Ladder
         Network Applied to Eddy Current Problems},
  journal={IEEE Trans. Magn.}, volume={57}, number={6}, year={2021},
  doi={10.1109/TMAG.2021.3059800}}
```

## Bloch equation (NMR-MRI) -- no Sugahara-Lab CLN-NMR paper

The 2019 design doc at public-safe curated corpus
sketches the analogy; no publication.  See `radia_mcp.nmr_mri` for
the actually-implemented NMR-MRI knowledge layer.
"""


_TOPICS = {
    "overview": OVERVIEW,
    "high_frequency": HIGH_FREQUENCY,
    "high_frequency_code": HIGH_FREQUENCY_CODE,
    "nonlinear": NONLINEAR,
    "nonlinear_code": NONLINEAR_CODE,
    "electrodynamic_circuit": ELECTRODYNAMIC_CIRCUIT,
    "bloch_equation": BLOCH_EQUATION,
    "fp_cln": FP_CLN,
    "fp_cln_code": FP_CLN_CODE,
    "lab_lessons": LAB_LESSONS,
    "key_files": KEY_FILES,
    "citations": CITATIONS,
    "all": "",
}
_TOPICS["all"] = "\n\n".join(v for k, v in _TOPICS.items() if k != "all")


def get_cln_advanced_documentation(topic: str = "all") -> str:
    """Return advanced-CLN documentation for the requested topic.

    Topics:
      overview               - which advanced extension fits which regime
      high_frequency         - HF / Helmholtz CLN (05 + 2020_12_07)
      high_frequency_code    - HF code excerpts (Helmholtz recursion)
      nonlinear              - nonlinear CLN lineage 2017-2023 (06)
      nonlinear_code         - nonlinear CLN code excerpts
      electrodynamic_circuit - Shindo electromagnet, FreeFEM++ (07)
      bloch_equation         - NMR/MRI Bloch tie-in (14, one doc only)
      fp_cln                 - Fixed-Point CLN, lab production (16)
      fp_cln_code            - FP-CLN MATLAB class + drivers
      lab_lessons            - 11 hard-won lessons across all regimes
      key_files              - top files referenced in the corpus
      citations              - BibTeX citations for the advanced layer
      all                    - everything concatenated
    """
    topic = topic.lower().strip()
    if topic in _TOPICS:
        return _TOPICS[topic]
    return ("Unknown topic '" + topic + "'. Available: "
            + ", ".join(sorted(_TOPICS.keys())))
