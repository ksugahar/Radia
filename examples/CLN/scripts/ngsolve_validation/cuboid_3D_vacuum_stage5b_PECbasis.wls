(* Stage 5b: Correct H(curl)-conforming PEC TE_z basis for cuboid.

   This corrects Stage 5a's incorrect formulation. The new basis:
     phi^TE_{p,q,r} = cos(p pi (x+a/2)/a) cos(q pi (y+b/2)/b) sin(r pi (z+c/2)/c)
     A^TE = curl(z_hat phi^TE)  with proper PEC tangential vanishing on all 6 faces.

   Parity (uniform B_z + D_2h): p even >= 0, q even >= 0, r odd >= 1.

   Source coupling discovery:
     b^{(0,q,r)} = -V/(r pi)   (q-independent, q even >= 2)
     b^{(p,0,r)} = -V/(r pi)   (p-independent, p even >= 2)
     b^{(p>=2,q>=2,r)} = 0     (cosine vanishes over full period)

   Closed-form alpha_2 from infinite series of edge modes.
*)
$HistoryLength = 0;
prec = 60;

Print["=== Stage 5b: Correct PEC TE_z basis (closed PEC, all closed form) ==="];
Print[""];

(* Constants *)
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;
ReffSq = 5 (aN^2 + bN^2)/24;
alphaLFRef = -mu0N sigmaN (aN^2 + bN^2)/48;

(* Sphere reference *)
alphaSphere[R_, k_] := Module[{x = mu0N sigmaN R^2, ck},
  ck = Switch[k, 1, -1/10, 2, +1/105, 3, -1/1050, 4, +1/10395, 5, -1/135135,
    _, (-1)^k * 9 * SetPrecision[Zeta[2 k + 2], prec]/PiN^(2 k + 2)];
  ck * x^k];

Print["Sphere reference (R_eff^2 = 5(a^2+b^2)/24):"];
Print["  alpha_1^sph = ", N[alphaSphere[Sqrt[ReffSq], 1], 16]];
Print["  alpha_2^sph = ", N[alphaSphere[Sqrt[ReffSq], 2], 16]];
Print["  alpha_3^sph = ", N[alphaSphere[Sqrt[ReffSq], 3], 16]];
Print[""];

(* === Closed-form alpha_2 from infinite TE_z mode series ===

   For (0, q, r) modes:
     A_x = -(q pi/b) sin(q pi v) sin(r pi w)   (cos(0)=1 in x)
     A_y = 0
     |A|^2 = (q pi/b)^2 sin^2(q pi v) sin^2(r pi w)
     M_ii = (V/4) (q pi/b)^2
     lambda_(0,q,r) = pi^2((q/b)^2 + (r/c)^2)
     K_ii = lambda * M_ii = (V/4)(q pi/b)^2 pi^2((q/b)^2 + (r/c)^2)
     b^{(0,q,r)} = -V/(r pi)
     b^2/K = (V^2/(r^2 pi^2)) / [(V/4)(q pi/b)^2 pi^2((q/b)^2 + (r/c)^2)]
           = 4 V b^2 / [q^2 r^2 pi^6 ((q/b)^2 + (r/c)^2)]
           = 4 V b^4 c^2 / [pi^6 q^2 r^2 (q^2 c^2 + r^2 b^2)]

   Sum over q even >= 2, r odd >= 1.  Substitute q = 2 q', q' >= 1:
     S0qr = sum (q'>=1, r odd) b^4 c^2 / [pi^6 (2q')^2 r^2 ((2q')^2 c^2 + r^2 b^2)]
          = b^4 c^2 / (4 pi^6) sum (q'>=1, r odd) 1/[q'^2 r^2 (4 q'^2 c^2 + r^2 b^2)]
*)

(* alpha_2 from (0, q, r) family *)
NMaxQ = 60;  (* truncation *)
NMaxR = 121;  (* odd up to 121 *)

S0qr = Sum[1/(qp^2 r^2 (4 qp^2 cN^2 + r^2 bN^2)),
  {qp, 1, NMaxQ}, {r, 1, NMaxR, 2}];
contrib0qr = (mu0N sigmaN)^2/VN * 4 VN bN^4 cN^2/(4 PiN^6) * S0qr;

(* Symmetric (p, 0, r) family with a <-> b: *)
Sp0r = Sum[1/(pp^2 r^2 (4 pp^2 cN^2 + r^2 aN^2)),
  {pp, 1, NMaxQ}, {r, 1, NMaxR, 2}];
contribp0r = (mu0N sigmaN)^2/VN * 4 VN aN^4 cN^2/(4 PiN^6) * Sp0r;

alpha2PEC = contrib0qr + contribp0r;

Print["alpha_2 closed PEC TE_z (truncated NMaxQ=", NMaxQ, ", NMaxR=", NMaxR, "):"];
Print["  contribution from (0, q, r): ", N[contrib0qr, 16]];
Print["  contribution from (p, 0, r): ", N[contribp0r, 16]];
Print["  total alpha_2^PEC = ", N[alpha2PEC, 16]];
Print[""];

(* Convergence study *)
Print["=== Convergence with truncation ==="];
Do[Module[{S0, Sp0, alpha2, NQ, NR},
  {NQ, NR} = pair;
  S0 = Sum[1/(qp^2 r^2 (4 qp^2 cN^2 + r^2 bN^2)),
    {qp, 1, NQ}, {r, 1, 2 NR + 1, 2}];
  Sp0 = Sum[1/(pp^2 r^2 (4 pp^2 cN^2 + r^2 aN^2)),
    {pp, 1, NQ}, {r, 1, 2 NR + 1, 2}];
  alpha2 = (mu0N sigmaN)^2 * (bN^4 cN^2 S0 + aN^4 cN^2 Sp0)/PiN^6;
  Print["  NMaxQ=", NQ, ", NMaxR=", 2 NR + 1, ": alpha_2^PEC = ",
    N[alpha2, 16]]],
  {pair, {{10, 20}, {20, 40}, {40, 60}, {60, 90}, {80, 120}}}];
Print[""];

(* Compare to sphere *)
sphAlpha2 = alphaSphere[Sqrt[ReffSq], 2];
Print["Sphere reference alpha_2^sph(R_eff) = ", N[sphAlpha2, 16]];
Print["Ratio alpha_2^PEC / alpha_2^sph = ", N[alpha2PEC/sphAlpha2, 8]];
Print[""];
Print["Note: PEC closed cuboid CANNOT match vacuum sphere because"];
Print["it lacks the eddy current bleed-out into vacuum.  Need exterior"];
Print["multipole basis enrichment (Step 5c)."];
Print[""];

(* === Try also: alpha_2 from sin sin sin scalar Dirichlet basis (Step 1+2 baseline) ===
   For comparison: that was a different (incorrect for vector PEC) basis.
   Just show numerical value to compare basis effects. *)
Print["=== Comparison to incorrect Step 1+2 sin-basis ==="];
Print["  Step 1+2 (sin sin sin scalar Dirichlet): alpha_2 = 8.349e-11 (Pmax=11)"];
Print["  Step 5b (correct PEC TE_z, infinite basis): alpha_2 = ", N[alpha2PEC, 12]];
Print["  Ratio = ", N[8.349*10^-11 / alpha2PEC, 6]];
Print[""];

Print["=== Closed PEC alpha_n series (TE_z complete basis) ==="];
(* General closed form for alpha_n from (0, q, r) edge modes, n >= 2:
   alpha_n contribution from (0,q,r) = (-1)^(n-1) (mu0 sigma)^n / V * b_i^2 M_ii^(n-2) / K_ii^(n-1)
   b^2 = V^2/(r^2 pi^2)
   M_ii = (V/4)(q pi/b)^2
   K_ii = lambda M_ii = pi^2((q/b)^2+(r/c)^2) (V/4)(q pi/b)^2

   M_ii^(n-2) / K_ii^(n-1) = M_ii^(n-2)/[lambda^(n-1) M_ii^(n-1)] = 1/[lambda^(n-1) M_ii]
   = 1/[pi^2((q/b)^2+(r/c)^2)]^(n-1) / [(V/4)(q pi/b)^2]
   = 4/[V (q pi/b)^2 pi^(2n-2) ((q/b)^2+(r/c)^2)^(n-1)]

   So contribution = (-1)^(n-1)(mu0 sigma)^n V/r^2 pi^2 * 4/[V (q pi/b)^2 pi^(2n-2) ((q/b)^2+(r/c)^2)^(n-1)]
                   = (-1)^(n-1) 4 (mu0 sigma)^n b^2 / [q^2 r^2 pi^(2n+2) ((q/b)^2+(r/c)^2)^(n-1)]
*)

(* Compute alpha_3, alpha_4, alpha_5 *)
computeAlphaN[n_, NMaxQ_, NMaxR_] := Module[{S0, Sp0},
  S0 = Sum[
    bN^(2 n - 2)/(qp^2 r^2 ((2 qp)^(2 n - 2) cN^(2 n - 2) +
        (2 qp)^(2 n - 4) (2 qp)^2 cN^(2 n - 4) bN^2 r^2 *
          0)) (* simpler: ((q/b)^2+(r/c)^2)^(n-1) = sum binomial *)
    , {qp, 1, NMaxQ}, {r, 1, 2 NMaxR + 1, 2}];
  (* Actually, just use direct computation: *)
  S0 = Sum[
    1/(qp^2 r^2 ((2 qp/bN)^2 + (r/cN)^2)^(n - 1)),
    {qp, 1, NMaxQ}, {r, 1, 2 NMaxR + 1, 2}];
  Sp0 = Sum[
    1/(pp^2 r^2 ((2 pp/aN)^2 + (r/cN)^2)^(n - 1)),
    {pp, 1, NMaxQ}, {r, 1, 2 NMaxR + 1, 2}];
  (-1)^(n - 1) * 4 (mu0N sigmaN)^n / (4 PiN^(2 n + 2)) * (bN^2 S0 + aN^2 Sp0)];

Do[Module[{val},
  val = computeAlphaN[k, 60, 60];
  Print["alpha_", k, "^PEC = ", N[val, 16],
    "  (sphere ref: ", N[alphaSphere[Sqrt[ReffSq], k], 12],
    ", ratio: ", N[val/alphaSphere[Sqrt[ReffSq], k], 8], ")"]],
  {k, 2, 5}];
Print[""];

Print["=== Stage 5b summary ==="];
Print["* Correct PEC TE_z basis derived (cos cos sin)"];
Print["* Source coupling: only edge modes (p=0 or q=0)"];
Print["* alpha_2^PEC (full TE basis) = ", N[alpha2PEC, 8]];
Print["  Closed PEC alpha_2 captured as DOUBLE infinite series"];
Print[""];
Print["Next: Step 5c -- enrich with exterior multipole for vacuum BC"];
