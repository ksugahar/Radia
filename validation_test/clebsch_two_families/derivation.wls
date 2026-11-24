#!/usr/bin/env wolframscript
(* ::Package:: *)
(*
  derivation.wls -- symbolic basis for two_families.py, built on the repository
  FEEC library (radia_mcp/mathematica/basis_functions/hcurl.wls).

  S1  analytic field: det DX = 1, inverse, B = dX/ds = grad u x grad v, div B = 0,
      labels constant along B; z = pi map is the identity; a whole field line
      lies in the element face x = c.
  S2  FEEC traces on the hex Nedelec / RT pair (HCurlHex[p], HDivHexT[p]):
      curl(HCurl_p) in HDiv_p, div curl = 0; zero tangential trace forces zero
      normal curl, but leaves A_n and tangential curl free (they jump).
      Degree of A_h.dl along a chord in a z = const plane (Gauss-6 exactness).
  S3  discrete Clebsch: u grad v in H(curl), grad u x grad v in HDivHexT[2p+1],
      not in HDivHexT[2p]; traces determined by face traces of u, v.
  S4  flux-form transport (Liouville): d/dz (B_z det J) = det J div B.
  S5  error orders: surface reconstruction leakage and geometry (series in the
      cell size on the analytic surfaces), polygon flux error per curve
      singularity (smooth, kink, fold, jump), and the fold produced by a
      tangential trace jump at a grazing crossing.
  Every check prints True/False or an explicit coefficient; nothing is fitted.
*)
here = DirectoryName[$InputFileName];
lib = FileNameJoin[{here, "..", "..", "packages", "radia-mcp", "src", "radia_mcp", "mathematica", "basis_functions"}];
SetDirectory[lib]; Block[{Print}, Get["hcurl.wls"]]; ResetDirectory[];
gr[f_] := D[f, {{x, y, z}}]; cu[F_] := Curl[F, {x, y, z}]; dvg[F_] := Div[F, {x, y, z}];
SeedRandom[20261120];
rk[fs_, pts_] := MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ fs)];
P[s__] := Print[s];

(* ------------------------------------------------------------------ S1 *)
P["== S1 analytic field (a, b symbolic)"];
Xm = {u + a Sin[v] Sin[s], v + b Sin[u + a Sin[v] Sin[s]] Sin[2 s], s};
P["det D(u,v,s)X == 1: ", Simplify[Det[D[Xm, {{u, v, s}}]]] === 1];
vI = y - b Sin[x] Sin[2 z]; uI = x - a Sin[vI] Sin[z];
P["inverse: X(uI, vI, z) == (x, y, z): ", Simplify[(Xm /. {u -> uI, v -> vI, s -> z}) - {x, y, z}] === {0, 0, 0}];
Bc = {a Sin[vI] Cos[z], b Cos[x] a Sin[vI] Cos[z] Sin[2 z] + 2 b Sin[x] Cos[2 z], 1};  (* code exact_b *)
P["code B == dX/ds at (uI, vI): ", Simplify[(D[Xm, s] /. {u -> uI, v -> vI, s -> z}) - Bc] === {0, 0, 0}];
P["code B == grad u x grad v: ", Simplify[Cross[gr[uI], gr[vI]] - Bc] === {0, 0, 0}];
P["div B == 0: ", Simplify[dvg[Bc]] === 0];
P["B.grad u == B.grad v == 0: ", Simplify[{Bc.gr[uI], Bc.gr[vI]}] === {0, 0}];
P["map at s = pi is the identity: ", Simplify[(Xm /. s -> Pi) - {u, v, Pi}] === {0, 0, 0},
  "  -> outlet patch at z = pi is the inlet square; use interior z for the patch-curvature test"];
P["field line (u = c, v = 0) lies in face x = c for all s: ", Simplify[(Xm[[1]] /. {u -> c, v -> 0})] === c,
  "; nearby lines: x - c = a sin v sin s (hug the face for |v| small)"];

(* ------------------------------------------------------------------ S2 *)
P["== S2 FEEC traces, hex Nedelec HCurlHex[p] / HDivHexT[p] (NGSolve order k = p)"];
Do[Module[{hc = HCurlHex[p], hd = HDivHexT[p], curls, ptsV, ptsF, trT, ker, nk, rows},
   curls = cu /@ hc;
   ptsV = Table[Thread[{x, y, z} -> RandomInteger[{1, 997}, 3]/1000], {Ceiling[Length[hd]/3] + 12}];
   ptsF = Table[Thread[{x, y, z} -> Prepend[RandomInteger[{1, 997}, 2]/1000, 0]], {(p + 2)^2 + 20}];
   (* tangential trace on x = 0 as sampled rows; kernel = zero-tangential-trace subspace *)
   trT = (Flatten[Table[{#[[2]], #[[3]]} /. pt, {pt, ptsF}]] & /@ hc);
   ker = NullSpace[Transpose[trT]]; nk = Length[ker];
   rows[f_] := MatrixRank[ker . (Flatten[Table[f[#] /. pt, {pt, ptsF}]] & /@ Range[Length[hc]])];
   P["p=", p, " dimHCurl=", Length[hc], " dimHDiv=", Length[hd],
    " | div curl = 0: ", Union[Expand /@ (dvg /@ curls)] === {0},
    " | curl(HCurl) in HDiv: ", rk[Join[hd, curls], ptsV] == rk[hd, ptsV],
    " | dim{tangential trace on x=0 vanishes} = ", nk,
    " | rank of (curl A)_x on x=0 over it = ", rows[(curls[[#, 1]]) &],
    " (0 => normal B continuous)",
    " | rank of A_x = ", rows[(hc[[#, 1]]) &],
    ", rank of (curl A)_{y,z} = ", rows[({curls[[#, 2]], curls[[#, 3]]}) &],
    " (>0 => A_n and tangential B jump across faces)"];
   Module[{t, d1, d2, x0, y0, z0, deg},
    deg = Max[Exponent[Expand[(#[[1]] d1 + #[[2]] d2) /. {x -> x0 + t d1, y -> y0 + t d2, z -> z0}], t] & /@ hc];
    P["   degree of A_h.dl along a chord in z = const: ", deg, " = 2p+1; Gauss-6 exact (<= 11): ", deg <= 11]]],
  {p, 2, 4}];

(* ------------------------------------------------------------------ S3 *)
P["== S3 discrete Clebsch with u, v in H1HexT[p]"];
Do[Module[{h1 = H1HexT[p], q = 2 p + 1, hd, hd2, uu, vv, A, B, pts},
   hd = HDivHexT[q]; hd2 = HDivHexT[q - 1];
   uu = RandomInteger[{-5, 5}, Length[h1]] . h1; vv = RandomInteger[{-5, 5}, Length[h1]] . h1;
   A = uu gr[vv]; B = Cross[gr[uu], gr[vv]];
   pts = Table[Thread[{x, y, z} -> RandomInteger[{1, 997}, 3]/1000], {Ceiling[Length[hd]/3] + 12}];
   P["p=", p, " curl(u grad v) == grad u x grad v: ", Expand[cu[A] - B] === {0, 0, 0},
    " | div == 0: ", Expand[dvg[B]] === 0,
    " | in HDivHexT[", q, "]: ", rk[Append[hd, B], pts] == rk[hd, pts],
    " | in HDivHexT[", q - 1, "]: ", rk[Append[hd2, B], pts] == rk[hd2, pts],
    " | B.n on x=0 from face traces: ",
    Expand[(B[[1]] /. x -> 0) - (D[uu /. x -> 0, y] D[vv /. x -> 0, z] - D[uu /. x -> 0, z] D[vv /. x -> 0, y])] === 0,
    " | A_t on x=0 from face traces: ",
    Expand[({A[[2]], A[[3]]} /. x -> 0) - (uu /. x -> 0) {D[vv /. x -> 0, y], D[vv /. x -> 0, z]}] === {0, 0}]],
  {p, 0, 1}];
P["   => grad u_h x grad v_h is H(div)-conforming and exactly div-free, but of degree 2p+1:",
  " a given B_h in HDiv_k is not generally one wedge of H1_k labels (nonlinear, degree doubling)."];

(* ------------------------------------------------------------------ S4 *)
P["== S4 flux-form transport"];
Module[{bx = Bx[x, y, z], by = By[x, y, z], bz = Bz[x, y, z], lhs},
  (* along dX/dz = (bx/bz, by/bz): d/dz(bz detJ) = detJ (D_z bz + bz div_xy F) *)
  lhs = (D[bz, z] + bx/bz D[bz, x] + by/bz D[bz, y]) + bz (D[bx/bz, x] + D[by/bz, y]);
  P["d/dz(B_z det J)/det J - div B == 0 (generic B): ", Simplify[lhs - (D[bx, x] + D[by, y] + D[bz, z])] === 0]];
P["   => B_z dx^dy is transported exactly by the exact flow of any div-free B; for piecewise B_h",
  " this needs only continuous B.n at faces (S2) and exact cell divergence."];

(* ------------------------------------------------------------------ S5 *)
P["== S5a surface reconstruction: linear across curves (dt), cubic Hermite along them (dz)"];
Module[{av = N[1/2, 80], bv = N[1/2, 80], u0 = N[3/10, 80], t0 = N[1/5, 80], z0 = N[7/10, 80], al = 1/3, be = 2/5, Xs, Bf, C, Dz, H, Sh, St, Sz, lk, gm, h00, h10, h01, h11, ser},
  Xs[t_, zz_] := {u0 + av Sin[t] Sin[zz], t + bv Sin[u0 + av Sin[t] Sin[zz]] Sin[2 zz], zz};
  Bf[{px_, py_, pz_}] := Module[{vv = py - bv Sin[px] Sin[2 pz]}, {av Sin[vv] Cos[pz], bv Cos[px] av Sin[vv] Cos[pz] Sin[2 pz] + 2 bv Sin[px] Cos[2 pz], 1}];
  h00 = 2 be^3 - 3 be^2 + 1; h10 = be^3 - 2 be^2 + be; h01 = -2 be^3 + 3 be^2; h11 = be^3 - be^2;
  H[t_] := h00 Xs[t, z0] + h10 zt ep Bf[Xs[t, z0]] + h01 Xs[t, z0 + zt ep] + h11 zt ep Bf[Xs[t, z0 + zt ep]];
  Sh = (1 - al) H[t0] + al H[t0 + ta ep];
  St = (H[t0 + ta ep] - H[t0])/(ta ep);
  (* Hermite derivative along z, written with the shape-function derivatives *)
  Sz = ((1 - al) (D[{2 bb^3 - 3 bb^2 + 1, bb^3 - 2 bb^2 + bb, -2 bb^3 + 3 bb^2, bb^3 - bb^2}, bb] /. bb -> be) .
        {Xs[t0, z0], zt ep Bf[Xs[t0, z0]], Xs[t0, z0 + zt ep], zt ep Bf[Xs[t0, z0 + zt ep]]} +
       al (D[{2 bb^3 - 3 bb^2 + 1, bb^3 - 2 bb^2 + bb, -2 bb^3 + 3 bb^2, bb^3 - bb^2}, bb] /. bb -> be) .
        {Xs[t0 + ta ep, z0], zt ep Bf[Xs[t0 + ta ep, z0]], Xs[t0 + ta ep, z0 + zt ep], zt ep Bf[Xs[t0 + ta ep, z0 + zt ep]]})/(zt ep);
  lk = Bf[Sh] . Cross[St, Sz];
  ser = Normal[Series[lk, {ep, 0, 3}]];
  P["leakage B(S_h).(S_t x S_z), coefficients of ep^0..3 (dt = ta ep, dz = zt ep):"];
  Do[P["  ep^", j, ": ", Chop[N[Expand[Coefficient[ser, ep, j]], 12], 10^-40]], {j, 0, 3}];
  gm = Sh - Xs[t0 + al ta ep, z0 + be zt ep];
  ser = Normal[Series[gm, {ep, 0, 4}]];
  P["geometry S_h - S, coefficients of ep^0..4:"];
  Do[P["  ep^", j, ": ", Chop[N[Expand[Coefficient[ser, ep, j]], 12], 10^-40]], {j, 0, 4}]];
P["   expected: leakage O(dt^2) + O(dz^3), geometry O(dt^2) + O(dz^4)"];

P["== S5b polygon flux error of one chord over a curve singularity at fraction lam of the chord"];
Module[{f, lin, area, lamv},
  area[f_] := Simplify[Integrate[f[sv] - (f[-lam h] + (f[(1 - lam) h] - f[-lam h]) (sv + lam h)/h), {sv, -lam h, (1 - lam) h},
       Assumptions -> {h > 0, 0 < lam < 1}], Assumptions -> {h > 0, 0 < lam < 1}];
  P["smooth  f = k s^2      : ", area[k #^2 &], "   -> per chord h^3, M chords: M^-2"];
  P["kink    f = th max(s,0): ", area[th Max[#, 0] &], "   -> h^2 per kink: M^-2"];
  P["fold    f = J sqrt(max(s,0)): ", area[J Sqrt[Max[#, 0]] &], "   -> h^(3/2) per fold: M^-3/2"];
  P["jump    f = g step(s)  : ", area[g UnitStep[#] &], "   -> h per jump: M^-1, and the flux itself is undefined (plateau)"]];

P["== S5c grazing crossing with a tangential trace jump J (normal component continuous)"];
Module[{xz, tin},
  (* path near a face x = 0 with B_x continuous and vanishing at the touch point: x(z) = d - c (z - z1)^2 *)
  xz = d - cc (z - z1)^2;
  tin = Simplify[(z /. Solve[xz == 0, z][[2]]) - (z /. Solve[xz == 0, z][[1]]), Assumptions -> {d > 0, cc > 0}];
  P["length inside the neighbour cell: ", Abs[tin] // Simplify[#, Assumptions -> {d > 0, cc > 0}] &,
   "  => outlet offset J*2 sqrt(d/c): a fold in the flow map of B_h, strength proportional to J"]];
P["   a field line lying in a face (S1) makes the fold a layer of width ~ (B_h error)/a in v,",
  " across which the outlet jumps by ~ J * (hugging length): M^-1 until the chord resolves the layer."];
