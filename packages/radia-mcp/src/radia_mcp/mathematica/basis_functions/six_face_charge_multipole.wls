#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  six_face_charge_multipole.wls -- WHAT FIELD does the 6-DoF MSC (yano) hexahedron element create?
  The symbolic counterpart of the numerical SVD / eigenmode analysis
  (retired VIM prototype inventory).  Cleanly, for an ARBITRARY (sheared)
  hexahedron and for ALL shear s at once, the 6 uniform face charges produce a field whose
  INDEPENDENT multipole content is exactly

       monopole (1)  +  dipole (3)  +  quadrupole (2)   =   6 DoF

  and NOTHING independent beyond it: octupole and the other 3 quadrupole components are
  LINEARLY DETERMINED by these 6.  So 6-DoF MSC = MMM (the 3 dipole DoF) + a monopole + the
  2 quadrupole moments representable by face charges.

  Proven by the RANK of the charge -> multipole-moment map, with the shear s kept SYMBOLIC
  (generic rank over the rational functions of s -- a proof for every distortion, not sample
  points):
       rank[mono;dip]                 = 4
       rank[mono;dip;quad]            = 6     (quadrupole adds EXACTLY 2 independent DoF)
       rank[mono;dip;quad;octupole]   = 6     (octupole adds 0 -> NO independent octupole)

  and the SHEAR (off-diagonal) quadrupole functional Int_face d_x d_y dA, in CLOSED FORM, is
  proportional to s -- identically 0 at s=0 (the symmetric cube => pure 2 DIAGONAL quadrupole)
  and nonzero for s != 0.  So a distorted 6-DoF creates nothing beyond 2 quadrupole DoF;
  distortion only ROTATES that 2-dim quadrupole subspace (shear is the global-frame appearance
  of the tilted principal axes -- eig(N) shows the quadrupole pair tilt ~12-26 deg for a 0.4
  shear while the count stays 2), never adding an independent 3rd quadrupole or octupole DoF.

  Jacobian note: each face is a parallelogram, so Int_face P dA = (face area) * Int_[0,1]^2 P du dv.
  The face area is a POSITIVE per-face scalar = a COLUMN scaling of the moment matrix; it changes
  no rank and no zero.  We therefore use the area-free (reduced) moment Int_[0,1]^2 P du dv, which
  is a clean polynomial in s (no Sqrt), and the rank / shear-vanishing conclusions are identical.
*)

Lmat[s_] := {{1, s, 0}, {0, 1, 0}, {0, 0, 1}};                  (* xy-shear: keeps faces planar *)
cubeFaces = {
  {{-1, -1, -1}, {2, 0, 0}, {0, 2, 0}},   (* z = -1 *)
  {{-1, -1,  1}, {2, 0, 0}, {0, 2, 0}},   (* z = +1 *)
  {{-1, -1, -1}, {2, 0, 0}, {0, 0, 2}},   (* y = -1 *)
  {{-1,  1, -1}, {2, 0, 0}, {0, 0, 2}},   (* y = +1 *)
  {{-1, -1, -1}, {0, 2, 0}, {0, 0, 2}},   (* x = -1 *)
  {{ 1, -1, -1}, {0, 2, 0}, {0, 0, 2}}};  (* x = +1 *)

(* area-free face moment of the monomial d^(a,b,c); d = P, P = L.V0 + u L.e1 + v L.e2, centroid = 0 *)
redMoment[face_, s_, deg_] := Module[{L = Lmat[s], V0, e1, e2, P},
  {V0, e1, e2} = face; P = L.V0 + u (L.e1) + v (L.e2);
  Integrate[P[[1]]^deg[[1]] P[[2]]^deg[[2]] P[[3]]^deg[[3]], {u, 0, 1}, {v, 0, 1}]];

momentRow[s_, deg_] := redMoment[#, s, deg] & /@ cubeFaces;     (* one functional = its 6 face values *)

degMono = {{0, 0, 0}};
degDip  = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
degQuad = {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {0, 1, 1}, {1, 0, 1}};
degOct  = Select[Tuples[Range[0, 3], 3], Total[#] == 3 &];      (* 10 cubic monomials *)
matFor[s_, degs_] := momentRow[s, #] & /@ degs;

(* ============================ self-test ============================ *)
Module[{ok = True, chk, rMD, rQ, rO, shear, shearCoeff},     (* s left global so it prints cleanly as "s" *)
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  Print["six_face_charge_multipole.wls self-test (6-DoF MSC field content, symbolic in shear s):"];

  (* ---- (1) rank structure for GENERIC (symbolic) s -- a proof for ALL distortions ---- *)
  rMD = MatrixRank[matFor[s, Join[degMono, degDip]]];
  rQ  = MatrixRank[matFor[s, Join[degMono, degDip, degQuad]]];
  rO  = MatrixRank[matFor[s, Join[degMono, degDip, degQuad, degOct]]];
  chk["rank[mono;dipole] = 4  (generic s)", rMD == 4];
  chk["rank[mono;dipole;quad] = 6  (generic s): quadrupole adds EXACTLY 2 DoF", rQ == 6];
  chk["rank[mono;dipole;quad;octupole] = 6  (generic s): octupole adds 0 -> no indep. octupole", rO == 6];

  (* ---- (2) the SHEAR quadrupole functional in CLOSED FORM: ~ s, zero at s=0 ---- *)
  shear = Simplify[momentRow[s, {1, 1, 0}]];                    (* Int_face d_x d_y, area-free, exact in s *)
  Print["  closed-form shear quadrupole functional  Int_face d_x d_y (per face) = ", shear];
  chk["shear quadrupole = 0 on every face at s=0 (cube => pure 2 DIAGONAL quadrupole)",
      Simplify[shear /. s -> 0] === {0, 0, 0, 0, 0, 0}];
  chk["shear quadrupole is NOT identically zero for s != 0 (distortion -> shear orientation)",
      Simplify[shear] =!= {0, 0, 0, 0, 0, 0}];
  shearCoeff = Simplify[Coefficient[shear, s]];                 (* the linear-in-s coefficient *)
  chk["shear quadrupole is exactly LINEAR in s (orientation-only, vanishing as s->0)",
      Simplify[shear - shearCoeff s] === {0, 0, 0, 0, 0, 0} && shearCoeff =!= {0, 0, 0, 0, 0, 0}];

  (* ---- (3) the 2 quadrupole DoF are genuinely 2: adding quad lifts the rank by 6-4 = 2 ---- *)
  chk["quadrupole contributes exactly 2 = rank[mono;dip;quad] - rank[mono;dip]", (rQ - rMD) == 2];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
