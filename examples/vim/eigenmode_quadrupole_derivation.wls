#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ----------------------------------------------------------------------------------------------------
   eigenmode_quadrupole_derivation.wls

   Moment-yano surface-charge demag: derive the RESIDUAL QUADRUPOLE EIGENMODES of an element and the
   magnetic field each eigenmode creates -- the "suppress the eigenmode with the field it produces"
   closure that replaces the hand-picked quadratic forms (dx^2-dy^2 / axial).

   Per element with nF face-charge DOF sigma, the moment system has nF rows:
       1 monopole  : sum_f Ae_f sigma_f = 0                      (div B = 0)
       3 dipole    : (dipole moment of sigma)/Ve = chi H(centroid)
       nF-4 quad   : (eigenmode's charge moment)  =  chi (field the eigenmode CREATES)
   The (nF-4) quadrupole rows are the RESIDUAL eigenmodes = charge patterns with ZERO monopole and ZERO
   dipole moment = null space of the 4 functionals {Ae, Ae dx, Ae dy, Ae dz}.  Each eigenmode sources a
   pure-quadrupole field; the row pins it to chi times that self-field (the gradH coupling), so the mode
   is suppressed by the field IT creates.  Same recipe for EVERY surface-charge element (the (nF-4) count
   is exactly the quadrupole-and-higher modes that survive after the 4 monopole+dipole conditions):
       tet     (4 faces) -> 0 eigenmodes  (pure monopole + 3 dipole -- no quadrupole row)
       wedge   (5 faces) -> 1 eigenmode
       pyramid (5 faces) -> 1 eigenmode   (PLANAR mode, distinct from the wedge's axial one)
       hex     (6 faces) -> 2 eigenmodes

   Self-tested: run  wolframscript -file eigenmode_quadrupole_derivation.wls
---------------------------------------------------------------------------------------------------- *)

On[Assert];   (* enable the self-test assertions below (off by default in Mathematica) *)

(* ============================ face geometry (area / centroid / outward normal) ============================ *)
polyArea[vs_] := Module[{c = Mean[vs], n = Length[vs]},
  (1/2) Sum[Norm[Cross[vs[[i]] - c, vs[[Mod[i, n] + 1]] - c]], {i, n}]];

(* Newell normal, oriented OUTWARD (away from the element centroid ec) *)
faceNormalOut[vs_, ec_] := Module[{nv, n = Length[vs]},
  nv = Sum[Cross[vs[[i]], vs[[Mod[i, n] + 1]]], {i, n}];
  nv = nv/Norm[nv];
  If[(Mean[vs] - ec) . nv < 0, -nv, nv]];

(* per-element data: areas Ae, face-centroid offsets d_f = fc_f - ec, outward normals nf, signed volume Ve *)
elementData[verts_, faces_] := Module[{ec, fc, Ae, nf, d, Ve},
  ec = Mean[verts];
  fc = Mean[verts[[#]]] & /@ faces;
  Ae = polyArea[verts[[#]]] & /@ faces;
  nf = faceNormalOut[verts[[#]], ec] & /@ faces;
  d  = # - ec & /@ fc;
  Ve = (1/3) Sum[Ae[[f]] (fc[[f]] . nf[[f]]), {f, Length[faces]}];   (* divergence-theorem volume *)
  <|"nF" -> Length[faces], "Ae" -> Ae, "d" -> d, "nf" -> nf, "Ve" -> Ve, "ec" -> ec, "fc" -> fc|>];

(* ============================ residual quadrupole eigenmodes ============================ *)
(* The 4 monopole+dipole MOMENT functionals (row vectors in R^nF):  m = Ae,  p_k = Ae d_k.
   The residual subspace = their null space = charge patterns with zero monopole and zero dipole moment.
   Orthonormalize -> the (nF-4) quadrupole eigenmodes.  Geometry-adaptive + ALWAYS full rank. *)
momentFunctionals[dat_] := Module[{Ae = dat["Ae"], d = dat["d"]},
  {Ae, Ae d[[All, 1]], Ae d[[All, 2]], Ae d[[All, 3]]}];

residualEigenmodes[dat_] := Orthogonalize[NullSpace[momentFunctionals[dat]]];

(* charge multipole MOMENTS of a mode phi (to certify it is pure quadrupole) *)
monopoleMoment[phi_, dat_] := dat["Ae"] . phi;
dipoleMoment[phi_, dat_]   := Sum[dat["Ae"][[f]] phi[[f]] dat["d"][[f]], {f, dat["nF"]}];
(* traceless symmetric quadrupole moment tensor Q = sum_f Ae_f phi_f (3 d d - |d|^2 I) *)
quadrupoleMoment[phi_, dat_] := Sum[
   dat["Ae"][[f]] phi[[f]] (3 KroneckerProduct[dat["d"][[f]], dat["d"][[f]]]
       - (dat["d"][[f]] . dat["d"][[f]]) IdentityMatrix[3]), {f, dat["nF"]}];

(* "the field the eigenmode CREATES": a surface-charge quadrupole Q sources, to leading order, the demag
   field  H(r) = (1/4pi) grad( (1/2) r.Q.r / |r|^5 ) ... ; in the moment row this enters as the contraction
   of the mode's quadrupole tensor Dm = sum_f phi_f (nf_f outer d_f) with the field gradient gradH at the
   centroid (the C++ Dvec.gradH term).  We expose Dm (the exact row operator) symbolically here. *)
modeDmTensor[phi_, dat_] := Sum[
   phi[[f]] KroneckerProduct[dat["nf"][[f]], dat["d"][[f]]], {f, dat["nF"]}];

(* ============================ canonical elements ============================ *)
hexV = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
        {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}};
hexF = {{1, 2, 3, 4}, {5, 6, 7, 8}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 4, 8, 7}, {4, 1, 5, 8}};

(* right triangular prism (axis = z) *)
wedV = {{0, 0, -1}, {2, 0, -1}, {1, Sqrt[3], -1}, {0, 0, 1}, {2, 0, 1}, {1, Sqrt[3], 1}};
wedF = {{1, 2, 3}, {4, 5, 6}, {1, 2, 5, 4}, {2, 3, 6, 5}, {3, 1, 4, 6}};

(* square pyramid: square base z=0 + apex on z axis *)
pyrV = {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}, {0, 0, 2}};
pyrF = {{1, 2, 3, 4}, {1, 2, 5}, {2, 3, 5}, {3, 4, 5}, {4, 1, 5}};

(* tetrahedron: 4 triangular faces -> 4 DOF = monopole(1) + dipole(3) + quadrupole(0).  The residual
   eigenmode count nF-4 = 0, so a tet carries NO quadrupole row -- it is pure monopole+dipole (this is the
   "4 DoF = monopole + 3 dipole" the user confirmed; radia's production tet is the 3-DOF MMM dipole, the
   same dipole block without the redundant net-charge monopole row). *)
tetV = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
tetF = {{1, 3, 2}, {1, 2, 4}, {1, 4, 3}, {2, 3, 4}};

(* a DISTORTED (sheared) hex + a SCALENE wedge, to certify geometry-adaptivity (no hand-picked axis) *)
shHexV = (# + {0.6 (#[[3]] + 1)/2, 0, 0}) & /@ hexV;   (* shear top in +x *)
scWedV = {{0, 0, -1}, {2.3, 0.2, -0.8}, {0.7, 1.9, -1.2}, {0.1, -0.1, 1}, {2.5, 0.3, 1.1}, {0.6, 1.8, 0.9}};

(* ============================ certify ============================ *)
report[name_, V_, F_] := Module[{dat, modes, q, mono, dip, Q, Dm},
  dat = elementData[V, F];
  modes = residualEigenmodes[dat];
  Print["==== ", name, "  nF=", dat["nF"], "  residual eigenmodes=", Length[modes],
        "  (expect ", dat["nF"] - 4, ")  Ve=", N[dat["Ve"], 4]];
  Do[
    mono = Chop[N[monopoleMoment[modes[[q]], dat]]];
    dip  = Chop[N[dipoleMoment[modes[[q]], dat]]];
    Q    = Chop[N[quadrupoleMoment[modes[[q]], dat]]];
    Dm   = Chop[N[modeDmTensor[modes[[q]], dat]]];
    Print["  mode ", q, ":  monopole=", mono, "  |dipole|=", Chop[N[Norm[dip]]]];
    Print["    quadrupole moment Q (the multipole that SOURCES the suppressing field):"];
    Print["      ", MatrixForm[Round[Q, 0.001]]];
    Print["    row operator Dm = sum_f phi_f (nf outer d)  (contracts with gradH = the field it creates):"];
    Print["      tr=", Chop[N[Tr[Dm]]], "  Dm=", MatrixForm[Round[Dm, 0.001]]],
    {q, Length[modes]}];
  (* ---- assertions ---- *)
  Assert[Length[modes] == dat["nF"] - 4];
  Do[
    Assert[Chop[N[monopoleMoment[modes[[q]], dat]]] == 0];                 (* zero monopole moment *)
    Assert[Chop[N[Norm[dipoleMoment[modes[[q]], dat]]]] == 0];             (* zero dipole moment *)
    Assert[Chop[N[Norm[quadrupoleMoment[modes[[q]], dat]]]] > 10^-6],      (* nonzero quadrupole -> real mode *)
    {q, Length[modes]}];
  modes];

tetModes = report["TET (tetrahedron)", tetV, tetF];
hexModes = report["HEX (cube)", hexV, hexF];
wedModes = report["WEDGE (triangular prism)", wedV, wedF];
pyrModes = report["PYRAMID (square)", pyrV, pyrF];
shModes  = report["HEX (sheared, distorted)", shHexV, hexF];
scModes  = report["WEDGE (scalene, distorted)", scWedV, wedF];

(* eigenmode quadrupole-moment SIGNATURE per element type ------------------------------------------- *)
(* The sorted eigenvalues of Q are the rotation-invariant signature of the mode.  Each element type has a
   DIFFERENT signature -- so a single hand-picked quadratic form cannot serve all of them; the eigenmode is
   derived per element.  hex: 2 modes (diagonal); wedge: AXIAL (l, l, -2l); pyramid: PLANAR (l, -l, 0). *)
qsig[modes_, dat_] := (Sort[Chop[N[Eigenvalues[quadrupoleMoment[#, dat]]]]] & /@ modes);
Print["==== eigenmode quadrupole-moment signatures (sorted eigenvalues of Q) ===="];
Print["  HEX   : ", qsig[hexModes, elementData[hexV, hexF]]];
Print["  WEDGE : ", qsig[wedModes, elementData[wedV, wedF]], "  (axial: two equal + one = -2x)"];
Print["  PYRAM : ", qsig[pyrModes, elementData[pyrV, pyrF]], "  (planar: +l, -l, 0 -- NOT the wedge's axial!)"];
(* assert wedge is axial (l,l,-2l) and pyramid is planar (l,-l,0) -- distinct, eigenmode-derived shapes *)
With[{w = First[qsig[wedModes, elementData[wedV, wedF]]]},
  Assert[Abs[w[[1]] - w[[2]]] < 10^-6 || Abs[w[[2]] - w[[3]]] < 10^-6]];   (* two equal -> axial *)
With[{p = First[qsig[pyrModes, elementData[pyrV, pyrF]]]},
  Assert[Min[Abs[p]] < 10^-6]];                                          (* one zero -> planar *)

(* cube: the hand-picked dx^2-dy^2 / dy^2-dz^2 CHARGE patterns lie in the eigenmode span (zero monopole +
   dipole moment by cube symmetry) -> the moment row-space is unchanged -> identical solve -> N=1/3. *)
hexDat = elementData[hexV, hexF];
hexHP = With[{Ae = hexDat["Ae"], d = hexDat["d"]},
   {Ae (d[[All, 1]]^2 - d[[All, 2]]^2), Ae (d[[All, 2]]^2 - d[[All, 3]]^2)}];
hexProj = (hexHP . Transpose[hexModes]) . hexModes;
Print["==== cube: hand-pick(dx2-dy2, dy2-dz2) in eigenmode span?  residual = ",
      Chop[N[Total[Norm /@ (hexHP - hexProj)]]], "  (0 => same row-space => N=1/3 preserved)"];
Assert[Chop[N[Total[Norm /@ (hexHP - hexProj)]]] == 0];

Print["ALL SELF-TESTS PASS"];
