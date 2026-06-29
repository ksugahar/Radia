#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  bem_cln_ih_workpiece.wls

  Phase 4 final-4: realistic induction-heating workpiece + multi-coil
  BEM-CLN demonstration.

  Geometry:
    Workpiece: cylindrical Cu billet, radius R_w = 15 mm, length L_w = 50 mm
    Coil: 3-turn helical coil around workpiece, R_c = 30 mm centre,
          axial spacing 10 mm between turns
    Total elements N = 1 workpiece + 3 coil turns = 4 elements

  Method:
    Each element treated as a 2D long cylinder (workpiece + each coil
    turn).  Per-element polarizability alpha_dip(s) from Phase 2.5
    rigorous Bessel framework.  Inter-element coupling via 2D
    dipole-dipole 1/D^2 (using projected distance in the cross-section
    plane perpendicular to the workpiece axis).

  This is a SIMPLIFIED 2D model of the 3D physics, suitable for
  demonstrating the BEM-CLN method on a real engineering geometry.
  Full 3D NGSolve FEM cross-check is a separate undertaking.

  Run: wolframscript -file bem_cln_ih_workpiece.wls
*)

ClearAll[a, sigma, mu, dParam, basisN];

(* Workpiece: 30 mm diameter (R_w = 15 mm) *)
Rw    = 15*^-3;
(* Coil turn: 5 mm diameter wire (a_c = 2.5 mm radius) *)
ac    = 2.5*^-3;
(* Coil radius (centre of turn) *)
Rc    = 30*^-3;
(* Materials *)
sigma = N[5.96*^7, 30];   (* Cu *)
mu    = N[4 Pi 1*^-7, 30];
dParam = N[2 Pi 1*^5, 30];
basisN = 4;

KsibcWp = 2 Pi Rw Sqrt[sigma/mu];   (* workpiece per unit length SIBC *)
KsibcCoil = 2 Pi ac Sqrt[sigma/mu]; (* coil wire per unit length SIBC *)

Print["=== bem_cln_ih_workpiece.wls ==="];
Print["IH workpiece + 3-turn coil BEM-CLN demonstration"];
Print["  Workpiece: R_w = ", N[Rw 1000, 4], " mm, K_SIBC^(wp) = ", N[KsibcWp, 4]];
Print["  Coil wire: a_c = ", N[ac 1000, 4], " mm, R_c = ", N[Rc 1000, 4], " mm"];
Print["  K_SIBC^(coil) per turn = ", N[KsibcCoil, 4]];
Print[""];

(* === Per-cylinder Bessel polarizability ================================ *)
YcylBessel[s_, aVal_] := Module[{gamma},
  gamma = Sqrt[s mu sigma];
  N[Pi aVal^2 sigma * 2 BesselI[1, gamma aVal] / (gamma aVal BesselI[0, gamma aVal]), 20]];

alphaCyl[s_, aVal_] := N[aVal^2 - YcylBessel[s, aVal] / (Pi sigma), 20];

(* === System geometry (projected onto axial cross-section) =============== *)
(*
  Workpiece at origin (axis along z).
  3 coil turns at angles 0, 120, 240 degrees in the xy plane around
  workpiece axis.  For simplicity, treat all 4 elements as parallel
  2D cylinders in the (x, y) plane.
  Position list: {workpiece, coil_turn_1, coil_turn_2, coil_turn_3}
*)

positions = {
  {0, 0},                                   (* workpiece *)
  {Rc, 0},                                  (* coil turn 1 *)
  {Rc Cos[2 Pi/3], Rc Sin[2 Pi/3]},         (* coil turn 2 *)
  {Rc Cos[4 Pi/3], Rc Sin[4 Pi/3]}          (* coil turn 3 *)
};

(* Radii per element: workpiece is large, coil turns are small *)
radiiList = {Rw, ac, ac, ac};
NList     = Length[positions];

Print["4-element configuration:"];
Print["  Element 1 (workpiece): pos=(0,0), R=", N[Rw 1000, 4], " mm"];
Print["  Element 2 (coil tn 1): pos=", N[positions[[2]] 1000, 4], " mm, R=", N[ac 1000, 4], " mm"];
Print["  Element 3 (coil tn 2): pos=", N[positions[[3]] 1000, 4], " mm, R=", N[ac 1000, 4], " mm"];
Print["  Element 4 (coil tn 3): pos=", N[positions[[4]] 1000, 4], " mm, R=", N[ac 1000, 4], " mm"];
Print[""];

(* === Coupling matrix G_ij = 1/|r_i - r_j|^2 for i != j ================= *)
couplingG = Table[
  If[i == j, 0, N[1/((positions[[i]] - positions[[j]]) . (positions[[i]] - positions[[j]])), 20]],
  {i, 1, NList}, {j, 1, NList}];

Print["Coupling magnitudes (1/D_ij^2):"];
Do[
  Do[
    If[i != j,
      Print["  G[", i, ",", j, "] = ", N[couplingG[[i, j]], 5],
            "  (D = ", N[Norm[positions[[i]] - positions[[j]]] 1000, 4], " mm)"]],
    {j, i + 1, NList}],
  {i, 1, NList}];
Print[""];

(* === Solve N x N BEM-CLN system ======================================== *)
solveBemCLNih[s_] := Module[
  {alphaList, alphaDiag, sysMat, rhs, cvec, B0 = 1.0},
  alphaList = Table[N[alphaCyl[s, radiiList[[k]]], 20], {k, 1, NList}];
  alphaDiag = DiagonalMatrix[alphaList];
  sysMat = IdentityMatrix[NList] - alphaDiag . couplingG;
  rhs = alphaList B0;
  cvec = LinearSolve[sysMat, rhs];
  cvec];

(* === Wall-band sweep ==================================================== *)
Print["BEM-CLN response across wall band:"];
Print[""];
Print[StringPadRight["f [Hz]", 12] <>
      StringPadRight["c_workpiece", 22] <>
      StringPadRight["c_coil_1", 22] <>
      "c_coil_2 (= c_coil_3)"];
Print[StringRepeat["-", 80]];

Do[
  With[{omega = 2 Pi (10.0^logf), s = I 2 Pi (10.0^logf)},
    With[{cvec = solveBemCLNih[s]},
      Print["f=1e", logf, "      ",
            StringPadRight[ToString[N[cvec[[1]], 5]], 22],
            StringPadRight[ToString[N[cvec[[2]], 5]], 22],
            N[cvec[[3]], 5]]
    ]],
  {logf, {3, 4, 5, 6, 7, 8}}];

Print[""];

(* === Workpiece response with vs without coil ============================ *)
Print["Effect of coil presence on workpiece response:"];
Print[""];
With[{f = 1*^5, omega = 2 Pi 1*^5, s = I 2 Pi 1*^5},
  With[{cvecWithCoil = solveBemCLNih[s]},
    With[{
      (* workpiece alone with same applied B_0 = 1 *)
      cAlone = alphaCyl[s, Rw] (* in absence of coil *)
    },
      Print["At f = 1e5 Hz:"];
      Print["  Workpiece alone (no coil)   c_wp = ", N[cAlone, 5]];
      Print["  Workpiece in 3-turn coil    c_wp = ", N[cvecWithCoil[[1]], 5]];
      Print["  Effect of coil on workpiece = ratio ",
            N[cvecWithCoil[[1]]/cAlone, 5]];
      Print["  Real part change            = ",
            N[100 (Re[cvecWithCoil[[1]]]/Re[cAlone] - 1), 4], " %"]
    ]]];

(* === Power-absorbed estimate ============================================ *)
Print[""];
Print["Power-absorbed estimate (engineering quantity):"];
Print[""];
Print["For applied B_0 = 1 T at f = 1e5 Hz, power dissipated in workpiece:"];
With[{f = 1*^5, omega = 2 Pi 1*^5, s = I 2 Pi 1*^5},
  With[{cvec = solveBemCLNih[s]},
    With[{Pwp = omega Im[cvec[[1]]] Re[1]},
      Print["  Im(c_wp) (resistive part)    = ", N[Im[cvec[[1]]], 5]];
      Print["  Pwp per unit z (= omega Im(c) B^2) = ", N[Pwp, 5], " W/m"]
    ]]];

(* === CSV output ========================================================= *)
Print[""];
Print["Writing IH workpiece CSV..."];
csvData = {{"f_Hz", "c_workpiece_re", "c_workpiece_im",
            "c_coil1_re", "c_coil1_im",
            "c_coil2_re", "c_coil2_im",
            "Pwp_per_z"}};
Do[
  With[{omega = 2 Pi (10.0^logf), s = I 2 Pi (10.0^logf)},
    With[{cvec = solveBemCLNih[s]},
      AppendTo[csvData,
        {10.0^logf,
         N[Re[cvec[[1]]], 6], N[Im[cvec[[1]]], 6],
         N[Re[cvec[[2]]], 6], N[Im[cvec[[2]]], 6],
         N[Re[cvec[[3]]], 6], N[Im[cvec[[3]]], 6],
         N[omega Im[cvec[[1]]], 6]}]
    ]],
  {logf, {2, 3, 4, 5, 6, 7, 8}}];

Export["bem_cln_ih_workpiece.csv", csvData, "CSV"];
Print["Wrote bem_cln_ih_workpiece.csv (", Length[csvData] - 1, " freq points)"];
Print[""];

Print["=== bem_cln_ih_workpiece.wls DONE ==="];
