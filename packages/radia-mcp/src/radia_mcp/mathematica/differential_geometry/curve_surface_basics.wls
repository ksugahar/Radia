#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  curve_surface_basics.wls -- curve and surface differential geometry distilled
  into executable Wolfram Language assertions for Radia.

  This is the Radia-facing, self-contained counterpart of the local background
  textbook "Mathematica 曲線と曲面の微分幾何.pdf".  The PDF is NOT a runtime
  dependency and is not copied into radia-mcp.  What ships is this executable
  knowledge:

    [1] curve geometry: speed, curvature, torsion;
    [2] surface geometry: parametrization, first fundamental form, area element,
        unit normal, second fundamental form, Gaussian / mean curvature;
    [3] stream-function geometry: grad_Gamma, K = n x grad_Gamma psi,
        tangentiality and div_Gamma K = 0;
    [4] surface Laplace-Beltrami as div_Gamma grad_Gamma.

  These are the foundations used by hodograph.wls, surface_derham.wls,
  dtn_geometry.wls, and the radia-streamfunction surface-current path.
*)

ClearAll["Global`*"];

$Assumptions =
  R > 0 && a > 0 && b > 0 && 0 < th < Pi && 0 < ph < 2 Pi &&
  Element[{R, a, b, t, th, ph, z}, Reals];

fs[x_] := FullSimplify[x];

curveSpeed[r_, u_] := fs[Sqrt[D[r, u].D[r, u]]];
curveCurvature[r_, u_] := fs[
  Sqrt[Cross[D[r, u], D[r, {u, 2}]].Cross[D[r, u], D[r, {u, 2}]]]/
    (D[r, u].D[r, u])^(3/2)
];
curveTorsion[r_, u_] := fs[
  Det[{D[r, u], D[r, {u, 2}], D[r, {u, 3}]}]/
    (Cross[D[r, u], D[r, {u, 2}]].Cross[D[r, u], D[r, {u, 2}]])
];

surfaceBasis[X_, vars_List] := D[X, #] & /@ vars;
surfaceMetric[X_, vars_List] := fs[
  Table[surfaceBasis[X, vars][[i]].surfaceBasis[X, vars][[j]],
    {i, Length[vars]}, {j, Length[vars]}]
];
areaDensity[X_, vars_List] := fs[Sqrt[Det[surfaceMetric[X, vars]]]];
unitNormal[X_, vars_List] := fs[
  Cross[surfaceBasis[X, vars][[1]], surfaceBasis[X, vars][[2]]]/
    areaDensity[X, vars]
];
secondFundamental[X_, vars_List] := Module[{n = unitNormal[X, vars], basis},
  basis = surfaceBasis[X, vars];
  fs[Table[n.D[basis[[i]], vars[[j]]], {i, Length[vars]}, {j, Length[vars]}]]
];
gaussianCurvature[X_, vars_List] := fs[
  Det[secondFundamental[X, vars]]/Det[surfaceMetric[X, vars]]
];
meanCurvature[X_, vars_List] := fs[
  Tr[Inverse[surfaceMetric[X, vars]].secondFundamental[X, vars]]/2
];

surfaceGradient[psi_, X_, vars_List] := Module[{g, gi, basis, nvar},
  g = surfaceMetric[X, vars]; gi = Inverse[g]; basis = surfaceBasis[X, vars];
  nvar = Length[vars];
  fs[Sum[gi[[i, j]] D[psi, vars[[j]]] basis[[i]], {i, nvar}, {j, nvar}]]
];
contravariantComponents[V_, X_, vars_List] := Module[{g, gi, basis, cov, nvar},
  g = surfaceMetric[X, vars]; gi = Inverse[g]; basis = surfaceBasis[X, vars];
  nvar = Length[vars];
  cov = Table[V.basis[[j]], {j, nvar}];
  fs[Table[Sum[gi[[i, j]] cov[[j]], {j, nvar}], {i, nvar}]]
];
surfaceDivergence[V_, X_, vars_List] := Module[{sqrtg, comp},
  sqrtg = areaDensity[X, vars];
  comp = contravariantComponents[V, X, vars];
  fs[(1/sqrtg) Sum[D[sqrtg comp[[i]], vars[[i]]], {i, Length[vars]}]]
];
laplaceBeltrami[psi_, X_, vars_List] := fs[
  surfaceDivergence[surfaceGradient[psi, X, vars], X, vars]
];

(* ============================ self-test ============================ *)
Module[{ok = True, chk, circle, helix, sphere, cylinder, plane, gSphere,
        nSphere, kSphere, gradSphere, kPlane, psi, xchart, jac},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  Print["curve_surface_basics.wls self-test (curves, surfaces, stream-function geometry):"];

  (* ---- [1] curve geometry ---- *)
  circle = {R Cos[t], R Sin[t], 0};
  chk["circle speed = R", fs[curveSpeed[circle, t] == R]];
  chk["circle curvature = 1/R", fs[curveCurvature[circle, t] == 1/R]];
  chk["circle torsion = 0", fs[curveTorsion[circle, t] == 0]];

  helix = {a Cos[t], a Sin[t], b t};
  chk["helix speed = sqrt(a^2+b^2)", fs[curveSpeed[helix, t] == Sqrt[a^2 + b^2]]];
  chk["helix curvature = a/(a^2+b^2)", fs[curveCurvature[helix, t] == a/(a^2 + b^2)]];
  chk["helix torsion = b/(a^2+b^2)", fs[curveTorsion[helix, t] == b/(a^2 + b^2)]];

  (* ---- [2] surface parametrization and curvature ---- *)
  sphere = {R Sin[th] Cos[ph], R Sin[th] Sin[ph], R Cos[th]};
  gSphere = {{R^2, 0}, {0, R^2 Sin[th]^2}};
  chk["sphere first fundamental form = diag(R^2, R^2 sin(th)^2)",
      fs[surfaceMetric[sphere, {th, ph}] == gSphere]];
  chk["sphere area density = R^2 sin(th)",
      fs[areaDensity[sphere, {th, ph}] == R^2 Sin[th]]];
  nSphere = unitNormal[sphere, {th, ph}];
  chk["sphere unit normal = x/R (outward orientation)",
      fs[nSphere == sphere/R]];
  chk["sphere Gaussian curvature = 1/R^2",
      fs[gaussianCurvature[sphere, {th, ph}] == 1/R^2]];
  chk["sphere mean curvature sign follows the chosen outward normal: H = -1/R",
      fs[meanCurvature[sphere, {th, ph}] == -1/R]];

  cylinder = {R Cos[ph], R Sin[ph], z};
  chk["cylinder first fundamental form = diag(R^2, 1)",
      fs[surfaceMetric[cylinder, {ph, z}] == {{R^2, 0}, {0, 1}}]];
  chk["cylinder area density = R",
      fs[areaDensity[cylinder, {ph, z}] == R]];
  chk["cylinder Gaussian curvature = 0",
      fs[gaussianCurvature[cylinder, {ph, z}] == 0]];
  chk["cylinder mean curvature = -1/(2R) for outward orientation",
      fs[meanCurvature[cylinder, {ph, z}] == -1/(2 R)]];

  (* metric = J^T J: this is the chart-level version of the Hodge-weight story. *)
  xchart = {u + v, u^2 - v, u v};
  jac = D[xchart, {{u, v}}];
  chk["first fundamental form is J^T J for an arbitrary surface chart",
      fs[surfaceMetric[xchart, {u, v}] == Transpose[jac].jac]];

  (* ---- [3] stream-function surface current geometry ---- *)
  gradSphere = surfaceGradient[Cos[th], sphere, {th, ph}];
  chk["sphere surface-gradient norm: |grad_Gamma cos(th)|^2 = sin(th)^2/R^2",
      fs[gradSphere.gradSphere == Sin[th]^2/R^2]];

  psi = psiF[th, ph];
  kSphere = fs[Cross[unitNormal[sphere, {th, ph}],
                     surfaceGradient[psi, sphere, {th, ph}]]];
  chk["stream-function current is tangential on the sphere: K.n = 0",
      fs[kSphere.unitNormal[sphere, {th, ph}] == 0]];
  chk["stream-function current is surface-divergence-free on the sphere: div_Gamma(n x grad_Gamma psi)=0",
      fs[surfaceDivergence[kSphere, sphere, {th, ph}] == 0]];

  plane = {x, y, 0};
  kPlane = fs[Cross[unitNormal[plane, {x, y}], surfaceGradient[psiP[x, y], plane, {x, y}]]];
  chk["planar stream-function formula gives K = {-psi_y, psi_x, 0}",
      fs[kPlane == {-D[psiP[x, y], y], D[psiP[x, y], x], 0}]];
  chk["planar stream-function current is divergence-free",
      fs[surfaceDivergence[kPlane, plane, {x, y}] == 0]];

  (* ---- [4] Laplace-Beltrami = div_Gamma grad_Gamma ---- *)
  chk["sphere Laplace-Beltrami formula recovered from div_Gamma grad_Gamma",
      fs[laplaceBeltrami[psiF[th, ph], sphere, {th, ph}] ==
         (1/(R^2 Sin[th])) D[Sin[th] D[psiF[th, ph], th], th] +
         (1/(R^2 Sin[th]^2)) D[psiF[th, ph], {ph, 2}]]];
  chk["Laplace-Beltrami eigenvalue for cos(th): Delta_Gamma cos(th) = -2 cos(th)/R^2",
      fs[laplaceBeltrami[Cos[th], sphere, {th, ph}] == -2 Cos[th]/R^2]];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
