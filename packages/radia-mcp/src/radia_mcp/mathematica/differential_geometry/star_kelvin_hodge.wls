(* Star-shaped Kelvin inversion: closed forms against the Jacobian of the map, at exact rational points.
   y = T(x) = rho(n)^2 x/r^2, n = x/r, rho homogeneous of degree 0, t = r grad(rho)/rho, w = n - t.
   Claims checked for each boundary:
     J = (rho/r)^2 (I - 2 n w^T), n.t = 0, A^2 = I, det A = -1;
     nu' = J^T J/|det J| = (r/rho)^2 (I - 2 (n t^T + t n^T) + 4 t t^T);
     mu' = J^-1 J^-T |det J| = (rho/r)^2 (I + 2 (n t^T + t n^T) + 4 |t|^2 n n^T), mu'.nu' = I;
     grad(Omega o T) = J^T grad(Omega) o T for a test potential;
     cond(mu') = (tan a + Sqrt[1 + tan^2 a])^4 with tan a = |t|.
   Pointwise evaluation to 60 digits; not a symbolic proof over all directions.
   Run: wolframscript -file star_kelvin_hodge.wls
   Prints one line per boundary with the largest residual, and exits nonzero above 10^-40. *)

digits = 60;
xv = {x1, x2, x3};
r = Sqrt[xv . xv];
nv = xv/r;
axial = Sqrt[x1^2 + x2^2]/r;
height = Sqrt[x3^2]/r;

(* Each boundary is rho as a function of position, degree 0; the points lie on a branch where it is smooth. *)
boundaries = {
  {"sphere", 3/10, {{1/10, 1/5, -1/7}, {-1/3, 1/11, 1/13}}},
  {"flat face", (1/10)/height, {{1/50, 1/40, 3/10}, {-1/30, 1/25, -1/5}}},
  {"rounded edge", (1/5)*axial + Sqrt[(1/25)*axial^2 - 1/25 + 1/100], {{3/5, 1/5, 1/10}, {-1/4, 1/2, -1/20}}},
  {"ellipsoid", 1/Sqrt[(nv[[1]]/(3/10))^2 + (nv[[2]]/(1/5))^2 + (nv[[3]]/(1/10))^2], {{1/7, 1/9, 1/11}, {-1/5, 1/6, -1/8}}}
};

residual[rho_, point_] := Module[{rules, tv, wv, ratio, jac, formJ, mat, nuClosed, muClosed, nuForm, muForm, numT, sigma, tanA, testO, lhs, rhs, values},
  rules = Thread[xv -> point];
  tv = r D[rho, {xv}]/rho;
  wv = nv - tv;
  ratio = rho/r;
  mat = IdentityMatrix[3] - 2 Outer[Times, nv, wv];
  formJ = ratio^2 mat;
  jac = D[rho^2 xv/r^2, {xv}];
  numT = rho^2 xv/r^2;
  nuClosed = Transpose[jac] . jac/Abs[Det[jac]];
  muClosed = Inverse[jac] . Transpose[Inverse[jac]] Abs[Det[jac]];
  nuForm = ratio^-2 (IdentityMatrix[3] - 2 (Outer[Times, nv, tv] + Outer[Times, tv, nv]) + 4 Outer[Times, tv, tv]);
  muForm = ratio^2 (IdentityMatrix[3] + 2 (Outer[Times, nv, tv] + Outer[Times, tv, nv]) + 4 (tv . tv) Outer[Times, nv, nv]);
  testO = {y1, y2, y3} . {1/3, -2/7, 5/11} + y1 y2 - y3^2/2;
  lhs = D[testO /. Thread[{y1, y2, y3} -> numT], {xv}];
  rhs = Transpose[jac] . (D[testO, {{y1, y2, y3}}] /. Thread[{y1, y2, y3} -> numT]);
  values = N[{jac - formJ, nv . tv, mat . mat - IdentityMatrix[3], Det[mat] + 1, nuClosed - nuForm, muClosed - muForm,
              muForm . nuForm - IdentityMatrix[3], lhs - rhs} /. rules, digits];
  sigma = N[Sqrt[tv . tv] /. rules, digits];
  tanA = sigma;
  Max[Max[Abs[Flatten[values]]],
      Abs[(Max[#]/Min[#] &)[Eigenvalues[N[muForm /. rules, digits]]] - (tanA + Sqrt[1 + tanA^2])^4]]];

worst = 0;
Do[Module[{value},
    value = Quiet[Max[residual[b[[2]], #] & /@ b[[3]]]];
    worst = Max[worst, value];
    Print[b[[1]], ": largest residual below 10^", Floor[Log10[Max[value, 10^-(digits - 5)]]] + 1]],
  {b, boundaries}];
If[worst > 10^-40, Print["FAILED"]; Exit[1], Print["PASSED"]; Exit[0]];
