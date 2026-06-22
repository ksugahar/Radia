#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  canonical.wls -- the HAMILTONIAN / CANONICAL reading of the weak form, made
  symbolic.  Companion to weakform_hodge.wls (the action + Hodge half) and
  hodograph.wls (the canonical-transformation half).  Answers the user's
  insight "弱形式ってハミルトニアンなの? 正準方程式だと思うと理解が進むかも":

  PRECISE STATEMENT.  The weak form itself is  delta(action) = 0  (Lagrangian /
  variational): for magnetostatics S[A] = INT (1/2) nu |curl A|^2 - J.A, and
  delta S = 0 IS the weak form INT nu (curl A).(curl w) = INT J.w.  The
  HAMILTONIAN / canonical structure is the LEGENDRE SHADOW of that action, and it
  shows up in four exact, verifiable ways:

    [1] FLUX LINES are literally Hamilton's canonical equations with A_z the
        Hamiltonian:  in 2-D B = (A_z,y, -A_z,x), so a field line dx/ds = B obeys
            dx/ds = +dA_z/dy ,   dy/ds = -dA_z/dx   (H = A_z, (x,y) = (q,p)).
        A_z is conserved along the flow (energy), the flow is area-preserving
        (Liouville / symplectic).  [Repo: flux_line_closure_symplectic.py,
        Noguchi/Sugahara 2020 symplectic flux-line tracking.]

    [2] The CONSTITUTIVE LAW is a LEGENDRE TRANSFORM:  H = dw/dB (energy density
        w(B)), coenergy w*(H) = H.B - w(B), B = dw*/dH -- exactly p = dL/dq',
        H_ham = p q' - L, q' = dH_ham/dp.  (B,H) are canonical conjugates; the
        tangent dH/dB is the Hessian of w (symmetric => derives from an energy),
        SPD <=> w convex <=> Legendre invertible <=> the ⋆_nu "no fold" of
        weakform_hodge.wls sec.6.

    [3] The HODOGRAPH is the LEGENDRE TRANSFORM of the potential (gas-dynamics
        fact) = a CANONICAL transformation: H = grad Phi <=> the (x, H) graph is
        LAGRANGIAN (the canonical 2-form sum dH_i ^ dx_i vanishes on it), and
        chi(H) = x.H - Phi gives x = dchi/dH (the hodograph inverse).  This is
        the canonical-transformation view of hodograph.wls; Chaplygin is then the
        action-angle / Hamilton-Jacobi step (the transform that trivialises it).

    [4] POISSON BRACKET:  {f,g} = f_x g_y - f_y g_x makes (x,y) a canonical pair
        ({x,y} = 1), and the flow is df/ds = {f, A_z}: {x,A_z} = dx/ds,
        {y,A_z} = dy/ds, {A_z,A_z} = 0 (conservation).

  Lagrangian (primal, one field A, MINIMUM, 2nd order) vs Hamiltonian (mixed,
  (A,B)/(H,B), SADDLE, 1st order) is the FEM face of the same duality (FEEC).
*)

(* ---- 2-D canonical helpers (q,p) = (x,y) ---- *)
hamVF[H_]   := {D[H, y], -D[H, x]};                 (* Hamiltonian vector field (dH/dp, -dH/dq) *)
pb[f_, g_]  := D[f, x] D[g, y] - D[f, y] D[g, x];   (* Poisson bracket {f,g} *)
div2[V_]    := D[V[[1]], x] + D[V[[2]], y];

(* ============================ self-test ============================ *)
Module[{ok = True, chk, zero, Bphys, w2, Hv, tang, r2},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  zero[m_] := Simplify[m] === ConstantArray[0, Dimensions[m]];
  Print["canonical.wls self-test (the Hamiltonian/canonical reading of the weak form):"];

  (* ---- [1] FLUX LINE = HAMILTON'S CANONICAL EQUATIONS, H = A_z ---- *)
  (* physical 2-D field-line velocity B = curl(A_z zhat) = (A_z,y, -A_z,x) *)
  Bphys = {D[az[x, y], y], -D[az[x, y], x]};
  chk["[1a] flux-line velocity B = the Hamiltonian vector field of A_z (H=A_z, (x,y)=(q,p))",
      zero[Simplify[Bphys - hamVF[az[x, y]]]]];
  chk["[1b] A_z is CONSERVED along the flow: dA_z/ds = grad(A_z).B = 0 (energy conservation)",
      Simplify[{D[az[x, y], x], D[az[x, y], y]}.Bphys] === 0];
  chk["[1c] the flow is AREA-PRESERVING (Liouville/symplectic): div B = 0",
      Simplify[div2[Bphys]] === 0];

  (* ---- [2] CONSTITUTIVE LAW = LEGENDRE TRANSFORM  B <-> H ---- *)
  (* 1-D: energy w(B), H = w'(B), coenergy w*(H) = B w'(B) - w(B); verify dw*/dH = B *)
  chk["[2a] Legendre B<->H: with H=w'(B), coenergy w*=B w'(B)-w(B) gives dw*/dH = B",
      Simplify[D[bb w'[bb] - w[bb], bb]/D[w'[bb], bb] - bb] === 0];
  (* vector: w = W(|B|), H = grad_B w; the tangent dH/dB is the HESSIAN => SYMMETRIC
     (the constitutive law DERIVES FROM AN ENERGY <=> tangent symmetric <=> canonical) *)
  w2 = W[Sqrt[bx^2 + by^2 + bz^2]];
  Hv = D[w2, {{bx, by, bz}}];
  tang = D[Hv, {{bx, by, bz}}];
  chk["[2b] tangent dH/dB is the Hessian of the energy => SYMMETRIC (derives from a potential)",
      zero[Simplify[tang - Transpose[tang], bx^2 + by^2 + bz^2 > 0]]];
  (* eigenstructure = weakform_hodge sec.6 (with nu = W'(r)/r): longitudinal W''(r), transverse W'(r)/r *)
  r2 = Sqrt[bx^2 + by^2 + bz^2];
  chk["[2c] longitudinal eigenvalue = W''(|B|)  (tangent.B = W''(|B|) B)  -- SPD iff W convex",
      Simplify[tang.{bx, by, bz} - W''[r2] {bx, by, bz}, bx^2 + by^2 + bz^2 > 0] === {0, 0, 0}];
  chk["[2d] transverse eigenvalue = W'(|B|)/|B| = nu  (tangent.t = (W'/r) t, t perp B)",
      Simplify[tang.{-by, bx, 0} - (W'[r2]/r2) {-by, bx, 0}, bx^2 + by^2 + bz^2 > 0] === {0, 0, 0}];

  (* ---- [3] HODOGRAPH = LEGENDRE TRANSFORM OF THE POTENTIAL = CANONICAL TRANSFORMATION ---- *)
  (* H = grad Phi  <=>  the (x,H) graph is LAGRANGIAN: the canonical 2-form
     dH_x ^ dx + dH_y ^ dy restricted to H=grad Phi has coefficient (curl H)_z = 0 *)
  chk["[3a] H=grad Phi => (x,H) graph is LAGRANGIAN (canonical 2-form sum dH_i^dx_i = 0 <=> curl H = 0)",
      Simplify[D[D[ph[x, y], y], x] - D[D[ph[x, y], x], y]] === 0];
  (* the hodograph IS the Legendre transform of the potential: chi(H) = x H - Phi(x),
     H = Phi'(x)  =>  x = dchi/dH  (same Legendre round-trip as [2a], on the POTENTIAL) *)
  chk["[3b] hodograph = Legendre transform of Phi: chi=x Phi'(x)-Phi(x) gives dchi/dH = x  (canonical)",
      Simplify[D[xx ph'[xx] - ph[xx], xx]/D[ph'[xx], xx] - xx] === 0];

  (* ---- [4] POISSON BRACKET: (x,y) canonical pair, flow = {., A_z} ---- *)
  chk["[4a] (x,y) is a CANONICAL pair: {x,y} = 1",
      Simplify[pb[x, y]] === 1 && Simplify[pb[x, x]] === 0 && Simplify[pb[y, y]] === 0];
  chk["[4b] Hamilton flow via bracket: {x,A_z} = dx/ds  and  {y,A_z} = dy/ds  (= B)",
      Simplify[pb[x, az[x, y]] - Bphys[[1]]] === 0 &&
      Simplify[pb[y, az[x, y]] - Bphys[[2]]] === 0];
  chk["[4c] conservation as a bracket: {A_z, A_z} = 0",
      Simplify[pb[az[x, y], az[x, y]]] === 0];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
