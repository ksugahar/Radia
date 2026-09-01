#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  qd_cauer_ladder.wls -- Rutishauser's quotient-difference (QD) algorithm: moments
  of an admittance in, Cauer ladder coefficients out, with a rigorous account of
  how many digits survive.

  A Cauer ladder network IS a continued fraction.  Given the moment sequence of
  Y(z) at an expansion point,

        Y(z) = c0 + c1 z + c2 z^2 + ...,

  the QD scheme fills a lozenge table of e and q entries by two rhombus rules and
  reads the S-fraction (= ladder) coefficients off its top edge:

        Y(z) = c0 / (1 - q1 z / (1 - e1 z / (1 - q2 z / (1 - e2 z / ...)))).

  Two facts make this file worth having.

  1. TERMINATION IS STRUCTURE.  If Y is rational with exactly m poles, the e_m
     column is identically zero and the fraction closes at q_m.  The zero column
     is not a numerical accident -- it is the rank of the underlying Hankel
     matrix -- so an exact-arithmetic QD detects the true order of a network.

  2. THE RECURSION, NOT THE INPUT, DESTROYS THE DIGITS.  e_k^(n) is built from a
     DIFFERENCE of two nearby q's, so each rung cancels leading digits, and
     q_{k+1} then DIVIDES by that thinned quantity.  The moments can be known to
     full precision and the ladder coefficients still decay to noise.  This is
     the reason a Cauer ladder cannot simply be extended: the ceiling is
     arithmetic, not modelling.  The self-test below proves the asymmetry
     directly -- moment widths stay put while ladder widths grow monotonically.

  Interval arithmetic is the honest instrument here: Mathematica's Interval[]
  carries the uncertainty through the rhombus rules, and when a divisor's
  interval straddles zero the result is an unbounded interval rather than a
  crash, so the breakdown is reported instead of thrown.

  Index convention (N moments, c_0 ... c_{N-1}):
        q_k^(n)  defined for  n = 0 .. N-2k,     k = 1, 2, ...
        e_k^(n)  defined for  n = 0 .. N-2k-1,   k = 1, 2, ...
  with the seeds  e_0^(n) = 0  and  q_1^(n) = c_{n+1}/c_n.

  Rhombus rules:
        e_k^(n)     = q_k^(n+1) - q_k^(n) + e_{k-1}^(n+1)
        q_{k+1}^(n) = q_k^(n+1) * e_k^(n+1) / e_k^(n)

  Reference: H. Rutishauser, "Der Quotienten-Differenzen-Algorithmus",
  Z. angew. Math. Phys. 5 (1954) 233-251.
*)

BeginPackage["QDCauerLadder`"];

QDTable::usage =
  "QDTable[c] returns an Association keyed by {\"q\",k,n} and {\"e\",k,n} holding \
the quotient-difference table of the moment list c (c[[1]] is c_0).";
LadderCoefficients::usage =
  "LadderCoefficients[c] returns {c0, {q1, e1, q2, e2, ...}}, the S-fraction \
(Cauer ladder) coefficients read off the top edge of the QD table.";
LadderValue::usage =
  "LadderValue[c0, a, z] evaluates the S-fraction with head c0 and coefficient \
list a at z, by backward recurrence.";
LadderOrder::usage =
  "LadderOrder[c] returns the number of rungs before the first vanishing e \
column, i.e. the pole count QD detects, or Missing[] when no column vanishes.";
SignificantDigits::usage =
  "SignificantDigits[x] returns the number of correct decimal digits implied by \
an Interval, as Log10[|mid|/radius]; a non-Interval is treated as exact.";
MomentsFromPoles::usage =
  "MomentsFromPoles[res, poles, n] returns the first n moments of \
Sum[res_i/(1 - poles_i z)].";
WidenMoments::usage =
  "WidenMoments[c, eps] wraps each moment in a relative-width Interval, the \
input uncertainty a finite-precision moment computation would carry.";

Begin["`Private`"];

(* --------------------------------------------------------------- QD table *)

QDTable[c_List] := Module[{n = Length[c], t = <||>, kMax, qn, en},
  If[n < 2, Return[<||>]];
  (* seeds *)
  Do[t[{"e", 0, j}] = 0, {j, 0, n - 1}];
  Do[t[{"q", 1, j}] = c[[j + 2]]/c[[j + 1]], {j, 0, n - 2}];
  kMax = Floor[(n - 1)/2];
  Do[
   (* e_k from the difference rule; one shorter than the q row above it *)
   en = n - 2 k - 1;
   Do[t[{"e", k, j}] =
      t[{"q", k, j + 1}] - t[{"q", k, j}] + t[{"e", k - 1, j + 1}],
    {j, 0, en}];
   (* q_{k+1} from the quotient rule; stop cleanly on a zero divisor *)
   qn = n - 2 (k + 1);
   Do[
    If[TrueQ[t[{"e", k, j}] == 0],
     t[{"q", k + 1, j}] = Indeterminate,
     t[{"q", k + 1, j}] =
      t[{"q", k, j + 1}] t[{"e", k, j + 1}]/t[{"e", k, j}]],
    {j, 0, qn}],
   {k, 1, kMax}];
  t];

(* ------------------------------------------------- ladder coefficient edge *)

LadderCoefficients[c_List] := Module[{t = QDTable[c], a = {}, k = 1},
  While[KeyExistsQ[t, {"q", k, 0}] && ! MatchQ[t[{"q", k, 0}], Indeterminate],
   AppendTo[a, t[{"q", k, 0}]];
   If[! KeyExistsQ[t, {"e", k, 0}], Break[]];
   If[TrueQ[t[{"e", k, 0}] == 0], Break[]];
   AppendTo[a, t[{"e", k, 0}]];
   k++];
  {First[c], a}];

LadderValue[c0_, a_List, z_] := Module[{t = 1},
  Do[t = 1 - a[[i]] z/t, {i, Length[a], 1, -1}];
  c0/t];

LadderOrder[c_List] := Module[{t = QDTable[c], k = 1},
  While[KeyExistsQ[t, {"e", k, 0}],
   If[TrueQ[t[{"e", k, 0}] == 0], Return[k]];
   k++];
  Missing["NoVanishingColumn"]];

(* ------------------------------------------------------- digit bookkeeping *)

SignificantDigits[x_Interval] := Module[{lo, hi, mid, rad},
  {lo, hi} = First[x];
  If[! (NumericQ[lo] && NumericQ[hi]), Return[0.]];
  mid = (lo + hi)/2; rad = (hi - lo)/2;
  If[rad == 0, Return[Infinity]];
  If[mid == 0, Return[0.]];
  Max[0., N[Log10[Abs[mid]/rad]]]];
SignificantDigits[x_] := Infinity;

MomentsFromPoles[res_List, poles_List, n_Integer] :=
  Table[res . (poles^j), {j, 0, n - 1}];

WidenMoments[c_List, eps_] :=
  Interval[{# (1 - eps), # (1 + eps)}] & /@ c;

End[];
EndPackage[];

(* ============================================================ self-test *)

If[! TrueQ[$QDNoSelfTest],
 Module[{ok = True, chk, res, poles, c, order, c0, a, z, probe, exact,
   cw, tw, momDigits, ladDigits, growing},

  chk[name_, cond_] := (Print["  ", If[cond, "ok   ", "FAIL "], name];
    ok = ok && TrueQ[cond]);

  Print["qd_cauer_ladder.wls self-test"];

  (* --- a rational admittance: QD must find its order and reproduce it ---
     Exact rationals throughout, so these are equalities, not tolerances. *)
  res = {3/2, 1/5, 7/3};
  poles = {1/2, -1/3, 1/7};
  c = QDCauerLadder`MomentsFromPoles[res, poles, 9];

  order = QDCauerLadder`LadderOrder[c];
  chk["3-pole admittance: the e_3 column vanishes exactly (QD finds order 3)",
   order === 3];

  {c0, a} = QDCauerLadder`LadderCoefficients[c];
  chk["ladder closes at q_3, giving 5 coefficients q1 e1 q2 e2 q3",
   Length[a] == 5];

  z = 1/11;
  exact = Total[res/(1 - poles z)];
  chk["truncating at the zero column reproduces Y(z) EXACTLY (rational check)",
   Simplify[QDCauerLadder`LadderValue[c0, a, z] - exact] === 0];

  probe = 2/23;
  chk["and again at a second point (not a coincidence of one z)",
   Simplify[QDCauerLadder`LadderValue[c0, a, probe] -
      Total[res/(1 - poles probe)]] === 0];

  (* --- the digit-loss asymmetry: input clean, recursion lossy ---
     Give every moment the SAME tiny relative width, then compare how many
     digits survive in the moments against the ladder coefficients. *)
  cw = QDCauerLadder`WidenMoments[N[c, 40], 10^-16];
  tw = QDCauerLadder`QDTable[cw];

  momDigits = QDCauerLadder`SignificantDigits /@ cw;
  chk["every input moment carries ~16 digits, and they all still do",
   Min[momDigits] > 15 && Max[momDigits] - Min[momDigits] < 1];

  ladDigits =
   Table[QDCauerLadder`SignificantDigits[tw[{"q", k, 0}]], {k, 1, 3}];
  chk["ladder coefficient digits decrease rung by rung (q1 > q2 > q3)",
   ladDigits[[1]] > ladDigits[[2]] > ladDigits[[3]]];

  chk["the loss is real, not rounding: q3 has lost digits the moments kept",
   ladDigits[[3]] < Min[momDigits] - 1];

  growing =
   Table[QDCauerLadder`SignificantDigits[tw[{"e", k, 0}]], {k, 1, 2}];
  chk["the e column, a difference of neighbouring q's, loses digits too",
   growing[[2]] < growing[[1]]];

  (* --- breakdown is reported, never thrown ---
     A divisor whose interval contains zero is exactly where a floating or
     affine implementation aborts.  Interval arithmetic returns an unbounded
     interval instead, so the caller can see the ladder has ended. *)
  chk["division by an interval straddling zero yields an unbounded interval",
   MatchQ[Interval[{1, 2}]/Interval[{-1, 1}], _Interval] &&
    ! FreeQ[Interval[{1, 2}]/Interval[{-1, 1}], DirectedInfinity]];

  chk["a repeated pole (rank-deficient Hankel) does not raise, it terminates",
   Head[QDCauerLadder`LadderCoefficients[
      QDCauerLadder`MomentsFromPoles[{1, 1}, {1/2, 1/2}, 7]]] === List];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
  If[! ok, Exit[1]];
 ];
];
