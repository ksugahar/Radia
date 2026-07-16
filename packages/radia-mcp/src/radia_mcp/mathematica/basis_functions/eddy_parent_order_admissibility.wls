#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  eddy_parent_order_admissibility.wls -- symbolic order ledger for deciding the
  parent HCurl order p in the topology-aware HCurl-VIM reduction.

  This file deliberately does NOT prove that "p=6 is good".  It proves the
  better statement:

      p is a parent-space admissibility order.

  The reduced model order is not p.  The reduced model order is the retained
  EVRS/CLN rank plus topology and surface branches:

      n_reduced = n_EVRS + n_cycle * dim(trace-polys) + n_surface.

  The parent HCurl order p must only be large enough to expose the spatial
  response classes before reduction:

      bulk current class:        curl(HCurl_p) is represented in HDiv/RT_p,
      bridge trace class:        conductor-conductor face traces up to r_bridge,
      SIBC/surface trace class:  conductor-air/exterior traces up to r_surface.

  Therefore:

      p >= max(r_bulk, r_bridge, r_surface).

  A surface-Omega branch with K = n x grad_Gamma Omega has one additional
  scalar-potential degree:

      deg(Omega) >= r_surface + 1.

  If the bridge class is enriched beyond one constant cycle per graph cycle,
  the number of bridge modes becomes

      n_bridge = cycle_rank * dim P_r(face).

  For the current coarse bridge-cycle basis r_bridge=0, so n_bridge=cycle_rank.
*)

ClearAll[
  DimP, DimFaceP, BridgeModeCount, ReducedModeCount,
  ParentOrderRequirement, AdmissibleParentOrderQ,
  SurfaceOmegaScalarDegreeRequirement, LedgerCompleteQ
];

DimP[dim_Integer, r_Integer] := Binomial[r + dim, dim];
DimFaceP[r_Integer] := DimP[2, r];

BridgeModeCount[cycleRank_Integer, bridgeTraceDegree_Integer] :=
  cycleRank * DimFaceP[bridgeTraceDegree];

ReducedModeCount[evrsRank_Integer, cycleRank_Integer, bridgeTraceDegree_Integer, surfaceModes_Integer] :=
  evrsRank + BridgeModeCount[cycleRank, bridgeTraceDegree] + surfaceModes;

ParentOrderRequirement[req_Association] := Max[
  req["bulkCurrentDegree"],
  req["bridgeTraceDegree"],
  req["surfaceTraceDegree"]
];

AdmissibleParentOrderQ[p_Integer, req_Association] := p >= ParentOrderRequirement[req];

SurfaceOmegaScalarDegreeRequirement[surfaceTraceDegree_Integer] := surfaceTraceDegree + 1;

LedgerCompleteQ[req_Association] := And @@ (
  KeyExistsQ[req, #] & /@
    {"bulkCurrentDegree", "bridgeTraceDegree", "surfaceTraceDegree"}
);

Module[
  {ok = True, chk, req, reqHighBridge, reqVeryHighBridge, pReq,
   nEVRS, cycleRank, surfaceModes,
   rBridge, gamma, bridgeModes, reducedModes},

  chk[name_, cond_] := (
    Print["  ", If[TrueQ[cond], "ok  ", "FAIL"], name];
    ok = ok && TrueQ[cond]
  );

  Print["eddy_parent_order_admissibility.wls self-test (p-order ledger):"];

  chk["dim P_r on a face is Binomial[r+2,2]",
      And @@ Table[DimFaceP[r] == Binomial[r + 2, 2], {r, 0, 5}]];
  chk["constant bridge trace keeps one mode per conductor graph cycle",
      BridgeModeCount[7, 0] == 7];
  chk["quadratic bridge trace would tensor each cycle by dim P_2(face)=6",
      BridgeModeCount[7, 2] == 42];

  req = <|
    "bulkCurrentDegree" -> 4,
    "bridgeTraceDegree" -> 0,
    "surfaceTraceDegree" -> 1
  |>;

  chk["ledger contains exactly the three spatial requirements used by p",
      LedgerCompleteQ[req]];
  chk["parent p requirement is max(bulk, bridge, surface), not CLN rank",
      ParentOrderRequirement[req] == 4];
  chk["p=4 is admissible for this ledger; p=3 is not",
      AdmissibleParentOrderQ[4, req] && ! AdmissibleParentOrderQ[3, req]];
  chk["p=6 can be conservative/admissible without being proven optimal",
      AdmissibleParentOrderQ[6, req] && ParentOrderRequirement[req] < 6];
  chk["surface-Omega scalar potential needs one degree above K trace",
      SurfaceOmegaScalarDegreeRequirement[req["surfaceTraceDegree"]] == 2];

  reqHighBridge = <|
    "bulkCurrentDegree" -> 4,
    "bridgeTraceDegree" -> 3,
    "surfaceTraceDegree" -> 1
  |>;
  chk["raising bridge trace degree raises the parent p requirement only if it is dominant",
      ParentOrderRequirement[reqHighBridge] == 4];
  reqVeryHighBridge = <|
    "bulkCurrentDegree" -> 4,
    "bridgeTraceDegree" -> 6,
    "surfaceTraceDegree" -> 1
  |>;
  chk["a high-order bridge trace can force p=6 symbolically",
      ParentOrderRequirement[reqVeryHighBridge] == 6];

  nEVRS = 24;
  cycleRank = 7;
  surfaceModes = 3;
  rBridge = 0;
  bridgeModes = BridgeModeCount[cycleRank, rBridge];
  reducedModes = ReducedModeCount[nEVRS, cycleRank, rBridge, surfaceModes];
  chk["current coarse p=6 smoke count is 24 + 7*dimP2(0) + 3 = 34",
      bridgeModes == 7 && reducedModes == 34];

  rBridge = 2;
  chk["if bridge faces need quadratic traces, mode count is no longer 34",
      ReducedModeCount[nEVRS, cycleRank, rBridge, surfaceModes] == 24 + 42 + 3];

  (* Symbolic statement: CLN/Krylov depth nEVRS is independent of p unless the
     spatial moment degree requested by that CLN construction changes. *)
  pReq = ParentOrderRequirement[<|
    "bulkCurrentDegree" -> rBulk,
    "bridgeTraceDegree" -> rBridgeSym,
    "surfaceTraceDegree" -> rSurface
  |>];
  chk["symbolic parent order is Max[rBulk, rBridge, rSurface]",
      pReq === Max[rBulk, rBridgeSym, rSurface]];

  gamma = Array[gCycle, cycleRank];
  chk["cycle coordinates are topology coefficients, not parent HCurl order coefficients",
      Length[gamma] == cycleRank && FreeQ[gamma, pReq]];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
  If[! ok, Exit[1]];
];
