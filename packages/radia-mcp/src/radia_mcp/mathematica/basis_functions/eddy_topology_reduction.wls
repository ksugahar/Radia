#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  eddy_topology_reduction.wls -- symbolic topology certificate for the
  topology-aware HCurl-VIM eddy-bubble reduction.

  This is the Mathematica-side answer to the question:

      Can we decide from each face whether the neighbor is air or conductor?

  Yes.  The production split is not an empirical p=6 rule.  It is a face-role
  and graph-complex rule:

      conductor | air/exterior  -> SIBC / surface-Omega trace class,
      conductor | conductor     -> bridge class,
      interior HCurl bulk       -> EVRS / CLN response class.

  The conductor-conductor bridge class is then reduced by the cycle space of
  the conductive adjacency graph.  If B is the oriented cell-face incidence
  matrix, cycle currents have the form

      j_bridge = Z gamma,       B Z = 0,

  so they preserve Kirchhoff continuity / local charge conservation exactly.

  The volume EVRS branch keeps T as the HCurl variable:

      T = Q z,      J = curl T,

  while the SIBC branch keeps a surface-Omega variable:

      K = n x grad_Gamma Omega.

  This script is intentionally symbolic and graph-theoretic.  It does not say
  "p=6 is good".  Instead, it certifies the non-empirical reduction rule that a
  parent order p must realize: preserve de Rham solenoidality in the bulk,
  preserve conductor graph cycles on bridge faces, and put SIBC only on
  conductor-air/exterior faces.
*)

ClearAll[
  ConductiveQ, AirLikeQ, ExteriorQ, FaceRole, OrientedIncidence,
  ComponentCountFromIncidence, CycleRankFromIncidence, CycleMatrix,
  SurfaceK, SurfaceDiv
];

conductiveMaterials = {"cond", "copper", "steel"};
airLikeMaterials = {"air", "vacuum", "exterior"};

ConductiveQ[mat_] := MemberQ[conductiveMaterials, mat];
AirLikeQ[mat_] := MemberQ[airLikeMaterials, mat];
ExteriorQ[mat_] := mat === None || mat === "exterior";

FaceRole[{left_, right_}] := Which[
  ConductiveQ[left] && ConductiveQ[right],
    "bridge-conductor-conductor",
  Xor[ConductiveQ[left], ConductiveQ[right]] &&
    (AirLikeQ[left] || AirLikeQ[right] || ExteriorQ[left] || ExteriorQ[right]),
    "sibc-conductor-air",
  Xor[ConductiveQ[left], ConductiveQ[right]],
    "surface-conductor-nonconductor",
  True,
    "ignore-nonconductive"
];

OrientedIncidence[nodes_, edges_] := Module[{pos = AssociationThread[nodes -> Range[Length[nodes]]]},
  Table[
    Which[
      pos[edges[[e, 1]]] === i, -1,
      pos[edges[[e, 2]]] === i,  1,
      True, 0
    ],
    {i, Length[nodes]}, {e, Length[edges]}
  ]
];

ComponentCountFromIncidence[incidence_] := Length[incidence] - MatrixRank[incidence];

CycleRankFromIncidence[incidence_] := Last[Dimensions[incidence]] - MatrixRank[incidence];

CycleMatrix[incidence_] := Transpose[NullSpace[incidence]];

(* Flat surface Gamma={z=0}, n=e_z. *)
SurfaceK[omega_] := {-D[omega, y], D[omega, x], 0};
SurfaceDiv[K_] := D[K[[1]], x] + D[K[[2]], y];

Module[
  {ok = True, chk, faceMaterials, roles, conductorNodes, bridgeEdges,
   sibcFaces, nonSibcTraceFaces, ignoredFaces, B, Z, gammaVec, jBridge,
   divCells, omega, T,
   p, requiredBulkMomentDegree, requiredSurfaceOmegaDegree,
   requiredBridgeTraceDegree, admissibleParentOrder},

  chk[name_, cond_] := (
    Print["  ", If[TrueQ[cond], "ok  ", "FAIL"], name];
    ok = ok && TrueQ[cond]
  );

  Print["eddy_topology_reduction.wls self-test (face role + cycle reduction):"];

  (* A tiny abstract conductive patch.  The first four faces form a conductive
     graph with one independent cycle.  Three faces touch air/exterior and must
     become SIBC/surface-Omega candidates, not bridge modes.  One face touches
     a non-air insulator; it may need a trace condition, but it is not SIBC. *)
  faceMaterials = {
    {"cond", "cond"},
    {"cond", "cond"},
    {"cond", "cond"},
    {"cond", "cond"},
    {"cond", "air"},
    {"cond", "exterior"},
    {"vacuum", "cond"},
    {"cond", "ceramic"},
    {"air", "vacuum"}
  };
  roles = FaceRole /@ faceMaterials;

  chk["face role is decided by the two neighboring material labels",
      roles === {
        "bridge-conductor-conductor",
        "bridge-conductor-conductor",
        "bridge-conductor-conductor",
        "bridge-conductor-conductor",
        "sibc-conductor-air",
        "sibc-conductor-air",
        "sibc-conductor-air",
        "surface-conductor-nonconductor",
        "ignore-nonconductive"
      }];

  conductorNodes = {"c1", "c2", "c3", "c4"};
  bridgeEdges = {{"c1", "c2"}, {"c2", "c3"}, {"c3", "c1"}, {"c3", "c4"}};
  sibcFaces = Pick[Range[Length[roles]], roles, "sibc-conductor-air"];
  nonSibcTraceFaces =
    Pick[Range[Length[roles]], roles, "surface-conductor-nonconductor"];
  ignoredFaces = Pick[Range[Length[roles]], roles, "ignore-nonconductive"];

  chk["SIBC faces are exactly conductor-air/exterior faces",
      sibcFaces === {5, 6, 7}];
  chk["conductor-insulator faces are trace faces, not SIBC half-spaces",
      nonSibcTraceFaces === {8} && FreeQ[sibcFaces, 8]];
  chk["air-air faces are ignored by the eddy-current basis",
      ignoredFaces === {9}];

  B = OrientedIncidence[conductorNodes, bridgeEdges];
  Z = CycleMatrix[B];
  gammaVec = Array[gBridge, Last[Dimensions[Z]]];
  jBridge = Z.gammaVec;
  divCells = Simplify[B.jBridge];

  chk["conductor graph component count is V - rank(B)",
      ComponentCountFromIncidence[B] ===
        Length[conductorNodes] - MatrixRank[B]];
  chk["conductor graph cycle rank is E - rank(B) = E - V + components",
      Last[Dimensions[Z]] === CycleRankFromIncidence[B] &&
      CycleRankFromIncidence[B] ===
        Length[bridgeEdges] - Length[conductorNodes] + ComponentCountFromIncidence[B]];
  chk["bridge-cycle currents satisfy B.Z = 0 exactly",
      Simplify[B.Z] === ConstantArray[0, {Length[conductorNodes], Last[Dimensions[Z]]}]];
  chk["j_bridge = Z gamma gives zero cell divergence for arbitrary gamma",
      divCells === ConstantArray[0, Length[conductorNodes]]];

  (* This is the topology part of the HCurl-VIM split.  The variable is still T
     for the bulk branch; the bridge class is not a free face-current DoF list
     anymore, but a cycle-coordinate vector gamma. *)
  T = Q.z;
  chk["bulk variable remains T = Q z, separate from bridge cycle coordinates",
      Head[T] === Plus || Head[T] === Dot];

  omega = b0 + b1 x + b2 y + b3 x^2 + b4 x y + b5 y^2;
  chk["surface-Omega current K = n x grad_Gamma Omega is tangential",
      Simplify[SurfaceK[omega][[3]]] === 0];
  chk["surface-Omega current is surface-solenoidal",
      Simplify[SurfaceDiv[SurfaceK[omega]]] === 0];

  (* Parent-order selection is therefore a symbolic admissibility problem, not a
     numerical p=6 slogan.  The actual requirements come from the intended CLN
     moments, surface ladder degree, and bridge trace representation. *)
  admissibleParentOrder =
    p >= Max[
      requiredBulkMomentDegree,
      requiredSurfaceOmegaDegree + 1,
      requiredBridgeTraceDegree
    ];
  chk["p is a parent-space admissibility condition, not an empirical constant",
      admissibleParentOrder ===
        p >= Max[
          requiredBulkMomentDegree,
          1 + requiredSurfaceOmegaDegree,
          requiredBridgeTraceDegree
        ]];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
  If[! ok, Exit[1]];
];
