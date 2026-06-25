#!/usr/bin/env wolframscript
(* ============================================================================
   MMMM analytic moment kernel -- closed form verified for ALL element face types:
   tet (4 triangular), pyramid (4 tri + 1 quad), wedge (2 tri + 3 quad), hex (6 quad).

   The kernel is per FLAT TRIANGLE; a quad face is fan-triangulated (2 triangles).
   So ONE closed form covers every element type.  Per source face at a target point:
       H_i  = Integral_T (obs-r')_i / |obs-r'|^3 dS'              (dipole field, = -grad I0)
       gH_ij= Integral_T [ delta_ij/r^3 - 3 (obs-r')_i (obs-r')_j / r^5 ] dS'  (quadrupole grad)
   with gH = grad_obs(H)  ->  closed-form gH = symbolic gradient of the closed-form H.

   Element-type usage of the kernel:
     tet     : 3-DOF dipole (MMM) -- needs H only (no quadrupole); faces are triangles.
     pyramid : 1 mono + 3 dipole + 1 quad eigenmode -> needs H AND gH.
     wedge   : 1 mono + 3 dipole + 1 quad eigenmode -> needs H AND gH.
     hex     : 1 mono + 3 dipole + 2 quad eigenmodes -> needs H AND gH.

   This .wls derives gH (symbolic gradient of H) and verifies analytic == 64pt Gauss for
   EACH element type's actual face set (triangular and/or quad faces).  Self-testing.
   ============================================================================ *)

Off[General::stop];
SetOptions[$Output, FormatType -> OutputForm];
gp = {0.01985507175123185, 0.10166676129318665, 0.23723379504183550, 0.40828267875217510,
      0.59171732124782490, 0.76276620495816450, 0.89833323870681340, 0.98014492824876820};
gw = {0.05061426814518815, 0.11119051722668725, 0.15685332293894365, 0.18134189168918100,
      0.18134189168918100, 0.15685332293894365, 0.11119051722668725, 0.05061426814518815};

(* production reference: bilinear-quad 8x8 Gauss over a quad {V0,V1,V2,V3} (a tri = degenerate V3=V2) *)
gaussFaceHG[obs_, V0_,V1_,V2_,V3_] := Module[
  {H={0.,0.,0.}, gH={0.,0.,0.,0.,0.,0.}, u,v,wuv,a0,a1,a2,a3,P,Tu,Tv,jv,dA,d,dx,dy,dz,r2,ir,ir3,ir5,c3,c5},
  Do[ u=gp[[iu]]; v=gp[[iv]]; wuv=gw[[iu]]gw[[iv]];
      a0=(1-u)(1-v); a1=u(1-v); a2=u v; a3=(1-u)v;
      P=a0 V0+a1 V1+a2 V2+a3 V3; Tu=(1-v)(V1-V0)+v(V2-V3); Tv=(1-u)(V3-V0)+u(V2-V1);
      jv=Cross[Tu,Tv]; dA=Sqrt[jv.jv]wuv; d=obs-P; {dx,dy,dz}=d; r2=d.d;
      ir=1/Sqrt[r2]; ir3=ir/r2; ir5=ir3/r2; c3=ir3 dA; c5=ir5 dA; H+=d c3;
      gH[[1]]+=c3-3dx dx c5; gH[[2]]+=c3-3dy dy c5; gH[[3]]+=c3-3dz dz c5;
      gH[[4]]+=-3dx dy c5;   gH[[5]]+=-3dx dz c5;   gH[[6]]+=-3dy dz c5;
    ,{iu,1,8},{iv,1,8}];
  {H,gH}];

(* analytic H (van Oosterom-Strackee, sigma=1, no 4pi); symbolic-obs capable *)
analyticH[ox_,oy_,oz_, V0_,V1_,V2_] := Module[
  {obs,e1,e2,e1len,ba,cc,clen,bc,bb,XYx,XYy,DS,AM,XD,YD,dd,E1,E2,E3,Xv,Yv,Hv,Ev,Rv,RM,RP,AL,HH1,HH2,HH3,AT,BT,l},
  obs={ox,oy,oz}; e1=V1-V0; e2=V2-V0; e1len=Sqrt[e1.e1]; ba=e1/e1len;
  cc=Cross[e1,e2]; clen=Sqrt[cc.cc]; bc=cc/clen; bb=Cross[bc,ba];
  XYx={0,e1len,e2.ba}; XYy={0,0,e2.bb}; DS=AM=XD=YD=ConstantArray[0,3];
  Do[ l=Mod[j,3]+1; Module[{dxj=XYx[[l]]-XYx[[j]], dyj=XYy[[l]]-XYy[[j]]},
        If[Abs[dxj] < 1.*^-12, dxj = If[dxj >= 0, 1.*^-12, -1.*^-12]];   (* clamp vertical local edge, matches C *)
        DS[[j]]=Sqrt[dxj^2+dyj^2]; AM[[j]]=dyj/dxj; XD[[j]]=-dxj/DS[[j]]; YD[[j]]=dyj/DS[[j]]];,{j,3}];
  dd=obs-V0; E1=dd.ba; E2=dd.bb; E3=dd.bc;
  Xv=Table[E1-XYx[[j]],{j,3}]; Yv=Table[E2-XYy[[j]],{j,3}];
  Hv=Table[Yv[[j]]Xv[[j]],{j,3}]; Ev=Table[E3^2+Xv[[j]]^2,{j,3}];
  Rv=Table[Sqrt[Xv[[j]]^2+Yv[[j]]^2+E3^2],{j,3}];
  RM=Table[Rv[[j]]+Rv[[Mod[j,3]+1]]-DS[[j]],{j,3}]; RP=Table[Rv[[j]]+Rv[[Mod[j,3]+1]]+DS[[j]],{j,3}];
  AL=Table[Log[RM[[j]]/RP[[j]]],{j,3}];
  HH1=-Sum[YD[[j]]AL[[j]],{j,3}]; HH2=-Sum[XD[[j]]AL[[j]],{j,3}];
  AT=Table[(AM[[j]]Ev[[j]]-Hv[[j]])/(E3 Rv[[j]]),{j,3}];
  BT=Table[(AM[[j]]Ev[[Mod[j,3]+1]]-Hv[[Mod[j,3]+1]])/(E3 Rv[[Mod[j,3]+1]]),{j,3}];
  HH3=-Sum[ArcTan[AT[[j]]],{j,3}]+Sum[ArcTan[BT[[j]]],{j,3}];
  HH1 ba+HH2 bb+HH3 bc];

(* analytic face = fan {V0,V1,V2}+{V0,V2,V3}; gH = Jacobian of the H-sum.  Degenerate 2nd tri (V3=V2) -> 0. *)
analyticFaceHG[obs_, V0_,V1_,V2_,V3_] := Module[{ox,oy,oz,Hs,M,Hval,Mval},
  Hs = analyticH[ox,oy,oz,V0,V1,V2] + If[Norm[Cross[V2-V0,V3-V0]] < 1.*^-12, {0,0,0}, analyticH[ox,oy,oz,V0,V2,V3]];
  M = D[Hs,{{ox,oy,oz}}];
  Hval = Hs /. {ox->obs[[1]],oy->obs[[2]],oz->obs[[3]]};
  Mval = M /. {ox->obs[[1]],oy->obs[[2]],oz->obs[[3]]};
  {Hval, {Mval[[1,1]],Mval[[2,2]],Mval[[3,3]],Mval[[1,2]],Mval[[1,3]],Mval[[2,3]]}}];

tri[a_,b_,c_] := {a,b,c,c};   (* triangular face as a degenerate quad *)

(* reference elements (unit-ish, general position) and their face vertex lists *)
tetV = {{0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
tetFaces = {tri[tetV[[1]],tetV[[2]],tetV[[3]]], tri[tetV[[1]],tetV[[2]],tetV[[4]]],
            tri[tetV[[1]],tetV[[3]],tetV[[4]]], tri[tetV[[2]],tetV[[3]],tetV[[4]]]};

pyrV = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.},{0.5,0.5,1.}};  (* 4 base + apex *)
pyrFaces = {{pyrV[[1]],pyrV[[2]],pyrV[[3]],pyrV[[4]]},               (* quad base *)
            tri[pyrV[[1]],pyrV[[2]],pyrV[[5]]], tri[pyrV[[2]],pyrV[[3]],pyrV[[5]]],
            tri[pyrV[[3]],pyrV[[4]],pyrV[[5]]], tri[pyrV[[4]],pyrV[[1]],pyrV[[5]]]};

wedV = {{0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.},{1.,0.,1.},{0.,1.,1.}};
wedFaces = {tri[wedV[[1]],wedV[[2]],wedV[[3]]], tri[wedV[[4]],wedV[[5]],wedV[[6]]],   (* 2 tri caps *)
            {wedV[[1]],wedV[[2]],wedV[[5]],wedV[[4]]}, {wedV[[2]],wedV[[3]],wedV[[6]],wedV[[5]]},
            {wedV[[3]],wedV[[1]],wedV[[4]],wedV[[6]]}};                                (* 3 quad sides *)

hexV = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.},{0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.}};
hexFaces = {{hexV[[1]],hexV[[2]],hexV[[3]],hexV[[4]]}, {hexV[[5]],hexV[[6]],hexV[[7]],hexV[[8]]},
            {hexV[[1]],hexV[[2]],hexV[[6]],hexV[[5]]}, {hexV[[2]],hexV[[3]],hexV[[7]],hexV[[6]]},
            {hexV[[3]],hexV[[4]],hexV[[8]],hexV[[7]]}, {hexV[[4]],hexV[[1]],hexV[[5]],hexV[[8]]}};

obs = {1.7, 1.3, 0.9};   (* external target point, off every face plane *)

Print["================================================================"];
Print["MMMM analytic kernel: closed-form H + gH vs 64pt Gauss, per element type"];
Print["================================================================"];

(* ---- gH = grad_obs(H) is the Hessian of the single-layer potential: symmetric + traceless ---- *)
Module[{ox,oy,oz, M, Mn},
  M  = D[analyticH[ox,oy,oz, {0.1,0.2,0.3},{1.,0.3,0.1},{0.2,1.1,0.4}], {{ox,oy,oz}}];
  Mn = M /. {ox->1.7, oy->1.3, oz->0.9};
  Print["  gH = grad_obs(H), the Hessian -grad grad I0 (I0 = single-layer potential):"];
  Print["    symmetry  |M - M^T| = ", ScientificForm[Norm[Mn - Transpose[Mn]], 3], "   (Hessian => exactly symmetric)"];
  Print["    traceless |Tr M|    = ", ScientificForm[Abs[Tr[Mn]], 3], "   (Laplace 1/r => div H = 0, auto)"];
];
Print["----------------------------------------------------------------"];
check[name_, faces_] := Module[{hrel=0., grel=0.},
  Do[ Module[{a=analyticFaceHG[obs, Sequence@@f], g=gaussFaceHG[obs, Sequence@@f]},
        hrel = Max[hrel, Norm[a[[1]]-g[[1]]]/Norm[g[[1]]]];
        grel = Max[grel, Norm[a[[2]]-g[[2]]]/Norm[g[[2]]]]; ], {f, faces}];
  Print["  ", name, ": ", Length[faces], " faces   max H rel=", ScientificForm[hrel,3],
        "   max gH rel=", ScientificForm[grel,3], "   -> ", If[hrel<1.*^-4 && grel<1.*^-4, "PASS", "FAIL"]];
  {hrel, grel}];

rt = check["tet    ", tetFaces];
rp = check["pyramid", pyrFaces];
rw = check["wedge  ", wedFaces];
rh = check["hex    ", hexFaces];
Print["----------------------------------------------------------------"];
Print["All element-type faces (tri + quad) use ONE per-triangle closed form."];
Print["Expected: H rel & gH rel ~1e-6..1e-9 (= 8x8 Gauss error; the closed form is EXACT)."];
allok = And @@ (#[[1]]<1.*^-4 && #[[2]]<1.*^-4 & /@ {rt,rp,rw,rh});
Print["SELF-TEST: ", If[allok, "PASS (all element types)", "FAIL"]];
Print["================================================================"];
