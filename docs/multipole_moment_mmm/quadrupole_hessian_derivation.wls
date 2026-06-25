#!/usr/bin/env wolframscript
(* ============================================================================
   MMMM quadrupole field-gradient closed-form derivation + verification.

   Production kernel CentroidFieldGradFromFace (rad_interaction.cpp:2210) computes,
   per unit-charge source face, at target centroid obs:
       H_i  = Integral_T (obs-r')_i / |obs-r'|^3 dS'              (= -grad I0)
       gH_ij= Integral_T [ delta_ij/r^3 - 3 (obs-r')_i (obs-r')_j / r^5 ] dS'
   BOTH via an 8x8 (=64 pt) Gauss bilinear-quad quadrature.

   KEY IDENTITY:  gH_ij = d/dobs_j [ Integral_T (obs-r')_i/r^3 dS' ] = dH_i/dobs_j.
   So gH = Jacobian of H w.r.t. obs.  H is ALREADY closed form
   (FieldFromChargedTriangleLocal, van Oosterom-Strackee).  Therefore the
   closed-form gH = symbolic gradient of that H.  This .wls:
     (A) ports the analytic H and the production 64pt-Gauss to WL;
     (B) verifies analytic H == Gauss H (the existing kernel is the right integral);
     (C) derives gH_closed = D[analyticH, obs] and verifies == Gauss gH,
         + symmetry + tracelessness (free physical checks);
     (D) emits a local-frame closed form (CForm) for the C implementation.
   ============================================================================ *)

Off[General::stop];
SetOptions[$Output, FormatType -> OutputForm];
pr[args__] := Print[args];

(* ---- 8-pt Gauss-Legendre on [0,1], matching MOM_GP/MOM_GW in rad_interaction.cpp ---- *)
gp = {0.01985507175123185, 0.10166676129318665, 0.23723379504183550, 0.40828267875217510,
      0.59171732124782490, 0.76276620495816450, 0.89833323870681340, 0.98014492824876820};
gw = {0.05061426814518815, 0.11119051722668725, 0.15685332293894365, 0.18134189168918100,
      0.18134189168918100, 0.15685332293894365, 0.11119051722668725, 0.05061426814518815};

(* ---- (A1) production reference: degenerate-quad bilinear 8x8 Gauss, isSelf (no center charge),
        IMA off.  Mirrors accumQG(V4={V0,V1,V2,V2}).  Returns {H(3), gH(6)} gH=(xx,yy,zz,xy,xz,yz). *)
gaussTriHG[obs_, V0_, V1_, V2_] := Module[
  {V3 = V2, H = {0.,0.,0.}, gH = {0.,0.,0.,0.,0.,0.}, u,v,wuv,a0,a1,a2,a3,P,Tu,Tv,jv,dA,
   dvec,dx,dy,dz,r2,invr,invr3,invr5,c3,c5},
  Do[
    u = gp[[iu]]; v = gp[[iv]]; wuv = gw[[iu]]*gw[[iv]];
    a0 = (1-u)(1-v); a1 = u (1-v); a2 = u v; a3 = (1-u) v;
    P  = a0 V0 + a1 V1 + a2 V2 + a3 V3;
    Tu = (1-v)(V1-V0) + v (V2-V3);
    Tv = (1-u)(V3-V0) + u (V2-V1);
    jv = Cross[Tu, Tv];
    dA = Sqrt[jv.jv]*wuv;
    dvec = obs - P; {dx,dy,dz} = dvec;
    r2 = dvec.dvec; invr = 1/Sqrt[r2]; invr3 = invr/r2; invr5 = invr3/r2;
    c3 = invr3 dA; c5 = invr5 dA;
    H += dvec c3;
    gH[[1]] += c3 - 3 dx dx c5; gH[[2]] += c3 - 3 dy dy c5; gH[[3]] += c3 - 3 dz dz c5;
    gH[[4]] += -3 dx dy c5;     gH[[5]] += -3 dx dz c5;     gH[[6]] += -3 dy dz c5;
    , {iu,1,8}, {iv,1,8}];
  {H, gH}];

(* ---- (A2) analytic H: WL port of FieldFromChargedTriangleLocal (sigma=1, no 4pi).
        Symbolic-obs capable.  Triangle {V0,V1,V2}; obs={ox,oy,oz}. ---- *)
analyticH[ox_, oy_, oz_, V0_, V1_, V2_] := Module[
  {obs, e1, e2, e1len, ba, cc, clen, bc, bb, XYx, XYy, dxj, dyj, DS, AM, XD, YD,
   dd, EE1, EE2, EE3, Xv, Yv, Hv, Ev, Rv, RM, RP, AL, HH1, HH2, HH3, AT, BT, l},
  obs = {ox, oy, oz};
  e1 = V1 - V0; e2 = V2 - V0;
  e1len = Sqrt[e1.e1];   ba = e1/e1len;
  cc = Cross[e1, e2];    clen = Sqrt[cc.cc];   bc = cc/clen;
  bb = Cross[bc, ba];
  XYx = {0, e1len, e2.ba};  XYy = {0, 0, e2.bb};
  DS = AM = XD = YD = ConstantArray[0, 3];
  Do[ l = Mod[j,3]+1;
      dxj = XYx[[l]] - XYx[[j]];  dyj = XYy[[l]] - XYy[[j]];
      DS[[j]] = Sqrt[dxj^2 + dyj^2];  AM[[j]] = dyj/dxj;
      XD[[j]] = -dxj/DS[[j]];  YD[[j]] = dyj/DS[[j]];
      , {j,1,3}];
  dd = obs - V0;  EE1 = dd.ba;  EE2 = dd.bb;  EE3 = dd.bc;
  Xv = Table[EE1 - XYx[[j]], {j,3}];
  Yv = Table[EE2 - XYy[[j]], {j,3}];
  Hv = Table[Yv[[j]] Xv[[j]], {j,3}];
  Ev = Table[EE3^2 + Xv[[j]]^2, {j,3}];
  Rv = Table[Sqrt[Xv[[j]]^2 + Yv[[j]]^2 + EE3^2], {j,3}];
  RM = Table[Rv[[j]] + Rv[[Mod[j,3]+1]] - DS[[j]], {j,3}];
  RP = Table[Rv[[j]] + Rv[[Mod[j,3]+1]] + DS[[j]], {j,3}];
  AL = Table[Log[RM[[j]]/RP[[j]]], {j,3}];
  HH1 = -Sum[YD[[j]] AL[[j]], {j,3}];
  HH2 = -Sum[XD[[j]] AL[[j]], {j,3}];
  AT = Table[(AM[[j]] Ev[[j]]            - Hv[[j]])           /(EE3 Rv[[j]]), {j,3}];
  BT = Table[(AM[[j]] Ev[[Mod[j,3]+1]]   - Hv[[Mod[j,3]+1]])  /(EE3 Rv[[Mod[j,3]+1]]), {j,3}];
  HH3 = -Sum[ArcTan[AT[[j]]],{j,3}] + Sum[ArcTan[BT[[j]]],{j,3}];
  HH1 ba + HH2 bb + HH3 bc
];

(* ---- test geometries (general position; obs well off the source plane) ---- *)
tris = {
  {{0.,0.,0.}, {1.0,0.,0.}, {0.3,0.9,0.}},                    (* axis-aligned flat *)
  {{0.1,0.2,0.3}, {1.0,0.3,0.1}, {0.2,1.1,0.4}}               (* tilted general position *)
};
obss = {
  {0.40, 0.20, 0.70}, {1.50,-0.60, 0.90}, {-0.50,0.40,-1.10}, {0.60,0.50,1.30}
};

pr["==================================================================="];
pr["MMMM Hessian closed-form: derive gH = grad_obs(analytic H), verify vs 64pt Gauss"];
pr["==================================================================="];

ClearAll[ox, oy, oz];
maxHrel = 0.; maxGrel = 0.; maxSym = 0.; maxTr = 0.;
Do[
  Module[{V0,V1,V2, Msym},
    {V0,V1,V2} = tris[[ti]];
    (* symbolic Jacobian of analytic H w.r.t. obs -> closed-form gH tensor *)
    Msym = D[analyticH[ox,oy,oz, V0,V1,V2], {{ox,oy,oz}}];
    Do[
      Module[{obs, gauss, Hg, gHg, Hc, Mc, gHc, hrel, grel, symerr, trerr},
        obs = obss[[oi]];
        gauss = gaussTriHG[obs, V0, V1, V2];  Hg = gauss[[1]];  gHg = gauss[[2]];
        Hc = analyticH[obs[[1]],obs[[2]],obs[[3]], V0,V1,V2];   (* analytic H value *)
        Mc = Msym /. {ox->obs[[1]], oy->obs[[2]], oz->obs[[3]]};
        (* pack closed-form tensor Mc[[i,j]] into (xx,yy,zz,xy,xz,yz) *)
        gHc = {Mc[[1,1]], Mc[[2,2]], Mc[[3,3]], Mc[[1,2]], Mc[[1,3]], Mc[[2,3]]};
        hrel = Norm[Hc - Hg]/Norm[Hg];
        grel = Norm[gHc - gHg]/Norm[gHg];
        symerr = Norm[Mc - Transpose[Mc]];        (* exact symmetry expected *)
        trerr  = Abs[Tr[Mc]];                       (* exact tracelessness expected (off source) *)
        maxHrel = Max[maxHrel, hrel]; maxGrel = Max[maxGrel, grel];
        maxSym = Max[maxSym, symerr]; maxTr = Max[maxTr, trerr];
        pr[StringForm["tri`1` obs=`2`  |H rel=`3`  gH rel=`4`  sym=`5`  trace=`6`",
            ti, NumberForm[obs,3],
            ScientificForm[hrel,3], ScientificForm[grel,3],
            ScientificForm[symerr,3], ScientificForm[trerr,3]]];
      ], {oi, Length[obss]}];
  ], {ti, Length[tris]}];

pr["-------------------------------------------------------------------"];
pr[StringForm["MAX over all cases:  H rel=`1`   gH rel=`2`   sym=`3`   trace=`4`",
    ScientificForm[maxHrel,3], ScientificForm[maxGrel,3],
    ScientificForm[maxSym,3], ScientificForm[maxTr,3]]];
pr["Expected: H rel & gH rel ~1e-7..1e-9 (8x8 Gauss accuracy); sym & trace ~machine 0."];
pr["==================================================================="];

(* ---- (D) local-frame closed form for the C implementation.
        Local triangle V0=(0,0,0), V1=(p,0,0), V2=(q,s,0); obs=(e1,e2,e3).
        Hloc=(HH1,HH2,HH3); Gloc = D[Hloc,{{e1,e2,e3}}].  CForm in (p,q,s,e1,e2,e3). ---- *)
ClearAll[p,q,s,e1,e2,e3];
locXYx = {0, p, q};  locXYy = {0, 0, s};
locDS = locAM = locXD = locYD = ConstantArray[0,3];
Do[ With[{l = Mod[j,3]+1},
      With[{dxj = locXYx[[l]]-locXYx[[j]], dyj = locXYy[[l]]-locXYy[[j]]},
        locDS[[j]] = Sqrt[dxj^2+dyj^2]; locAM[[j]] = dyj/dxj;
        locXD[[j]] = -dxj/locDS[[j]]; locYD[[j]] = dyj/locDS[[j]]; ]], {j,1,3}];
locX = Table[e1 - locXYx[[j]], {j,3}];
locY = Table[e2 - locXYy[[j]], {j,3}];
locHh= Table[locY[[j]] locX[[j]], {j,3}];
locE = Table[e3^2 + locX[[j]]^2, {j,3}];
locR = Table[Sqrt[locX[[j]]^2 + locY[[j]]^2 + e3^2], {j,3}];
locRM= Table[locR[[j]] + locR[[Mod[j,3]+1]] - locDS[[j]], {j,3}];
locRP= Table[locR[[j]] + locR[[Mod[j,3]+1]] + locDS[[j]], {j,3}];
locAL= Table[Log[locRM[[j]]/locRP[[j]]], {j,3}];
HH1loc = -Sum[locYD[[j]] locAL[[j]], {j,3}];
HH2loc = -Sum[locXD[[j]] locAL[[j]], {j,3}];
locAT = Table[(locAM[[j]] locE[[j]]          - locHh[[j]])         /(e3 locR[[j]]), {j,3}];
locBT = Table[(locAM[[j]] locE[[Mod[j,3]+1]] - locHh[[Mod[j,3]+1]])/(e3 locR[[Mod[j,3]+1]]), {j,3}];
HH3loc = -Sum[ArcTan[locAT[[j]]],{j,3}] + Sum[ArcTan[locBT[[j]]],{j,3}];
Hloc = {HH1loc, HH2loc, HH3loc};
Gloc = D[Hloc, {{e1,e2,e3}}];

(* numeric cross-check that the local-frame Gloc, transformed to global, equals Msym *)
pr["Local-frame Gloc derived (3x3). Spot-check local vs global path on tilted tri/obs:"];
Module[{V0,V1,V2,obs, e1g,e2g,e1len,ba,cc,clen,bc,bb,dd,pv,qv,sv,ev1,ev2,ev3,
        GlocN, B, Gglob, Mglob, rel},
  {V0,V1,V2} = tris[[2]];  obs = obss[[4]];
  e1g = V1-V0; e2g = V2-V0; e1len = Sqrt[e1g.e1g]; ba = e1g/e1len;
  cc = Cross[e1g,e2g]; clen = Sqrt[cc.cc]; bc = cc/clen; bb = Cross[bc,ba];
  pv = e1len; qv = e2g.ba; sv = e2g.bb;
  dd = obs - V0; ev1 = dd.ba; ev2 = dd.bb; ev3 = dd.bc;
  GlocN = Gloc /. {p->pv,q->qv,s->sv,e1->ev1,e2->ev2,e3->ev3};
  B = Transpose[{ba,bb,bc}];                 (* columns = local basis vectors *)
  Gglob = B . GlocN . Transpose[B];          (* congruence transform of rank-2 tensor *)
  Mglob = D[analyticH[ox,oy,oz,V0,V1,V2],{{ox,oy,oz}}] /. {ox->obs[[1]],oy->obs[[2]],oz->obs[[3]]};
  rel = Norm[Gglob - Mglob]/Norm[Mglob];
  pr[StringForm["  local->global vs direct-global gH:  rel=`1`  (expect ~machine 0)",
      ScientificForm[rel,3]]];
];

pr["==================================================================="];
pr["CForm of local-frame Gloc components (symmetric; emit xx,yy,zz,xy,xz,yz):"];
pr["--- Gxx ---"]; pr[CForm[Gloc[[1,1]]]];
pr["--- Gyy ---"]; pr[CForm[Gloc[[2,2]]]];
pr["--- Gzz ---"]; pr[CForm[Gloc[[3,3]]]];
pr["--- Gxy ---"]; pr[CForm[Gloc[[1,2]]]];
pr["--- Gxz ---"]; pr[CForm[Gloc[[1,3]]]];
pr["--- Gyz ---"]; pr[CForm[Gloc[[2,3]]]];
pr["=== DONE ==="];
