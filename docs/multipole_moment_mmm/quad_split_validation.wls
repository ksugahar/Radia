#!/usr/bin/env wolframscript
(* Validate the hex-face handling: analytic 2-triangle-fan sum vs production bilinear-quad 8x8 Gauss.
   Planar quad -> 2 flat triangles tile it EXACTLY -> match to Gauss accuracy.
   Warped quad -> flat triangulation differs from bilinear surface by O(warp); quantify it.
   This is the SAME modeling choice the existing analytic FIELD eval already makes (Compute3x3BlockFast
   triangulates faces), so switching the moment kernel to triangulation is a consistency improvement. *)

Off[General::stop];
gp = {0.01985507175123185, 0.10166676129318665, 0.23723379504183550, 0.40828267875217510,
      0.59171732124782490, 0.76276620495816450, 0.89833323870681340, 0.98014492824876820};
gw = {0.05061426814518815, 0.11119051722668725, 0.15685332293894365, 0.18134189168918100,
      0.18134189168918100, 0.15685332293894365, 0.11119051722668725, 0.05061426814518815};

(* production bilinear-quad 8x8 Gauss for (H,gH) over a GENUINE quad {V0,V1,V2,V3} (accumQG surface part) *)
gaussQuadHG[obs_, V0_,V1_,V2_,V3_] := Module[
  {H={0.,0.,0.}, gH={0.,0.,0.,0.,0.,0.}, u,v,wuv,a0,a1,a2,a3,P,Tu,Tv,jv,dA,d,dx,dy,dz,r2,ir,ir3,ir5,c3,c5},
  Do[ u=gp[[iu]]; v=gp[[iv]]; wuv=gw[[iu]]gw[[iv]];
      a0=(1-u)(1-v); a1=u(1-v); a2=u v; a3=(1-u)v;
      P=a0 V0+a1 V1+a2 V2+a3 V3;
      Tu=(1-v)(V1-V0)+v(V2-V3); Tv=(1-u)(V3-V0)+u(V2-V1);
      jv=Cross[Tu,Tv]; dA=Sqrt[jv.jv]wuv;
      d=obs-P; {dx,dy,dz}=d; r2=d.d; ir=1/Sqrt[r2]; ir3=ir/r2; ir5=ir3/r2; c3=ir3 dA; c5=ir5 dA;
      H+=d c3;
      gH[[1]]+=c3-3dx dx c5; gH[[2]]+=c3-3dy dy c5; gH[[3]]+=c3-3dz dz c5;
      gH[[4]]+=-3dx dy c5;   gH[[5]]+=-3dx dz c5;   gH[[6]]+=-3dy dz c5;
    ,{iu,1,8},{iv,1,8}];
  {H,gH}];

(* analytic H (port of FieldFromChargedTriangleLocal, sigma=1, no 4pi), symbolic-obs capable *)
analyticH[ox_,oy_,oz_, V0_,V1_,V2_] := Module[
  {obs,e1,e2,e1len,ba,cc,clen,bc,bb,XYx,XYy,DS,AM,XD,YD,dd,EE1,EE2,EE3,Xv,Yv,Hv,Ev,Rv,RM,RP,AL,HH1,HH2,HH3,AT,BT,l},
  obs={ox,oy,oz}; e1=V1-V0; e2=V2-V0; e1len=Sqrt[e1.e1]; ba=e1/e1len;
  cc=Cross[e1,e2]; clen=Sqrt[cc.cc]; bc=cc/clen; bb=Cross[bc,ba];
  XYx={0,e1len,e2.ba}; XYy={0,0,e2.bb};
  DS=AM=XD=YD=ConstantArray[0,3];
  Do[ l=Mod[j,3]+1; With[{dxj=XYx[[l]]-XYx[[j]], dyj=XYy[[l]]-XYy[[j]]},
        DS[[j]]=Sqrt[dxj^2+dyj^2]; AM[[j]]=dyj/dxj; XD[[j]]=-dxj/DS[[j]]; YD[[j]]=dyj/DS[[j]]];,{j,3}];
  dd=obs-V0; EE1=dd.ba; EE2=dd.bb; EE3=dd.bc;
  Xv=Table[EE1-XYx[[j]],{j,3}]; Yv=Table[EE2-XYy[[j]],{j,3}];
  Hv=Table[Yv[[j]]Xv[[j]],{j,3}]; Ev=Table[EE3^2+Xv[[j]]^2,{j,3}];
  Rv=Table[Sqrt[Xv[[j]]^2+Yv[[j]]^2+EE3^2],{j,3}];
  RM=Table[Rv[[j]]+Rv[[Mod[j,3]+1]]-DS[[j]],{j,3}]; RP=Table[Rv[[j]]+Rv[[Mod[j,3]+1]]+DS[[j]],{j,3}];
  AL=Table[Log[RM[[j]]/RP[[j]]],{j,3}];
  HH1=-Sum[YD[[j]]AL[[j]],{j,3}]; HH2=-Sum[XD[[j]]AL[[j]],{j,3}];
  AT=Table[(AM[[j]]Ev[[j]]-Hv[[j]])/(EE3 Rv[[j]]),{j,3}];
  BT=Table[(AM[[j]]Ev[[Mod[j,3]+1]]-Hv[[Mod[j,3]+1]])/(EE3 Rv[[Mod[j,3]+1]]),{j,3}];
  HH3=-Sum[ArcTan[AT[[j]]],{j,3}]+Sum[ArcTan[BT[[j]]],{j,3}];
  HH1 ba+HH2 bb+HH3 bc];

(* analytic quad = fan {V0,V1,V2}+{V0,V2,V3}; G = Jacobian of the H-sum *)
analyticQuadHG[obs_, V0_,V1_,V2_,V3_] := Module[{ox,oy,oz, Hs, M, Hval, Mval},
  Hs = analyticH[ox,oy,oz,V0,V1,V2] + analyticH[ox,oy,oz,V0,V2,V3];
  M  = D[Hs, {{ox,oy,oz}}];
  Hval = Hs /. {ox->obs[[1]],oy->obs[[2]],oz->obs[[3]]};
  Mval = M  /. {ox->obs[[1]],oy->obs[[2]],oz->obs[[3]]};
  {Hval, {Mval[[1,1]],Mval[[2,2]],Mval[[3,3]],Mval[[1,2]],Mval[[1,3]],Mval[[2,3]]}}];

planar = {{0.,0.,0.},{1.0,0.,0.},{1.1,0.9,0.},{0.2,1.0,0.}};               (* z=0 planar quad *)
warped = {{0.,0.,0.},{1.0,0.,0.05},{1.1,0.9,-0.04},{0.2,1.0,0.06}};        (* non-planar *)
obss = {{0.5,0.4,0.8},{1.6,-0.5,1.1},{-0.4,0.5,-1.2}};

Print["=== planar quad: analytic 2-tri fan vs production bilinear 8x8 Gauss ==="];
Do[Module[{q=gaussQuadHG[ob,Sequence@@planar], a=analyticQuadHG[ob,Sequence@@planar], hr,gr},
   hr=Norm[a[[1]]-q[[1]]]/Norm[q[[1]]]; gr=Norm[a[[2]]-q[[2]]]/Norm[q[[2]]];
   Print[StringForm["  obs=`1`  H rel=`2`  gH rel=`3`  (expect ~Gauss accuracy 1e-6..1e-9)",
      NumberForm[ob,2], ScientificForm[hr,3], ScientificForm[gr,3]]];],{ob,obss}];

Print["=== warped quad: flat-triangulation vs bilinear surface (quantify modeling diff) ==="];
Do[Module[{q=gaussQuadHG[ob,Sequence@@warped], a=analyticQuadHG[ob,Sequence@@warped], hr,gr},
   hr=Norm[a[[1]]-q[[1]]]/Norm[q[[1]]]; gr=Norm[a[[2]]-q[[2]]]/Norm[q[[2]]];
   Print[StringForm["  obs=`1`  H rel=`2`  gH rel=`3`  (O(warp) modeling diff, NOT an error)",
      NumberForm[ob,2], ScientificForm[hr,3], ScientificForm[gr,3]]];],{ob,obss}];
Print["=== DONE ==="];
