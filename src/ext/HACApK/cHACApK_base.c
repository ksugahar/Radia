/*
!=====================================================================*
!                                                                     *
!   Software Name : HACApK                                            *
!         Version : 1.3.0                                             *
!                                                                     *
!   License                                                           *
!     This file is part of HACApK.                                    *
!     HACApK is a free software, you can use it under the terms       *
!     of The MIT License (MIT). See LICENSE file and User's guide     *
!     for more details.                                               *
!                                                                     *
!   ppOpen-HPC project:                                               *
!     Open Source Infrastructure for Development and Execution of     *
!     Large-Scale Scientific Applications on Post-Peta-Scale          *
!     Supercomputers with Automatic Tuning (AT).                      *
!                                                                     *
!   Sponsorship:                                                      *
!     Japan Science and Technology Agency (JST), Basic Research       *
!     Programs: CREST, Development of System Software Technologies    *
!     for post-Peta Scale High Performance Computing.                 *
!                                                                     *
!   Copyright (c) 2015 <Akihiro Ida and Takeshi Iwashita>             *
!                                                                     *
!=====================================================================*
!C**************************************************************************
!C  This file includes basic routines for H-matrices
!C  created by Akihiro Ida at Kyoto University on May 2012
!C  added functions related to ACA+ to HACApK1.0.0 on Nov. 2016
!C  corrected the allocation for st_ctl%lthr on Nov. 2016
!C  added a function related to HACApK_view to HACApK1.1.0 on May 2017
!C  added a function related to writing H-matrix to HACApK1.2.0 on May 2017
!C  added functions related to Block clustering to HACApK1.2.0 on May 2017
!C  translated to C language by Akihiro Ida and Kazuya Goto
!C**************************************************************************
*/
#include "cHACApK_base.h"
#include "cHACApK_calc_entry_ij.h"
#include "cHACApK_lib.h"
#include "cHACApK_pca_cluster.h"   /* PCA cluster split for flat meshes */
#include "hacapk_log.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi_stub.h"

#include "rad_hacapk_parallel.h"

/* Windows defines min/max as macros - undefine to avoid conflicts */
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

/* Symmetric-fill mode (opt-in, default OFF): when set, the leaf fill SKIPS every strictly-lower leaf
 * (nstrtl > nstrtt) -- HACApK_matvec_sym_wrapper never reads them (it mirrors the upper triangle), so for
 * a manager whose applies all route through the symmetric matvec (the HDiv charge Gram) their ACA/dense
 * fill is pure waste (~half the build time and leaf memory).  A matrix built this way MUST only be applied
 * through matvec_sym: the owning manager is responsible for routing plain matvec/matvec_transpose to it
 * (RadHACApKChargeGram does).  Set around ONE build and reset after; not safe across CONCURRENT builds
 * (builds are serial at the Radia API level -- BuildHMatrix stands up its own TaskManager region). */
static volatile int g_sym_fill = 0;
void cHACApK_set_sym_fill(int flag) { g_sym_fill = flag; }

//***cHACApK_generate_frame_blrleaf
void cHACApK_generate_frame_blrleaf(
  st_cHACApK_leafmtxp st_leafmtxp,
  int i_bemv,
  st_cHACApK_lcontrol st_ctl,
  double **gmid_t,  // 2D array [ndim+1][nofc+1]
  int lnmtx[4+1],
  int nofc,
  int nffc,
  int ndim)
{
  st_cHACApK_cluster st_clt;
  st_cHACApK_leafmtx *st_leafmtx, *st_leafmtx_lcl;
  int64_t mem8, nlfall;
  int *lhp, *lnp;
  double *param;
  int *lpmd, *lod, *lthr, *lodfc;
  int mpinr, mpilog, nrank, irank, icomm, nthr, nd;
  int nsrt,ndf,nclst,ndpth,ndscd,nblall,nlfalt,ill,itt;
  int npgl,npgt,ierr,il,ig,ip,in,it,is,ikey, iclr, icommn;
  int nlft, nlfth, mpinrth, nlfl, nlflh, mpinrlh, nbl;
  int nrank_t,nrank_l,irank_t,irank_l,inml,inmt;
  int ilh,ith,nlf,ipgclr,ilf,itf,ndlfs, ndtfs,iw,isnlf;
  int ltmtx,ndl,ndt,ns;
  double zzz,ktp;
  MPI_Comm comm,commn;
  char fname[32];
  FILE *fmpilog;

  param = st_ctl->param;
  lpmd = st_ctl->lpmd; lod = st_ctl->lod; lthr = &(st_ctl->lthr[1]);
  mpinr=lpmd[3]; mpilog=lpmd[4]; nrank=lpmd[2]; icomm=lpmd[1]; nthr=lpmd[20];

  comm=MPI_Comm_f2c(icomm);

  ierr = MPI_Comm_rank ( comm, &irank );
  snprintf(fname,sizeof(fname),"log%04d.txt",irank);
  fmpilog=fopen(fname,"a");
  if(fmpilog==NULL) {
    fprintf(stderr, "Error: cHACApK_generarte_frame_blrleaf: fopen %s\n",fname);
    goto error;
  }

  nd=nofc*nffc;
  lodfc = (int *) malloc(sizeof(int)*(nofc+1));
  if(lodfc==NULL) {
    fprintf(stderr, "Error: cHACApK_generate_frame_blrleaf: malloc lodfc\n");
    goto error;
  }
  for(il=1; il<=nofc; il++) {
    lodfc[il]=il;
  }
  //!!!!!!!!!!!!!!!! start clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nsrt=1; ndf=nofc; nclst=0; ndpth=0; ndscd=0;
  cHACApK_generate_cbitree(&st_clt,gmid_t,param,lpmd,lodfc,&ndpth,ndscd,nsrt,ndf,nofc,ndim,&nclst);
  if(st_ctl->param[1]>0 && mpinr==0) printf("No. of cluster=%12d\n",nclst);
  if(st_ctl->param[1]>1)  printf("No. of cluster=%12d\n",nclst);

  cHACApK_bndbox(st_clt,gmid_t,lodfc,nofc);
  for(il=1; il<=nofc; il++) {
    for(ig=1; ig<=nffc; ig++) {
      is=ig+(il-1)*nffc;
      lod[is]=lodfc[il];
    }
  }
  //!!!!!!!!!!!!!!!! end clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  free(lodfc);

  
  //for(ill=1; ill<10; ill++) {
    //ill=1; itt=1; zzz=cHACApK_entry_ij(lod[ill],lod[itt],i_bemv);
    //itt=1; zzz=cHACApK_entry_ij(lod[ill],lod[itt],i_bemv);
    //printf("ill=%12d; itt=%12d; zzz=%21.6lf\n",ill,itt,zzz);
  //}
  //!!!!!!!!!!!!!!!! start construction of H-matrix  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  npgl=param[41]; if(npgl==0) npgl=sqrt((double)nrank);
  // if(param[42]==0) param[42]=sqrt((double)nd);
  if(param[42]==0) param[42]=nd/param[43]/npgl;
  if(param[42]<param[21]) {
    if(mpinr==0) printf("sub HACApK_generate_frame_blrleaf; param[42]=%21.6lf param[21]=%21.6lf\n",param[42],param[21]);
    if(mpinr==0) printf("Error: sub cHACApK_generate_frame_blrleaf; param[42](block size) must be larger than param[21](leaf size) !!!\n"); goto error;
  }
  ndpth=0;
  for(il=1; il<=4; il++) lnmtx[il]=0;
  cHACApK_count_blrnmb(st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,&ndpth);
  nblall=lnmtx[4];
  st_leafmtx = (st_cHACApK_leafmtx *) malloc(sizeof(st_cHACApK_leafmtx)*(nblall+1));
  if(st_leafmtx==NULL) {
    fprintf(stderr, "Error: cHACApK_generate_frame_blrleaf: malloc st_leafmtx\n");
    goto error;
  }
  for(il=1; il<=nblall; il++) {
    st_leafmtx[il] = (st_cHACApK_leafmtx) calloc(1,sizeof(st_cHACApK_leafmtx_t));
    if(st_leafmtx[il]==NULL) {
      fprintf(stderr, "Error: cHACApK_generate_frame_blrleaf: malloc st_leafmtx[%d]\n",il);
      goto error;
    }
  }
  nlfalt=sqrt((double)nblall); st_leafmtxp->nlfalt=nlfalt;
  if(st_ctl->param[1]>0 && mpinr==0) printf("Number of MPI_Blocks=%12d; sqrt(nblall)=%12d\n",nblall,nlfalt);
  ndpth=0;
  for(il=1; il<=4; il++) lnmtx[il]=0;
  cHACApK_count_blrleaf(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,&ndpth);
  if(st_ctl->param[1]>0 && mpinr==0) printf("No. of nsmtx %12d %12d %12d %12d\n",lnmtx[1],lnmtx[2],lnmtx[3],lnmtx[4]);
  if(st_ctl->param[1]>0 && mpinr==0) printf("   1:Rk-matrix 2: dense-mat 3:H-matrix 4:MPI_Block\n");
  st_leafmtxp->nlfkt=lnmtx[1];
  nlfall=lnmtx[1]+lnmtx[2];
  if(st_ctl->param[1]>0 && mpinr==0) printf("nlf global=%12ld\n",nlfall);
  if(nlfall<nthr) {
    printf("Error; HACApK_generate_frame_blrleaf; # of leaves must be larger than # of threads.\n");
    exit(EXIT_FAILURE);
  }

  nblall=0; ndpth=0;
  for(il=1; il<=4; il++) lnmtx[il]=0;
  cHACApK_generate_blrleaf(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,&nblall,&ndpth);
  if(st_ctl->param[1]>1 && mpinr==0) printf("HACApK_generate_frame_blrleaf; HACApK_generate_leafmtx end\n");
  cHACApK_sort_leafmtx(st_leafmtx,nblall);
  for(ip=1; ip<=nblall; ip++) {
    if(st_leafmtx[ip]->ltmtx==4) cHACApK_sort_leafmtx(st_leafmtx[ip]->st_lf,st_leafmtx[ip]->nlf);
  }
  if(st_ctl->param[1]>1 && mpinr==0) printf("HACApK_generate_frame_blrleaf; HACApK_sort_leafmtx end\n");
  cHACApK_free_st_clt(st_clt);

  //!!!!!!!!!!!!!!!! start MPI load balance  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npgt=nrank/npgl;
  if(st_ctl->param[1]>0 && mpinr==0) printf(" npgl=%12d npgt=%12d\n",npgl,npgt);
  if(st_ctl->param[1]>1) printf(" npgl=%12d npgt=%12d\n",npgl,npgt);
  if(npgt>nlfalt || npgl>nlfalt) {
    ierr = MPI_Barrier( comm );
    if(mpinr==0) printf("Error: HACApK_generate_frame_blrleaf; Too few blocks compared with #MPI !!!\n"); goto error;
  }
  if(npgt*npgl!=nrank) {
    ierr = MPI_Barrier( comm );
    if(mpinr==0) printf("Error: HACApK_generate_frame_blrleaf; Invalid processor grid!!!\n"); goto error;
  }

// Split MPI communicator
  ikey=0; iclr=mpinr/npgt;
  ierr = MPI_Comm_split(comm, iclr, ikey, &commn); icommn=MPI_Comm_c2f(commn); st_ctl->lpmd[31]=icommn;
  if(ierr!=0) {
    if(mpinr==0) printf("Error: sub HACApK_generate_frame_blrleaf; MPI_COMM_SPLIT failed !!!\n"); goto error;
  }
  ierr = MPI_Comm_size ( commn, &nrank ); st_ctl->lpmd[32]=nrank;
  if(ierr!=0) {
    if(mpinr==0) printf("Error: sub HACApK_generate_frame_blrleaf; MPI_Comm_size failed !!!\n"); goto error;
  }
  ierr = MPI_Comm_rank ( commn, &irank ); st_ctl->lpmd[33]=irank;
  if(ierr!=0) {
    if(mpinr==0) printf("Error: sub HACApK_generate_frame_blrleaf; MPI_Comm_rank failed !!!\n"); goto error;
  }
  ikey=0; iclr=mpinr%npgt;
  ierr = MPI_Comm_split(comm, iclr, ikey, &commn); icommn=MPI_Comm_c2f(commn); st_ctl->lpmd[35]=icommn;
  if(ierr!=0) {
    if(mpinr==0) printf("Error: sub HACApK_generate_frame_blrleaf; MPI_COMM_SPLIT failed !!!\n"); goto error;
  }
  ierr = MPI_Comm_size ( commn, &nrank ); st_ctl->lpmd[36]=nrank;
  if(ierr!=0) {
    if(mpinr==0) printf("Error: sub HACApK_generate_frame_blrleaf; MPI_Comm_size failed !!!\n"); goto error;
  }
  ierr = MPI_Comm_rank ( commn, &irank ); st_ctl->lpmd[37]=irank;
  if(ierr!=0) {
    if(mpinr==0) printf("Error: sub HACApK_generate_frame_blrleaf; MPI_Comm_rank failed !!!\n"); goto error;
  }

  if(st_ctl->param[1]>1) printf("irank=%12d; irank_t=%12d; irank_l=%12d\n",mpinr,st_ctl->lpmd[33],st_ctl->lpmd[37]);
  if(st_ctl->param[1]>1) fprintf(fmpilog,"irank_t=%12d; nrank_t=%12d\n",st_ctl->lpmd[33],st_ctl->lpmd[32]);
  if(st_ctl->param[1]>1) fprintf(fmpilog,"irank_l=%12d; nrank_l=%12d\n",st_ctl->lpmd[37],st_ctl->lpmd[36]);

  // stop

  nlft=nlfalt/npgt; nlfth=nlfalt%npgt; mpinrth=mpinr%npgt;
  if(mpinrth<nlfth) nlft=nlft+1;
  nlfl=nlfalt/npgl; nlflh=nlfalt%npgl; mpinrlh=mpinr/npgt;
  if(mpinrlh<nlflh) nlfl=nlfl+1;
  nbl=nlfl*nlft; st_leafmtxp->nbl=nbl; st_leafmtxp->nlfl=nlfl; st_leafmtxp->nlft=nlft;
  if(st_ctl->param[1]>1) printf("irank=%12d; nbl=%12d; nblall=%12d\n",mpinr,nbl,nblall);
  if(st_ctl->param[1]>1) fprintf(fmpilog,"No. of blocks; nbl=%12d; row=%12d; column=%12d; global nbl=%12d\n",nbl,nlfl,nlft,nblall);


  ierr = MPI_Barrier( comm );
  // stop

  nrank_t=st_ctl->lpmd[32]; nrank_l=st_ctl->lpmd[36];

  st_leafmtxp->lbl2t = (int *) calloc(npgl,sizeof(int));
  if(st_leafmtxp->lbl2t==NULL) {
    fprintf(stderr, "Error: cHACApK_generate_frame_blrleaf: malloc st_leafmtxp->lbl2t\n");
    goto error;
  }
  irank_t=st_ctl->lpmd[33]; irank_l=st_ctl->lpmd[37];
  for(in=0; in<nlfalt; in++) {
    inml=in%npgl; inmt=in%npgt;
    if(inmt==irank_t) st_leafmtxp->lbl2t[inml]=1;
  }
  if(st_ctl->param[1]>1) fprintf(fmpilog,"st_leafmtxp->lbl2t\n");
  if(st_ctl->param[1]>1) {
    for(il=0; il<npgl; il++) fprintf(fmpilog,"%12d\n",st_leafmtxp->lbl2t[il]);
  }

  st_leafmtxp->lbstrtl = (int *) malloc(sizeof(int)*(nlfalt+1+1));
  st_leafmtxp->lbstrtt = (int *) malloc(sizeof(int)*(nlfalt+1+1));
  st_leafmtxp->lbndl = (int *) malloc(sizeof(int)*(nlfalt+1));
  st_leafmtxp->lbndt = (int *) malloc(sizeof(int)*(nlfalt+1));
  st_leafmtxp->lbndlfs = (int *) calloc(nrank_l,sizeof(int));
  st_leafmtxp->lbndtfs = (int *) calloc(nrank_t,sizeof(int));
  if(st_leafmtxp->lbstrtl==NULL ||
     st_leafmtxp->lbstrtt==NULL ||
     st_leafmtxp->lbndl  ==NULL ||
     st_leafmtxp->lbndt  ==NULL ||
     st_leafmtxp->lbndlfs==NULL ||
     st_leafmtxp->lbndtfs==NULL) {
    fprintf(stderr, "Error: cHACApK_generate_frame_blrleaf: malloc st_leafmtxp->lb***\n");
    goto error;
  }
  for(il=0; il<nlfalt; il++) {
    is=nlfalt*il+1;
    ilh=il%npgl;
    st_leafmtxp->lbndlfs[ilh]=st_leafmtxp->lbndlfs[ilh]+st_leafmtx[is]->ndl;
    st_leafmtxp->lbstrtl[il+1]=st_leafmtx[is]->nstrtl;
    st_leafmtxp->lbndl[il+1]=st_leafmtx[is]->ndl;
  }
  st_leafmtxp->lbstrtl[nlfalt+1]=nd+1;
  for(it=0; it<nlfalt; it++) {
    is=it+1;
    ith=it%npgt;
    st_leafmtxp->lbndtfs[ith]=st_leafmtxp->lbndtfs[ith]+st_leafmtx[is]->ndt;
    st_leafmtxp->lbstrtt[it+1]=st_leafmtx[is]->nstrtt;
    st_leafmtxp->lbndt[it+1]=st_leafmtx[is]->ndt;
  }
  st_leafmtxp->lbstrtt[nlfalt+1]=nd+1;
  if(0) {
// if(mpinr==0) {
    printf("lbstrtl=\n");
    for(il=1; il<=nlfalt+1; il++) printf("%12d\n",st_leafmtxp->lbstrtl[il]);
    printf("lbstrtt=\n");
    for(il=1; il<=nlfalt+1; il++) printf("%12d\n",st_leafmtxp->lbstrtt[il]);
    printf("lbndlfs=\n");
    for(il=0; il<nrank_l; il++) printf("%12d\n",st_leafmtxp->lbndlfs[il]);
    printf("lbndtfs=\n");
    for(il=0; il<nrank_t; il++) printf("%12d\n",st_leafmtxp->lbndtfs[il]);
  }

  st_leafmtx_lcl = (st_cHACApK_leafmtx *) malloc(sizeof(st_cHACApK_leafmtx)*(nbl+1));
  st_leafmtxp->lnlfl2g_t = (int64_t **) malloc(sizeof(int64_t *)*(nlfl+1));
  if(st_leafmtxp->lnlfl2g_t==NULL) {
    fprintf(stderr, "Error: cHACApK_generate_frame_blrleaf: malloc st_leafmtxp->lnlfl2g\n");
    goto error;
  }
  for(il=1; il<=nlfl; il++) {
    st_leafmtxp->lnlfl2g_t[il] = (int64_t *) malloc(sizeof(int64_t)*(nlft+1));
    if(st_leafmtxp->lnlfl2g_t[il]==NULL) {
      fprintf(stderr, "Error: cHACApK_generate_frame_blrleaf: malloc st_leafmtxp->lnlfl2g[%d]\n",il);
      goto error;
    }
  }
  ip=0; nlf=0;
  for(il=0; il<nlfalt; il++) {
    for(it=0; it<nlfalt; it++) {
      is=it+nlfalt*il+1;
      ilh=il%npgl; ith=it%npgt;
      ipgclr=ith+ilh*npgt;
      if(ipgclr==mpinr) {
        ilf=ip/nlft+1; itf=ip%nlft+1; ip=ip+1;
        st_leafmtxp->lnlfl2g_t[ilf][itf]=is;
        st_leafmtx_lcl[ip]=st_leafmtx[is];
        if(st_leafmtx[is]->ltmtx==1) {
          nlf=nlf+1;
        } else if(st_leafmtx[is]->ltmtx==4) {
          nlf=nlf+st_leafmtx[is]->nlf;
        }
      }
    }
  }

  ndlfs=0;
  for(ilf=1; ilf<=nlfl; ilf++) {
    ip=(ilf-1)*nlft+1;
    ndlfs=ndlfs+st_leafmtx_lcl[ip]->ndl;
  }
  st_leafmtxp->ndlfs=ndlfs;

  ndtfs=0;
  for(ip=1; ip<=nlft; ip++) {
    ndtfs=ndtfs+st_leafmtx_lcl[ip]->ndt;
  }
  st_leafmtxp->ndtfs=ndtfs;

//  print*,'mpinr=',mpinr,'; nlf=',nlf
  st_leafmtxp->nlf=nlf;
  st_leafmtxp->st_lf = (st_cHACApK_leafmtx *) malloc(sizeof(st_cHACApK_leafmtx)*(nlf+1));
  if(st_leafmtxp->st_lf==NULL) {
    fprintf(stderr, "Error: cHACApK_generate_frame_blrleaf: malloc st_leafmtxp->st_lf\n");
    goto error;
  }
  ip=0; ndlfs=0; ndtfs=0;
  for(il=0; il<nlfalt; il++) {
    for(it=0; it<nlfalt; it++) {
      is=it+nlfalt*il+1;
      ilh=il%npgl; ith=it%npgt;
      ipgclr=ith+ilh*npgt;
      if(ipgclr==mpinr) {
        if(st_leafmtx[is]->ltmtx==1) {
          ip=ip+1;
          st_leafmtxp->st_lf[ip]=st_leafmtx[is];
        } else if(st_leafmtx[is]->ltmtx==4) {
          isnlf=st_leafmtx[is]->nlf;
          for(iw=1; iw<=isnlf; iw++) {
            st_leafmtxp->st_lf[ip+iw]=st_leafmtx[is]->st_lf[iw];
          }
          ip=ip+isnlf;
          // if(st_leafmtx[is]->nstrtt==1) ndlfs=ndlfs+st_leafmtx[is]->ndl:
          // if(st_leafmtx[is]->nstrtl==1) ndtfs=ndtfs+st_leafmtx[is]->ndt;
        }
      }
    }
  }

  free(st_leafmtx);
  free(st_leafmtx_lcl);

//! print*,'mpinr=',mpinr,'; ndlfs=',st_leafmtxp%ndlfs,'; ndtfs=',st_leafmtxp%ndtfs
//! call MPI_Barrier( icomm, ierr )
// stop

  if(st_ctl->param[1]>1) printf("irank=%12d; ndlfs=%12d; ndtfs=%12d\n",mpinr,st_leafmtxp->ndlfs,st_leafmtxp->ndtfs);
  if(st_ctl->param[1]>1) fprintf(fmpilog,"Vector sizes; nd=%12d; ndlfs=%12d; ndtfs=%12d\n",nd,st_leafmtxp->ndlfs,st_leafmtxp->ndtfs);

  if(st_ctl->param[1]>1) {
    fprintf(fmpilog,"lnlfl2g=\n");
    for(il=1; il<=nlfl; il++) {
      for(it=1; it<=nlft; it++) {
        fprintf(fmpilog,"%9ld",st_leafmtxp->lnlfl2g_t[il][it]);
      }
    }
  }

  ierr = MPI_Barrier( comm );

  for(il=1; il<=4; il++) lnmtx[il]=0;
  mem8=0; ktp=param[62];
  for(ip=1; ip<=nlf; ip++) {
    ltmtx=st_leafmtxp->st_lf[ip]->ltmtx; ndl=st_leafmtxp->st_lf[ip]->ndl; ndt=st_leafmtxp->st_lf[ip]->ndt; ns=ndl*ndt;
    if(ltmtx==1) {
      lnmtx[1]=lnmtx[1]+1; mem8=mem8+(ndt+ndl)*ktp;
    } else {
      lnmtx[2]=lnmtx[2]+1; mem8=mem8+ns;
    }
  }
  if(st_ctl->param[1]>1)  fprintf(fmpilog,"No. of nsmtx %12d %12d\n",lnmtx[1],lnmtx[2]);
  cHACApK_setcutthread(lthr,st_leafmtxp,st_ctl,mem8,nthr,ktp);

  fclose(fmpilog);
  return;

error:
  exit(EXIT_FAILURE);
}

//***cHACApK_setcutthread
void cHACApK_setcutthread(
  int *lthr,
  st_cHACApK_leafmtxp st_leafmtxp,
  st_cHACApK_lcontrol st_ctl,
  int64_t mem8,
  int nthr,
  int ktp)
{
  int nlf,ith,kt,il,ltmtx,ndl,ndt;
  int64_t nth1_mem,imem;

  nlf=st_leafmtxp->nlf;
  nth1_mem=mem8/nthr;
  // if(st_ctl->param[1]>1) printf("HACApK_setcutthread; nlf=%12d mem8=%12ld nthr=%12d\n",nlf,mem8,nthr);
  lthr[0]=1; lthr[nthr]=nlf+1;
  imem=0; ith=1; kt=ktp;
  for(il=1; il<=nlf; il++) {
    ltmtx=st_leafmtxp->st_lf[il]->ltmtx;
    ndl=st_leafmtxp->st_lf[il]->ndl; ndt=st_leafmtxp->st_lf[il]->ndt;
    if(ltmtx==1) {
      if(ktp==0) kt=st_leafmtxp->st_lf[il]->kt;
      imem=imem+(ndl+ndt)*kt;
    } else {
      imem=imem+ndl*ndt;
    }
    if(imem>nth1_mem*ith) {
      lthr[ith]=il;
      ith=ith+1;
      if(ith==nthr) break;
    }
  }
  // if(st_ctl->param[1]>1) printf("HACApK_setcutthread; lthr=%12d\n",lthr[0:nthr]);
}

//***min
static int min(
  int i1,
  int i2)
{
  return i1<i2 ? i1 : i2;
}

//***max
static int max(
  int i1,
  int i2)
{
  return i1>i2 ? i1 : i2;
}

extern void dgeqp3_(int*,int*,double*,int*,int*,double*,double*,int*,int*);
extern void dorgqr_(int*,int*,int*,double*,int*,double*,double*,int*,int*);

//***cHACApK_RRQR
int cHACApK_RRQR(
  double *zaa, // zaa(ndl,kmax)
  double *zab, // zab(ndt,kmax)
  double *param,
  int ndl,
  int ndt,
  int nstrtl,
  int nstrtt,
  int *lod,
  int i_bemv,
  int kmax,
  double eps,
  double znrmmat,
  double pRRQR_EPS)
{
  double *tau,*work,work1[1];
  double *waa;
  int *jpvt;
  int kmin,kRRQR,nn,lwork,lda,il,it,ill,itt,info;
  double zzz;

  kmin=param[64];
  kRRQR=0;

  nn=min(ndl,ndt);
  lda=ndl;
  waa = (double *) calloc(lda*ndt,sizeof(double));
  jpvt = (int *) calloc(ndt,sizeof(int));
  tau = (double *) calloc(nn,sizeof(double));
  if(waa==NULL || jpvt==NULL || tau==NULL) {
    fprintf(stderr, "Error: cHACApK_RRQR: malloc waa jpvt tau\n");
    goto error;
  }

  for (il=0; il<ndl; il++) {
    for (it=0; it<ndt; it++) {
      ill=il+nstrtl; itt=it+nstrtt;
      waa[il+lda*it]=cHACApK_entry_ij(lod[ill],lod[itt],i_bemv);
    }
  }

  lwork=-1;
  dgeqp3_(&ndl,&ndt,waa,&lda,jpvt,tau,work1,&lwork,&info);
  lwork=work1[0];
  work = (double *) malloc(lwork*sizeof(double));
  if(work==NULL) {
    fprintf(stderr, "Error: cHACApK_RRQR: malloc work\n");
    goto error;
  }
  dgeqp3_(&ndl,&ndt,waa,&lda,jpvt,tau,work,&lwork,&info);

  for (il=1; il<nn; il++) {
    zzz=fabs(waa[il+lda*il]/waa[0]);
    if(zzz<eps) break;
  }

  kRRQR=min(il,nn);
  if(kRRQR<kmin) {
    kRRQR=kmin;
  } else if(kRRQR>kmax) {
    kRRQR=kmax;
  }

  for (it=0; it<kRRQR; it++) {
    for (il=0; il<it; il++) {
      zab[(jpvt[il]-1)+ndt*it]=0.0;
    }
    for (il=it; il<ndt; il++) {
      zab[(jpvt[il]-1)+ndt*it]=waa[it+lda*il];
    }
  }

  // dorgqr_(&ndl,&nn,&nn,waa,&lda,tau,work,&lwork,&info);
  // dorgqr_(&ndl,&nn,&kRRQR,waa,&lda,tau,work,&lwork,&info);
  dorgqr_(&ndl,&kRRQR,&kRRQR,waa,&lda,tau,work,&lwork,&info);
  for (it=0; it<kRRQR; it++) {
    for (il=0; il<ndl; il++) {
      zaa[il+ndl*it]=waa[il+lda*it];
    }
  }

  free(work); free(waa); free(jpvt); free(tau);
  return kRRQR;
error:
  exit(EXIT_FAILURE);
}

extern void dgesvd_(char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);

//***cHACApK_SVD
int cHACApK_SVD(
  double *zaa, // zaa(ndl,kmax)
  double *zab, // zab(ndt,kmax)
  double *param,
  int ndl,
  int ndt,
  int nstrtl,
  int nstrtt,
  int *lod,
  int i_bemv,
  int kmax,
  double eps,
  double znrmmat,
  double pSVD_EPS)
{
  double *prow,*pcol;
  int *lrow_msk,*lcol_msk;
  double *w,*work;
  double *u,*vt,*waa;
  int kmin,krank,kSVD,k,lstop_aca,nn,lwork,lda,ldu,ldvt,il,it,ill,itt,info;
  double znrm,zzz;
  char jobu,jobvt;

  kmin=param[64];
  // printf("nstrtl=%12d nstrtt=%12d ndl=%12d ndt=%12d kmax=%12d\n",nstrtl,nstrtt,ndl,ndt,kmax);
  krank=min(ndl,ndt);
  znrm=znrmmat*sqrt((double)ndl*(double)ndt);
  // allocate(lrow_msk(ndl),lcol_msk(ndt)); lrow_msk(:)=0; lcol_msk(:)=0; nrow_done=0; ncol_done=0
  kSVD=0; k=1; lstop_aca=0;

  nn=min(ndl,ndt);
  //print*,'nn=',nn,' eps=',eps
  lwork=10*nn; lda=ndl; ldu=ndl; ldvt=ndt;
  //allocate(w(nn),work(lwork),u(ldu,nn),vt(ldvt,nn),waa(ndl,ndt),waa2(ndl,ndt))
  w = (double *) calloc(nn,sizeof(double));
  work = (double *) calloc(lwork,sizeof(double));
  u = (double *) calloc(ldu*nn,sizeof(double));
  vt = (double *) calloc(nn*ldvt,sizeof(double));
  waa = (double *) calloc(ndl*ndt,sizeof(double));
  if(w==NULL || work==NULL || u==NULL || vt==NULL || waa==NULL) {
    fprintf(stderr, "Error: cHACApK_SVD: malloc w work u vt waa\n");
    goto error;
  }

  for (il=0; il<ndl; il++) {
    for (it=0; it<ndt; it++) {
      ill=il+nstrtl; itt=it+nstrtt;
      waa[il+ndl*it]=cHACApK_entry_ij(lod[ill],lod[itt],i_bemv);
    }
  }

  // HACApK_SVD=kmax; zaa=waa; zab=0.0d0; do il=1,nn; zab(il,il)=1.0d0; enddo; return
  //   waa2=waa

  //call dgesvd ( 'A', 'A', ndl, ndt, waa, lda, w, u, ldu, vt, ldvt, work, lwork, info)
  //call dgesvd ( 'S', 'S', ndl, ndt, waa, lda, w, u, ldu, vt, ldvt, work, lwork, info)
  jobu = 'S'; jobvt = 'S';
  dgesvd_(&jobu, &jobvt, &ndl, &ndt, waa, &lda, w, u, &ldu, vt, &nn, work, &lwork, &info);
  // call HACApK_gesvd (waa,w,u,vt)

  //    print*, 'info_dgesvd=',info
  //    print*, 'eigenvalues_large=',w(1:10)
  //    print*, 'eigenvalues_small=',w(nn-10:nn)
  for (il=1; il<nn; il++) {
    zzz=w[il]/w[0];
    if(zzz<eps) {
      //!        print*,'HACApK_SVD; rank_e4=',il-1,'/',nn
      //        print*,'eigen_max=',w(1)
      //        print*,'eigen_il=',w(il)
      break;
    }
  }
  //!    if(il==nn+1) then
  //!      print*,'HACApK_SVD; rank_e4=',nn,'/',nn
  //!    endif

  kSVD=min(il,nn);
  if(kSVD<kmin) {
    kSVD=kmin;
    //   print*,'HACApK_SVD; rank_e4 is changed to',kmin
  } else if(kSVD>kmax) {
    kSVD=kmax;
    //   print*,'HACApK_SVD; rank_e4 is changed to',kmax
  }

  //za2=0.0; za1=0.0;
  for (it=0; it<kSVD; it++) {
    for (il=0; il<ndl; il++) {
      zaa[il+ndl*it]=u[il+it*ndl]*w[it];
    }
    for (il=0; il<ndt; il++) {
      zab[il+ndt*it]=vt[it+il*nn];
    }
  }

  free(w); free(work); free(u); free(vt); free(waa);
  return kSVD;
error:
  exit(EXIT_FAILURE);
}

//***cHACApK_calc_vec
// ld==0: row direction, ld==1: column direction
// Optimized with BLAS dcopy for strided column extraction
// workspace: pre-allocated buffer of size >= kmax to avoid malloc in hot path
void cHACApK_calc_vec(
  double *zaa,
  double *zab,
  int ndp,
  int ndt,
  int k,
  int ip,
  double *vec,
  int nstrtl,
  int nstrtt,
  int *lod,
  int i_bemv,
  int *lmsk,
  int ld,
  double *workspace)  /* Pre-allocated workspace of size >= k */
{
  int ii,ill,itt,il;

  for (ii=0; ii<ndp; ii++) {
    if(lmsk[ii]==0) {
      if(ld==0) {
        ill=ip+nstrtl; itt=ii+nstrtt;
      } else {
        ill=ii+nstrtl; itt=ip+nstrtt;
      }
      vec[ii]=cHACApK_entry_ij(lod[ill],lod[itt],i_bemv);
    }
  }
  if(k==0) return;
  /* Use pre-allocated workspace instead of malloc/free */
  /* Extract column with stride using BLAS dcopy (optimized strided access) */
  cHACApK_extract_col(workspace, zab, ip, ndt, k);
  cHACApK_adotsub_dsm(vec,zaa,workspace,ndp,k,ndp);
  for (il=0; il<ndp; il++) {
    if(lmsk[il]==1) vec[il]=0.0;
  }
}

//***cHACApK_aca
// Basic Adaptive Cross Approximation (ACA) - simpler and faster than ACA+
// Based on the Fortran implementation from HACApK
int cHACApK_aca(
  double *zaa, // zaa(ndl,kmax) - output left factor
  double *zab, // zab(ndt,kmax) - output right factor
  double *param,
  int ndl,
  int ndt,
  int nstrtl,
  int nstrtt,
  int *lod,
  int i_bemv,
  int kmax,
  double eps,
  double znrmmat,
  double pACA_EPS)
{
  int *lrow_msk, *lcol_msk;
  double *prow, *pcol;
  double *workspace;

  double znrm, ACA_EPS, col_maxval, row_maxval, zdltinv, zeps;
  int krank, kstop, k, ist, jst, istn, lstop_aca;
  int il, it;

  krank = (ndl < ndt) ? ndl : ndt;  /* min(ndl, ndt) */
  znrm = znrmmat * sqrt((double)ndl * (double)ndt);

  if ((int)param[61] == 1) ACA_EPS = pACA_EPS;
  else if ((int)param[61] == 2 || (int)param[61] == 3) ACA_EPS = pACA_EPS * znrm;
  else ACA_EPS = pACA_EPS;

  lrow_msk = (int *) calloc(ndl, sizeof(int));
  lcol_msk = (int *) calloc(ndt, sizeof(int));
  if (lrow_msk == NULL || lcol_msk == NULL) {
    fprintf(stderr, "Error: cHACApK_aca: malloc lrow_msk lcol_msk\n");
    goto error;
  }

  /* Pre-allocate workspace for cHACApK_calc_vec */
  workspace = (double *) calloc(kmax, sizeof(double));
  if (workspace == NULL) {
    fprintf(stderr, "Error: cHACApK_aca: malloc workspace\n");
    goto error;
  }

  k = 0;  /* 0-indexed (Fortran uses 1-indexed) */
  lstop_aca = 0;
  kstop = (kmax < krank) ? kmax : krank;

  /* Initial row index selection */
  if (nstrtl > nstrtt) {
    ist = 0;  /* 0-indexed */
  } else {
    ist = ndl - 1;  /* 0-indexed */
  }

  while (k < kstop && lstop_aca == 0) {
    pcol = zaa + ndl * k;  /* k-th column of zaa */
    prow = zab + ndt * k;  /* k-th column of zab */

    /* Compute row: prow = A(ist, :) - approximation */
    cHACApK_calc_vec(zab, zaa, ndt, ndl, k, ist, prow, nstrtl, nstrtt, lod, i_bemv, lcol_msk, 0, workspace);

    /* Find max abs value in prow (with mask) */
    cHACApK_maxabsvallocm_d(prow, &row_maxval, &jst, ndt, lcol_msk);

    /* Scale row by pivot */
    if (fabs(prow[jst]) > 1.0e-20) {
      zdltinv = 1.0 / prow[jst];
      for (it = 0; it < ndt; it++) prow[it] *= zdltinv;
    } else {
      /* Pivot too small, stop */
      break;
    }

    /* Compute column: pcol = A(:, jst) - approximation */
    cHACApK_calc_vec(zaa, zab, ndl, ndt, k, jst, pcol, nstrtl, nstrtt, lod, i_bemv, lrow_msk, 1, workspace);

    /* Mark used indices */
    lrow_msk[ist] = 1;
    lcol_msk[jst] = 1;

    /* Find next row index (max abs in pcol with mask) */
    cHACApK_maxabsvallocm_d(pcol, &col_maxval, &istn, ndl, lrow_msk);

    /* Check stopping criterion */
    if (fabs(row_maxval) < ACA_EPS && fabs(col_maxval) < ACA_EPS && k >= (int)param[64]) {
      lstop_aca = 1;
      break;
    }

    /* Compute approximation norm */
    zeps = cHACApK_unrm_d(ndl, pcol) * cHACApK_unrm_d(ndt, prow);
    if (k == 0 && (int)param[61] == 1) znrm = zeps;
    zeps = zeps / znrm;

    if (zeps < eps || k == kstop - 1) lstop_aca = 1;
    if (lstop_aca == 1 && k >= (int)param[64]) {
      k++;
      break;
    }

    ist = istn;  /* Next row index */
    k++;
  }

  free(lrow_msk);
  free(lcol_msk);
  free(workspace);
  return k;

error:
  exit(EXIT_FAILURE);
}

//***cHACApK_acaplus
int cHACApK_acaplus(
  double *zaa, // zaa(ndl,kmax)
  double *zab, // zab(ndt,kmax)
  double *param,
  int ndl,
  int ndt,
  int nstrtl,
  int nstrtt,
  int *lod,
  int i_bemv,
  int kmax,
  double eps,
  double znrmmat,
  double pACA_EPS)
{
  int *lrow_msk,*lcol_msk;
  double *pa_ref,*pb_ref;
  double *prow,*pcol;
  double *workspace;  /* Pre-allocated workspace for cHACApK_calc_vec */

  const double za_ACA_EPS=1.0e-30;
  double znrm,ACA_EPS,colnorm,rownorm,apxnorm,col_maxval,row_maxval,zinvmax,blknorm;
  int kacaplus,ntries,ntries_row,ntries_col,k,j_ref,i_ref,lstop_aca,i,j,il,it;
  // write(6,1000) 'nstrtl=',nstrtl,' nstrtt=',nstrtt,' ndl=',ndl,' ndt=',ndt
  znrm=znrmmat*sqrt((double)ndl*(double)ndt);
  if(param[61]==2 || param[61]==1) ACA_EPS=pACA_EPS;
  if(param[61]==3) ACA_EPS=pACA_EPS*znrm;

  kacaplus=0; ntries = max(ndl,ndt)+1; ntries_row = 6; ntries_col = 6;
  lrow_msk = (int *) calloc(ndl,sizeof(int));
  lcol_msk = (int *) calloc(ndt,sizeof(int));
  if(lrow_msk==NULL || lcol_msk==NULL) {
    fprintf(stderr, "Error: cHACApK_acaplus: malloc lrow_msk lcol_msk\n");
    goto error;
  }

  /* Pre-allocate workspace for cHACApK_calc_vec to avoid malloc/free in hot path */
  workspace = (double *) calloc(kmax, sizeof(double));
  if(workspace==NULL) {
    fprintf(stderr, "Error: cHACApK_acaplus: malloc workspace\n");
    goto error;
  }

  k = 0;

  j_ref=0; // arbitrary j_ref
  pa_ref = (double *) calloc(ndl,sizeof(double));
  if(pa_ref==NULL) {
    fprintf(stderr, "Error: cHACApK_acaplus: malloc pa_ref\n");
    goto error;
  }

  cHACApK_calc_vec(zaa,zab,ndl,ndt,k,j_ref,pa_ref,nstrtl,nstrtt,lod,i_bemv,lrow_msk,1,workspace);
  //  print*,'pa_ref=',pa_ref
  colnorm = cHACApK_unrm_d(ndl,pa_ref);

  cHACApK_minabsvalloc_d(pa_ref,&rownorm,&i_ref,ndl); // determine i_ref:=argmin ||pa_ref(1:ndl)||
  //    print*,'i_ref=',i_ref
  pb_ref = (double *) calloc(ndt,sizeof(double));
  if(pb_ref==NULL) {
    fprintf(stderr, "Error: cHACApK_acaplus: malloc pb_ref\n");
    goto error;
  }
  cHACApK_calc_vec(zab,zaa,ndt,ndl,k,i_ref,pb_ref,nstrtl,nstrtt,lod,i_bemv,lcol_msk,0,workspace);
  //  print*,'pb_ref=',pb_ref
  rownorm=cHACApK_unrm_d(ndt,pb_ref);

  apxnorm = 0.0; lstop_aca = 0;

  while((k<kmax) && (ntries_row>0 || ntries_col>0) && (ntries>0)) {
    ntries=ntries-1;
    pcol = zaa + ndl*k; // (k+1)th column of zaa
    prow = zab + ndt*k; // (k+1)th column of zab
    col_maxval = 0.0; cHACApK_maxabsvalloc_d(pa_ref,&col_maxval,&i,ndl);
    row_maxval = 0.0; cHACApK_maxabsvalloc_d(pb_ref,&row_maxval,&j,ndt);

    //    write(6,1000) 'i=',i,' i_ref=',i_ref,' j=',j,' j_ref=',j_ref

    if(row_maxval>col_maxval) {
      if(j!=j_ref) {
        cHACApK_calc_vec(zaa,zab,ndl,ndt,k,j,pcol,nstrtl,nstrtt,lod,i_bemv,lrow_msk,1,workspace);
      } else {
        for (il=0; il<ndl; il++) pcol[il]=pa_ref[il];
      }
      cHACApK_maxabsvalloc_d(pcol,&col_maxval,&i,ndl);

      if(col_maxval < ACA_EPS && k>=param[64]) {
        lstop_aca = 1;
        //         print*,'2***************lstop_aca==1***********************2'
      } else {
        cHACApK_calc_vec(zab,zaa,ndt,ndl,k,i,prow,nstrtl,nstrtt,lod,i_bemv,lcol_msk,0,workspace);
        if(fabs(pcol[i])>1.0e-20) {
          zinvmax=1.0/pcol[i];
        } else {
          k=max(k-1,0); break;
        }
        //        if(isnan(zinvmax))then
        //          print*,'1.0/pcol(i)=NaN',' k=',k
        //          exit
        //          stop
        //        endif
        for (il=0; il<ndl; il++) pcol[il]*=zinvmax;
      }
    } else {
      if(i!=i_ref) {
        cHACApK_calc_vec(zab,zaa,ndt,ndl,k,i,prow,nstrtl,nstrtt,lod,i_bemv,lcol_msk,0,workspace);
      } else {
        for (it=0; it<ndt; it++) prow[it]=pb_ref[it];
      }
      cHACApK_maxabsvalloc_d(prow,&row_maxval,&j,ndt);

      if(row_maxval < ACA_EPS && k>=param[64]) {
        lstop_aca = 1;
        //         print*,'3***************lstop_aca==1***********************3'
      } else {
        cHACApK_calc_vec(zaa,zab,ndl,ndt,k,j,pcol,nstrtl,nstrtt,lod,i_bemv,lrow_msk,1,workspace);
        if(fabs(prow[j])>1.0e-20) {
          zinvmax=1.0/prow[j];
        } else {
          k=max(k-1,0); break;
        }
        //        if(isnan(zinvmax))then
        //          print*,'1.0/prow(j)=NaN',' k=',k
        //          exit
        //          stop
        //        endif
        for (it=0; it<ndt; it++) prow[it]*=zinvmax;
      }
    }
    lrow_msk[i] = 1; lcol_msk[j] = 1;
    //    write(6,1000) 'i=',i,' i_ref=',i_ref,' j=',j,' j_ref=',j_ref

    if(i!=i_ref) {
      zinvmax = -pcol[i_ref];
      for (it=0; it<ndt; it++) pb_ref[it]+=prow[it]*zinvmax;
      rownorm = cHACApK_unrm_d(ndt,pb_ref);
    }
    if(i==i_ref || rownorm<ACA_EPS) {
      if(i==i_ref) ntries_row++;
      if(ntries_row>0) {
        rownorm = 0.0; i=i_ref;
        //        print*,'lrow_msk',lrow_msk
        while(i!=(i_ref+ndl-1)%ndl && rownorm<za_ACA_EPS && ntries_row>0) {
          //          print*,'i=',i,' ii=',mod((i_ref+ndl-2),ndl)+1
          if(lrow_msk[i]==0) {
            //            write(6,1000) 'i=',i
            cHACApK_calc_vec(zab,zaa,ndt,ndl,k+1,i,pb_ref,nstrtl,nstrtt,lod,i_bemv,lcol_msk,0,workspace);
            rownorm = cHACApK_unrm_d(ndt,pb_ref);
            if(rownorm<ACA_EPS) lrow_msk[i] = 1;
            ntries_row--;
          } else {
            rownorm = 0.0;
          }
          i=(i+1)%ndl;
        }
        i_ref=(i+ndl-1)%ndl;
      }
    }
    //    print*,'i_ref=',i_ref

    if(j!=j_ref) {
      zinvmax = -prow[j_ref];
      for (il=0; il<ndl; il++) pa_ref[il]+=pcol[il]*zinvmax;
      colnorm = cHACApK_unrm_d(ndl,pa_ref);
    }
    if(j==j_ref || colnorm<ACA_EPS) {
      if(j==j_ref) ntries_col++;
      if(ntries_col>0) {
        colnorm = 0.0; j=j_ref;
        //        print*,'lcol_msk',lcol_msk
        while(j!=(j_ref+ndt-1)%ndt && colnorm<za_ACA_EPS && ntries_col>0) {
          if(lcol_msk[j]==0) {
            cHACApK_calc_vec(zaa,zab,ndl,ndt,k+1,j,pa_ref,nstrtl,nstrtt,lod,i_bemv,lrow_msk,1,workspace);
            colnorm = cHACApK_unrm_d(ndl,pa_ref);
            if(colnorm<ACA_EPS) lcol_msk[j]=1;
            ntries_col--;
          } else {
            colnorm = 0.0;
          }
          j=(j+1)%ndt;
        }
        j_ref=(j+ndt-1)%ndt;
      }
    }

    //    write(6,2000) 'colnorm=',colnorm,' rownorm=',rownorm
    if(colnorm<ACA_EPS && rownorm<ACA_EPS && k>=param[64]) {
      lstop_aca=1; k++;
      //       print*,'1***************lstop_aca==1***********************1'
    }

    if(lstop_aca==0) {
      blknorm = (cHACApK_unrm_d(ndl,pcol)*cHACApK_unrm_d(ndt,prow));
      if(k == 0) {
        if(param[61]==1) {
          apxnorm = blknorm;
        } else if(param[61]==2 || param[61]==3) {
          apxnorm =znrm;
        } else {
          fprintf(stderr, "ERROR!:: invalid param[61]=%lf\n",param[61]);
          exit(EXIT_FAILURE);
        }
      } else {
        if(blknorm < apxnorm * eps &&
           rownorm < apxnorm * eps &&
           colnorm < apxnorm * eps &&
           k>=param[64]) lstop_aca = 1;
      }
    }
    if(0) {
      {
        printf("pcol\n");
        for (il=0; il<ndl; il++) printf("%lf\n",pcol[il]);
        printf("prow\n");
        for (it=0; it<ndt; it++) printf("%lf\n",prow[it]);
      }
    }
    if(lstop_aca==1 && k>=param[64]) break;
    k++;
  }
  //  if(k==kmax .or. ntries_row==0 .or. ntries_col==0 .or. ntries==0)then
  //    k=k-1
  //  endif

  if(k<param[64]) {
    {
      fprintf(stderr, "colnorm=%lf rownorm=%lf ACA_EPS=%lf\n",colnorm,rownorm,ACA_EPS);
      fprintf(stderr, "col_maxval=%lf, row_maxval=%lf\n",col_maxval,row_maxval);
      fprintf(stderr, "ntries_row=%d ntries_col=%d ntries=%d\n",ntries_row,ntries_col,ntries);
      fprintf(stderr, "k=%d\n",k);
      //    k=k-1; if(k<1) stop
      //    stop
    }
  }
  free(lrow_msk); free(lcol_msk); free(pa_ref); free(pb_ref); free(workspace);
  kacaplus=k;
  //  print*,'HACApK_acaplus=',HACApK_acaplus
  //  write(6,2000) 'blknorm=',blknorm/apxnorm,' colnorm=',colnorm/apxnorm,' rownorm=',rownorm/apxnorm
  //  if(nstrtt==         113) stop
  return kacaplus;
error:
  exit(EXIT_FAILURE);
}

/* Context and callback for parallel leaf fill via TaskManager */
typedef struct {
  st_cHACApK_leafmtx *st_lf;
  int kparam;
  double eps, ACA_EPS, znrmmat;
  double *param;
  int *lodl, *lodt;
  int i_bemv;
} fill_leaf_ctx;

static void fill_one_leaf_block(int idx, void *data) {
  fill_leaf_ctx *ctx = (fill_leaf_ctx*)data;
  int ip = idx + 1;  /* 0-based idx to 1-based ip */
  st_cHACApK_leafmtx *st_lf = ctx->st_lf;
  int kparam = ctx->kparam;
  double eps = ctx->eps;
  double ACA_EPS = ctx->ACA_EPS;
  double znrmmat = ctx->znrmmat;
  double *param = ctx->param;
  int *lodl = ctx->lodl;
  int *lodt = ctx->lodt;
  int i_bemv = ctx->i_bemv;

  int ndl   = st_lf[ip]->ndl;
  int ndt   = st_lf[ip]->ndt;
  int nstrtl= st_lf[ip]->nstrtl;
  int nstrtt= st_lf[ip]->nstrtt;
  int ltmtx = st_lf[ip]->ltmtx;

  /* Symmetric fill: strictly-lower leaves stay EMPTY (a1=a2=NULL, kt=0) -- matvec_sym skips
   * nstrtl > nstrtt and mirrors the upper leaf, so their fill would never be read. */
  if (g_sym_fill && nstrtl > nstrtt) { st_lf[ip]->kt = 0; return; }

  if(ltmtx==1) {
    /* Low-rank block: use ACA+ */
    double *zab = (double *) calloc(ndt*kparam, sizeof(double));
    double *zaa = (double *) calloc(ndl*kparam, sizeof(double));
    if(zab == NULL || zaa == NULL) {
      fprintf(stderr, "sub cHACApK_fill_leafmtx_hyp; zab,zaa Memory allocation failed !\n");
      fprintf(stderr, "ip=%d ndt=%d ndl=%d kparam=%d\n",ip,ndt,ndl,kparam);
      if(zab) free(zab);
      if(zaa) free(zaa);
      return;
    }

    int kt = 0;
    if(param[60]==1) {
      kt=cHACApK_aca(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,i_bemv,kparam,eps,znrmmat,ACA_EPS);
    } else if(param[60]==2) {
      kt=cHACApK_acaplus(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,i_bemv,kparam,eps,znrmmat,ACA_EPS);
    } else if(param[60]==3) {
      kt=cHACApK_SVD(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,i_bemv,kparam,eps,znrmmat,ACA_EPS);
    } else if(param[60]==5) {
      kt=cHACApK_RRQR(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,i_bemv,kparam,eps,znrmmat,ACA_EPS);
    } else {
      fprintf(stderr, "Only ACA and ACA+ is available! Set param[60]=1-5.\n");
      free(zab); free(zaa);
      return;
    }

    st_lf[ip]->kt=kt;
    st_lf[ip]->a1 = (double *) calloc(ndt*kt,sizeof(double));
    st_lf[ip]->a2 = (double *) calloc(ndl*kt,sizeof(double));
    if(st_lf[ip]->a1 == NULL || st_lf[ip]->a2 == NULL) {
      fprintf(stderr, "sub cHACApK_fill_leafmtx_hyp; a1,a2 Memory allocation failed !\n");
      fprintf(stderr, "ip=%d ndt=%d ndl=%d kt=%d\n",ip,ndt,ndl,kt);
      free(zab); free(zaa);
      return;
    }
    for (int il=0; il<ndt*kt; il++) st_lf[ip]->a1[il]=zab[il];
    for (int il=0; il<ndl*kt; il++) st_lf[ip]->a2[il]=zaa[il];
    free(zab); free(zaa);

  } else if(ltmtx==2) {
    /* Dense block */
    st_lf[ip]->a1 = (double *) calloc(ndt*ndl,sizeof(double));
    if(st_lf[ip]->a1 == NULL) {
      fprintf(stderr, "sub cHACApK_fill_leafmtx_hyp; a1 Memory allocation failed !\n");
      fprintf(stderr, "ip=%d ndt=%d ndl=%d\n",ip,ndt,ndl);
      return;
    }
    for (int il=0; il<ndl; il++) {
      int ill=il+nstrtl;
      for (int it=0; it<ndt; it++) {
        int itt=it+nstrtt;
        double val = cHACApK_entry_ij(lodl[ill],lodt[itt],i_bemv);
        st_lf[ip]->a1[it+ndt*il] = val;
      }
    }
  }
}

//***cHACApK_fill_leafmtx_hyp
// TaskManager parallelization over leaf blocks (replaces OpenMP dynamic scheduling)
void cHACApK_fill_leafmtx_hyp(
  st_cHACApK_leafmtx *st_lf,
  int i_bemv,
  double *param,
  double znrmmat,
  int *lpmd,
  int *lnmtx,
  int *lodl, // [nd]
  int *lodt, // [nd]
  int nd,
  int nlf,
  int *lnps,
  int *lnpe,
  int *lthr) // [0:nthr] - thread-to-block assignment array (unused with dynamic)
{
  int mpinr,mpilog,nrank,icomm,kparam;
  int ip;
  double eps,ACA_EPS;

  mpinr=lpmd[3]; mpilog=lpmd[4]; nrank=lpmd[2]; icomm=lpmd[1];
  eps=param[71]; ACA_EPS=param[72]*eps; kparam=(int)param[63];

  /* Parallel leaf fill using TaskManager (via C wrapper)
   * Each leaf block is filled independently */
  {
    fill_leaf_ctx fctx;
    fctx.st_lf = st_lf;
    fctx.kparam = kparam;
    fctx.eps = eps;
    fctx.ACA_EPS = ACA_EPS;
    fctx.znrmmat = znrmmat;
    fctx.param = param;
    fctx.lodl = lodl;
    fctx.lodt = lodt;
    fctx.i_bemv = i_bemv;

    hacapk_parallel_for(nlf, fill_one_leaf_block, &fctx);
  }

  /* Update lnps/lnpe (sequential post-processing) */
  for (int ip=1; ip<=nlf; ip++) {
    int ndl=st_lf[ip]->ndl;
    int nstrtl=st_lf[ip]->nstrtl;
    if(nstrtl < *lnps) *lnps=nstrtl;
    if(nstrtl+ndl > *lnpe) *lnpe=nstrtl+ndl;
  }
}

//***cHACApK_count_blrnmb
void cHACApK_count_blrnmb(
  st_cHACApK_cluster st_cltl,
  st_cHACApK_cluster st_cltt,
  double *param,
  int *lpmd,
  int *lnmtx,
  int nofc,
  int nffc,
  int *p_ndpth)
{
  int ndpth,ndl,ndt,nstrtl,nstrtt,nnsonl,nnsont,nleaf,nlmax,mdpth;
  int it,il;

  ndpth=*p_ndpth;

  ndl=st_cltl->nsize*nffc; ndt=st_cltt->nsize*nffc;
  nstrtl=st_cltl->nstrt; nstrtt=st_cltt->nstrt;
  nnsonl=st_cltl->nnson; nnsont=st_cltt->nnson;
  nleaf=param[42]+1; nlmax=param[22]*nofc;

  ndpth=ndpth+1;
  mdpth=param[53];

  if(ndpth==mdpth || (ndl<nleaf && ndt<nleaf)) {
    lnmtx[4]=lnmtx[4]+1;
    return;
  }

  lnmtx[3]=lnmtx[3]+1;

  if(ndl<nleaf) {
    for(it=1; it<=nnsont; it++) {
      cHACApK_count_blrnmb(st_cltl,st_cltt->pc_sons[it],param,lpmd,lnmtx,nofc,nffc,&ndpth);
      ndpth=ndpth-1;
    }
  } else if(ndt<nleaf) {
    for(il=1; il<=nnsonl; il++) {
      cHACApK_count_blrnmb(st_cltl->pc_sons[il],st_cltt,param,lpmd,lnmtx,nofc,nffc,&ndpth);
      ndpth=ndpth-1;
    }
  } else {
    for(il=1; il<=nnsonl; il++) {
      for(it=1; it<=nnsont; it++) {
        cHACApK_count_blrnmb(st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lpmd,lnmtx,nofc,nffc,&ndpth);
        ndpth=ndpth-1;
      }
    }
  }
  *p_ndpth=ndpth;
}

//***cHACApK_count_blrleaf
void cHACApK_count_blrleaf(
  st_cHACApK_leafmtx *st_leafmtx,
  st_cHACApK_cluster st_cltl,
  st_cHACApK_cluster st_cltt,
  double *param,
  int *lpmd,
  int *lnmtx,
  int nofc,
  int nffc,
  int *p_ndpth)
{
  int ndpth,ndl,ndt,nstrtl,nstrtt,nnsonl,nnsont,nleaf,nlmax,mdpth,iblnlf;
  int lnmtx2[3+1];
  int id,il,it,ibl;
  double zs,zdistlt,zeta;

  ndpth=*p_ndpth;

  ndl=st_cltl->nsize*nffc; ndt=st_cltt->nsize*nffc;
  nstrtl=st_cltl->nstrt; nstrtt=st_cltt->nstrt;
  nnsonl=st_cltl->nnson; nnsont=st_cltt->nnson;
  nleaf=param[42]+1; nlmax=param[22]*nofc;

  ndpth=ndpth+1;
  mdpth=param[53];

  // printf("ndl=%12d; ndt=%12d; nleaf=%12d\n",ndl,ndt,nleaf);

  if(ndpth==mdpth || (ndl<nleaf && ndt<nleaf)) {
    lnmtx[4]=lnmtx[4]+1;
    ibl=lnmtx[4];
    zs=0.0;
    for(id=1; id<=st_cltl->ndim; id++) {
      if(st_cltl->bmax[id]<st_cltt->bmin[id]) {
        zs=zs+(st_cltt->bmin[id]-st_cltl->bmax[id])*(st_cltt->bmin[id]-st_cltl->bmax[id]);
      } else if(st_cltt->bmax[id]<st_cltl->bmin[id]) {
        zs=zs+(st_cltl->bmin[id]-st_cltt->bmax[id])*(st_cltl->bmin[id]-st_cltt->bmax[id]);
      } else {
      }
    }
    zdistlt=sqrt(zs);
    zeta=param[51];

    if(st_cltl->zwdth<=zeta*zdistlt || st_cltt->zwdth<=zeta*zdistlt) {
      st_leafmtx[ibl]->nlf=1;
      lnmtx[1]=lnmtx[1]+1;
      *p_ndpth=ndpth;
      return;
    }

    lnmtx[3]=lnmtx[3]+1;
    iblnlf=0;
    for(il=1; il<=nnsonl; il++) {
      for(it=1; it<=nnsont; it++) {
        lnmtx2[1]=0; lnmtx2[2]=0; lnmtx2[3]=0;
        cHACApK_count_lntmx(st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lpmd,lnmtx2,nofc,nffc,&ndpth);
        lnmtx[1]=lnmtx[1]+lnmtx2[1]; lnmtx[2]=lnmtx[2]+lnmtx2[2]; lnmtx[3]=lnmtx[3]+lnmtx2[3];
        iblnlf=iblnlf+lnmtx2[1]+lnmtx2[2];
        ndpth=ndpth-1;
      }
    }
    st_leafmtx[ibl]->nlf=iblnlf;
    st_leafmtx[ibl]->st_lf = (st_cHACApK_leafmtx *) malloc(sizeof(st_cHACApK_leafmtx)*(iblnlf+1));
    if(st_leafmtx[ibl]->st_lf==NULL) {
      fprintf(stderr, "Error: cHACApK_count_blrleaf: malloc st_leafmtx[%d]->st_lf\n",ibl);
      goto error;
    }
    for(il=1; il<=iblnlf; il++) {
      st_leafmtx[ibl]->st_lf[il] = (st_cHACApK_leafmtx) calloc(1,sizeof(st_cHACApK_leafmtx_t));
      if(st_leafmtx[ibl]->st_lf[il]==NULL) {
        fprintf(stderr, "Error: cHACApK_count_blrleaf: malloc st_leafmtx[%d]->st_lf[il]\n",ibl,il);
        goto error;
      }
    }
    *p_ndpth=ndpth;
    return;
  }

  lnmtx[3]=lnmtx[3]+1;

  if(ndl<nleaf) {
    for(it=1; it<=nnsont; it++) {
      cHACApK_count_blrleaf(st_leafmtx,st_cltl,st_cltt->pc_sons[it],param,lpmd,lnmtx,nofc,nffc,&ndpth);
      ndpth=ndpth-1;
    }
  } else if(ndt<nleaf) {
    for(il=1; il<=nnsonl; il++) {
      cHACApK_count_blrleaf(st_leafmtx,st_cltl->pc_sons[il],st_cltt,param,lpmd,lnmtx,nofc,nffc,&ndpth);
      ndpth=ndpth-1;
    }
  } else {
    for(il=1; il<=nnsonl; il++) {
      for(it=1; it<=nnsont; it++) {
        cHACApK_count_blrleaf(st_leafmtx,st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lpmd,lnmtx,nofc,nffc,&ndpth);
        ndpth=ndpth-1;
      }
    }
  }
  *p_ndpth=ndpth;
  return;
error:
  exit(EXIT_FAILURE);
}

//***cHACApK_generate_blrleaf
void cHACApK_generate_blrleaf(
  st_cHACApK_leafmtx *st_leafmtx,
  st_cHACApK_cluster st_cltl,
  st_cHACApK_cluster st_cltt,
  double *param,
  int *lpmd,
  int *lnmtx,
  int nofc,
  int nffc,
  int *p_nlf,
  int *p_ndpth)
{
  int nlf,ndpth,ndl,ndt,nstrtl,nstrtt,nnsonl,nnsont,nleaf,nlmax,mdpth;
  int id,il,it,ibl,iblnlf;
  double zs,zdistlt,zeta;

  nlf=*p_nlf;
  ndpth=*p_ndpth;

  ndl=st_cltl->nsize*nffc; ndt=st_cltt->nsize*nffc;
  /* Convert element-based cluster start indices to DOF-based indices */
  nstrtl=(st_cltl->nstrt-1)*nffc+1; nstrtt=(st_cltt->nstrt-1)*nffc+1;
  nnsonl=st_cltl->nnson; nnsont=st_cltt->nnson;
  // printf("%12d %12d %12d %12d\n",nnsonl,ndl,nnsont,ndt);
  nleaf=param[42]+1; nlmax=param[22]*nofc;
  // printf("nleaf=%12d\n",nleaf); stop

  ndpth=ndpth+1;
  mdpth=param[53];

  if(ndpth==mdpth || (ndl<nleaf && ndt<nleaf)) {
    lnmtx[4]=lnmtx[4]+1;
    ibl=lnmtx[4];
    zs=0.0;
    for(id=1; id<=st_cltl->ndim; id++) {
      if(st_cltl->bmax[id]<st_cltt->bmin[id]) {
        zs=zs+(st_cltt->bmin[id]-st_cltl->bmax[id])*(st_cltt->bmin[id]-st_cltl->bmax[id]);
      } else if(st_cltt->bmax[id]<st_cltl->bmin[id]) {
        zs=zs+(st_cltl->bmin[id]-st_cltt->bmax[id])*(st_cltl->bmin[id]-st_cltt->bmax[id]);
      } else {
      }
    }
    // zdistlt=max(sqrt(zs)-st_cltl->zwdth/ndl-st_cltt->zwdth/ndt,0.0);
    zdistlt=sqrt(zs);
    zeta=param[51];

    nlf=nlf+1;
    st_leafmtx[nlf]->nstrtl=nstrtl; st_leafmtx[nlf]->ndl=ndl;
    st_leafmtx[nlf]->nstrtt=nstrtt; st_leafmtx[nlf]->ndt=ndt;
    if(st_cltl->zwdth<=zeta*zdistlt || st_cltt->zwdth<=zeta*zdistlt) {
      st_leafmtx[nlf]->kt=0;
      st_leafmtx[nlf]->ltmtx=1;
      // printf("ibl=%12d; iblnlf=%12d; ltmtx=%12d\n",ibl,1,1);
      *p_nlf=nlf;
      *p_ndpth=ndpth;
      return;
    } else {
      // st_leafmtx[ibl]->st_lf = (st_cHACApK_leafmtx *) malloc(sizeof(st_cHACApK_leafmtx)*(nnsonl*nnsont));
      iblnlf=0;
      for(il=1; il<=nnsonl; il++) {
        for(it=1; it<=nnsont; it++) {
          cHACApK_generate_leafmtx(st_leafmtx[ibl]->st_lf,st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lpmd,lnmtx,nofc,nffc,&iblnlf,&ndpth);
          ndpth=ndpth-1;
        }
      }
      st_leafmtx[nlf]->kt=0;
      st_leafmtx[nlf]->ltmtx=4;
      // printf("ibl=%12d; iblnlf=%12d; ltmtx=%12d\n",ibl,st_leafmtx[ibl]->nlf,4);
      *p_nlf=nlf;
      *p_ndpth=ndpth;
      return;

    }
  }

  if(ndl<nleaf) {
    for(it=1; it<=nnsont; it++) {
      cHACApK_generate_blrleaf(st_leafmtx,st_cltl,st_cltt->pc_sons[it],param,lpmd,lnmtx,nofc,nffc,&nlf,&ndpth);
      ndpth=ndpth-1;
    }
  } else if(ndt<nleaf) {
    for(il=1; il<=nnsonl; il++) {
      cHACApK_generate_blrleaf(st_leafmtx,st_cltl->pc_sons[il],st_cltt,param,lpmd,lnmtx,nofc,nffc,&nlf,&ndpth);
      ndpth=ndpth-1;
    }
  } else {
    for(il=1; il<=nnsonl; il++) {
      for(it=1; it<=nnsont; it++) {
        cHACApK_generate_blrleaf(st_leafmtx,st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lpmd,lnmtx,nofc,nffc,&nlf,&ndpth);
        ndpth=ndpth-1;
      }
    }
  }
  *p_nlf=nlf;
  *p_ndpth=ndpth;
}

//***cHACApK_count_lntmx
void cHACApK_count_lntmx(
  st_cHACApK_cluster st_cltl,
  st_cHACApK_cluster st_cltt,
  double *param,
  int *lpmd,
  int *lnmtx,
  int nofc,
  int nffc,
  int *p_ndpth)
{
  int ndpth,ndl,ndt,nstrtl,nstrtt,nnsonl,nnsont,nleaf,nlmax,mdpth;
  double zs,zdistlt,zeta;
  int id,il,it;

  ndpth=*p_ndpth;

  ndl=st_cltl->nsize*nffc; ndt=st_cltt->nsize*nffc;
  /* Convert element-based cluster start indices to DOF-based indices */
  nstrtl=(st_cltl->nstrt-1)*nffc+1; nstrtt=(st_cltt->nstrt-1)*nffc+1;
  nnsonl=st_cltl->nnson; nnsont=st_cltt->nnson;
  /* ELF-compatible: convert element-based leaf_size to DOF-based threshold */
  nleaf=((int)param[21]+1)*nffc; nlmax=(int)(param[22]*nofc)*nffc;

  // printf("ndl=%12d; nleaf=%12d",ndl,nleaf); stop

  ndpth=ndpth+1;
  mdpth=param[53];

  zs=0.0;
  for(id=1; id<=st_cltl->ndim; id++) {
    if(st_cltl->bmax[id]<st_cltt->bmin[id]) {
      zs=zs+(st_cltt->bmin[id]-st_cltl->bmax[id])*(st_cltt->bmin[id]-st_cltl->bmax[id]);
    } else if(st_cltt->bmax[id]<st_cltl->bmin[id]) {
      zs=zs+(st_cltl->bmin[id]-st_cltt->bmax[id])*(st_cltl->bmin[id]-st_cltt->bmax[id]);
    } else {
    }
  }
  zdistlt=sqrt(zs);
  zeta=param[51];

 if((st_cltl->zwdth<=zeta*zdistlt || st_cltt->zwdth<=zeta*zdistlt)
    && (ndl>=nleaf && ndt>=nleaf && ndl<=nlmax && ndt<=nlmax)
   ) {
    if(param[52]==0
       || param[52]==1 &&((nstrtl+ndl)!=nstrtt && (nstrtt+ndt)!=nstrtl)) {
      lnmtx[1]=lnmtx[1]+1;
      *p_ndpth=ndpth;
      return;
    } else if(param[52]/=1) {
      printf("Invalid admissiblity!; Set param[52]=0 or 1.\n");
      exit(EXIT_FAILURE);
    }
 }
 if(ndpth==mdpth || (nnsonl==0 || nnsont==0 || ndl<=nleaf || ndt<=nleaf)) {
   lnmtx[2]=lnmtx[2]+1;
   *p_ndpth=ndpth;
   return;
 }
 lnmtx[3]=lnmtx[3]+1;
 for(il=1; il<=nnsonl; il++) {
   for(it=1; it<=nnsont; it++) {
     cHACApK_count_lntmx(st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lpmd,lnmtx,nofc,nffc,&ndpth);
     ndpth=ndpth-1;
   }
 }
 *p_ndpth=ndpth;
}

//***cHACApK_generate_leafmtx
void cHACApK_generate_leafmtx(
  st_cHACApK_leafmtx *st_leafmtx,
  st_cHACApK_cluster st_cltl,
  st_cHACApK_cluster st_cltt,
  double *param,
  int *lpmd,
  int *lnmtx,
  int nofc,
  int nffc,
  int *p_nlf,
  int *p_ndpth)
{
  int nlf,ndpth,ndl,ndt,nstrtl,nstrtt,nnsonl,nnsont,nleaf,nlmax,mdpth;
  double zs,zdistlt,zeta;
  int id,il,it;

  nlf=*p_nlf;
  ndpth=*p_ndpth;

  ndl=st_cltl->nsize*nffc; ndt=st_cltt->nsize*nffc;
  /* Convert element-based cluster start indices to DOF-based indices */
  nstrtl=(st_cltl->nstrt-1)*nffc+1; nstrtt=(st_cltt->nstrt-1)*nffc+1;
  nnsonl=st_cltl->nnson; nnsont=st_cltt->nnson;
  // printf("%12d %12d %12d %12d\n",nnsonl,ndl,nnsont,ndt);
  /* ELF-compatible: convert element-based leaf_size to DOF-based threshold */
  nleaf=((int)param[21]+1)*nffc; nlmax=(int)(param[22]*nofc)*nffc;
  // printf("nlmax=%12d\n",nlmax); stop

  ndpth=ndpth+1;
  mdpth=param[53];

  zs=0.0;
  for(id=1; id<=st_cltl->ndim; id++) {
    if(st_cltl->bmax[id]<st_cltt->bmin[id]) {
      zs=zs+(st_cltt->bmin[id]-st_cltl->bmax[id])*(st_cltt->bmin[id]-st_cltl->bmax[id]);
    } else if(st_cltt->bmax[id]<st_cltl->bmin[id]) {
      zs=zs+(st_cltl->bmin[id]-st_cltt->bmax[id])*(st_cltl->bmin[id]-st_cltt->bmax[id]);
    } else {
    }
  }
  // zdistlt=max(sqrt(zs)-st_cltl->zwdth/ndl-st_cltt->zwdth/ndt,0.0);
  zdistlt=sqrt(zs);
  zeta=param[51];

 if((st_cltl->zwdth<=zeta*zdistlt || st_cltt->zwdth<=zeta*zdistlt)
    && (ndl>=nleaf && ndt>=nleaf && ndl<=nlmax && ndt<=nlmax)
   ) {
    if(param[52]==0
       || param[52]==1 &&((nstrtl+ndl)!=nstrtt && (nstrtt+ndt)!=nstrtl)) {
      nlf=nlf+1;
      st_leafmtx[nlf]->nstrtl=nstrtl; st_leafmtx[nlf]->ndl=ndl;
      st_leafmtx[nlf]->nstrtt=nstrtt; st_leafmtx[nlf]->ndt=ndt;
      st_leafmtx[nlf]->kt=0;
      st_leafmtx[nlf]->ltmtx=1;
      *p_nlf=nlf;
      *p_ndpth=ndpth;
      return;
    }
 }
 // if(ndpth==mdpth || (nnsonl==0 || nnsont==0 || (ndl<=nleaf && ndt<=nleaf))) {
 if(ndpth==mdpth || (nnsonl==0 || nnsont==0 || ndl<=nleaf || ndt<=nleaf)) {
 // if((nnsonl==0 || nnsont==0 || ndl<=nleaf || ndt<=nleaf)) {
   nlf=nlf+1;
   st_leafmtx[nlf]->nstrtl=nstrtl; st_leafmtx[nlf]->ndl=ndl;
   st_leafmtx[nlf]->nstrtt=nstrtt; st_leafmtx[nlf]->ndt=ndt;
   st_leafmtx[nlf]->ltmtx=2;
   // st_leafmtx[nlf]->a1 = (double *) malloc(sizeof(double)*(ndt*ndl));
   *p_nlf=nlf;
   *p_ndpth=ndpth;
   return;
 }
 for(il=1; il<=nnsonl; il++) {
   for(it=1; it<=nnsont; it++) {
     cHACApK_generate_leafmtx(st_leafmtx,st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lpmd,lnmtx,nofc,nffc,&nlf,&ndpth);
     ndpth=ndpth-1;
   }
 }
 *p_nlf=nlf;
 *p_ndpth=ndpth;
}

//***cHACApK_sort_leafmtx
void cHACApK_sort_leafmtx(
  st_cHACApK_leafmtx *st_leafmtx,
  int nlf)
{
  int ilp,ips,ip,il;
  cHACApK_qsort_row_leafmtx(st_leafmtx,1,nlf);
  ilp=1; ips=1;
  for(ip=1; ip<=nlf; ip++) {
    il=st_leafmtx[ip]->nstrtl;
    if(il<ilp) {    printf("Error!; HACApK_sort_leafmtx row_sort\n");
    } else if(il>ilp) {
      cHACApK_qsort_col_leafmtx(st_leafmtx,ips,ip-1);
      ilp=il; ips=ip;
    }
  }
  cHACApK_qsort_col_leafmtx(st_leafmtx,ips,nlf);
}

//***cHACApK_qsort_col_leafmtx
void cHACApK_qsort_col_leafmtx(
  st_cHACApK_leafmtx *st_leafmtx,
  int nlf_s,
  int nlf_e)
{
  st_cHACApK_leafmtx st_www;
  int nl,nr,nlr2,nlt,nrt,nlrt,nmid;

  if(nlf_s>=nlf_e) return;
  nl = nlf_s; nr = nlf_e; nlr2=nl+(nr-nl)/2;
  nlt=st_leafmtx[nl]->nstrtt; nrt=st_leafmtx[nr]->nstrtt; nlrt=st_leafmtx[nlr2]->nstrtt;
  nmid=cHACApK_med3(nlt,nrt,nlrt);
  // printf("nlf_s=%12d nlf_e=%12d nlr2=%12d nmid=%12d\n",nlf_s,nlf_e,nlr2,nmid);
  for(;;) {
    while(st_leafmtx[nl]->nstrtt < nmid) { nl=nl+1; }
    while(st_leafmtx[nr]->nstrtt > nmid) { nr=nr-1; }
    if(nl >= nr) break;
    st_www = st_leafmtx[nl]; st_leafmtx[nl] = st_leafmtx[nr]; st_leafmtx[nr] = st_www;
    nl=nl+1; nr=nr-1;
  }
  cHACApK_qsort_col_leafmtx(st_leafmtx,nlf_s,nl-1);
  cHACApK_qsort_col_leafmtx(st_leafmtx,nr+1 ,nlf_e);
}

//***cHACApK_qsort_row_leafmtx
void cHACApK_qsort_row_leafmtx(
  st_cHACApK_leafmtx *st_leafmtx,
  int nlf_s,
  int nlf_e)
{
  st_cHACApK_leafmtx st_www;
  int nl,nr,nlr2,nmid;

  if(nlf_s>=nlf_e) return;
  nl = nlf_s; nr = nlf_e; nlr2=nl+(nr-nl)/2;
  nmid=cHACApK_med3(st_leafmtx[nl]->nstrtl,st_leafmtx[nr]->nstrtl,st_leafmtx[nlr2]->nstrtl);
  // nmid=st_leafmtx[nlr]->nstrtl;
  // printf("nlf_s=%12d nlf_e=%12d nlr2=%12d nmid=%12d\n",nlf_s,nlf_e,nlr2,nmid);
  for(;;) {
    while(st_leafmtx[nl]->nstrtl < nmid) { nl=nl+1; }
    while(st_leafmtx[nr]->nstrtl > nmid) { nr=nr-1; }
    if(nl >= nr) break;
    st_www = st_leafmtx[nl]; st_leafmtx[nl] = st_leafmtx[nr]; st_leafmtx[nr] = st_www;
    nl=nl+1; nr=nr-1;
  }
  cHACApK_qsort_row_leafmtx(st_leafmtx,nlf_s,nl-1);
  cHACApK_qsort_row_leafmtx(st_leafmtx,nr+1 ,nlf_e);
}

//***cHACApK_free_st_clt
void cHACApK_free_st_clt(
  st_cHACApK_cluster st_clt)
{
  int nnson,ic;
  nnson=st_clt->nnson;
  for(ic=1; ic<=nnson; ic++) {
    cHACApK_free_st_clt(st_clt->pc_sons[ic]);
  }
  free(st_clt->bmin);
  free(st_clt->bmax);
  free(st_clt->pc_sons);
  free(st_clt);
}

//***cHACApK_generate_cluster
st_cHACApK_cluster cHACApK_generate_cluster(
  int *p_nmbr,
  int ndpth,
  int nstrt,
  int nsize,
  int ndim,
  int nson)
{
  st_cHACApK_cluster st_clt;
  int nmbr;
  nmbr=*p_nmbr;

  st_clt = (st_cHACApK_cluster) calloc(1,sizeof(st_cHACApK_cluster_t));
  if(st_clt==NULL) {
    fprintf(stderr, "Error: cHACApK_generate_cluster: malloc st_clt\n");
    goto error;
  }
  st_clt->bmin = NULL;
  st_clt->bmax = NULL;

  nmbr=nmbr+1;
  st_clt->nstrt=nstrt; st_clt->nsize=nsize; st_clt->ndim=ndim; st_clt->nnson=nson;
  st_clt->nmbr=nmbr; st_clt->ndpth=ndpth;
  st_clt->pc_sons = (st_cHACApK_cluster *) malloc(sizeof(st_cHACApK_cluster)*(nson+1));
  if(st_clt->pc_sons==NULL) {
    fprintf(stderr, "Error: cHACApK_generate_cluster: malloc st_clt->pc_sons\n");
    goto error;
  }

  *p_nmbr=nmbr;
  return st_clt;
error:
  exit(EXIT_FAILURE);
}

//***cHACApK_bndbox
void cHACApK_bndbox(
  st_cHACApK_cluster st_clt,
  double **zgmid_t, // 2D array [st_clt->ndim+1][nofc+1]
  int *lod,
  int nofc)
{
  int ic,l,ndim,id,il;
  double zwdth;

  for(ic=1; ic<=st_clt->nnson; ic++) {
    if(ic==1) { l=1; }
    else { l=l+st_clt->pc_sons[ic-1]->nsize; }
    cHACApK_bndbox(st_clt->pc_sons[ic],zgmid_t,&(lod[l-1]),nofc);
  }
  ndim=st_clt->ndim;
  st_clt->bmin = (double *) malloc(sizeof(double)*(ndim+1));
  st_clt->bmax = (double *) malloc(sizeof(double)*(ndim+1));
  if(st_clt->bmin==NULL || st_clt->bmax==NULL) {
    fprintf(stderr, "Error: cHACApK_bndbox: malloc st_clt->bmin st_clt->bmax\n");
    goto error;
  }
  if(st_clt->nnson == 0) {
    for(id=1; id<=ndim; id++) {
      st_clt->bmin[id]=zgmid_t[id][lod[1]]; st_clt->bmax[id]=zgmid_t[id][lod[1]];
    }
    for(id=1; id<=ndim; id++) {
      for(il=2; il<=st_clt->nsize; il++) {
        if(zgmid_t[id][lod[il]] < st_clt->bmin[id]) st_clt->bmin[id] = zgmid_t[id][lod[il]];
        if(st_clt->bmax[id] < zgmid_t[id][lod[il]]) st_clt->bmax[id] = zgmid_t[id][lod[il]];
      }
    }
  } else {
    for(id=1; id<=ndim; id++) {
      st_clt->bmin[id]=st_clt->pc_sons[1]->bmin[id];
      st_clt->bmax[id]=st_clt->pc_sons[1]->bmax[id];
    }
    for(il=2; il<=st_clt->nnson; il++) {
      for(id=1; id<=ndim; id++) {
        if(st_clt->pc_sons[il]->bmin[id] < st_clt->bmin[id]) st_clt->bmin[id]=st_clt->pc_sons[il]->bmin[id];
        if(st_clt->bmax[id] < st_clt->pc_sons[il]->bmax[id]) st_clt->bmax[id]=st_clt->pc_sons[il]->bmax[id];
      }
    }
  }
  zwdth=(st_clt->bmax[1]-st_clt->bmin[1])*(st_clt->bmax[1]-st_clt->bmin[1]);
  for(id=2; id<=ndim; id++) {
    zwdth=zwdth+(st_clt->bmax[id]-st_clt->bmin[id])*(st_clt->bmax[id]-st_clt->bmin[id]);
  }
  st_clt->zwdth=sqrt(zwdth);
  return;
error:
  exit(EXIT_FAILURE);
}

//***cHACApK_generate_cbitree
// ITERATIVE VERSION: Uses explicit stack instead of recursion
// to avoid stack overflow for large problems (>30,000 elements)
void cHACApK_generate_cbitree(
  st_cHACApK_cluster *p_st_clt,
  double **zgmid_t, // 2D array [ndim+1][md+1]
  double *param,
  int *lpmd,
  int *lod,
  int *p_ndpth,
  int ndscd,
  int nsrt,
  int nd,
  int md,
  int ndim,
  int *p_nclst)
{
  // Stack frame structure for iterative traversal
  typedef struct {
    st_cHACApK_cluster *p_result;  // Where to store the cluster
    int *lod_ptr;                   // Pointer into lod array
    int nsrt;                       // Start index
    int nd;                         // Number of elements
    int ndpth;                      // Current depth
    int phase;                      // 0=process, 1=after left, 2=after right
    st_cHACApK_cluster parent_clt;  // Parent cluster for storing result
    int nl;                         // Split position (saved for phase 1->2)
  } StackFrame;

  int minsz = param[21];
  int nclst = *p_nclst;
  int ndpth_init = *p_ndpth;

  // Estimate max depth: log2(nd/minsz) + safety margin
  int max_depth = 64;  // Should be enough for billions of elements

  // Allocate stack
  StackFrame *stack = (StackFrame*)malloc(sizeof(StackFrame) * max_depth);
  if (!stack) {
    fprintf(stderr, "Error: cHACApK_generate_cbitree: malloc stack failed\n");
    exit(EXIT_FAILURE);
  }

  // Allocate reusable min/max arrays (avoid malloc per call)
  double *zlmin = (double*)malloc(sizeof(double) * (ndim + 1));
  double *zlmax = (double*)malloc(sizeof(double) * (ndim + 1));
  if (!zlmin || !zlmax) {
    fprintf(stderr, "Error: cHACApK_generate_cbitree: malloc zlmin/zlmax failed\n");
    free(stack);
    exit(EXIT_FAILURE);
  }

  // Initialize stack with root node
  int sp = 0;
  stack[sp].p_result = p_st_clt;
  stack[sp].lod_ptr = lod;
  stack[sp].nsrt = nsrt;
  stack[sp].nd = nd;
  stack[sp].ndpth = ndpth_init;
  stack[sp].phase = 0;
  stack[sp].parent_clt = NULL;
  stack[sp].nl = 0;
  sp++;

  while (sp > 0) {
    sp--;
    StackFrame *f = &stack[sp];

    if (f->phase == 0) {
      // Phase 0: Process this node
      int cur_ndpth = f->ndpth + 1;
      int cur_nd = f->nd;
      int cur_nsrt = f->nsrt;
      int *cur_lod = f->lod_ptr;

      if (cur_nd <= minsz) {
        // Leaf node
        st_cHACApK_cluster st_clt = cHACApK_generate_cluster(&nclst, cur_ndpth, cur_nsrt, cur_nd, ndim, 0);
        st_clt->ndscd = cur_nd;
        *(f->p_result) = st_clt;
        // Done with this node, continue to next
      } else {
        // Internal node: need to split
        // Find bounding box
        int id, il;
        for (id = 1; id <= ndim; id++) {
          zlmin[id] = zgmid_t[id][cur_lod[1]];
          zlmax[id] = zlmin[id];
          for (il = 2; il <= cur_nd; il++) {
            double zg = zgmid_t[id][cur_lod[il]];
            if (zg < zlmin[id]) zlmin[id] = zg;
            else if (zlmax[id] < zg) zlmax[id] = zg;
          }
        }

        // Find split dimension (largest extent) -- needed for the BBOX path AND
        // as fallback if PCA fails (collinear points / dsyev info != 0).
        double zdiff = zlmax[1] - zlmin[1];
        int ncut = 1;
        for (id = 1; id <= ndim; id++) {
          double zidiff = zlmax[id] - zlmin[id];
          if (zidiff > zdiff) {
            zdiff = zidiff;
            ncut = id;
          }
        }
        double zlmid = (zlmax[ncut] + zlmin[ncut]) / 2.0;

        // Strategy dispatch: PCA split (better on flat / elongated geometries)
        // or BBOX midpoint split (historical default, all existing tests rely on it).
        int nl = -1;  /* will hold the 1-based "first right-side" pointer */
        if (cHACApK_get_cluster_strategy() == CHACAPK_CLUSTER_PCA) {
          int n_left = cHACApK_pca_split(zgmid_t, cur_lod, cur_nd, ndim);
          if (n_left > 0) {
            nl = n_left + 1;  /* 1-based first right-side index */
          }
          /* n_left == -1 (covariance singular, dsyev fail, or degenerate split)
           * -> fall through to BBOX path below */
        }
        if (nl < 0) {
          // BBOX midpoint partition (default / fallback)
          int nl_bb = 1, nr = cur_nd;
          while (nl_bb < nr) {
            while (nl_bb < cur_nd && zgmid_t[ncut][cur_lod[nl_bb]] <= zlmid) nl_bb++;
            while (nr >= 0 && zgmid_t[ncut][cur_lod[nr]] > zlmid) nr--;
            if (nl_bb < nr) {
              int nh = cur_lod[nl_bb];
              cur_lod[nl_bb] = cur_lod[nr];
              cur_lod[nr] = nh;
            }
          }
          nl = nl_bb;
        }

        // Create parent cluster
        st_cHACApK_cluster st_clt = cHACApK_generate_cluster(&nclst, cur_ndpth, cur_nsrt, cur_nd, ndim, 2);
        st_clt->ndscd = cur_nd;
        *(f->p_result) = st_clt;

        // Save state for returning after children
        f->phase = 1;
        f->parent_clt = st_clt;
        f->nl = nl;
        f->ndpth = cur_ndpth;
        sp++;  // Push this frame back

        // Check stack overflow
        if (sp >= max_depth - 2) {
          fprintf(stderr, "Error: cHACApK_generate_cbitree: stack overflow (depth=%d)\n", sp);
          free(zlmin);
          free(zlmax);
          free(stack);
          exit(EXIT_FAILURE);
        }

        // Push left child
        stack[sp].p_result = &(st_clt->pc_sons[1]);
        stack[sp].lod_ptr = cur_lod;
        stack[sp].nsrt = cur_nsrt;
        stack[sp].nd = nl - 1;
        stack[sp].ndpth = cur_ndpth;
        stack[sp].phase = 0;
        stack[sp].parent_clt = NULL;
        stack[sp].nl = 0;
        sp++;
      }
    } else if (f->phase == 1) {
      // Phase 1: After left child, now push right child
      f->phase = 2;
      sp++;  // Push this frame back

      st_cHACApK_cluster parent = f->parent_clt;
      int nl = f->nl;
      int cur_ndpth = f->ndpth;
      int cur_nd = f->nd;
      int cur_nsrt = f->nsrt;
      int *cur_lod = f->lod_ptr;

      // Push right child
      stack[sp].p_result = &(parent->pc_sons[2]);
      stack[sp].lod_ptr = &(cur_lod[nl - 1]);
      stack[sp].nsrt = cur_nsrt + nl - 1;
      stack[sp].nd = cur_nd - nl + 1;
      stack[sp].ndpth = cur_ndpth;
      stack[sp].phase = 0;
      stack[sp].parent_clt = NULL;
      stack[sp].nl = 0;
      sp++;
    }
    // Phase 2: After right child, just continue (parent already set)
  }

  // Cleanup
  free(zlmin);
  free(zlmax);
  free(stack);

  *p_nclst = nclst;
  *p_ndpth = ndpth_init;  // Restore original depth (caller expects this)
}
