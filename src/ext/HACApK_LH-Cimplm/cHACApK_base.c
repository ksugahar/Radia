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
#include "hacapk_log.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi_stub.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* Windows defines min/max as macros - undefine to avoid conflicts */
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

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
    st_leafmtx[il] = (st_cHACApK_leafmtx) calloc(1,sizeof(struct st_cHACApK_leafmtx));
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
  int ld)
{
  double *zz;
  int ii,ill,itt,it,il;

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
  zz = (double *) calloc(k,sizeof(double));
  if(zz == NULL) {
#pragma omp critical
    printf("sub cHACApK_calc_vec; zz allocation failed !\n");
    exit(EXIT_FAILURE);
  }
  for (it=0; it<k; it++) zz[it]=zab[ip+ndt*it];
  cHACApK_adotsub_dsm(vec,zaa,zz,ndp,k,ndp);
  for (il=0; il<ndp; il++) {
    if(lmsk[il]==1) vec[il]=0.0;
  }
  free(zz);
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
  k = 0;

  j_ref=0; // arbitrary j_ref
  pa_ref = (double *) calloc(ndl,sizeof(double));
  if(pa_ref==NULL) {
    fprintf(stderr, "Error: cHACApK_acaplus: malloc pa_ref\n");
    goto error;
  }

  cHACApK_calc_vec(zaa,zab,ndl,ndt,k,j_ref,pa_ref,nstrtl,nstrtt,lod,i_bemv,lrow_msk,1);
  //  print*,'pa_ref=',pa_ref
  colnorm = cHACApK_unrm_d(ndl,pa_ref);

  cHACApK_minabsvalloc_d(pa_ref,&rownorm,&i_ref,ndl); // determine i_ref:=argmin ||pa_ref(1:ndl)||
  //    print*,'i_ref=',i_ref
  pb_ref = (double *) calloc(ndt,sizeof(double));
  if(pb_ref==NULL) {
    fprintf(stderr, "Error: cHACApK_acaplus: malloc pb_ref\n");
    goto error;
  }
  cHACApK_calc_vec(zab,zaa,ndt,ndl,k,i_ref,pb_ref,nstrtl,nstrtt,lod,i_bemv,lcol_msk,0);
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
        cHACApK_calc_vec(zaa,zab,ndl,ndt,k,j,pcol,nstrtl,nstrtt,lod,i_bemv,lrow_msk,1);
      } else {
        for (il=0; il<ndl; il++) pcol[il]=pa_ref[il];
      }
      cHACApK_maxabsvalloc_d(pcol,&col_maxval,&i,ndl);

      if(col_maxval < ACA_EPS && k>=param[64]) {
        lstop_aca = 1; 
        //         print*,'2***************lstop_aca==1***********************2'
      } else {
        cHACApK_calc_vec(zab,zaa,ndt,ndl,k,i,prow,nstrtl,nstrtt,lod,i_bemv,lcol_msk,0);
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
        cHACApK_calc_vec(zab,zaa,ndt,ndl,k,i,prow,nstrtl,nstrtt,lod,i_bemv,lcol_msk,0);
      } else {
        for (it=0; it<ndt; it++) prow[it]=pb_ref[it];  
      }
      cHACApK_maxabsvalloc_d(prow,&row_maxval,&j,ndt);

      if(row_maxval < ACA_EPS && k>=param[64]) {
        lstop_aca = 1;
        //         print*,'3***************lstop_aca==1***********************3'
      } else {
        cHACApK_calc_vec(zaa,zab,ndl,ndt,k,j,pcol,nstrtl,nstrtt,lod,i_bemv,lrow_msk,1);
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
            cHACApK_calc_vec(zab,zaa,ndt,ndl,k+1,i,pb_ref,nstrtl,nstrtt,lod,i_bemv,lcol_msk,0);
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
            cHACApK_calc_vec(zaa,zab,ndl,ndt,k+1,j,pa_ref,nstrtl,nstrtt,lod,i_bemv,lrow_msk,1);
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
#pragma omp critical
          printf("ERROR!:: invalid param[61]=%lf\n",param[61]);
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
#pragma omp critical
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
#pragma omp critical
    {
      printf("colnorm=%lf rownorm=%lf ACA_EPS=%lf\n",colnorm,rownorm,ACA_EPS);
      printf("col_maxval=%lf, row_maxval=%lf\n",col_maxval,row_maxval);
      printf("ntries_row=%d ntries_col=%d ntries=%d\n",ntries_row,ntries_col,ntries);
      printf("k=%d\n",k);
      //    k=k-1; if(k<1) stop
      //    stop
    }
  }
  free(lrow_msk); free(lcol_msk); free(pa_ref); free(pb_ref);
  kacaplus=k;
  //  print*,'HACApK_acaplus=',HACApK_acaplus
  //  write(6,2000) 'blknorm=',blknorm/apxnorm,' colnorm=',colnorm/apxnorm,' rownorm=',rownorm/apxnorm
  //  if(nstrtt==         113) stop
  return kacaplus;
error:
  exit(EXIT_FAILURE);
}

//***cHACApK_fill_leafmtx_hyp
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
  int *lthr) // [0:]
{
  // type(st_HACApK_leafmtxp) ::  st_leafmtxp
  double *zab,*zaa;
  int mpinr,mpilog,nrank,icomm,kparam;
  int ith,nthr,ith1,nths,nthe,ierr,ip,ndl,ndt,ns,nstrtl,nstrtt,ltmtx,kt;
  int il,it,ill,itt;
  double eps,ACA_EPS;

  mpinr=lpmd[3]; mpilog=lpmd[4]; nrank=lpmd[2]; icomm=lpmd[1];
  eps=param[71]; ACA_EPS=param[72]*eps; kparam=param[63];

  /* Force single thread - OpenMP completely disabled */
  ith = 0;
  nthr = 1;
  ith1 = ith + 1;
  nths = lthr[ith];
  nthe = lthr[ith1] - 1;

  ierr = 0;
  for (ip = nths; ip <= nthe; ip++) {
      ndl   =st_lf[ip]->ndl   ; ndt   =st_lf[ip]->ndt   ; ns=ndl*ndt;
      nstrtl=st_lf[ip]->nstrtl; nstrtt=st_lf[ip]->nstrtt; ltmtx=st_lf[ip]->ltmtx;
      // write(mpilog,1000) 'sub HACApK_fill_leafmtx_hyp; ip=',ip,'; ndl=',ndl,'; ndt=',ndt,'; mpinr',mpinr

      if(ltmtx==1) {
        zab = (double *) calloc(ndt*kparam,sizeof(double *));
        zaa = (double *) calloc(ndl*kparam,sizeof(double *));
        if(zab == NULL || zaa == NULL) {
          printf("sub cHACApK_fill_leafmtx_hyp; zab,zaa Memory allocation failed !\n");
          printf("ip=%d ndt=%d ndl=%d kparam=%d\n",ip,ndt,ndl,kparam);
          exit(10);
        }
        if(param[60]==1) {
          // kt=HACApK_aca(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,st_bemv,kparam,eps,znrmmat,ACA_EPS)
        } else if(param[60]==2) {
          kt=cHACApK_acaplus(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,i_bemv,kparam,eps,znrmmat,ACA_EPS);
        } else if(param[60]==3) {
          kt=cHACApK_SVD(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,i_bemv,kparam,eps,znrmmat,ACA_EPS);
        } else if(param[60]==4) {
          // kt=HACApK_rSVD(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,st_bemv,kparam,eps,znrmmat,ACA_EPS)
        } else if(param[60]==5) {
          kt=cHACApK_RRQR(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,i_bemv,kparam,eps,znrmmat,ACA_EPS);
        } else {
          printf("Only ACA and ACA+ is available! Set param[60]=1-5.\n");
          exit(EXIT_FAILURE);
        }
        if(kt>kparam-1) {
          //!!$omp critical
          //        write(*,1000) 'WARNING: Insufficient k: kt=',kt,', kparam=',kparam, &
          //                      ' nstrtl=',nstrtl,' nstrtt=',nstrtt,' ndl=',ndl,' ndt=',ndt
          //!!$omp end critical
        }
        st_lf[ip]->kt=kt;
        st_lf[ip]->a1 = (double *) calloc(ndt*kt,sizeof(double));
        st_lf[ip]->a2 = (double *) calloc(ndl*kt,sizeof(double));
        if(st_lf[ip]->a1 == NULL || st_lf[ip]->a2 == NULL) {
          printf("sub cHACApK_fill_leafmtx_hyp; a1,a2 Memory allocation failed !\n");
          printf("ip=%d ndt=%d ndl=%d kt=%d\n",ip,ndt,ndl,kt);
          exit(20);
        }
        for (il=0; il<ndt*kt; il++) st_lf[ip]->a1[il]=zab[il];
        for (il=0; il<ndl*kt; il++) st_lf[ip]->a2[il]=zaa[il];
        free(zab); free(zaa);
      } else if(ltmtx==2) {
        st_lf[ip]->a1 = (double *) calloc(ndt*ndl,sizeof(double));
        if(st_lf[ip]->a1 == NULL) {
          printf("sub cHACApK_fill_leafmtx_hyp; a1 Memory allocation failed !\n");
          printf("ip=%d ndt=%d ndl=%d\n",ip,ndt,ndl);
          exit(30);
        }
        for (il=0; il<ndl; il++) {
          ill=il+nstrtl;
          for (it=0; it<ndt; it++) {
            itt=it+nstrtt;
            double val = cHACApK_entry_ij(lodl[ill],lodt[itt],i_bemv);
            st_lf[ip]->a1[it+ndt*il] = val;
          }
        }
      } else {
        printf("HACApK_fill_leafmtx_hyp; ip=%d ltmtx=%d\n",ip,ltmtx);
      }
  }
  for (ip=1; ip<=nlf; ip++) {
    ndl=st_lf[ip]->ndl; nstrtl=st_lf[ip]->nstrtl;
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
      st_leafmtx[ibl]->st_lf[il] = (st_cHACApK_leafmtx) calloc(1,sizeof(struct st_cHACApK_leafmtx));
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
  nleaf=param[21]+1; nlmax=param[22]*nofc;

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
 // if(ndpth==mdpth || (nnsonl==0 || nnsont==0 || (ndl<=nleaf && ndt<=nleaf))) {
 if(ndpth==mdpth || (nnsonl==0 || nnsont==0 || ndl<=nleaf || ndt<=nleaf)) {
 // if((nnsonl==0 || nnsont==0 || ndl<=nleaf || ndt<=nleaf)) {
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
  nleaf=param[21]+1; nlmax=param[22]*nofc;
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

  st_clt = (st_cHACApK_cluster) calloc(1,sizeof(struct st_cHACApK_cluster));
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
  st_cHACApK_cluster st_clt;
  double *zlmin,*zlmax;
  int ndpth,nclst,minsz,nson,id,il,ncut,nl,nr,nh,nsrt1,nd1;
  double zdiff,zlmid,zg,zidiff;

  ndpth=*p_ndpth;
  nclst=*p_nclst;

  minsz=param[21];
  // minsz=param[21]/4+1:
  ndpth=ndpth+1;
  // ndscd=ndscd+1:
  // if(i>26) stop
  // printf("\n");
  // printf("nsrt=%12d nd=%12d\n",nsrt,nd);
  if(nd <= minsz) {
    nson=0;
    // nclst=nclst+1;
    st_clt=cHACApK_generate_cluster(&nclst,ndpth,nsrt,nd,ndim,nson);
  } else {
    zlmin = (double *) malloc(sizeof(double)*(ndim+1));
    zlmax = (double *) malloc(sizeof(double)*(ndim+1));
    if(zlmin==NULL || zlmax==NULL) {
      fprintf(stderr, "Error: cHACApK_generate_cbitree: malloc zlmin zlmax\n");
      goto error;
    }
    for(id=1; id<=ndim; id++) {
      zlmin[id]=zgmid_t[id][lod[1]]; zlmax[id]=zlmin[id];
      for(il=2; il<=nd; il++) {
        zg=zgmid_t[id][lod[il]];
        if     (zg<zlmin[id]) { zlmin[id]=zg; }
        else if(zlmax[id]<zg) { zlmax[id]=zg; }
      }
    }
    // printf("zlmin=%21.6lf\n",zlmin);
    // printf("zlmax=%21.6lf\n",zlmax);

    zdiff=zlmax[1]-zlmin[1]; ncut = 1;
    for(id=1; id<=ndim; id++) {
      zidiff=zlmax[id]-zlmin[id];
      if(zidiff>zdiff) {
        zdiff =zidiff; ncut=id;
      }
    }
    zlmid= (zlmax[ncut]+zlmin[ncut])/2;
    // printf("ncut=%12d; zlmid=%21.6lf\n",ncut,zlmid);

    nl = 1; nr = nd;
    while(nl < nr) {
      while(nl < nd && zgmid_t[ncut][lod[nl]] <= zlmid) { nl=nl+1; }
      while(nr >= 0 && zgmid_t[ncut][lod[nr]] > zlmid) { nr=nr-1; }
      if(nl < nr) { nh = lod[nl]; lod[nl] = lod[nr]; lod[nr] = nh; }
    }

    // printf("nd=%12d;ncut=%12d; nsrt=%12d; nl=%12d\n",nd,ncut,nsrt,nl);

    nson=2;
    st_clt=cHACApK_generate_cluster(&nclst,ndpth,nsrt,nd,ndim,nson);
    nsrt1=nsrt; nd1=nl-1;
    cHACApK_generate_cbitree(&(st_clt->pc_sons[1]),zgmid_t,param,lpmd,lod,&ndpth,ndscd,nsrt1,nd1,md,ndim,&nclst);
    ndpth=ndpth-1;
    // ndscd=ndscd+st_clt->pc_sons[1].ndscd;
    nsrt1=nsrt+nl-1; nd1=nd-nl+1;
    cHACApK_generate_cbitree(&(st_clt->pc_sons[2]),zgmid_t,param,lpmd,&(lod[nl-1]),&ndpth,ndscd,nsrt1,nd1,md,ndim,&nclst);
    ndpth=ndpth-1;
    // ndscd=ndscd+st_clt->pc_sons[2].ndscd;
  }
  st_clt->ndscd=nd;
  *p_st_clt=st_clt;
  *p_ndpth=ndpth;
  *p_nclst=nclst;
  return;
error:
  exit(EXIT_FAILURE);
}
