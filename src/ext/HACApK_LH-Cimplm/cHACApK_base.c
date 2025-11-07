#include "cHACApK_base.h"
#include "cHACApK_calc_entry_ij.h"
#include "cHACApK_lib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

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
  for(il=1; il<=nblall; il++) {
    st_leafmtx[il] = (st_cHACApK_leafmtx) calloc(1,sizeof(struct st_cHACApK_leafmtx));
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
  for(il=1; il<=nlfl; il++) {
    st_leafmtxp->lnlfl2g_t[il] = (int64_t *) malloc(sizeof(int64_t)*(nlft+1));
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
    for(il=1; il<=iblnlf; il++) {
      st_leafmtx[ibl]->st_lf[il] = (st_cHACApK_leafmtx) calloc(1,sizeof(struct st_cHACApK_leafmtx));
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
  nstrtl=st_cltl->nstrt; nstrtt=st_cltt->nstrt;
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
  nstrtl=st_cltl->nstrt; nstrtt=st_cltt->nstrt;
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
  nstrtl=st_cltl->nstrt; nstrtt=st_cltt->nstrt;
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
  st_clt->bmin = NULL;
  st_clt->bmax = NULL;

  nmbr=nmbr+1;
  st_clt->nstrt=nstrt; st_clt->nsize=nsize; st_clt->ndim=ndim; st_clt->nnson=nson;
  st_clt->nmbr=nmbr; st_clt->ndpth=ndpth;
  st_clt->pc_sons = (st_cHACApK_cluster *) malloc(sizeof(st_cHACApK_cluster)*(nson+1));

  *p_nmbr=nmbr;
  return st_clt;
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
}
