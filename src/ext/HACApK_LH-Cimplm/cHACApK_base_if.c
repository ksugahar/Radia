#include "cHACApK_base.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

static st_cHACApK_lcontrol st_ctl;
static st_cHACApK_leafmtxp st_leafmtxp;

//***cHACApK_copy_lcontrol_init
int cHACApK_copy_lcontrol_init()
{
  st_ctl=(st_cHACApK_lcontrol) calloc(1,sizeof(struct st_cHACApK_lcontrol));
  if(st_ctl==NULL) {
    printf("Error: cHACApK_copy_lcontrol_init: calloc st_ctl\n");
    goto error;
  }
  return 0;
error:
  return 1;
}

//***cHACApK_copy_lcontrol_finalize
int cHACApK_copy_lcontrol_finalize()
{
  free(st_ctl);
  return 0;
}

//***cHACApK_copy_lcontrol_f2c
int cHACApK_copy_lcontrol_f2c(
  int *lod,
  int *lsp,
  int *lnp,
  int *lthr,
  int *lpmd,
  double *param,
  double *time,
  int lf_umpi)
{
  st_ctl->lod=lod;
  st_ctl->lsp=lsp;
  st_ctl->lnp=lnp;
  st_ctl->lthr=lthr;
  st_ctl->lpmd=lpmd;
  st_ctl->param=param;
  st_ctl->time=time;
  st_ctl->lf_umpi=lf_umpi;
  return 0;
}

//***cHACApK_copy_lcontrol_c2f
int cHACApK_copy_lcontrol_c2f(
  int **lod,
  int **lsp,
  int **lnp,
  int **lthr,
  int **lpmd,
  double **param,
  double **time,
  int *lf_umpi)
{
  *lod=st_ctl->lod;
  *lsp=st_ctl->lsp;
  *lnp=st_ctl->lnp;
  *lthr=st_ctl->lthr;
  *lpmd=st_ctl->lpmd;
  *param=st_ctl->param;
  *time=st_ctl->time;
  *lf_umpi=st_ctl->lf_umpi;
  return 0;
}

//***cHACApK_copy_leafmtxp_init
int cHACApK_copy_leafmtxp_init()
{
  st_leafmtxp=(st_cHACApK_leafmtxp) calloc(1,sizeof(struct st_cHACApK_leafmtxp));
  if(st_leafmtxp==NULL) {
    printf("Error: cHACApK_copy_leafmtxp_init: calloc st_leafmtxp\n");
    goto error;
  }
  return 0;
error:
  return 1;
}

//***cHACApK_copy_leafmtxp_finalize
int cHACApK_copy_leafmtxp_finalize()
{
  for(int il=1; il<=st_leafmtxp->nlfl; il++) {
    free(st_leafmtxp->lnlfl2g_t[il]);
  }
  free(st_leafmtxp->lnlfl2g_t);
  free(st_leafmtxp->lbstrtl);
  free(st_leafmtxp->lbstrtt);
  free(st_leafmtxp->lbndl);
  free(st_leafmtxp->lbndt);
  free(st_leafmtxp->lbndlfs);
  free(st_leafmtxp->lbndtfs);
  free(st_leafmtxp->lbl2t);
  free(st_leafmtxp);
  return 0;
}

//***cHACApK_copy_leafmtxp_c2f_scalars
int cHACApK_copy_leafmtxp_c2f_scalars(
  int *nd,
  int *nlf,
  int *nlfkt,
  int *ktmax,
  int *nbl,
  int *nlfalt,
  int *nlfl,
  int *nlft,
  int *ndlfs,
  int *ndtfs,
  int *nrank_l,
  int *nrank_t,
  int *npgl)
{
  *nd=st_leafmtxp->nd;
  *nlf=st_leafmtxp->nlf;
  *nlfkt=st_leafmtxp->nlfkt;
  *ktmax=st_leafmtxp->ktmax;
  *nbl=st_leafmtxp->nbl;
  *nlfalt=st_leafmtxp->nlfalt;
  *nlfl=st_leafmtxp->nlfl;
  *nlft=st_leafmtxp->nlft;
  *ndlfs=st_leafmtxp->ndlfs;
  *ndtfs=st_leafmtxp->ndtfs;
  *nrank_l=st_ctl->lpmd[36];
  *nrank_t=st_ctl->lpmd[32];
  int nrank=st_ctl->lpmd[2];
  *npgl=(st_ctl->param[41]==0) ? sqrt((double)nrank) : st_ctl->param[41];
  return 0;
}

//***cHACApK_copy_leafmtxp_c2f_arrays
int cHACApK_copy_leafmtxp_c2f_arrays(
  int lnlfl2g[],
  int **lbstrtl,
  int **lbstrtt,
  int **lbndl,
  int **lbndt,
  int **lbndlfs,
  int **lbndtfs,
  int **lbl2t)
{
  for(int il=1; il<=st_leafmtxp->nlfl; il++) {
    for(int it=1; it<=st_leafmtxp->nlft; it++) {
      lnlfl2g[st_leafmtxp->nlft*(il-1)+(it-1)]=st_leafmtxp->lnlfl2g_t[il][it];
    }
  }
  *lbstrtl=st_leafmtxp->lbstrtl;
  *lbstrtt=st_leafmtxp->lbstrtt;
  *lbndl=st_leafmtxp->lbndl;
  *lbndt=st_leafmtxp->lbndt;
  *lbndlfs=st_leafmtxp->lbndlfs;
  *lbndtfs=st_leafmtxp->lbndtfs;
  *lbl2t=st_leafmtxp->lbl2t;
  return 0;
error:
  return 1;
}

//***cHACApK_copy_leafmtx_c2f_scalars
int cHACApK_copy_leafmtx_c2f_scalars(
  int ndpth,
  int idx[],
  int *ltmtx,
  int *kt,
  int *nstrtl,
  int *ndl,
  int *nstrtt,
  int *ndt,
  int *nlf)
{
  st_cHACApK_leafmtx st_lf;
  assert(ndpth>=1);
  st_lf=st_leafmtxp->st_lf[idx[0]];
  for(int il=1; il<ndpth; il++) {
    assert(0<idx[il] && idx[il]<=st_lf->nlf);
    if(st_lf->st_lf==NULL) goto error;
    st_lf=st_lf->st_lf[idx[il]];
  }
  *ltmtx=st_lf->ltmtx;
  *kt=st_lf->kt;
  *nstrtl=st_lf->nstrtl;
  *ndl=st_lf->ndl;
  *nstrtt=st_lf->nstrtt;
  *ndt=st_lf->ndt;
  *nlf=st_lf->nlf;
  return 0;
error:
  return 1;
}

//***cHACApK_generate_frame_blrleaf_if
int cHACApK_generate_frame_blrleaf_if(
  int i_bemv,
  double gmid[],
  int lnmtx[],
  int nofc,
  int nffc,
  int ndim)
{
  double **gmid_t;

  gmid_t=(double **) calloc(ndim+1,sizeof(double *));
  if(gmid_t==NULL) {
    printf("Error: cHACApK_generate_frame_blrleaf_c: calloc gmid_t\n");
    goto error;
  }
  for(int il=1; il<=ndim; il++) {
    gmid_t[il]=(double *) calloc(nofc+1,sizeof(double));
    if(gmid_t[il]==NULL) {
      printf("Error: cHACApK_generate_frame_blrleaf_c: calloc gmid_t[%d]\n",il);
      goto error;
    }
  }

  for(int il=1; il<=ndim; il++) {
    for(int it=1; it<=nofc; it++) {
      gmid_t[il][it]=gmid[nofc*(il-1)+(it-1)];
    }
  }
  cHACApK_generate_frame_blrleaf(st_leafmtxp,i_bemv,st_ctl,gmid_t,lnmtx,nofc,nffc,ndim);

  for(int il=1; il<=ndim; il++) {
    free(gmid_t[il]);
  }
  free(gmid_t);
  return 0;
error:
  return 1;
}
