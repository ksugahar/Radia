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

  cHACApK_active_entry_vector(ndp,ip,nstrtl,nstrtt,lod,i_bemv,lmsk,ld,vec);
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
  double *workspace, *row_bounds;

  double znrm, col_maxval, row_maxval, pivot, zeps;
  double approximation_norm = 0.0;
  long double approximation_norm2 = 0.0;
  int check_residual = 0;
  int verification_cursor = 0;
  int krank, kstop, k, ist, jst, istn, lstop_aca;
  int il, it;

  krank = (ndl < ndt) ? ndl : ndt;  /* min(ndl, ndt) */
  znrm = znrmmat * sqrt((double)ndl * (double)ndt);

  /* Absolute pivot cutoffs are not invariant under physical-unit scaling. */
  (void)pACA_EPS;

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
  row_bounds = (double *)malloc((size_t)ndl * sizeof(double));
  if (row_bounds == NULL) {
    fprintf(stderr, "Error: cHACApK_aca: malloc row_bounds\n");
    goto error;
  }
  for (il = 0; il < ndl; ++il) row_bounds[il] = HUGE_VAL;

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

    /* A small selected row/update does not certify the residual block.
     * Bounding every unused residual row by eps*||U V^T||_F/sqrt(ndl)
     * bounds the remaining Frobenius residual relative to the approximation,
     * modulo roundoff in pivoted rows/columns. Stop-time verification can
     * cost O(ndl*ndt) entries.
     */
    {
      double row_limit = eps * approximation_norm / sqrt((double)ndl);
      double row_norm = cHACApK_unrm_d(ndt, prow);
      row_bounds[ist] = nextafter(row_norm, HUGE_VAL);
      if (check_residual || row_norm <= row_limit) {
        int found = 0, offset;
        for (offset = 0; offset < ndl; ++offset) {
          int row = (verification_cursor + offset) % ndl;
          if (lrow_msk[row] || row_bounds[row] <= row_limit) continue;
          cHACApK_calc_vec(zab, zaa, ndt, ndl, k, row, prow,
                          nstrtl, nstrtt, lod, i_bemv, lcol_msk, 0, workspace);
          row_norm = cHACApK_unrm_d(ndt, prow);
          row_bounds[row] = nextafter(row_norm, HUGE_VAL);
          if (row_norm > row_limit) {
            ist = row;
            verification_cursor = (row + 1) % ndl;
            found = 1;
            break;
          }
        }
        if (!found) break;
        check_residual = 0;
      }
    }
    cHACApK_maxabsvallocm_d(prow, &row_maxval, &jst, ndt, lcol_msk);
    pivot = prow[jst];
    /* Divide directly: forming 1/pivot can overflow for tiny matrices. */
    for (it = 0; it < ndt; it++) prow[it] /= pivot;

    /* Compute column: pcol = A(:, jst) - approximation */
    cHACApK_calc_vec(zaa, zab, ndl, ndt, k, jst, pcol, nstrtl, nstrtt, lod, i_bemv, lrow_msk, 1, workspace);

    /* Mark used indices */
    lrow_msk[ist] = 1;
    lcol_msk[jst] = 1;

    /* Find next row index (max abs in pcol with mask) */
    cHACApK_maxabsvallocm_d(pcol, &col_maxval, &istn, ndl, lrow_msk);

    /* Compute approximation norm */
    {
      double right_norm = cHACApK_unrm_d(ndt, prow);
      zeps = cHACApK_unrm_d(ndl, pcol) * right_norm;
      /* Triangle inequality carries verified row bounds across updates.
       * No new entries are needed while the upper bound meets the gate.
       */
      for (il = 0; il < ndl; ++il) {
        if (!lrow_msk[il])
          row_bounds[il] = nextafter(row_bounds[il] + fabs(pcol[il]) * right_norm, HUGE_VAL);
      }
    }
    /* Update ||sum_j u_j v_j^T||_F, including cross terms. Using only the
     * first update over-tightens the residual gate as the rank increases.
     */
    {
      int previous;
      approximation_norm2 += (long double)zeps * zeps;
      for (previous = 0; previous < k; ++previous) {
        long double dot_left = 0.0, dot_right = 0.0;
        for (il = 0; il < ndl; ++il)
          dot_left += (long double)pcol[il] * zaa[il + ndl * previous];
        for (it = 0; it < ndt; ++it)
          dot_right += (long double)prow[it] * zab[it + ndt * previous];
        approximation_norm2 += 2.0L * dot_left * dot_right;
      }
      approximation_norm = sqrt((double)fmaxl(0.0L, approximation_norm2));
    }
    if ((int)param[61] == 1) znrm = approximation_norm;
    zeps = zeps / znrm;

    if (zeps < eps && k + 1 >= (int)param[64]) check_residual = 1;

    ist = istn;  /* Next row index */
    k++;
  }

  free(lrow_msk);
  free(lcol_msk);
  free(workspace);
  free(row_bounds);
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
        /* Inspect every row once, including the final cyclic candidate. */
        int scanned_rows = 0;
        while(scanned_rows++<ndl && rownorm<za_ACA_EPS && ntries_row>0) {
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
        /* Inspect every column once, including the final cyclic candidate. */
        int scanned_columns = 0;
        while(scanned_columns++<ndt && colnorm<za_ACA_EPS && ntries_col>0) {
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

