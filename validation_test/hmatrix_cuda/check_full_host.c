/* Full existing HACApK cluster/leaf/fill integration with serial host control.
 * A synthetic smooth kernel exercises infrastructure, not HDiv physics. */
#include "cHACApK_base.h"
#include "cHACApK_pca_cluster.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cHACApK_dense_leaf_sink.h"

int HACApK_build_hmatrix_wrapper(void*,void*,double*,int,int,int,double,int,double,int,int);
void HACApK_matvec_wrapper(void*,void*,const double*,double*,int);
void HACApK_free_hmatrix_wrapper(void*,void*);
void HACApK_free_leafmtxp(void*);
static double coordinates[132*3];
static size_t entry_calls=0;
double cHACApK_entry_ij(int i,int j,int kind) {
    double d=1.; int k; (void)kind;
    ++entry_calls;
    for(k=0;k<3;++k) { double v=coordinates[3*(i-1)+k]-coordinates[3*(j-1)+k]; d+=v*v; }
    return 1./sqrt(d)+(i==j ? 2. : 0.);
}
void cHACApK_entry_batch(int n,const int*r,const int*c,int kind,double*v) {
    int k; for(k=0;k<n;++k) v[k]=cHACApK_entry_ij(r[k],c[k],kind);
}
#ifdef RAD_HACAPK_FRAME_ADMISSION
#include "cHACApK_frame_admission.h"
static int admission_enabled=0,admission_calls=0,admission_invalid=0;
int cHACApK_frame_admission(st_cHACApK_leafmtx const *leaves,int count,int nd,int rank_cap) {
    int i;++admission_calls;
    if(!leaves||count<1||nd!=132||rank_cap!=32) admission_invalid=1;
    for(i=1;!admission_invalid&&i<=count;++i) {
        const st_cHACApK_leafmtx_t *leaf=leaves[i];
        if(!leaf||leaf->a1||leaf->a2||leaf->kt||leaf->ndl<1||leaf->ndt<1||leaf->nstrtl<1||leaf->nstrtt<1||leaf->nstrtl>nd||leaf->nstrtt>nd||leaf->ndl>nd-leaf->nstrtl+1||leaf->ndt>nd-leaf->nstrtt+1) admission_invalid=1;
    }
    return admission_enabled&&!admission_invalid;
}
#endif
#ifndef RAD_HACAPK_DENSE_LEAF_SINK
int main(void) {
    int strategy,i,j,lowrank=0; double worst=0.;
    for(i=0;i<132;++i) { coordinates[3*i]=(i<66 ? 0. : 20.)+.01*(i%66); coordinates[3*i+1]=.03*(i%7); coordinates[3*i+2]=.02*(i%11); }
#ifdef RAD_HACAPK_FRAME_ADMISSION
    /* Reject repeatedly before entry evaluation, then retry normal assembly. */
    for(strategy=0;strategy<3;++strategy) {
        st_cHACApK_leafmtxp lp=calloc(1,sizeof(*lp));st_cHACApK_lcontrol ctl=calloc(1,sizeof(*ctl));
        if(!lp||!ctl)return 20;entry_calls=0;
        if(HACApK_build_hmatrix_wrapper(lp,ctl,coordinates,132,1,3,1e-10,64,2.,32,0)!=-2)return 21;
        if(entry_calls||lp->st_lf||lp->nlf||lp->st_clt_root||admission_invalid)return 22;
        HACApK_free_hmatrix_wrapper(lp,ctl);HACApK_free_leafmtxp(lp);free(ctl);
    }
    admission_enabled=1;
#endif
    for(strategy=0;strategy<2;++strategy) {
        st_cHACApK_leafmtxp lp=calloc(1,sizeof(*lp));
        st_cHACApK_lcontrol ctl=calloc(1,sizeof(*ctl));
        double x[132],y[132],diff=0.,norm=0.;
        if(!lp || !ctl) return 2;
        cHACApK_set_cluster_strategy(strategy);
        if(HACApK_build_hmatrix_wrapper(lp,ctl,coordinates,132,1,3,1e-10,64,2.,32,0)) return 3;
        if((int)ctl->param[60]!=1) return 4; /* Preserve production ACA default. */
        for(i=0;i<132;++i) x[i]=sin(.17*i)+.3;
        HACApK_matvec_wrapper(lp,ctl,x,y,132);
        for(i=0;i<132;++i) { double ref=0.; for(j=0;j<132;++j) ref+=cHACApK_entry_ij(i+1,j+1,0)*x[j]; diff+=(y[i]-ref)*(y[i]-ref); norm+=ref*ref; }
        diff=sqrt(diff/norm); if(!isfinite(diff)||diff>1e-8) return 5;
        if(diff>worst) worst=diff; lowrank+=lp->nlfkt;
        HACApK_free_hmatrix_wrapper(lp,ctl); HACApK_free_leafmtxp(lp); free(ctl);
    }
    cHACApK_set_cluster_strategy(CHACAPK_CLUSTER_BBOX);
    if(!lowrank) return 6;
#ifdef RAD_HACAPK_FRAME_ADMISSION
    if(admission_invalid||admission_calls!=5)return 23;
    printf("{\"frame_rejections\":3,\"frame_acceptances\":2,\"entries_before_rejection\":0,\"output_ownership_unpublished_on_rejection\":true,\"maximum_relative_difference\":%.17g}\n",worst);
    return 0;
#endif
    printf("{\"cases\":2,\"maximum_relative_difference\":%.17g,\"lowrank_leaves\":%d,\"taskmanager_started\":false}\n",worst,lowrank);
    return 0;
}
#else
void HACApK_matvec_sym_wrapper(void*,void*,const double*,double*,int);
static int sink_enabled=0,sink_count=0,sink_invalid=0;
static size_t deferred_entries=0;
static int captured_rows[132*132],captured_columns[132*132];
int cHACApK_dense_leaf_sink(const st_cHACApK_leafmtx_t *leaf,const int *rows,const int *columns,int problem) {
    int i,j;if(!sink_enabled)return 0;
    if(!leaf||!rows||!columns||problem||leaf->ltmtx!=2||leaf->kt||leaf->a1||leaf->a2||
       leaf->ndl<=0||leaf->ndt<=0||leaf->nstrtl<1||leaf->nstrtt<1||
       leaf->ndl>133-leaf->nstrtl||leaf->ndt>133-leaf->nstrtt) {sink_invalid=1;return 1;}
    for(i=0;i<leaf->ndl;++i)for(j=0;j<leaf->ndt;++j) {
        int r=rows[leaf->nstrtl+i],c=columns[leaf->nstrtt+j];
        if(r<1||r>132||c<1||c>132||deferred_entries>=132*132){sink_invalid=1;return 1;}
        captured_rows[deferred_entries]=r;captured_columns[deferred_entries]=c;++deferred_entries;
    }
    ++sink_count;return 1;
}
static void release(st_cHACApK_leafmtxp lp,st_cHACApK_lcontrol ctl) {
    HACApK_free_hmatrix_wrapper(lp,ctl);HACApK_free_leafmtxp(lp);free(ctl);
}
int main(void) {
    int strategy,cap,symmetric,i,j,case_count=0,total_sinks=0,permuted=0;size_t removed=0;double worst=0.;
    for(i=0;i<132;++i){int point=(37*i)%132;coordinates[3*i]=(point<66?0.:20.)+.01*(point%66);coordinates[3*i+1]=.03*(point%7);coordinates[3*i+2]=.02*(point%11);}
    for(strategy=0;strategy<2;++strategy)for(cap=1;cap<=32;cap+=31)for(symmetric=0;symmetric<2;++symmetric) {
        st_cHACApK_leafmtxp ref=calloc(1,sizeof(*ref)),lp=calloc(1,sizeof(*lp));
        st_cHACApK_lcontrol cr=calloc(1,sizeof(*cr)),ctl=calloc(1,sizeof(*ctl));
        size_t reference_calls,without_dense,cursor=0;double x[132],y[132],diff=0.,norm=0.;
        if(!ref||!lp||!cr||!ctl)return 20;
        cHACApK_set_cluster_strategy(strategy);cHACApK_set_sym_fill(symmetric);
        sink_enabled=0;entry_calls=0;
        if(HACApK_build_hmatrix_wrapper(ref,cr,coordinates,132,1,3,1e-10,64,2.,cap,0))return 21;
        reference_calls=entry_calls;sink_enabled=1;sink_count=0;sink_invalid=0;deferred_entries=0;entry_calls=0;
        if(HACApK_build_hmatrix_wrapper(lp,ctl,coordinates,132,1,3,1e-10,64,2.,cap,0))return 22;
        without_dense=entry_calls;
        if(sink_invalid||!sink_count||ref->nlf!=lp->nlf||ref->nlfkt!=lp->nlfkt||ref->ktmax!=lp->ktmax||
           (int)ctl->param[60]!=1||reference_calls-without_dense!=deferred_entries)return 23;
        for(i=1;i<=132;++i){if(cr->lod[i]!=ctl->lod[i])return 24;permuted+=ctl->lod[i]!=i;}
        for(i=1;i<=lp->nlf;++i) {
            st_cHACApK_leafmtx a=ref->st_lf[i],leaf=lp->st_lf[i];
            if(a->ltmtx!=leaf->ltmtx||a->kt!=leaf->kt||a->ndl!=leaf->ndl||a->ndt!=leaf->ndt||a->nstrtl!=leaf->nstrtl||a->nstrtt!=leaf->nstrtt)return 25;
            if(leaf->ltmtx==1) {
                if(memcmp(a->a1,leaf->a1,sizeof(double)*leaf->ndt*leaf->kt)||memcmp(a->a2,leaf->a2,sizeof(double)*leaf->ndl*leaf->kt))return 26;
            }else if(a->a1) {
                size_t k,total=(size_t)leaf->ndl*leaf->ndt;
                if(leaf->a1||leaf->a2)return 27;
                leaf->a1=calloc(total,sizeof(double));if(!leaf->a1)return 28;
                for(k=0;k<total;++k) {
                    int r=ctl->lod[leaf->nstrtl+(int)(k/(size_t)leaf->ndt)],c=ctl->lod[leaf->nstrtt+(int)(k%(size_t)leaf->ndt)];
                    if(cursor>=deferred_entries||captured_rows[cursor]!=r||captured_columns[cursor]!=c)return 29;
                    leaf->a1[k]=cHACApK_entry_ij(r,c,0);++cursor;
                }
                if(memcmp(a->a1,leaf->a1,total*sizeof(double)))return 30;
            }
        }
        if(cursor!=deferred_entries||!permuted)return 31;
        for(i=0;i<132;++i)x[i]=sin(.17*i)+.3;
        if(symmetric)HACApK_matvec_sym_wrapper(lp,ctl,x,y,132);else HACApK_matvec_wrapper(lp,ctl,x,y,132);
        for(i=0;i<132;++i){double expected=0.;for(j=0;j<132;++j)expected+=cHACApK_entry_ij(i+1,j+1,0)*x[j];diff+=(y[i]-expected)*(y[i]-expected);norm+=expected*expected;}
        diff=sqrt(diff/norm);if(!isfinite(diff)||diff>1e-8)return 32;if(diff>worst)worst=diff;
        total_sinks+=sink_count;removed+=deferred_entries;++case_count;release(ref,cr);release(lp,ctl);
        /* Destruction must also tolerate never-materialized deferred leaves. */
        lp=calloc(1,sizeof(*lp));ctl=calloc(1,sizeof(*ctl));if(!lp||!ctl)return 33;
        sink_count=0;deferred_entries=0;entry_calls=0;
        if(HACApK_build_hmatrix_wrapper(lp,ctl,coordinates,132,1,3,1e-10,64,2.,cap,0)||sink_invalid)return 34;
        release(lp,ctl);
    }
    sink_enabled=0;cHACApK_set_sym_fill(0);cHACApK_set_cluster_strategy(CHACAPK_CLUSTER_BBOX);
    printf("{\"cases\":%d,\"maximum_relative_difference\":%.17g,\"deferred_leaves\":%d,\"removed_host_entry_queries\":%llu,\"permuted_indices\":%d,\"gpu_connected\":false,\"taskmanager_started\":false}\n",case_count,worst,total_sinks,(unsigned long long)removed,permuted);
    return 0;
}
#endif
