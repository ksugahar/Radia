/* Full existing HACApK cluster/leaf/fill integration with serial host control.
 * A synthetic smooth kernel exercises infrastructure, not HDiv physics. */
#include "cHACApK_base.h"
#include "cHACApK_pca_cluster.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int HACApK_build_hmatrix_wrapper(void*,void*,double*,int,int,int,double,int,double,int,int);
void HACApK_matvec_wrapper(void*,void*,const double*,double*,int);
void HACApK_free_hmatrix_wrapper(void*,void*);
void HACApK_free_leafmtxp(void*);
static double coordinates[132*3];
double cHACApK_entry_ij(int i,int j,int kind) {
    double d=1.; int k; (void)kind;
    for(k=0;k<3;++k) { double v=coordinates[3*(i-1)+k]-coordinates[3*(j-1)+k]; d+=v*v; }
    return 1./sqrt(d)+(i==j ? 2. : 0.);
}
void cHACApK_entry_batch(int n,const int*r,const int*c,int kind,double*v) {
    int k; for(k=0;k<n;++k) v[k]=cHACApK_entry_ij(r[k],c[k],kind);
}
int main(void) {
    int strategy,i,j,lowrank=0; double worst=0.;
    for(i=0;i<132;++i) { coordinates[3*i]=(i<66 ? 0. : 20.)+.01*(i%66); coordinates[3*i+1]=.03*(i%7); coordinates[3*i+2]=.02*(i%11); }
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
    printf("{\"cases\":2,\"maximum_relative_difference\":%.17g,\"lowrank_leaves\":%d,\"taskmanager_started\":false}\n",worst,lowrank);
    return 0;
}
