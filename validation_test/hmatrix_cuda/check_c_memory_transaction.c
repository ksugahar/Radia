/* Host-only transaction contract. Link with hmatrix_c_memory.c instead of
 * actual HACApK assembly. The MEX numerical/failure gate remains separate. */
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "cHACApK_base.h"
#define RAD_HACAPK_C_MEMORY
#include "hmatrix_c_memory.h"
static int forced_exit;
int HACApK_build_hmatrix_wrapper(void *matrix,void *control,double *xyz,
    int n,int nffc,int dim,double eps,int leaf,double eta,int rank,int print){
    char *p=(char*)malloc(17);
    char *q=(char*)calloc(2,19);
    (void)matrix;(void)control;(void)xyz;(void)n;(void)nffc;(void)dim;
    (void)eps;(void)leaf;(void)eta;(void)rank;(void)print;
    assert(p && q);memset(p,7,17);
    p=(char*)realloc(p,23);assert(p && p[0]==7 && p[16]==7);
    if(forced_exit)exit(EXIT_FAILURE);
    free(p);free(q);return 0;
}
int main(void){
    st_cHACApK_leafmtxp_t matrix;
    st_cHACApK_lcontrol_t control;
    size_t fail_at;int repeat;
    for(fail_at=0;fail_at<=7;++fail_at){
        for(repeat=0;repeat<10;++repeat){
            rad_hmatrix_c_memory_report report;int status;
            memset(&matrix,0,sizeof(matrix));memset(&control,0,sizeof(control));
            rad_hmatrix_c_memory_configure(fail_at);
            status=rad_hmatrix_guarded_build(&matrix,&control,NULL,1,1,3,1e-8,1,2.,1,0);
            report=rad_hmatrix_c_memory_last();
            assert((fail_at==0 && status==0) || (fail_at==1 && status==-92) || (fail_at>1 && status==-93));
            assert(report.live_bytes==0);
            if(fail_at){assert(report.failure_code==1);assert(report.allocations==fail_at);}
        }
    }
    rad_hmatrix_c_memory_configure(0);forced_exit=1;
    assert(rad_hmatrix_guarded_build(&matrix,&control,NULL,1,1,3,1e-8,1,2.,1,0)==-93);
    assert(rad_hmatrix_c_memory_last().failure_code==4);
    assert(rad_hmatrix_c_memory_last().live_bytes==0);
    forced_exit=0;
    assert(rad_hmatrix_guarded_build(&matrix,&control,NULL,1,1,3,1e-8,1,2.,1,0)==0);
    return 0;
}
