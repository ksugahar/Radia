/* Scalar/batched C differential contract of shared ACA+ cyclic search. */
#include "cHACApK_base.h"
#include <math.h>
#include <stdio.h>
static int rows,columns,kind;
double cHACApK_entry_ij(int row,int column,int problem){
    int i=row-1,j=column-rows-1;(void)problem;
    return kind?(i+1.)*(j+2.):1./(1.+fabs((double)i-j));
}
#ifdef RAD_HACAPK_BATCH_ENTRIES
void cHACApK_entry_batch(int count,const int* r,const int* c,int problem,double* values){
    int i;for(i=0;i<count;++i)values[i]=cHACApK_entry_ij(r[i],c[i],problem);
}
#endif
int main(void){
    int sizes[][2]={{1,1},{2,2},{3,5},{5,3},{7,7},{13,9}},test,i,j,k,rank,limit,checks=0;
    double worst=0.;
    for(test=0;test<6;++test)for(kind=0;kind<2;++kind){
        int lod[24];double a[256]={0},b[256]={0},parameters[100]={0},error=0.,norm=0.;
        rows=sizes[test][0];columns=sizes[test][1];limit=rows<columns?rows:columns;
        for(i=0;i<=rows+columns;++i)lod[i]=i;
        parameters[61]=1;parameters[62]=limit;parameters[64]=1;
        rank=cHACApK_acaplus(a,b,parameters,rows,columns,1,rows+1,lod,0,limit,1e-10,1.,1e-13);
        if(rank<1||rank>limit)return 1;
        for(i=0;i<rows;++i)for(j=0;j<columns;++j){
            double actual=0.,expected=cHACApK_entry_ij(i+1,rows+j+1,0);
            for(k=0;k<rank;++k)actual+=a[i+rows*k]*b[j+columns*k];
            error+=(actual-expected)*(actual-expected);norm+=expected*expected;
        }
        error=sqrt(error/norm);if(error>worst)worst=error;
        if(error>1e-10)return 2;++checks;
    }
    printf("{\"accepted\":true,\"cases\":%d,\"max_relative_error\":%.17g}\n",checks,worst);return 0;
}
