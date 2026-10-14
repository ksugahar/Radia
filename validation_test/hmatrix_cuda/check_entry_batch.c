/* Contract of the unchanged ACA active row/column extraction. Not GPU parity. */
#include "cHACApK_calc_entry_ij.h"
#include <stdio.h>
static int calls,requested;
double cHACApK_entry_ij(int row,int column,int problem){
    return row*100000.0+column*10.0+problem;
}
#ifdef RAD_HACAPK_BATCH_ENTRIES
void cHACApK_entry_batch(int count,const int* rows,const int* columns,int problem,double* values){
    int i;++calls;requested+=count;
    if(count<1||count>256){requested=-1000000;return;}
    for(i=0;i<count;++i)values[i]=cHACApK_entry_ij(rows[i],columns[i],problem);
}
#endif
int main(void){
    int sizes[]={0,1,255,256,257,1024},lod[1100],mask[1024];double vec[1024];
    int s,d,mode,i,checks=0;
    for(i=0;i<1100;++i)lod[i]=1100-i;
    for(s=0;s<6;++s)for(d=0;d<2;++d)for(mode=0;mode<3;++mode){
        int active=0;calls=0;requested=0;
        for(i=0;i<sizes[s];++i){mask[i]=mode==0?0:mode==1?(i%3==0):1;vec[i]=-7.;if(!mask[i])++active;}
        cHACApK_active_entry_vector(sizes[s],3,7,11,lod,2,mask,d,vec);
        for(i=0;i<sizes[s];++i){
            const double expected=mask[i]?-7.:cHACApK_entry_ij(lod[(d==0?3:i)+7],lod[(d==0?i:3)+11],2);
            if(vec[i]!=expected)return 1;
        }
#ifdef RAD_HACAPK_BATCH_ENTRIES
        if(requested!=active||calls!=(active+255)/256)return 2;
#endif
        ++checks;
    }
    printf("{\"accepted\":true,\"cases\":%d}\n",checks);return 0;
}
