/* Independent scalar C execution of the shared ACA+ on a 4x4 toy Gram input.
 * Diagnose compression separately from GPU evaluation; not a GPU fallback. */
#include "cHACApK_base.h"
#include <stdio.h>
#include <math.h>
static double gram[16],scale[4];
double cHACApK_entry_ij(int row,int column,int problem){
    (void)problem;return gram[(row-1)*4+column-1]/scale[row-1]/scale[column-1];
}
int main(int argc,char** argv){
    FILE* input;int i,j,k,rank,lod[5]={0,1,2,3,4};
    double a[4]={0},b[4]={0},parameters[100]={0},error=0,norm=0;
    if(argc!=2)return 1;input=fopen(argv[1],"r");if(!input)return 2;
    for(i=0;i<16;++i)if(fscanf(input,"%lf",&gram[i])!=1){fclose(input);return 3;}
    fclose(input);for(i=0;i<4;++i)scale[i]=sqrt(gram[i*4+i]);
    parameters[61]=1;parameters[62]=2;parameters[64]=1;
    rank=cHACApK_acaplus(a,b,parameters,2,2,1,3,lod,0,2,1e-10,1.,1e-13);
    for(i=0;i<2;++i)for(j=0;j<2;++j){
        double actual=0,expected=gram[i*4+j+2];
        for(k=0;k<rank;++k)actual+=a[i+2*k]*b[j+2*k]*scale[i]*scale[j+2];
        error+=(actual-expected)*(actual-expected);norm+=expected*expected;
    }
    printf("{\"rank\":%d,\"relative_error\":%.17g,\"route\":\"independent scalar C, no CUDA\"}\n",rank,sqrt(error/norm));
    return 0;
}
