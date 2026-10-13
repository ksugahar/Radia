#include "mex.h"
extern "C" void RadiaProjectionGateway(int,mxArray**,int,const mxArray**);
void mexFunction(int nl,mxArray* out[],int nr,const mxArray* in[]){RadiaProjectionGateway(nl,out,nr,in);}
