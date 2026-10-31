// Guards: source linear moment, target monomial and local-frame construction.
// Assumes: shared production potentials; algebraic controls, not independent accuracy validation.
#include "hdiv_p1_diagonal_input.hpp"
#include "rad_hdiv_vim.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#define RADIA_P1_KERNEL
#define RADIA_P1_FIELD rad_hdiv
#include "hdiv_p1_diagonal_math.inc"
#undef RADIA_P1_FIELD
#undef RADIA_P1_KERNEL

int main(){
    try{
        for(int axis=0;axis<2;++axis){
            radia_hdiv_python::P1Charge charge;
            charge.kind=1;charge.vertices={0,0,0,1,0,0,0,1,0,0,0,0};
            charge.exponent[axis]=1;charge.count=1;
            const auto source=radia_hdiv_python::p1ChargeSource(charge);
            radia_hdiv_python::PhysicalChargePoint point{{source.center[0],source.center[1],0.},1.};
            double localVertices[3][3],local[3]={0,0,0},moment[3];
            for(int i=0;i<3;++i)for(int k=0;k<3;++k)
                localVertices[i][k]=source.vertices[3*i+k]-source.center[k];
            const auto constant=source.constant*rad_hdiv::TriPotential(localVertices,local);
            rad_hdiv::TriMoment1(localVertices,local,moment);
            if(moment[axis]==0.)throw std::runtime_error("source moment control is unobservable");
            const double expected=source.constant*(constant+moment[axis])*0.07957747154594766788;
            const double actual=radia_hdiv_p1_diagonal::selfEntry(source,&point);
            if(actual!=expected){std::cerr<<"source moment actual="<<actual<<" expected="<<expected<<'\n';return 1;}
            // This reference vertex has monomial zero and a nonzero source potential.
            point.point[0]=point.point[1]=point.point[2]=0.;
            const double atVertex=radia_hdiv_p1_diagonal::selfEntry(source,&point);
            if(atVertex!=0.){std::cerr<<"target monomial actual="<<atVertex<<" expected=0\n";return 1;}
            std::swap(charge.vertices[3],charge.vertices[6]);
            std::swap(charge.vertices[4],charge.vertices[7]);
            const auto permuted=radia_hdiv_python::p1ChargeSource(charge);
            if(permuted.gradient[1-axis]!=1.||permuted.gradient[axis]!=0.)
                throw std::runtime_error("permuted reference gradient actual differs from expected swapped basis");
        }
        std::cout<<"both source moments, target zeros and reference permutations passed\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
