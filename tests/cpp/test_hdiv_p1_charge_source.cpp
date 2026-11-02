// Guards: the affine P1 charge frame shared by the CPU ChargeGram and the device
// flat Gram (gradient, oblique frame, reference permutation, source moment, every
// field written when built in place over garbage) and its loud rejection of
// singular or nonfinite faces.
// Assumes: shared production potentials; algebraic controls, not independent accuracy validation.
#include "hdiv_p1_diagonal_input.hpp"
#include "rad_hdiv_vim.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#define RADIA_P1_KERNEL
#define RADIA_P1_FIELD rad_hdiv
#include "hdiv_p1_diagonal_math.inc"
#undef RADIA_P1_FIELD
#undef RADIA_P1_KERNEL
using namespace radia_hdiv_python;

static P1Charge unitFace(int axis){
    P1Charge charge;charge.kind=1;charge.vertices={0,0,0,1,0,0,0,1,0,0,0,0};charge.exponent[axis]=1;
    return charge;
}
template<class Mutation>static void rejects(const char* name,Mutation mutate){
    auto charge=unitFace(0);mutate(charge);
    try{p1ChargeSource(charge);}
    catch(const std::invalid_argument&){return;}
    throw std::runtime_error(std::string(name)+": actual=accepted; expected=invalid_argument");
}
int main(){
    try{
        for(int axis=0;axis<2;++axis){
            auto charge=unitFace(axis);
            const auto source=p1ChargeSource(charge);
            if(!source.linear||source.gradient[axis]!=1.||source.gradient[1-axis]!=0.||source.gradient[2]!=0.||
               source.center[0]!=1./3.||source.center[1]!=1./3.||source.constant!=1./3.)
                throw std::runtime_error("actual P1 gradient differs from expected unit-triangle barycentric gradient");
            // Source potential at the centre is A*TriPotential + gradient . TriMoment1 in the centred frame.
            double localVertices[3][3],local[3]={0,0,0},moment[3];
            for(int i=0;i<3;++i)for(int k=0;k<3;++k)localVertices[i][k]=source.vertices[3*i+k]-source.center[k];
            rad_hdiv::TriMoment1(localVertices,local,moment);
            if(moment[axis]==0.)throw std::runtime_error("source moment control is unobservable");
            const double expected=source.constant*rad_hdiv::TriPotential(localVertices,local)+moment[axis];
            const double actual=radia_hdiv_p1_diagonal::potentialAt(source,source.center);
            if(actual!=expected)
                throw std::runtime_error("source moment actual="+std::to_string(actual)+" expected="+std::to_string(expected));
            // Oblique frame: x=xi+eta, y=eta, hence xi=x-y and eta=y.
            auto oblique=charge;oblique.vertices[6]=1.;
            const auto skew=p1ChargeSource(oblique);
            if(skew.gradient[0]!=(axis?0.:1.)||skew.gradient[1]!=(axis?1.:-1.))
                throw std::runtime_error("actual oblique P1 gradient differs from independent affine-coordinate identity");
            auto permutedCharge=charge;
            std::swap(permutedCharge.vertices[3],permutedCharge.vertices[6]);
            std::swap(permutedCharge.vertices[4],permutedCharge.vertices[7]);
            const auto permuted=p1ChargeSource(permutedCharge);
            if(permuted.gradient[1-axis]!=1.||permuted.gradient[axis]!=0.)
                throw std::runtime_error("permuted reference gradient actual differs from expected swapped basis");
        }
        // The device builds sources in place in uninitialized memory: every field must be written.
        const double nan=std::numeric_limits<double>::quiet_NaN();
        auto dirty=[&]{P1ChargeSource s;for(auto& v:s.vertices)v=nan;for(int k=0;k<3;++k)s.center[k]=s.gradient[k]=nan;
                       s.constant=nan;s.kind=7;s.linear=true;return s;};
        for(int kind=0;kind<2;++kind){
            P1Charge constant;constant.kind=kind;constant.vertices={0,0,0,1,0,0,0,1,0,0,0,1};
            auto inPlace=dirty();
            if(!p1ChargeSourceInto(kind,constant.vertices.data(),constant.exponent.data(),inPlace))
                throw std::runtime_error("constant charge source rejected");
            if(inPlace.linear||inPlace.kind!=kind||inPlace.constant!=1.||inPlace.center[0]!=0.||inPlace.center[1]!=0.||
               inPlace.center[2]!=0.||inPlace.gradient[0]!=0.||inPlace.gradient[1]!=0.||inPlace.gradient[2]!=0.)
                throw std::runtime_error("in-place constant source kind="+std::to_string(kind)+" actual constant="+
                    std::to_string(inPlace.constant)+" expected=1 with zero centre and gradient");
        }
        {
            const auto face=unitFace(1);auto inPlace=dirty();const auto fresh=p1ChargeSource(face);
            if(!p1ChargeSourceInto(1,face.vertices.data(),face.exponent.data(),inPlace)||inPlace.constant!=fresh.constant||
               inPlace.gradient[1]!=fresh.gradient[1]||inPlace.center[2]!=fresh.center[2]||inPlace.vertices[11]!=fresh.vertices[11])
                throw std::runtime_error("in-place linear face source differs from the fresh source");
        }
        rejects("collinear face",[](auto& c){c.vertices[7]=0.;});
        rejects("coincident face vertices",[](auto& c){c.vertices[6]=1.;c.vertices[7]=0.;});
        rejects("nonfinite face",[](auto& c){c.vertices[0]=std::numeric_limits<double>::infinity();});
        rejects("NaN face",[](auto& c){c.vertices[4]=std::numeric_limits<double>::quiet_NaN();});
        std::cout<<"both linear TRI frames, source moments, oblique and permuted frames, in-place constant and linear sources, 4 rejections passed\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
