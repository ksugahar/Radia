// Guards: the affine P1 self/touching outer rule (pairEntry) against regression to
// a fixed low-order outer rule, which was 0.5 % high on constant TRI self entries.
// Assumes: references from methods that import no Radia potential: TRI by the
// closed-form in-plane triangle potential with a GL120 Duffy outer rule, TET by
// the star-simplex reference (validation_test/hmatrix_cuda/p1_self_integral_reference.py,
// outer 24 / inner 64). DETECTOR_REL is a proposed regression detector pending user
// approval, not an adopted accuracy criterion.
#include "hdiv_p1_diagonal_input.hpp"
#include "rad_hdiv_vim.h"
#include <cmath>
#include <iostream>
#define RADIA_P1_KERNEL
#define RADIA_P1_FIELD rad_hdiv
#include "hdiv_p1_diagonal_math.inc"
#undef RADIA_P1_FIELD
#undef RADIA_P1_KERNEL

static constexpr double DETECTOR_REL = 1e-5;

static int check(const char* name, int kind, std::array<double, 12> vertices, double reference){
    radia_hdiv_python::P1Charge charge;
    charge.kind=kind;charge.vertices=vertices;charge.exponent={0,0,0};charge.count=1;
    const auto source=radia_hdiv_python::p1ChargeSource(charge);
    // Unit-density integral / (4 pi); pairEntry carries the normalized charge constant twice.
    const double actual=radia_hdiv_p1_diagonal::pairEntry(source,source)/(source.constant*source.constant);
    const double relative=std::fabs(actual-reference)/reference;
    if(!(relative<=DETECTOR_REL)){
        std::cerr<<name<<" actual="<<actual<<" expected="<<reference<<" relative="<<relative<<'\n';
        return 1;
    }
    return 0;
}

int main(){
    int failures=0;
    failures+=check("tri unit right self",1,{0,0,0,1,0,0,0,1,0,0,0,0},0.0798214472);
    failures+=check("tet unit corner self",0,{0,0,0,1,0,0,0,1,0,0,0,1},0.006846366);
    return failures;
}
