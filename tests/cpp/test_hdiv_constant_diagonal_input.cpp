// Guards: host admission before CUDA; no GPU/mathematical acceptance claim.
// Assumes: caller owns mapped physical quadrature and charge-space semantics.
#include "hdiv_constant_diagonal_input.hpp"
#include <iostream>
#include <limits>
#include <string>
using namespace radia_hdiv_python;

ConstantDiagonalInput fixture(){
    ConstantDiagonalInput input;
    ConstantCharge tet;tet.vertices={0,0,0,1,0,0,0,1,0,0,0,1};tet.first=0;tet.count=1;
    ConstantCharge tri;tri.vertices={0,0,0,1,0,0,0,1,0,0,0,0};tri.kind=1;tri.first=1;tri.count=1;
    input.charges={tet,tri};
    input.quadrature={{{.25,.25,.25},1./6.},{{1./3.,1./3.,0},.5}};
    return input;
}
template<class Mutation>void rejects(const char* name,Mutation mutate){
    auto input=fixture();mutate(input);
    try{input.validate();}
    catch(const std::invalid_argument&){return;}
    throw std::runtime_error(std::string(name)+": actual=accepted; expected=invalid_argument");
}
int main(){
    try{
        fixture().validate();
        rejects("empty charge set",[](auto& x){x.charges.clear();x.quadrature.clear();});
        rejects("charge budget",[](auto& x){
            auto charge=x.charges[0];auto point=x.quadrature[0];
            x.charges.clear();x.quadrature.clear();
            for(std::size_t i=0;i<4609;++i){charge.first=i;x.charges.push_back(charge);x.quadrature.push_back(point);}
        });
        rejects("unsupported support",[](auto& x){x.charges[0].kind=2;});
        rejects("negative support",[](auto& x){x.charges[0].kind=-1;});
        rejects("nonconstant monomial",[](auto& x){x.charges[0].exponent[0]=1;});
        rejects("negative exponent",[](auto& x){x.charges[0].exponent[2]=-1;});
        rejects("partition gap",[](auto& x){x.charges[1].first=0;});
        rejects("empty rule",[](auto& x){x.charges.resize(1);x.charges[0].count=0;x.quadrature.clear();});
        rejects("per-charge point budget",[](auto& x){
            x.quadrature.resize(258);x.charges[0].count=257;x.charges[1].first=257;
        });
        rejects("partition bounds",[](auto& x){x.charges[0].count=3;});
        rejects("unused quadrature",[](auto& x){x.quadrature.push_back(x.quadrature[0]);});
        rejects("nonfinite geometry",[](auto& x){x.charges[0].vertices[0]=std::numeric_limits<double>::infinity();});
        rejects("triangle padding",[](auto& x){x.charges[1].vertices[9]=1;});
        rejects("degenerate tetrahedron",[](auto& x){x.charges[0].vertices[11]=0;});
        rejects("degenerate triangle",[](auto& x){x.charges[1].vertices[7]=0;});
        rejects("nonfinite point",[](auto& x){x.quadrature[0].point[1]=std::numeric_limits<double>::quiet_NaN();});
        rejects("nonfinite weight",[](auto& x){x.quadrature[1].weight=std::numeric_limits<double>::quiet_NaN();});
        // A general quadrature rule may have signed weights; computed diagonal
        // positivity is checked on device, not guessed by the host contract.
        auto signedRule=fixture();signedRule.quadrature[0].weight=-.1;signedRule.validate();
        std::cout<<"host admission: 17 negative controls and valid TET/TRI passed\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
