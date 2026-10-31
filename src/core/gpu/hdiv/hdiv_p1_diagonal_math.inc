// Shared arithmetic instantiated with production CPU potentials in host controls
// and the same production device potentials in CUDA. No quadrature selection.
namespace radia_hdiv_p1_diagonal {
RADIA_P1_KERNEL inline double selfEntry(const radia_hdiv_python::P1ChargeSource& source,
        const radia_hdiv_python::PhysicalChargePoint* points){
    double total=0.;
    for(std::size_t q=0;q<source.count;++q){
        const auto& point=points[source.first+q];double potential,monomial=1.;
        if(!source.linear){
            potential=source.kind==0?
                RADIA_P1_FIELD::PhiTet(reinterpret_cast<const double(*)[3]>(source.vertices),point.point):
                RADIA_P1_FIELD::TriPotential(reinterpret_cast<const double(*)[3]>(source.vertices),point.point);
        }else{
            double vertices[3][3],local[3],moment[3];
            for(int k=0;k<3;++k){
                local[k]=point.point[k]-source.center[k];
                for(int v=0;v<3;++v)vertices[v][k]=source.vertices[3*v+k]-source.center[k];
            }
            potential=source.constant*RADIA_P1_FIELD::TriPotential(vertices,local);
            RADIA_P1_FIELD::TriMoment1(vertices,local,moment);
            monomial=source.constant;
            for(int k=0;k<3;++k){
                potential+=source.gradient[k]*moment[k];
                monomial+=source.gradient[k]*local[k];
            }
        }
        total+=point.weight*monomial*potential;
    }
    return total*0.07957747154594766788; // Same CPU QuadDot 1/(4*pi).
}
}
