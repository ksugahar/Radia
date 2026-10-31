// Shared arithmetic instantiated with production CPU potentials in host controls
// and the same production device potentials in CUDA. No quadrature selection.
namespace radia_hdiv_p1_diagonal {
// Bounded endpoint-graded outer rule for affine P1 self/touching pairs.
// The analytic inner potential is unchanged. No CPU/GPU bitwise claim.
RADIA_P1_KERNEL inline double potentialAt(const radia_hdiv_python::P1ChargeSource& s,
        const double* p){
    if(!s.linear)return s.kind==0?
        RADIA_P1_FIELD::PhiTet(reinterpret_cast<const double(*)[3]>(s.vertices),p):
        RADIA_P1_FIELD::TriPotential(reinterpret_cast<const double(*)[3]>(s.vertices),p);
    double vertices[3][3],local[3],moment[3];
    for(int k=0;k<3;++k){
        local[k]=p[k]-s.center[k];
        for(int v=0;v<3;++v)vertices[v][k]=s.vertices[3*v+k]-s.center[k];
    }
    double value=s.constant*RADIA_P1_FIELD::TriPotential(vertices,local);
    RADIA_P1_FIELD::TriMoment1(vertices,local,moment);
    for(int k=0;k<3;++k)value+=s.gradient[k]*moment[k];
    return value;
}
template<int N> RADIA_P1_KERNEL inline void gradedRule(double (&x)[N],double (&w)[N]){
    for(int i=0;i<(N+1)/2;++i){
        double z=cos(3.14159265358979323846*(i+.75)/(N+.5)),derivative=0.;
        for(int iteration=0;iteration<16;++iteration){
            double p=1.,previous=0.;
            for(int k=1;k<=N;++k){double old=p;p=((2*k-1)*z*p-(k-1)*previous)/k;previous=old;}
            derivative=N*(z*p-previous)/(z*z-1.);
            const double step=p/derivative;z-=step;
            if(fabs(step)<2.e-16)break;
        }
        // Re-evaluate derivative at the final root (also after iteration cap).
        double p=1.,previous=0.;
        for(int k=1;k<=N;++k){double old=p;p=((2*k-1)*z*p-(k-1)*previous)/k;previous=old;}
        derivative=N*(z*p-previous)/(z*z-1.);
        const double weight=1./((1-z*z)*derivative*derivative);
        const double a=(1-z)/2,b=(1+z)/2;
        x[i]=a*a*(3-2*a);x[N-1-i]=b*b*(3-2*b);
        w[i]=weight*6*a*(1-a);w[N-1-i]=weight*6*b*(1-b);
    }
}
template<int N=16> RADIA_P1_KERNEL inline double pairEntry(
        const radia_hdiv_python::P1ChargeSource& target,
        const radia_hdiv_python::P1ChargeSource& source){
    double x[N],w[N],edges[3][3]={};gradedRule<N>(x,w);
    for(int v=0;v<(target.kind?2:3);++v)
        for(int k=0;k<3;++k)edges[v][k]=target.vertices[3*(v+1)+k]-target.vertices[k];
    double cross[3]={edges[0][1]*edges[1][2]-edges[0][2]*edges[1][1],
        edges[0][2]*edges[1][0]-edges[0][0]*edges[1][2],
        edges[0][0]*edges[1][1]-edges[0][1]*edges[1][0]};
    const double measure=target.kind?sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]):
        fabs(cross[0]*edges[2][0]+cross[1]*edges[2][1]+cross[2]*edges[2][2]);
    double total=0.;
    for(int i=0;i<N;++i)for(int j=0;j<N;++j)for(int l=0;l<(target.kind?1:N);++l){
        const double u=x[i],v=(1-u)*x[j],r=target.kind?0.:(1-u)*(1-x[j])*x[l];
        double point[3],mono=target.constant;
        for(int k=0;k<3;++k){
            point[k]=target.vertices[k]+u*edges[0][k]+v*edges[1][k]+r*edges[2][k];
            if(target.linear)mono+=target.gradient[k]*(point[k]-target.center[k]);
        }
        const double weight=w[i]*w[j]*(1-u)*(target.kind?1.:w[l]*(1-u)*(1-x[j]));
        total+=weight*mono*potentialAt(source,point);
    }
    return total*measure*0.07957747154594766788;
}
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
