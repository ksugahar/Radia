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
// Order-independent placement of the target vertices for the graded rule:
// the rule collapses onto slot 1 and is not symmetric in its slots, so the
// slots follow isometry-invariant keys of the PAIR geometry (squared
// distance to the nearest, then sum of squared distances to all, source
// vertices), not the mesh vertex numbering.  A touching vertex (key 0)
// becomes the collapse apex; reflected or renumbered copies of the same
// pair get the same rule.  Exact key ties keep the input order.
// perm[slot] = input vertex index.
RADIA_P1_KERNEL inline void gradedSlots(const double* t,int nt,const double* s,int ns,int (&perm)[4]){
    double near[4],sum[4];int order[4];
    for(int i=0;i<nt;++i){
        near[i]=0.;sum[i]=0.;order[i]=i;
        for(int j=0;j<ns;++j){
            double d=0.;
            for(int k=0;k<3;++k){const double x=t[3*i+k]-s[3*j+k];d+=x*x;}
            if(j==0||d<near[i])near[i]=d;
            sum[i]+=d;
        }
    }
    for(int i=1;i<nt;++i)for(int j=i;j>0;--j){
        const int a=order[j-1],b=order[j];
        if(near[b]<near[a]||(near[b]==near[a]&&sum[b]<sum[a])){order[j-1]=b;order[j]=a;}
        else break;
    }
    perm[1]=order[0];perm[0]=order[1];
    for(int i=2;i<nt;++i)perm[i]=order[i];
    for(int i=nt;i<4;++i)perm[i]=i;
}
template<int N> RADIA_P1_KERNEL inline double pairEntryOrder(
        const radia_hdiv_python::P1ChargeSource& target,
        const radia_hdiv_python::P1ChargeSource& source){
    double x[N],w[N],edges[3][3]={},slot[4][3];gradedRule<N>(x,w);
    int perm[4];gradedSlots(target.vertices,target.kind?3:4,source.vertices,source.kind?3:4,perm);
    for(int v=0;v<(target.kind?3:4);++v)
        for(int k=0;k<3;++k)slot[v][k]=target.vertices[3*perm[v]+k];
    for(int v=0;v<(target.kind?2:3);++v)
        for(int k=0;k<3;++k)edges[v][k]=slot[v+1][k]-slot[0][k];
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
            point[k]=slot[0][k]+u*edges[0][k]+v*edges[1][k]+r*edges[2][k];
            if(target.linear)mono+=target.gradient[k]*(point[k]-target.center[k]);
        }
        const double weight=w[i]*w[j]*(1-u)*(target.kind?1.:w[l]*(1-u)*(1-x[j]));
        total+=weight*mono*potentialAt(source,point);
    }
    return total*measure*0.07957747154594766788;
}
// Order per TARGET kind: a TET target is a 3-D tensor rule (N^3 potentials),
// a TRI target 2-D (N^2). Measured host self entries vs GL48 on the oblique
// tet: TET N=8 1.3e-6 at 1.4 ms (N=16 4e-11 at 15 ms); TRI N=12 <=1.7e-7.
// Observed diagnostics, not an adopted tolerance.
template<int NTri=12,int NTet=8> RADIA_P1_KERNEL inline double pairEntry(
        const radia_hdiv_python::P1ChargeSource& target,
        const radia_hdiv_python::P1ChargeSource& source){
    return target.kind?pairEntryOrder<NTri>(target,source):pairEntryOrder<NTet>(target,source);
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
