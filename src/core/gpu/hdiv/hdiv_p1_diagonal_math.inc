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
template<int N> RADIA_P1_KERNEL inline double pairEntryOrder(
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
// All entries of one host pair at once, the same arithmetic as pairEntryOrder per
// entry (a /fp:fast build may round differently): the targets share the host's
// vertices, so the outer points are common,
// and every source potential is evaluated once per point instead of once per
// target mode.  The linear modes of one face share TriPotential/TriMoment1
// arguments (vertices and centroid), so those parts are reused too.  Each
// entry keeps pairEntryOrder's operation order.
constexpr int kP1BlockMaxModes=4;
template<int N> RADIA_P1_KERNEL inline bool pairBlockOrder(
        const radia_hdiv_python::P1ChargeSource* targets,int nt,
        const radia_hdiv_python::P1ChargeSource* sources,int ns,double* out){
    if(nt<1||ns<1||nt>kP1BlockMaxModes||ns>kP1BlockMaxModes)return false;
    const auto& target=targets[0];
    for(int t=1;t<nt;++t){
        if(targets[t].kind!=target.kind)return false;
        for(int k=0;k<12;++k)if(targets[t].vertices[k]!=target.vertices[k])return false;
    }
    double x[N],w[N],edges[3][3]={};gradedRule<N>(x,w);
    for(int v=0;v<(target.kind?2:3);++v)
        for(int k=0;k<3;++k)edges[v][k]=target.vertices[3*(v+1)+k]-target.vertices[k];
    double cross[3]={edges[0][1]*edges[1][2]-edges[0][2]*edges[1][1],
        edges[0][2]*edges[1][0]-edges[0][0]*edges[1][2],
        edges[0][0]*edges[1][1]-edges[0][1]*edges[1][0]};
    const double measure=target.kind?sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]):
        fabs(cross[0]*edges[2][0]+cross[1]*edges[2][1]+cross[2]*edges[2][2]);
    // share[s]: earlier linear source with the same vertices and centroid.
    int share[kP1BlockMaxModes];
    for(int s=0;s<ns;++s){
        share[s]=-1;
        if(!sources[s].linear)continue;
        for(int q=0;q<s&&share[s]<0;++q){
            if(!sources[q].linear||share[q]>=0||sources[q].kind!=sources[s].kind)continue;
            bool same=true;
            for(int k=0;k<12&&same;++k)same=sources[q].vertices[k]==sources[s].vertices[k];
            for(int k=0;k<3&&same;++k)same=sources[q].center[k]==sources[s].center[k];
            if(same)share[s]=q;
        }
    }
    double total[kP1BlockMaxModes][kP1BlockMaxModes]={};
    for(int i=0;i<N;++i)for(int j=0;j<N;++j)for(int l=0;l<(target.kind?1:N);++l){
        const double u=x[i],v=(1-u)*x[j],r=target.kind?0.:(1-u)*(1-x[j])*x[l];
        double point[3],mono[kP1BlockMaxModes];
        for(int t=0;t<nt;++t)mono[t]=targets[t].constant;
        for(int k=0;k<3;++k){
            point[k]=target.vertices[k]+u*edges[0][k]+v*edges[1][k]+r*edges[2][k];
            for(int t=0;t<nt;++t)
                if(targets[t].linear)mono[t]+=targets[t].gradient[k]*(point[k]-targets[t].center[k]);
        }
        const double weight=w[i]*w[j]*(1-u)*(target.kind?1.:w[l]*(1-u)*(1-x[j]));
        double potential[kP1BlockMaxModes],tri[kP1BlockMaxModes],moment[kP1BlockMaxModes][3];
        for(int s=0;s<ns;++s){
            const auto& src=sources[s];
            if(!src.linear){potential[s]=potentialAt(src,point);continue;}
            if(share[s]<0){
                double vertices[3][3],local[3];
                for(int k=0;k<3;++k){
                    local[k]=point[k]-src.center[k];
                    for(int c=0;c<3;++c)vertices[c][k]=src.vertices[3*c+k]-src.center[k];
                }
                tri[s]=RADIA_P1_FIELD::TriPotential(vertices,local);
                RADIA_P1_FIELD::TriMoment1(vertices,local,moment[s]);
            }else{
                tri[s]=tri[share[s]];
                for(int k=0;k<3;++k)moment[s][k]=moment[share[s]][k];
            }
            double value=src.constant*tri[s];
            for(int k=0;k<3;++k)value+=src.gradient[k]*moment[s][k];
            potential[s]=value;
        }
        for(int t=0;t<nt;++t)for(int s=0;s<ns;++s)total[t][s]+=weight*mono[t]*potential[s];
    }
    for(int t=0;t<nt;++t)for(int s=0;s<ns;++s)out[t*ns+s]=total[t][s]*measure*0.07957747154594766788;
    return true;
}
template<int NTri=24,int NTet=16> RADIA_P1_KERNEL inline bool pairBlock(
        const radia_hdiv_python::P1ChargeSource* targets,int nt,
        const radia_hdiv_python::P1ChargeSource* sources,int ns,double* out){
    return targets[0].kind?pairBlockOrder<NTri>(targets,nt,sources,ns,out):
        pairBlockOrder<NTet>(targets,nt,sources,ns,out);
}
// Order per TARGET kind: a TET target is a 3-D tensor rule (N^3 potentials),
// a TRI target 2-D (N^2). Measured host self entries vs GL48 on the oblique
// tet: TET N=8 1.3e-6 at 1.4 ms, N=16 4e-11 at 15 ms; TRI N=12 <=1.7e-7.
// The rule is anchored at local vertex 0, so its error also differs between
// renumbered or reflected copies of one pair; N=16/24 keeps that mirror
// discrepancy small while the rule stays smooth in the vertex positions.
// Observed diagnostics, not an adopted tolerance.
template<int NTri=24,int NTet=16> RADIA_P1_KERNEL inline double pairEntry(
        const radia_hdiv_python::P1ChargeSource& target,
        const radia_hdiv_python::P1ChargeSource& source){
    return target.kind?pairEntryOrder<NTri>(target,source):pairEntryOrder<NTet>(target,source);
}
}
