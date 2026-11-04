// Shared arithmetic of tetrahedral, pyramidal and triangular charge hosts, instantiated with the
// production CPU potentials in host controls and the same device potentials in CUDA.
// A host is stored as the Q2 lattice of its map from the unit cube, which is trilinear for every
// type: the cube corners collapse onto the tet (V3 nodes 0, V0 2, V1 6, V2 18), the pyramid
// (V0 0, V1 2, V2 8, V3 6, V4 18) and the triangle (V2 0, V0 2, V1 6). A charge density is a monomial
// of the NGSolve reference coordinates in the reference measure, whose Jacobian cancels the Piola one.
// The order 2 pyramid adds the rational div image x y / (1 - z), coded by a z exponent of -1.
namespace radia_hdiv_typed {
enum HostType {hexHost=0,quadHost=1,tetHost=2,triHost=3,pyramidHost=4};
// Polynomial coefficients of a typed source: a tet or pyramid charge has degree <= 1, a triangle charge <= 2.
constexpr int cellPolynomials=4,facePolynomials=10,polynomialSlots=10;
// Reference point y at cube point c; returns the Jacobian of the collapse.
RADIA_TYPED_KERNEL inline double collapse(int type,const double* c,double* y){
    if(type==tetHost){const double a=1.0-c[1],b=1.0-c[2];y[0]=c[0]*a*b;y[1]=c[1]*b;y[2]=c[2];return a*b*b;}
    if(type==pyramidHost){const double b=1.0-c[2];y[0]=c[0]*b;y[1]=c[1]*b;y[2]=c[2];return b*b;}
    if(type==triHost){const double a=1.0-c[1];y[0]=c[0]*a;y[1]=c[1];y[2]=0.0;return a;}
    y[0]=c[0];y[1]=c[1];y[2]=c[2];return 1.0;
}
RADIA_TYPED_KERNEL inline double monomial(const int* e,const double* y){
    double v=1.0;
    for(int k=0;k<2;++k)for(int i=0;i<e[k];++i)v*=y[k];
    if(e[2]<0){const double top=1.0-y[2];return top>1e-300?v/top:0.0;}
    for(int i=0;i<e[2];++i)v*=y[2];
    return v;
}
// Density of one charge per unit cube measure.
RADIA_TYPED_KERNEL inline double density(int type,const int* e,const double* c){
    double y[3];const double det=collapse(type,c,y);return det*monomial(e,y);
}
RADIA_TYPED_KERNEL inline bool isCell(int type){return type==hexHost||type==tetHost||type==pyramidHost;}
// A pair of two hosts integrates one direction: the cheapest inner source, otherwise the lower host (cells before
// faces) as the target, so that both orders of a pair take the same block.
RADIA_TYPED_KERNEL inline int sourceCost(int type){return type==triHost?0:type==tetHost?1:type==pyramidHost?2:type==quadHost?3:4;}
RADIA_TYPED_KERNEL inline bool lowerHostIsSource(int typeLower,int typeUpper){return sourceCost(typeLower)<sourceCost(typeUpper);}
RADIA_TYPED_KERNEL inline int polyIndex(int ax,int ay,int az){
    const int deg=ax+ay+az;int idx=0;
    for(int k=0;k<deg;++k)idx+=(k+1)*(k+2)/2;
    for(int a=0;a<ax;++a)idx+=deg-a+1;
    return idx+ay;
}
RADIA_TYPED_KERNEL inline void polyTimesLinear(double* poly,int& deg,const double lin[4],int ncoeff){
    double next[polynomialSlots];
    for(int i=0;i<polynomialSlots;++i)next[i]=0.0;
    for(int total=0;total<=deg;++total)for(int ax=0;ax<=total;++ax)for(int ay=0;ay<=total-ax;++ay){
        const int az=total-ax-ay;const double c=poly[polyIndex(ax,ay,az)];if(c==0.0)continue;
        next[polyIndex(ax,ay,az)]+=c*lin[0];next[polyIndex(ax+1,ay,az)]+=c*lin[1];
        next[polyIndex(ax,ay+1,az)]+=c*lin[2];next[polyIndex(ax,ay,az+1)]+=c*lin[3];
    }
    ++deg;for(int i=0;i<ncoeff;++i)poly[i]=next[i];
}
RADIA_TYPED_KERNEL inline double norm3(const double* a){return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);}
// Whether a tet, triangle or pyramid map is affine in the reference coordinates; if so the inverse
// linear forms about the origin, a reference point inside the host and the reference-measure scale.
RADIA_TYPED_KERNEL inline bool affineForms(int type,const double* nd,double lin[3][4],double& scaleInv){
    const double* o=nd;const double* px=nd+6;const double* py=nd+18;
    if(type==triHost){
        double a[3],b[3];for(int k=0;k<3;++k){a[k]=px[k]-o[k];b[k]=py[k]-o[k];}
        const double aa=a[0]*a[0]+a[1]*a[1]+a[2]*a[2],ab=a[0]*b[0]+a[1]*b[1]+a[2]*b[2],bb=b[0]*b[0]+b[1]*b[1]+b[2]*b[2],det=aa*bb-ab*ab;
        if(!(det>1e-300))return false;
        const double i00=bb/det,i01=-ab/det,i11=aa/det;
        for(int k=0;k<3;++k){lin[0][k+1]=i00*a[k]+i01*b[k];lin[1][k+1]=i01*a[k]+i11*b[k];lin[2][k+1]=0.0;}
        for(int f=0;f<3;++f)lin[f][0]=-(lin[f][1]*o[0]+lin[f][2]*o[1]+lin[f][3]*o[2]);
        scaleInv=1.0/sqrt(det);
        return true;
    }
    const double* pz=nd+54;double A[3][3];
    for(int k=0;k<3;++k){A[k][0]=px[k]-o[k];A[k][1]=py[k]-o[k];A[k][2]=pz[k]-o[k];}
    const double det=A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
    if(fabs(det)<1e-300)return false;
    if(type==pyramidHost){
        const double* top=nd+24;double scale=0.0,gap=0.0;
        for(int k=0;k<3;++k){const double v=top[k]-px[k]-py[k]+o[k];gap+=v*v;}
        for(int j=0;j<3;++j){const double column[3]={A[0][j],A[1][j],A[2][j]};scale=fmax(scale,norm3(column));}
        if(sqrt(gap)>1e-10*scale+1e-12)return false;
    }
    const double id=1.0/det;double B[3][3];
    B[0][0]=(A[1][1]*A[2][2]-A[1][2]*A[2][1])*id;B[0][1]=-(A[0][1]*A[2][2]-A[0][2]*A[2][1])*id;B[0][2]=(A[0][1]*A[1][2]-A[0][2]*A[1][1])*id;
    B[1][0]=-(A[1][0]*A[2][2]-A[1][2]*A[2][0])*id;B[1][1]=(A[0][0]*A[2][2]-A[0][2]*A[2][0])*id;B[1][2]=-(A[0][0]*A[1][2]-A[0][2]*A[1][0])*id;
    B[2][0]=(A[1][0]*A[2][1]-A[1][1]*A[2][0])*id;B[2][1]=-(A[0][0]*A[2][1]-A[0][1]*A[2][0])*id;B[2][2]=(A[0][0]*A[1][1]-A[0][1]*A[1][0])*id;
    for(int q=0;q<3;++q){lin[q][0]=-(B[q][0]*o[0]+B[q][1]*o[1]+B[q][2]*o[2]);lin[q][1]=B[q][0];lin[q][2]=B[q][1];lin[q][3]=B[q][2];}
    scaleInv=1.0/fabs(det);
    return true;
}
// Physical polynomials of nS charges about center. A rational charge gets zero coefficients: its
// potential comes from the radial rule.
RADIA_TYPED_KERNEL inline bool chargePolynomials(const double lin[3][4],double scaleInv,const double* center,bool cell,const int* exponents,int nS,double* out){
    const int axes=cell?3:2,ncoeff=cell?cellPolynomials:facePolynomials;double shifted[3][4];
    for(int a=0;a<axes;++a){shifted[a][0]=lin[a][0]+lin[a][1]*center[0]+lin[a][2]*center[1]+lin[a][3]*center[2];
        for(int k=1;k<4;++k)shifted[a][k]=lin[a][k];}
    bool finite=true;
    for(int l=0;l<nS;++l){const int* e=exponents+3*l;double* dst=out+std::size_t(l)*ncoeff;
        for(int i=0;i<ncoeff;++i)dst[i]=0.0;
        if(e[2]<0)continue;
        if(e[0]+e[1]+e[2]>(cell?1:2))return false;
        double poly[polynomialSlots];for(int i=0;i<polynomialSlots;++i)poly[i]=0.0;
        poly[0]=1.0;int deg=0;
        for(int a=0;a<axes;++a)for(int r=0;r<e[a];++r)polyTimesLinear(poly,deg,shifted[a],ncoeff);
        for(int i=0;i<ncoeff;++i){dst[i]=scaleInv*poly[i];finite=finite&&std::isfinite(dst[i]);}
    }
    return finite;
}
// Exact potentials of the polynomial charges of an affine tet, pyramid or triangle.
RADIA_TYPED_KERNEL inline void affinePotentials(int type,const double* nd,const double* center,const double* p,int nS,const double* coeff,double* inner){
    const double lp[3]={p[0]-center[0],p[1]-center[1],p[2]-center[2]};
    if(type==triHost){
        const int corner[3]={0,2,6};double V[3][3];
        for(int i=0;i<3;++i)for(int k=0;k<3;++k)V[i][k]=nd[3*corner[i]+k]-center[k];
        double moments[facePolynomials];RADIA_TYPED_FIELD::TriPotentialMomentsUpTo2(V,lp,moments);
        for(int l=0;l<nS;++l){const double* c=coeff+std::size_t(l)*facePolynomials;double v=0.0;
            for(int i=0;i<facePolynomials;++i)v+=c[i]*moments[i];inner[l]+=v;}
        return;
    }
    const int tets[2][4]={{0,2,8,18},{0,8,6,18}},tet[1][4]={{0,2,6,18}};const int count=type==tetHost?1:2;
    for(int sub=0;sub<count;++sub){const int* c4=type==tetHost?tet[0]:tets[sub];double V[4][3];
        for(int i=0;i<4;++i)for(int k=0;k<3;++k)V[i][k]=nd[3*c4[i]+k]-center[k];
        double moments[cellPolynomials];RADIA_TYPED_FIELD::TetPotentialMomentsUpTo1(V,lp,moments);
        for(int l=0;l<nS;++l){const double* c=coeff+std::size_t(l)*cellPolynomials;double v=0.0;
            for(int i=0;i<cellPolynomials;++i)v+=c[i]*moments[i];inner[l]+=v;}
    }
}
RADIA_TYPED_KERNEL inline void pyramidVertexMap(const double* v0,const double* v1,const double* v2,const double* v3,const double* v4,const double* y,double X[3],double J[3][3]){
    const double top=fmax(1.0-y[2],1e-300),u=y[0]/top,v=y[1]/top;
    for(int k=0;k<3;++k){const double e1=v1[k]-v0[k],e3=v3[k]-v0[k],b=v2[k]-v1[k]-v3[k]+v0[k],ez=v4[k]-v0[k];
        X[k]=v0[k]+y[0]*e1+y[1]*e3+y[2]*ez+y[0]*v*b;
        J[k][0]=e1+v*b;J[k][1]=e3+u*b;J[k][2]=ez+u*v*b;}
}
RADIA_TYPED_KERNEL inline void pyramidMap(const double* nd,const double* y,double X[3],double J[3][3]){
    pyramidVertexMap(nd,nd+6,nd+24,nd+18,nd+54,y,X,J);
}
RADIA_TYPED_KERNEL inline void clampPyramid(double* y){
    y[2]=fmin(fmax(y[2],0.0),1.0);const double top=1.0-y[2];
    y[0]=fmin(fmax(y[0],0.0),top);y[1]=fmin(fmax(y[1],0.0),top);
}
// Reference point of the pyramid nearest p: Newton on the map, staying inside the pyramid.
RADIA_TYPED_KERNEL inline void pyramidAnchor(const double* nd,const double* p,double y[3]){
    y[0]=y[1]=0.375;y[2]=0.25;
    for(int it=0;it<8;++it){
        double X[3],J[3][3];pyramidMap(nd,y,X,J);
        const double r[3]={p[0]-X[0],p[1]-X[1],p[2]-X[2]};
        const double det=J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
        if(!std::isfinite(det)||std::fabs(det)<1e-300)return;
        const double id=1.0/det;
        y[0]+=id*(r[0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(r[1]*J[2][2]-J[1][2]*r[2])+J[0][2]*(r[1]*J[2][1]-J[1][1]*r[2]));
        y[1]+=id*(J[0][0]*(r[1]*J[2][2]-J[1][2]*r[2])-r[0]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*r[2]-r[1]*J[2][0]));
        y[2]+=id*(J[0][0]*(J[1][1]*r[2]-r[1]*J[2][1])-J[0][1]*(J[1][0]*r[2]-r[1]*J[2][0])+r[0]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]));
        clampPyramid(y);
    }
}
// Radial cones from the anchor over the two reference sub-tets of a pyramid, for charges first..nS-1.
RADIA_TYPED_KERNEL inline void pyramidRadial(const double* nd,const double* p,const double* anchor,const int* exponents,int first,int nS,const double* gx,const double* gw,int nq,double* inner){
    const double corners[5][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1}};
    const int subs[2][4]={{0,1,2,4},{0,2,3,4}},faces[4][3]={{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
    for(int sub=0;sub<2;++sub){
        double V[4][3];for(int i=0;i<4;++i)for(int k=0;k<3;++k)V[i][k]=corners[subs[sub][i]][k];
        double x0[3];RADIA_TYPED_FIELD::ClosestPointTet(V,anchor,x0);
        for(int f=0;f<4;++f){
            const double* b1=V[faces[f][0]];const double* b2=V[faces[f][1]];const double* b3=V[faces[f][2]];
            double d1[3],d2[3],d3[3],e21[3],e32[3];
            for(int k=0;k<3;++k){d1[k]=b1[k]-x0[k];d2[k]=b2[k]-x0[k];d3[k]=b3[k]-x0[k];e21[k]=b2[k]-b1[k];e32[k]=b3[k]-b2[k];}
            const double D=fabs(d1[0]*(d2[1]*d3[2]-d2[2]*d3[1])-d1[1]*(d2[0]*d3[2]-d2[2]*d3[0])+d1[2]*(d2[0]*d3[1]-d2[1]*d3[0]));
            if(D<1e-300)continue;
            for(int a=0;a<nq;++a)for(int b=0;b<nq;++b)for(int c=0;c<nq;++c){
                const double u=gx[a],v=gx[b],w=gx[c];double y[3],X[3],J[3][3];
                for(int k=0;k<3;++k)y[k]=x0[k]+u*(d1[k]+v*(e21[k]+w*e32[k]));
                pyramidMap(nd,y,X,J);
                const double dr[3]={p[0]-X[0],p[1]-X[1],p[2]-X[2]};const double r=norm3(dr);if(r<1e-300)continue;
                const double weight=gw[a]*gw[b]*gw[c]*u*u*v*D/r;
                for(int l=first;l<nS;++l)inner[l]+=weight*monomial(exponents+3*l,y);
            }
        }
    }
}
// Potential of the rational charge x y/(1-z) of an affine pyramid. Each level z is a parallelogram carrying the
// polynomial x y, whose potential is exact; the rule in z is split at the level of p, where the kernel peaks, and graded.
RADIA_TYPED_KERNEL inline double pyramidRationalPotential(const double* nd,const double* center,const double* p,const double* gx,const double* gw,int nq){
    double lin[3][4],scaleInv,e1[3],e3[3],ez[3];
    if(!affineForms(pyramidHost,nd,lin,scaleInv))return NAN;
    for(int k=0;k<3;++k){e1[k]=nd[6+k]-nd[k];e3[k]=nd[18+k]-nd[k];ez[k]=nd[54+k]-nd[k];}
    const double normal[3]={e1[1]*e3[2]-e1[2]*e3[1],e1[2]*e3[0]-e1[0]*e3[2],e1[0]*e3[1]-e1[1]*e3[0]};
    const int product[3]={1,1,0};double coeff[facePolynomials];
    if(!chargePolynomials(lin,1.0/norm3(normal),center,false,product,1,coeff))return NAN;
    const double lp[3]={p[0]-center[0],p[1]-center[1],p[2]-center[2]};
    const double level=fmin(fmax(lin[2][0]+lin[2][1]*p[0]+lin[2][2]*p[1]+lin[2][3]*p[2],0.0),1.0);
    double total=0.0;
    for(int part=0;part<2;++part){
        const double a=part?level:0.0,b=part?1.0:level;
        if(b<=a)continue;
        for(int i=0;i<nq;++i){
            const double s=gx[i],t=s*s*s*(10.0-15.0*s+6.0*s*s),dt=30.0*s*s*(1.0-s)*(1.0-s);
            const double z=a+(b-a)*t,top=1.0-z;
            if(top<=0.0)continue;
            double corners[4][3];
            for(int k=0;k<3;++k){corners[0][k]=nd[k]+z*ez[k]-center[k];corners[1][k]=corners[0][k]+top*e1[k];
                corners[2][k]=corners[1][k]+top*e3[k];corners[3][k]=corners[0][k]+top*e3[k];}
            double value=0.0;
            for(int tri=0;tri<2;++tri){
                double V[3][3],moments[facePolynomials];
                for(int k=0;k<3;++k){V[0][k]=corners[0][k];V[1][k]=corners[1+tri][k];V[2][k]=corners[2+tri][k];}
                RADIA_TYPED_FIELD::TriPotentialMomentsUpTo2(V,lp,moments);
                for(int m=0;m<facePolynomials;++m)value+=coeff[m]*moments[m];
            }
            total+=gw[i]*dt*(b-a)*value/top;
        }
    }
    return total;
}
} // namespace radia_hdiv_typed
