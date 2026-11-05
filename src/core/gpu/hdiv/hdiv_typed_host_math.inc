// Shared arithmetic of tetrahedral, pyramidal and triangular charge hosts, instantiated with the
// production CPU potentials in host controls and the same device potentials in CUDA.
// A host is stored as the Q2 lattice of its map from the unit cube, which is trilinear for every
// type: the cube corners collapse onto the tet (V3 nodes 0, V0 2, V1 6, V2 18), the pyramid
// (V0 0, V1 2, V2 8, V3 6, V4 18) and the triangle (V2 0, V0 2, V1 6). A charge density is a monomial
// of the NGSolve reference coordinates in the reference measure, whose Jacobian cancels the Piola one.
// The order 2 pyramid adds the rational div image x y / (1 - z), coded by a z exponent of -1.
namespace radia_hdiv_typed {
enum HostType {hexHost=0,quadHost=1,tetHost=2,triHost=3,pyramidHost=4};
// The math headers of the three compilers of this source declare isfinite differently; false for NaN and infinities.
RADIA_TYPED_KERNEL inline bool finiteValue(double v){return fabs(v)<=1.7976931348623157e308;}
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
        for(int i=0;i<ncoeff;++i){dst[i]=scaleInv*poly[i];finite=finite&&finiteValue(dst[i]);}
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
        if(!finiteValue(det)||fabs(det)<1e-300)return;
        const double id=1.0/det;
        y[0]+=id*(r[0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(r[1]*J[2][2]-J[1][2]*r[2])+J[0][2]*(r[1]*J[2][1]-J[1][1]*r[2]));
        y[1]+=id*(J[0][0]*(r[1]*J[2][2]-J[1][2]*r[2])-r[0]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*r[2]-r[1]*J[2][0]));
        y[2]+=id*(J[0][0]*(J[1][1]*r[2]-r[1]*J[2][1])-J[0][1]*(J[1][0]*r[2]-r[1]*J[2][0])+r[0]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]));
        clampPyramid(y);
    }
}
// Radial cones from the anchor over the two reference sub-tets of a pyramid, for charges first..nS-1. The source, the
// exponents and the sums stay in locals, so that nothing is reloaded past a store at each of the many points of the rule.
RADIA_TYPED_KERNEL inline void pyramidRadial(const double* nd,const double* p,const double* anchor,const int* exponents,int first,int nS,const double* gx,const double* gw,int nq,double* inner){
    constexpr int pyramidCharges=5;
    if(nS>pyramidCharges){for(int l=first;l<nS;++l)inner[l]=NAN;return;}
    double v0[3],e1[3],e3[3],ez[3],bilinear[3],sum[pyramidCharges];int powers[pyramidCharges][3];
    for(int k=0;k<3;++k){v0[k]=nd[k];e1[k]=nd[6+k]-nd[k];e3[k]=nd[18+k]-nd[k];bilinear[k]=nd[24+k]-nd[6+k]-nd[18+k]+nd[k];ez[k]=nd[54+k]-nd[k];}
    for(int l=first;l<nS;++l){sum[l]=inner[l];for(int k=0;k<3;++k)powers[l][k]=exponents[3*l+k];}
    const double px=p[0],py=p[1],pz=p[2];
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
            for(int a=0;a<nq;++a){const double u=gx[a],wa=gw[a];
                for(int b=0;b<nq;++b){const double v=gx[b],wab=wa*gw[b];
                    for(int c=0;c<nq;++c){const double w=gx[c];double y[3],X[3];
                        for(int k=0;k<3;++k)y[k]=x0[k]+u*(d1[k]+v*(e21[k]+w*e32[k]));
                        const double top=fmax(1.0-y[2],1e-300),share=y[1]/top;
                        for(int k=0;k<3;++k)X[k]=v0[k]+y[0]*e1[k]+y[1]*e3[k]+y[2]*ez[k]+y[0]*share*bilinear[k];
                        const double dr[3]={px-X[0],py-X[1],pz-X[2]};const double r=norm3(dr);if(r<1e-300)continue;
                        const double weight=wab*gw[c]*u*u*v*D/r;
                        for(int l=first;l<nS;++l)sum[l]+=weight*monomial(powers[l],y);
                    }
                }
            }
        }
    }
    for(int l=first;l<nS;++l)inner[l]=sum[l];
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
RADIA_TYPED_KERNEL inline void hexLag3(double t,double v[3],double d[3]){
    v[0]=2.0*(t-0.5)*(t-1.0);v[1]=4.0*t*(1.0-t);v[2]=2.0*t*(t-0.5);
    d[0]=4.0*t-3.0;d[1]=4.0-8.0*t;d[2]=4.0*t-1.0;
}
RADIA_TYPED_KERNEL inline void hexMapX(const double* nd,const double* xi,double X[3]){
    double vx[3],dx[3],vy[3],dy[3],vz[3],dz[3];hexLag3(xi[0],vx,dx);hexLag3(xi[1],vy,dy);hexLag3(xi[2],vz,dz);
    X[0]=X[1]=X[2]=0.0;
    for(int iz=0;iz<3;++iz)for(int iy=0;iy<3;++iy){const double syz=vy[iy]*vz[iz];
        for(int ix=0;ix<3;++ix){const double* n=nd+3*(ix+3*iy+9*iz);const double s=vx[ix]*syz;X[0]+=s*n[0];X[1]+=s*n[1];X[2]+=s*n[2];}}
}
RADIA_TYPED_KERNEL inline void hexMapJ(const double* nd,const double* xi,double X[3],double J[3][3]){
    double vx[3],dx[3],vy[3],dy[3],vz[3],dz[3];hexLag3(xi[0],vx,dx);hexLag3(xi[1],vy,dy);hexLag3(xi[2],vz,dz);
    for(int k=0;k<3;++k){X[k]=0.0;J[k][0]=J[k][1]=J[k][2]=0.0;}
    for(int iz=0;iz<3;++iz)for(int iy=0;iy<3;++iy)for(int ix=0;ix<3;++ix){const double* n=nd+3*(ix+3*iy+9*iz);
        const double s=vx[ix]*vy[iy]*vz[iz],gx=dx[ix]*vy[iy]*vz[iz],gy=vx[ix]*dy[iy]*vz[iz],gz=vx[ix]*vy[iy]*dz[iz];
        for(int k=0;k<3;++k){X[k]+=s*n[k];J[k][0]+=gx*n[k];J[k][1]+=gy*n[k];J[k][2]+=gz*n[k];}}
}
RADIA_TYPED_KERNEL inline void quadMapX(const double* nd,const double* uv,double X[3]){
    double vu[3],du[3],vv[3],dv[3];hexLag3(uv[0],vu,du);hexLag3(uv[1],vv,dv);X[0]=X[1]=X[2]=0.0;
    for(int iv=0;iv<3;++iv)for(int iu=0;iu<3;++iu){const double* n=nd+3*(iu+3*iv);const double s=vu[iu]*vv[iv];X[0]+=s*n[0];X[1]+=s*n[1];X[2]+=s*n[2];}
}
RADIA_TYPED_KERNEL inline void quadMapT(const double* nd,const double* uv,double X[3],double T[3][2]){
    double vu[3],du[3],vv[3],dv[3];hexLag3(uv[0],vu,du);hexLag3(uv[1],vv,dv);
    for(int k=0;k<3;++k){X[k]=0.0;T[k][0]=T[k][1]=0.0;}
    for(int iv=0;iv<3;++iv)for(int iu=0;iu<3;++iu){const double* n=nd+3*(iu+3*iv);const double s=vu[iu]*vv[iv];
        for(int k=0;k<3;++k){X[k]+=s*n[k];T[k][0]+=du[iu]*vv[iv]*n[k];T[k][1]+=vu[iu]*dv[iv]*n[k];}}
}
// Tangent rule of a warped source. With F linearized at a base point to T, the potential of a charge at p is the closed
// form on the straight image T(host) plus the integral of density*(1/|p-F| - 1/|p-T|), which stays bounded at the base
// point because F - T is of second order there. A density is an exact polynomial in the reference coordinates, so the
// closed form has no fit error, and the remainder takes a plain rule on the unit cube.
// Reference point of a hex nearest p: Newton on the map, staying inside the cube.
RADIA_TYPED_KERNEL inline void hexAnchor(const double* nd,const double* p,double y[3]){
    y[0]=y[1]=y[2]=0.5;
    for(int it=0;it<8;++it){
        double X[3],J[3][3];hexMapJ(nd,y,X,J);
        const double r[3]={p[0]-X[0],p[1]-X[1],p[2]-X[2]};
        const double det=J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
        if(!finiteValue(det)||fabs(det)<1e-300)return;
        const double id=1.0/det;
        y[0]+=id*(r[0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(r[1]*J[2][2]-J[1][2]*r[2])+J[0][2]*(r[1]*J[2][1]-J[1][1]*r[2]));
        y[1]+=id*(J[0][0]*(r[1]*J[2][2]-J[1][2]*r[2])-r[0]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*r[2]-r[1]*J[2][0]));
        y[2]+=id*(J[0][0]*(J[1][1]*r[2]-r[1]*J[2][1])-J[0][1]*(J[1][0]*r[2]-r[1]*J[2][0])+r[0]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]));
        for(int k=0;k<3;++k)y[k]=fmin(fmax(y[k],0.0),1.0);
    }
}
// Reference point of a quad nearest p: projected Gauss-Newton on the unit square.
RADIA_TYPED_KERNEL inline void quadAnchor(const double* nd,const double* p,double y[3]){
    y[0]=y[1]=0.5;y[2]=0.0;
    for(int it=0;it<8;++it){
        double X[3],T[3][2];quadMapT(nd,y,X,T);
        const double r[3]={p[0]-X[0],p[1]-X[1],p[2]-X[2]};
        const double aa=T[0][0]*T[0][0]+T[1][0]*T[1][0]+T[2][0]*T[2][0],ab=T[0][0]*T[0][1]+T[1][0]*T[1][1]+T[2][0]*T[2][1],bb=T[0][1]*T[0][1]+T[1][1]*T[1][1]+T[2][1]*T[2][1];
        const double ra=T[0][0]*r[0]+T[1][0]*r[1]+T[2][0]*r[2],rb=T[0][1]*r[0]+T[1][1]*r[1]+T[2][1]*r[2],det=aa*bb-ab*ab;
        if(!finiteValue(det)||det<1e-300)return;
        y[0]=fmin(fmax(y[0]+(bb*ra-ab*rb)/det,0.0),1.0);y[1]=fmin(fmax(y[1]+(aa*rb-ab*ra)/det,0.0),1.0);
    }
}
// One linear form on a moment functional: out[k] = l(x) x^k, for the moments of degree <= deg.
RADIA_TYPED_KERNEL inline void momentsTimesLinear(const double* in,const double l[4],int deg,double* out){
    for(int total=0;total<=deg;++total)for(int ax=0;ax<=total;++ax)for(int ay=0;ay<=total-ax;++ay){const int az=total-ax-ay;
        out[polyIndex(ax,ay,az)]=l[0]*in[polyIndex(ax,ay,az)]+l[1]*in[polyIndex(ax+1,ay,az)]+l[2]*in[polyIndex(ax,ay+1,az)]+l[3]*in[polyIndex(ax,ay,az+1)];}
}
// The first and, at order 2, second power of one form on the moments, leaving degree <= keep after the last power.
RADIA_TYPED_KERNEL inline void axisPowers(const double* in,const double l[4],int order,int keep,double* first,double* second){
    momentsTimesLinear(in,l,keep+order-1,first);
    if(order==2)momentsTimesLinear(first,l,keep,second);
}
// Potentials of the tensor monomials of exponents 0..order on a straight parallelogram (index a+3b) or parallelepiped
// (index a+3b+9c), from the physical moments: the inverse forms act on the moments, one axis at a time.
RADIA_TYPED_KERNEL inline void faceTensorPotentials(const double* M,const double forms[][4],int order,int offset,double* N){
    double y1[20],y2[10];axisPowers(M,forms[1],order,order,y1,y2);
    for(int b=0;b<=order;++b){const double* Mb=b==0?M:b==1?y1:y2;double x1[4],x2[1];axisPowers(Mb,forms[0],order,0,x1,x2);
        N[offset+3*b]=Mb[0];N[offset+1+3*b]=x1[0];if(order==2)N[offset+2+3*b]=x2[0];}
}
RADIA_TYPED_KERNEL inline void cellTensorPotentials(const double* M,const double forms[3][4],int order,double* N){
    double z1[56],z2[35];axisPowers(M,forms[2],order,2*order,z1,z2);
    for(int c=0;c<=order;++c)faceTensorPotentials(c==0?M:c==1?z1:z2,forms,order,9*c,N);
}
// Difference of the kernels at p; a rule point on p itself takes the bounded integrand as zero.
RADIA_TYPED_KERNEL inline double kernelDifference(const double* p,const double* F,const double* T,double guard){
    const double a[3]={p[0]-F[0],p[1]-F[1],p[2]-F[2]},b[3]={p[0]-T[0],p[1]-T[1],p[2]-T[2]};
    const double ra=norm3(a),rb=norm3(b);
    return ra>guard&&rb>guard?1.0/ra-1.0/rb:0.0;
}
RADIA_TYPED_KERNEL inline bool hexTangent(const double* nd,const double* base,const double* p,const int* exponents,int nS,int order,const double* gx,const double* gw,int nq,double* inner){
    double X0[3],J[3][3];hexMapJ(nd,base,X0,J);
    const double m00=J[1][1]*J[2][2]-J[1][2]*J[2][1],m01=J[1][0]*J[2][2]-J[1][2]*J[2][0],m02=J[1][0]*J[2][1]-J[1][1]*J[2][0];
    const double det=J[0][0]*m00-J[0][1]*m01+J[0][2]*m02;
    const double col[3][3]={{J[0][0],J[1][0],J[2][0]},{J[0][1],J[1][1],J[2][1]},{J[0][2],J[1][2],J[2][2]}};
    const double n0=norm3(col[0]),n1=norm3(col[1]),n2=norm3(col[2]);
    if(!finiteValue(det)||!(fabs(det)>1e-10*n0*n1*n2))return false;
    const double id=1.0/det,mid[3]={0.5-base[0],0.5-base[1],0.5-base[2]};
    double forms[3][4]={{0.5,m00*id,-(J[0][1]*J[2][2]-J[0][2]*J[2][1])*id,(J[0][1]*J[1][2]-J[0][2]*J[1][1])*id},
                        {0.5,-m01*id,(J[0][0]*J[2][2]-J[0][2]*J[2][0])*id,-(J[0][0]*J[1][2]-J[0][2]*J[1][0])*id},
                        {0.5,m02*id,-(J[0][0]*J[2][1]-J[0][1]*J[2][0])*id,(J[0][0]*J[1][1]-J[0][1]*J[1][0])*id}};
    double center[3],lp[3],C[8][3];
    for(int k=0;k<3;++k){center[k]=X0[k]+J[k][0]*mid[0]+J[k][1]*mid[1]+J[k][2]*mid[2];lp[k]=p[k]-center[k];}
    // One face in order and then the opposite one.
    const double hx[8]={-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5},hy[8]={-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5},hz[8]={-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5};
    for(int i=0;i<8;++i)for(int k=0;k<3;++k)C[i][k]=J[k][0]*hx[i]+J[k][1]*hy[i]+J[k][2]*hz[i];
    double moments[84],N[27];
    if(order==2)RADIA_TYPED_FIELD::HexPotentialMomentsUpTo6(C,lp,moments);else RADIA_TYPED_FIELD::HexPotentialMomentsUpTo3(C,lp,moments);
    cellTensorPotentials(moments,forms,order,N);
    const double scale=1.0/fabs(det),guard=1e-6*fmax(n0,fmax(n1,n2));
    const int modes=(order+1)*(order+1)*(order+1);
    for(int m=0;m<modes;++m)N[m]*=scale;
    for(int iz=0;iz<nq;++iz)for(int iy=0;iy<nq;++iy)for(int ix=0;ix<nq;++ix){
        const double eta[3]={gx[ix],gx[iy],gx[iz]},d[3]={eta[0]-base[0],eta[1]-base[1],eta[2]-base[2]};
        double F[3],T[3];hexMapX(nd,eta,F);
        for(int k=0;k<3;++k)T[k]=X0[k]+J[k][0]*d[0]+J[k][1]*d[1]+J[k][2]*d[2];
        const double delta=gw[ix]*gw[iy]*gw[iz]*kernelDifference(p,F,T,guard);
        double pz=delta;
        for(int c=0;c<=order;++c,pz*=eta[2]){double py=pz;
            for(int b=0;b<=order;++b,py*=eta[1]){double px=py;
                for(int a=0;a<=order;++a,px*=eta[0])N[a+3*b+9*c]+=px;}}
    }
    for(int l=0;l<nS;++l){const int* e=exponents+3*l;if(e[0]>order||e[1]>order||e[2]>order)return false;inner[l]+=N[e[0]+3*e[1]+9*e[2]];}
    return true;
}
RADIA_TYPED_KERNEL inline bool quadTangent(const double* nd,const double* base,const double* p,const int* exponents,int nS,int order,const double* gx,const double* gw,int nq,double* inner){
    double X0[3],T0[3][2];quadMapT(nd,base,X0,T0);
    const double aa=T0[0][0]*T0[0][0]+T0[1][0]*T0[1][0]+T0[2][0]*T0[2][0],ab=T0[0][0]*T0[0][1]+T0[1][0]*T0[1][1]+T0[2][0]*T0[2][1],bb=T0[0][1]*T0[0][1]+T0[1][1]*T0[1][1]+T0[2][1]*T0[2][1];
    const double det=aa*bb-ab*ab;
    if(!finiteValue(det)||!(det>1e-20*aa*bb))return false;
    const double i00=bb/det,i01=-ab/det,i11=aa/det,mid[2]={0.5-base[0],0.5-base[1]};
    double forms[2][4]={{0.5,0.0,0.0,0.0},{0.5,0.0,0.0,0.0}},center[3],lp[3],V[4][3];
    for(int k=0;k<3;++k){forms[0][k+1]=i00*T0[k][0]+i01*T0[k][1];forms[1][k+1]=i01*T0[k][0]+i11*T0[k][1];
        center[k]=X0[k]+T0[k][0]*mid[0]+T0[k][1]*mid[1];lp[k]=p[k]-center[k];}
    const double hu[4]={-0.5,0.5,0.5,-0.5},hv[4]={-0.5,-0.5,0.5,0.5};
    for(int i=0;i<4;++i)for(int k=0;k<3;++k)V[i][k]=T0[k][0]*hu[i]+T0[k][1]*hv[i];
    double moments[35],N[9];
    if(order==2)RADIA_TYPED_FIELD::QuadPotentialMomentsUpTo4(V,lp,moments);else RADIA_TYPED_FIELD::QuadPotentialMomentsUpTo2(V,lp,moments);
    faceTensorPotentials(moments,forms,order,0,N);
    const double scale=1.0/sqrt(det),guard=1e-6*sqrt(fmax(aa,bb));
    const int modes=(order+1)*(order+1);
    for(int m=0;m<modes;++m)N[m]*=scale;
    for(int iy=0;iy<nq;++iy)for(int ix=0;ix<nq;++ix){
        const double eta[2]={gx[ix],gx[iy]},d[2]={eta[0]-base[0],eta[1]-base[1]};
        double F[3],T[3];quadMapX(nd,eta,F);
        for(int k=0;k<3;++k)T[k]=X0[k]+T0[k][0]*d[0]+T0[k][1]*d[1];
        double py=gw[ix]*gw[iy]*kernelDifference(p,F,T,guard);
        for(int b=0;b<=order;++b,py*=eta[1]){double px=py;
            for(int a=0;a<=order;++a,px*=eta[0])N[a+3*b]+=px;}
    }
    for(int l=0;l<nS;++l){const int* e=exponents+3*l;if(e[0]>order||e[1]>order)return false;inner[l]+=N[e[0]+3*e[1]];}
    return true;
}
// At the reference point (x0,y0,z0), with u = x0/(1-z0) and v = y0/(1-z0), the first order Taylor polynomial of the rational
// charge x y/(1-z) is q = v x + u y - u v (1-z), and x y/(1-z) - q = (1-z)(x/(1-z) - u)(y/(1-z) - v). So the charge is a
// combination of the polynomial charges plus a difference that vanishes at the point, where a plain rule then converges.
RADIA_TYPED_KERNEL inline void rationalBase(const double* y,double& u,double& v){
    const double top=1.0-y[2];u=top>1e-300?y[0]/top:0.0;v=top>1e-300?y[1]/top:0.0;
}
RADIA_TYPED_KERNEL inline double rationalTaylorPotential(const int* exponents,int nS,const double* potentials,double u,double v){
    double sum=0.0;
    for(int l=0;l<nS;++l){const int* e=exponents+3*l;
        if(e[2]<0)continue;
        sum+=(e[0]==1?v:e[1]==1?u:e[2]==1?u*v:-u*v)*potentials[l];}
    return sum;
}
// Rational charge of a straight pyramid at p; inner holds the potentials of its polynomial charges and nothing else.
RADIA_TYPED_KERNEL inline void straightPyramidRational(const double* nd,const double* p,const int* exponents,int nS,const double* gx,const double* gw,int nq,double* inner){
    int rational=-1;for(int l=0;l<nS;++l)if(exponents[3*l+2]<0)rational=l;
    if(rational<0)return;
    double lin[3][4],scaleInv=0.0;
    if(!affineForms(pyramidHost,nd,lin,scaleInv)){inner[rational]=NAN;return;}
    double y[3],u,v;for(int a=0;a<3;++a)y[a]=lin[a][0]+lin[a][1]*p[0]+lin[a][2]*p[1]+lin[a][3]*p[2];
    clampPyramid(y);rationalBase(y,u,v);
    double sum=rationalTaylorPotential(exponents,nS,inner,u,v);
    for(int iz=0;iz<nq;++iz)for(int iy=0;iy<nq;++iy)for(int ix=0;ix<nq;++ix){
        const double c[3]={gx[ix],gx[iy],gx[iz]},top=1.0-c[2];double point[3],X[3],J[3][3];
        collapse(pyramidHost,c,point);pyramidMap(nd,point,X,J);
        const double d[3]={p[0]-X[0],p[1]-X[1],p[2]-X[2]},r=norm3(d);
        if(r>1e-300)sum+=gw[ix]*gw[iy]*gw[iz]*top*top*top*(c[0]-u)*(c[1]-v)/r;
    }
    inner[rational]=sum;
}
// The rational charge takes its Taylor polynomial on the tangent pyramid, whose base is a parallelogram.
RADIA_TYPED_KERNEL inline bool pyramidTangent(const double* nd,const double* base,const double* p,const int* exponents,int nS,const double* gx,const double* gw,int nq,double* inner){
    constexpr int pyramidCharges=5;
    if(nS>pyramidCharges)return false;
    double X0[3],J[3][3];pyramidMap(nd,base,X0,J);
    const double col[3][3]={{J[0][0],J[1][0],J[2][0]},{J[0][1],J[1][1],J[2][1]},{J[0][2],J[1][2],J[2][2]}};
    const double guard=1e-6*fmax(norm3(col[0]),fmax(norm3(col[1]),norm3(col[2])));
    // The five corners at their lattice nodes; the straight routines read no other node.
    double tangent[57],lin[3][4],scaleInv=0.0,center[3],coeff[pyramidCharges*cellPolynomials];
    const int node[5]={0,2,8,6,18};const double corner[5][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1}},middle[3]={0.25,0.25,0.5};
    for(int i=0;i<5;++i)for(int k=0;k<3;++k)tangent[3*node[i]+k]=X0[k]+J[k][0]*(corner[i][0]-base[0])+J[k][1]*(corner[i][1]-base[1])+J[k][2]*(corner[i][2]-base[2]);
    for(int k=0;k<3;++k)center[k]=X0[k]+J[k][0]*(middle[0]-base[0])+J[k][1]*(middle[1]-base[1])+J[k][2]*(middle[2]-base[2]);
    if(!affineForms(pyramidHost,tangent,lin,scaleInv)||!chargePolynomials(lin,scaleInv,center,true,exponents,nS,coeff))return false;
    double potentials[pyramidCharges],u,v;
    for(int l=0;l<nS;++l)potentials[l]=0.0;
    affinePotentials(pyramidHost,tangent,center,p,nS,coeff,potentials);
    rationalBase(base,u,v);
    for(int l=0;l<nS;++l)inner[l]+=exponents[3*l+2]<0?rationalTaylorPotential(exponents,nS,potentials,u,v):potentials[l];
    for(int iz=0;iz<nq;++iz)for(int iy=0;iy<nq;++iy)for(int ix=0;ix<nq;++ix){
        const double c[3]={gx[ix],gx[iy],gx[iz]};double y[3],F[3],T[3],JF[3][3];
        const double det=collapse(pyramidHost,c,y);pyramidMap(nd,y,F,JF);
        for(int k=0;k<3;++k)T[k]=X0[k]+J[k][0]*(y[0]-base[0])+J[k][1]*(y[1]-base[1])+J[k][2]*(y[2]-base[2]);
        const double a[3]={p[0]-F[0],p[1]-F[1],p[2]-F[2]},b[3]={p[0]-T[0],p[1]-T[1],p[2]-T[2]},ra=norm3(a),rb=norm3(b);
        if(!(ra>guard&&rb>guard))continue;
        const double w=gw[ix]*gw[iy]*gw[iz]*det,polynomial=v*y[0]+u*y[1]-u*v*(1.0-y[2]);
        for(int l=0;l<nS;++l){const int* e=exponents+3*l;const double density=monomial(e,y);
            inner[l]+=w*(e[2]<0?density/ra-polynomial/rb:density*(1.0/ra-1.0/rb));}
    }
    return true;
}
// Potentials of a warped hex, quad or pyramid source at p. The base is the target's own point on a self pair, whose unit
// cube point comes in selfCube, otherwise the clamped preimage of p.
RADIA_TYPED_KERNEL inline bool warpedPotentials(int type,const double* nd,const double* p,const double* selfCube,const int* exponents,int nS,int order,const double* gx,const double* gw,int nq,double* inner){
    if(order<1||order>2)return false;
    double base[3];
    if(type==pyramidHost){
        if(selfCube)collapse(type,selfCube,base);else pyramidAnchor(nd,p,base);
        return pyramidTangent(nd,base,p,exponents,nS,gx,gw,nq,inner);
    }
    if(type==hexHost){
        if(selfCube){base[0]=selfCube[0];base[1]=selfCube[1];base[2]=selfCube[2];}else hexAnchor(nd,p,base);
        return hexTangent(nd,base,p,exponents,nS,order,gx,gw,nq,inner);
    }
    if(type==quadHost){
        if(selfCube){base[0]=selfCube[0];base[1]=selfCube[1];base[2]=0.0;}else quadAnchor(nd,p,base);
        return quadTangent(nd,base,p,exponents,nS,order,gx,gw,nq,inner);
    }
    return false;
}
} // namespace radia_hdiv_typed
