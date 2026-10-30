// CUDA loop fusion of existing CurvedTri/TetPotential Duffy rules.
// Geometry is shared across monomial components; each sum keeps its original order.
#pragma once
struct MomentSourceGroup {std::size_t member[6];std::size_t count=0;};
__device__ void curvedTriMomentPotentials(const CurvedCharge* charges,const MomentSourceGroup& group,const double* p,
    const double* gl,const double* gw,int nq,double* acc){
    const auto* nodes=reinterpret_cast<const double(*)[3]>(charges[group.member[0]].nodes);
    double xi0[2];rad_hdiv_cuda::ClosestRefTri(nodes,p,xi0);
    const double C[3][2]={{0,0},{1,0},{0,1}};
    for(int k=0;k<3;++k){const double* A=C[k];const double* B=C[(k+1)%3];
        const double e1x=A[0]-xi0[0],e1y=A[1]-xi0[1],e2x=B[0]-xi0[0],e2y=B[1]-xi0[1];
        const double sgn2=e1x*e2y-e1y*e2x;
        for(int a=0;a<nq;++a){const double u=gl[a];for(int b=0;b<nq;++b){const double v=gl[b];
            const double xi=xi0[0]+u*e1x+u*v*(e2x-e1x),eta=xi0[1]+u*e1y+u*v*(e2y-e1y);
            double X[3],Xu[3],Xv[3];rad_hdiv_cuda::CurvedTriEval(nodes,xi,eta,X,Xu,Xv);
            const double cr[3]={Xu[1]*Xv[2]-Xu[2]*Xv[1],Xu[2]*Xv[0]-Xu[0]*Xv[2],Xu[0]*Xv[1]-Xu[1]*Xv[0]};
            const double J=sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
            const double dx=p[0]-X[0],dy=p[1]-X[1],dz=p[2]-X[2],r=sqrt(dx*dx+dy*dy+dz*dz);if(r<1e-300)continue;
            for(std::size_t m=0;m<group.count;++m){const auto& charge=charges[group.member[m]];
                acc[m]+=gw[a]*gw[b]*(u*sgn2)*J*rad_hdiv_cuda::_ipow(xi,charge.expo[0])*rad_hdiv_cuda::_ipow(eta,charge.expo[1])/r;
            }
        }}
    }
}
__device__ void curvedTetMomentPotentials(const CurvedCharge* charges,const MomentSourceGroup& group,const double* p,
    const double* gl,const double* gw,int nq,double* acc){
    const auto* nodes=reinterpret_cast<const double(*)[3]>(charges[group.member[0]].nodes);
    double xi0[3];rad_hdiv_cuda::ClosestRefTet(nodes,p,xi0);
    const double C[4][3]={{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
    const int FC[4][3]={{1,2,3},{0,3,2},{0,1,3},{2,1,0}};
    for(int f=0;f<4;++f)for(int lead=0;lead<3;++lead){
        const double* b1=C[FC[f][lead]];const double* b2=C[FC[f][(lead+1)%3]];const double* b3=C[FC[f][(lead+2)%3]];
        double d1[3],d2[3],d3[3],e21[3],e32[3];
        for(int k=0;k<3;++k){d1[k]=b1[k]-xi0[k];d2[k]=b2[k]-xi0[k];d3[k]=b3[k]-xi0[k];e21[k]=b2[k]-b1[k];e32[k]=b3[k]-b2[k];}
        const double cr[3]={d2[1]*d3[2]-d2[2]*d3[1],d2[2]*d3[0]-d2[0]*d3[2],d2[0]*d3[1]-d2[1]*d3[0]};
        const double D=d1[0]*cr[0]+d1[1]*cr[1]+d1[2]*cr[2];if(fabs(D)<1e-300)continue;
        for(int a=0;a<nq;++a){const double u=gl[a];for(int b=0;b<nq;++b){const double v=gl[b];for(int c=0;c<nq;++c){const double w=gl[c];
            double z[3];for(int k=0;k<3;++k)z[k]=xi0[k]+u*(d1[k]+v*(e21[k]+w*e32[k]));
            double X[3],Jac[3][3];rad_hdiv_cuda::CurvedTetEval(nodes,z[0],z[1],z[2],X,Jac);
            const double Jv=fabs(rad_hdiv_cuda::det3(Jac));const double dx=p[0]-X[0],dy=p[1]-X[1],dz=p[2]-X[2],r=sqrt(dx*dx+dy*dy+dz*dz);if(r<1e-300)continue;
            for(std::size_t m=0;m<group.count;++m){const auto& charge=charges[group.member[m]];
                acc[m]+=(gw[a]*gw[b]*gw[c]/3.0)*(u*u*v*D)*Jv*rad_hdiv_cuda::_ipow(z[0],charge.expo[0])*rad_hdiv_cuda::_ipow(z[1],charge.expo[1])*rad_hdiv_cuda::_ipow(z[2],charge.expo[2])/r;
            }
        }}}
    }
}
