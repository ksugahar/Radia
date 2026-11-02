// Device construction of the TET HDiv charge map B and B^T.
// NGSolve tabulates charge-projection tables (affine) or reference moments
// (curved P2) per vertex-rank class on the host; only those tables, the Gram
// geometry, DOF numbers, classes and the curved mass rules cross.
#pragma once
#include <climits>
#include <cstdint>
#include <cub/device/device_radix_sort.cuh>
#include "hdiv_resident_python_api.hpp"
#include "hdiv_flat_evaluation.cuh"
namespace radia_hdiv_bh {
// Curved P2 hosts: x solves G x = s r for one table column r, with
// G_ab = sum_q w_q m_a m_b |dX| on the declared mass rule (cells |det J| in xi,
// faces |X_e0 x X_e1| in the boundary element's eta) and s = sign det J on
// cells. One thread per (host, column); every column recomputes the same G
// and pivots identically (the host solveSmall elimination).
__device__ inline double curvedChargeMonomial(const int* e,const double* x,bool cell){
    double m=1.;for(int c=0;c<(cell?3:2);++c)for(int k=0;k<e[c];++k)m*=x[c];return m;
}
__global__ void curvedChargeTables(std::size_t cells,std::size_t faces,std::size_t cm,std::size_t fm,std::size_t cw,std::size_t fw,
    std::size_t cellClasses,std::size_t faceClasses,const double* cellNodes,const double* faceNodes,const int* cellClass,const int* faceClass,
    const double* cellMoments,const double* faceMoments,const double* cellRule,int cellPoints,const double* faceRule,int facePoints,
    const int* cellExponents,const int* faceExponents,double* cellOut,double* faceOut,int* invalid){
    const auto t=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;const auto cellColumns=cells*cw;if(t>=cellColumns+faces*fw)return;
    const bool cell=t<cellColumns;const auto n=cell?cm:fm,width=cell?cw:fw;
    const auto host=cell?t/cw:(t-cellColumns)/fw,column=cell?t%cw:(t-cellColumns)%fw;
    const int cls=cell?cellClass[host]:faceClass[host];
    if(cls<0||std::size_t(cls)>=(cell?cellClasses:faceClasses)){atomicExch(invalid,1);return;}
    double G[36]={},x[6];int sign=0;
    const double* rule=cell?cellRule:faceRule;const int points=cell?cellPoints:facePoints;const int* exponents=cell?cellExponents:faceExponents;
    for(int q=0;q<points;++q){const double* row=rule+4*q;double X[3],measure;
        if(cell){double Jac[3][3];rad_hdiv_cuda::CurvedTetEval(reinterpret_cast<const double(*)[3]>(cellNodes+30*host),row[0],row[1],row[2],X,Jac);
            const double det=rad_hdiv_cuda::det3(Jac);const int s=det>0.?1:det<0.?-1:0;
            if(!isfinite(det)||!s||(sign&&s!=sign)){atomicExch(invalid,1);return;}sign=s;measure=fabs(det);}
        else rad_hdiv_cuda::CurvedTriMapMeasure(reinterpret_cast<const double(*)[3]>(faceNodes+18*host),row[0],row[1],X,measure);
        if(!isfinite(measure)||!(measure>0.)){atomicExch(invalid,1);return;}
        double m[6];for(std::size_t a=0;a<n;++a)m[a]=curvedChargeMonomial(exponents+(cell?3:2)*a,row,cell);
        for(std::size_t a=0;a<n;++a)for(std::size_t b=0;b<n;++b)G[a*n+b]+=row[3]*measure*m[a]*m[b];
    }
    const double* R=(cell?cellMoments:faceMoments)+std::size_t(cls)*n*width;
    for(std::size_t a=0;a<n;++a)x[a]=(cell?double(sign):1.)*R[a*width+column];
    for(std::size_t c=0;c<n;++c){std::size_t p=c;for(std::size_t r=c+1;r<n;++r)if(fabs(G[r*n+c])>fabs(G[p*n+c]))p=r;
        if(!(fabs(G[p*n+c])>0.)){atomicExch(invalid,1);return;}
        if(p!=c){for(std::size_t j=0;j<n;++j){const double v=G[c*n+j];G[c*n+j]=G[p*n+j];G[p*n+j]=v;}const double v=x[c];x[c]=x[p];x[p]=v;}
        for(std::size_t r=0;r<n;++r){if(r==c)continue;const double f=G[r*n+c]/G[c*n+c];if(f==0.)continue;
            for(std::size_t j=c;j<n;++j)G[r*n+j]-=f*G[c*n+j];x[r]-=f*x[c];}}
    double* out=(cell?cellOut:faceOut)+host*n*width;
    for(std::size_t a=0;a<n;++a){const double v=x[a]/G[a*n+a];if(!isfinite(v))atomicExch(invalid,1);out[a*width+column]=v;}
}
// Rows: cell monomials (cell-major), then exterior-face monomials; columns
// ascend per row. Affine: class tables, cell rows scaled by 1/det J, face rows
// by 1/|(F1-F0)x(F2-F0)|. Curved: per-host tables from curvedChargeTables.
__global__ void flatChargeEntries(std::size_t cells,std::size_t faces,std::size_t cm,std::size_t fm,std::size_t cw,std::size_t fw,
    std::size_t cellClasses,std::size_t faceClasses,const double* cellVertices,const double* faceVertices,const int* cellClass,
    const int* cellDofs,const int* faceClass,const int* faceDofs,const double* cellTables,const double* faceTables,bool perHost,
    std::size_t* indices,double* values,int* invalid){
    const auto r=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;const auto cellRows=cells*cm;if(r>=cellRows+faces*fm)return;
    const bool cell=r<cellRows;const auto host=cell?r/cm:(r-cellRows)/fm,a=cell?r%cm:(r-cellRows)%fm;
    const auto width=cell?cw:fw,base=cell?r*cw:cellRows*cw+(r-cellRows)*fw;
    double scale;int cls;
    if(perHost){scale=1.;cls=int(host);}
    else if(cell){const double* V=cellVertices+12*host;double J[3][3];
        for(int c=0;c<3;++c)for(int d=0;d<3;++d)J[c][d]=V[3*d+c]-V[9+c];
        scale=J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
        cls=cellClass[host];if(!isfinite(scale)||scale==0.||cls<0||std::size_t(cls)>=cellClasses){atomicExch(invalid,1);return;}}
    else{const double* F=faceVertices+9*host;const double u[3]={F[3]-F[0],F[4]-F[1],F[5]-F[2]},v[3]={F[6]-F[0],F[7]-F[1],F[8]-F[2]};
        const double n[3]={u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
        scale=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        cls=faceClass[host];if(!isfinite(scale)||!(scale>0.)||cls<0||std::size_t(cls)>=faceClasses){atomicExch(invalid,1);return;}}
    const int* dofs=cell?cellDofs+host*cw:faceDofs+host*fw;
    const double* T=cell?cellTables+(std::size_t(cls)*cm+a)*cw:faceTables+(std::size_t(cls)*fm+a)*fw;
    // Widths are bounded by 64 on the host.
    int order[64];
    for(std::size_t i=0;i<width;++i){std::size_t j=i;while(j>0&&dofs[order[j-1]]>dofs[i]){order[j]=order[j-1];--j;}order[j]=int(i);}
    for(std::size_t k=0;k<width;++k){const int i=order[k];const double value=T[i]/scale;
        if(!isfinite(value))atomicExch(invalid,1);indices[base+k]=std::size_t(dofs[i]);values[base+k]=value;}
}
__global__ void flatChargeOffsets(std::size_t rows,std::size_t cellRows,std::size_t cw,std::size_t fw,std::size_t* offsets){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<=rows)offsets[i]=i<=cellRows?i*cw:cellRows*cw+(i-cellRows)*fw;
}
// Stable radix order keeps rows ascending within each column (deterministic B^T).
__global__ void flatChargeTransposeEntries(std::size_t nnz,std::size_t cellEntries,std::size_t cellRows,std::size_t cw,std::size_t fw,
    const unsigned* entries,const double* values,std::size_t* indices,double* transposed){
    const auto j=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(j>=nnz)return;
    const std::size_t e=entries[j];indices[j]=e<cellEntries?e/cw:cellRows+(e-cellEntries)/fw;transposed[j]=values[e];
}
struct FlatChargeMapBuilder {
    template<class T>using Buffer=radia_hmatrix_cuda::Buffer<T>;
    radia_hmatrix_cuda::HostVector<FlatEvaluationTransfer> transfers;std::size_t launches=0,scratchBytes=0;
    static std::size_t nonzeros(const radia_hdiv_python::ChargeMapInput& in){
        using radia_hmatrix_cuda::checkedProduct;
        return checkedProduct(checkedProduct(in.cellClass.size(),in.cellMonomials),in.cellWidth)+
               checkedProduct(checkedProduct(in.faceClass.size(),in.faceMonomials),in.faceWidth);
    }
    // Every transient buffer build() allocates; the CUB size query allocates nothing.
    static std::size_t scratchBound(const radia_hdiv_python::FlatInput& flat,const radia_hdiv_python::ChargeMapInput& in){
        using radia_hmatrix_cuda::checkedProduct;
        const auto nnz=nonzeros(in);int bits=1;while((std::size_t(1)<<bits)<in.dofs)++bits;std::size_t scratch=0;
        radia_hmatrix_cuda::check(cub::DeviceRadixSort::SortPairs(nullptr,scratch,static_cast<const unsigned*>(nullptr),static_cast<unsigned*>(nullptr),
            static_cast<const unsigned*>(nullptr),static_cast<unsigned*>(nullptr),int(nnz),0,bits));
        // Curved: rules, exponents and the per-host tables (cells x cm x cw plus faces x fm x fw).
        const auto hostTables=in.curved?nonzeros(in):std::size_t(0);
        const auto doubles=checkedProduct(flat.cells.size()+flat.faces.size()+in.cellTables.size()+in.faceTables.size()+
                                          in.cellRule.size()+in.faceRule.size()+hostTables,sizeof(double));
        const auto ints=checkedProduct(in.cellClass.size()+in.cellDofs.size()+in.faceClass.size()+in.faceDofs.size()+1+
                                       (in.curved?in.cellExponents.size()+in.faceExponents.size():0),sizeof(int));
        const auto keys=checkedProduct(4*nnz,sizeof(unsigned));
        if(doubles>SIZE_MAX-ints||doubles+ints>SIZE_MAX-keys||doubles+ints+keys>SIZE_MAX-scratch)throw std::overflow_error("Charge-map scratch overflow");
        return doubles+ints+keys+scratch;
    }
    // Fills B (charges x dofs) and B^T (dofs x charges) in caller-owned device maps.
    void build(const radia_hdiv_python::FlatInput& flat,const radia_hdiv_python::ChargeMapInput& in,
        radia_hmatrix_cuda::DeviceMap& B,radia_hmatrix_cuda::DeviceMap& BT){
        using namespace radia_hmatrix_cuda;
        const auto cells=in.cellClass.size(),faces=in.faceClass.size(),cm=in.cellMonomials,fm=in.faceMonomials,cw=in.cellWidth,fw=in.faceWidth;
        const auto cellRows=cells*cm,rows=cellRows+faces*fm,cellEntries=cellRows*cw,nnz=nonzeros(in);
        scratchBytes=scratchBound(flat,in);
        Buffer<double> cellVertices,faceVertices,cellTables,faceTables;Buffer<int> cellClass,cellDofs,faceClass,faceDofs,invalid;
        auto upload=[&](auto& buffer,const auto& values,const char* payload){buffer.upload(values);
            transfers.push_back({payload,"H2D",checkedProduct(values.size(),sizeof(values[0]))});};
        upload(cellVertices,flat.cells,"charge_map_cell_vertices_m");upload(faceVertices,flat.faces,"charge_map_face_vertices_m");
        upload(cellTables,in.cellTables,"charge_cell_tables");upload(faceTables,in.faceTables,"charge_face_tables");
        upload(cellClass,in.cellClass,"charge_cell_class");upload(cellDofs,in.cellDofs,"charge_cell_dofs");
        upload(faceClass,in.faceClass,"charge_face_class");upload(faceDofs,in.faceDofs,"charge_face_dofs");
        upload(invalid,HostVector<int>{0},"charge_map_invalid_flag");
        const auto cellClasses=in.cellTables.size()/(cm*cw),faceClasses=in.faceTables.size()/(fm*fw);
        Buffer<double> cellRule,faceRule,hostCell,hostFace;Buffer<int> cellExponents,faceExponents;
        if(in.curved){
            upload(cellRule,in.cellRule,"charge_cell_mass_rule");upload(faceRule,in.faceRule,"charge_face_mass_rule");
            upload(cellExponents,in.cellExponents,"charge_cell_exponents");upload(faceExponents,in.faceExponents,"charge_face_exponents");
            hostCell.allocate(cells*cm*cw);hostFace.allocate(faces*fm*fw);const auto columns=cells*cw+faces*fw;
            curvedChargeTables<<<unsigned((columns+127)/128),128>>>(cells,faces,cm,fm,cw,fw,cellClasses,faceClasses,cellVertices.data,faceVertices.data,
                cellClass.data,faceClass.data,cellTables.data,faceTables.data,cellRule.data,int(in.cellRule.size()/4),faceRule.data,int(in.faceRule.size()/4),
                cellExponents.data,faceExponents.data,hostCell.data,hostFace.data,invalid.data);check(cudaGetLastError());++launches;
        }
        B.rows=rows;B.offsets.allocate(rows+1);B.indices.allocate(nnz);B.values.allocate(nnz);
        flatChargeEntries<<<unsigned((rows+127)/128),128>>>(cells,faces,cm,fm,cw,fw,in.curved?cells:cellClasses,in.curved?faces:faceClasses,
            cellVertices.data,faceVertices.data,cellClass.data,cellDofs.data,faceClass.data,faceDofs.data,
            in.curved?hostCell.data:cellTables.data,in.curved?hostFace.data:faceTables.data,in.curved,
            B.indices.data,B.values.data,invalid.data);check(cudaGetLastError());++launches;
        flatChargeOffsets<<<unsigned((rows+256)/256),256>>>(rows,cellRows,cw,fw,B.offsets.data);check(cudaGetLastError());++launches;
        int failed=0;check(cudaMemcpy(&failed,invalid.data,sizeof(int),cudaMemcpyDeviceToHost));
        transfers.push_back({"charge_map_flag","D2H",sizeof(int)});
        if(failed)throw std::runtime_error("Device HDiv charge map has a degenerate or inverted cell or face, invalid class, singular curved Gram or nonfinite value");
        Buffer<unsigned> keys,sortedKeys,entries,sortedEntries;keys.allocate(nnz);sortedKeys.allocate(nnz);entries.allocate(nnz);sortedEntries.allocate(nnz);
        columnKeys<<<unsigned((nnz+255)/256),256>>>(nnz,B.indices.data,keys.data,entries.data);check(cudaGetLastError());++launches;
        int bits=1;while((std::size_t(1)<<bits)<in.dofs)++bits;std::size_t scratch=0;
        check(cub::DeviceRadixSort::SortPairs(nullptr,scratch,keys.data,sortedKeys.data,entries.data,sortedEntries.data,int(nnz),0,bits));
        Buffer<unsigned char> temporary;temporary.allocate(scratch);
        check(cub::DeviceRadixSort::SortPairs(temporary.data,scratch,keys.data,sortedKeys.data,entries.data,sortedEntries.data,int(nnz),0,bits));++launches;
        BT.rows=in.dofs;BT.offsets.allocate(in.dofs+1);BT.indices.allocate(nnz);BT.values.allocate(nnz);
        flatChargeTransposeEntries<<<unsigned((nnz+255)/256),256>>>(nnz,cellEntries,cellRows,cw,fw,sortedEntries.data,B.values.data,
            BT.indices.data,BT.values.data);check(cudaGetLastError());++launches;
        columnOffsets<<<unsigned((in.dofs+256)/256),256>>>(in.dofs,nnz,sortedKeys.data,BT.offsets.data);check(cudaGetLastError());++launches;
        check(cudaDeviceSynchronize());
    }
};
// Host bounds before any upload: the B rows must be the FlatInput charges in
// order, every index is range-checked and every table value is finite. The
// curved map requires the curved Gram geometry (10/6 P2 nodes per host).
inline void validateFlatChargeMap(const radia_hdiv_python::FlatInput& flat,const radia_hdiv_python::ChargeMapInput& in){
    const auto cells=in.cellClass.size(),faces=in.faceClass.size(),cm=in.cellMonomials,fm=in.faceMonomials,cw=in.cellWidth,fw=in.faceWidth;
    // Structural limits only: int32 cells, faces, DOFs and charge rows; the CUB sort counts nonzeros as int.
    constexpr std::size_t top=INT_MAX;
    if(!cells||cells>top||!faces||faces>top||!in.dofs||in.dofs>top||!cm||cm>4||!fm||fm>6||!cw||cw>64||!fw||fw>cw)
        throw std::invalid_argument("Charge tables outside the bounded TET scope");
    if(in.curved!=(flat.family=="curved-tet"))throw std::invalid_argument("Charge map and Gram geometry families disagree");
    const std::size_t cellStride=in.curved?30:12,faceStride=in.curved?18:9;
    if(in.curved){
        auto rule=[](const auto& rows,int dimension){if(rows.empty()||rows.size()%4||rows.size()/4>4096)return false;
            for(std::size_t k=0;k<rows.size();k+=4){if(!std::isfinite(rows[k])||!std::isfinite(rows[k+1])||!std::isfinite(rows[k+2])||!std::isfinite(rows[k+3])||!(rows[k+3]>0))return false;
                if(dimension==2&&rows[k+2]!=0.)return false;}return true;};
        if(!rule(in.cellRule,3)||!rule(in.faceRule,2))throw std::invalid_argument("Curved charge mass rules require finite (xi,weight>0) rows");
    }else if(!in.cellRule.empty()||!in.faceRule.empty())throw std::invalid_argument("Affine charge tables take no curved mass rules");
    if(faces>top/fm||cells>(top-faces*fm)/cm||FlatChargeMapBuilder::nonzeros(in)>top)
        throw std::invalid_argument("Flat charge rows or map nonzeros exceed the checked int range");
    if(in.cellTables.empty()||in.cellTables.size()%(cm*cw)||in.cellTables.size()/(cm*cw)>24||
       in.faceTables.empty()||in.faceTables.size()%(fm*fw)||in.faceTables.size()/(fm*fw)>6)
        throw std::invalid_argument("Flat charge class tables outside the vertex-rank class bound");
    const auto rows=cells*cm+faces*fm;
    if(rows!=flat.host.size()||flat.kind.size()!=rows||flat.expo.size()!=3*rows||flat.cells.size()!=cellStride*cells||flat.faces.size()!=faceStride*faces||
       in.cellDofs.size()!=cells*cw||in.faceDofs.size()!=faces*fw||in.cellExponents.size()!=3*cm||in.faceExponents.size()!=2*fm)
        throw std::invalid_argument("Flat charge table shape mismatch");
    for(std::size_t r=0;r<rows;++r){const bool cell=r<cells*cm;const auto host=cell?r/cm:(r-cells*cm)/fm,a=cell?r%cm:(r-cells*cm)%fm;
        const int* e=cell?&in.cellExponents[3*a]:&in.faceExponents[2*a];
        if(flat.kind[r]!=(cell?0:1)||std::size_t(flat.host[r])!=host||flat.expo[3*r]!=e[0]||flat.expo[3*r+1]!=e[1]||flat.expo[3*r+2]!=(cell?e[2]:0))
            throw std::invalid_argument("Charge map rows disagree with the Gram charge order");}
    for(double v:in.cellTables)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite charge cell table");
    for(double v:in.faceTables)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite charge face table");
    for(int c:in.cellClass)if(c<0||std::size_t(c)>=in.cellTables.size()/(cm*cw))throw std::invalid_argument("Charge cell class out of range");
    for(int c:in.faceClass)if(c<0||std::size_t(c)>=in.faceTables.size()/(fm*fw))throw std::invalid_argument("Charge face class out of range");
    auto dofRows=[&](const auto& dofs,std::size_t count,std::size_t width){
        for(int d:dofs)if(d<0||std::size_t(d)>=in.dofs)throw std::invalid_argument("Charge DOF index out of range");
        for(std::size_t e=0;e<count;++e)for(std::size_t i=0;i<width;++i)for(std::size_t j=0;j<i;++j)
            if(dofs[e*width+i]==dofs[e*width+j])throw std::invalid_argument("Duplicate charge local DOF");};
    dofRows(in.cellDofs,cells,cw);dofRows(in.faceDofs,faces,fw);
}
}
