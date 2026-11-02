// Device construction of the affine-TET HDiv charge map B and B^T.
// NGSolve tabulates charge-projection tables per vertex-rank class on the
// host; only those tables, the Gram vertices, DOF numbers and classes cross.
#pragma once
#include <cub/device/device_radix_sort.cuh>
#include "hdiv_resident_python_api.hpp"
#include "hdiv_flat_evaluation.cuh"
namespace radia_hdiv_bh {
// Rows: cell monomials (cell-major), then exterior-face monomials. Cell rows
// scale by 1/det J, face rows by 1/|(F1-F0)x(F2-F0)|; columns ascend per row.
__global__ void flatChargeEntries(std::size_t cells,std::size_t faces,std::size_t cm,std::size_t fm,std::size_t cw,std::size_t fw,
    std::size_t cellClasses,std::size_t faceClasses,const double* cellVertices,const double* faceVertices,const int* cellClass,
    const int* cellDofs,const int* faceClass,const int* faceDofs,const double* cellTables,const double* faceTables,
    std::size_t* indices,double* values,int* invalid){
    const auto r=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;const auto cellRows=cells*cm;if(r>=cellRows+faces*fm)return;
    const bool cell=r<cellRows;const auto host=cell?r/cm:(r-cellRows)/fm,a=cell?r%cm:(r-cellRows)%fm;
    const auto width=cell?cw:fw,base=cell?r*cw:cellRows*cw+(r-cellRows)*fw;
    double scale;int cls;
    if(cell){const double* V=cellVertices+12*host;double J[3][3];
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
    // Fills B (charges x dofs) and B^T (dofs x charges) in caller-owned device maps.
    void build(const radia_hdiv_python::FlatInput& flat,const radia_hdiv_python::ChargeMapInput& in,
        radia_hmatrix_cuda::DeviceMap& B,radia_hmatrix_cuda::DeviceMap& BT){
        using namespace radia_hmatrix_cuda;
        const auto cells=in.cellClass.size(),faces=in.faceClass.size(),cm=in.cellMonomials,fm=in.faceMonomials,cw=in.cellWidth,fw=in.faceWidth;
        const auto cellRows=cells*cm,rows=cellRows+faces*fm,cellEntries=cellRows*cw,nnz=nonzeros(in);
        Buffer<double> cellVertices,faceVertices,cellTables,faceTables;Buffer<int> cellClass,cellDofs,faceClass,faceDofs,invalid;
        auto upload=[&](auto& buffer,const auto& values,const char* payload){buffer.upload(values);
            transfers.push_back({payload,"H2D",checkedProduct(values.size(),sizeof(values[0]))});};
        upload(cellVertices,flat.cells,"charge_map_cell_vertices_m");upload(faceVertices,flat.faces,"charge_map_face_vertices_m");
        upload(cellTables,in.cellTables,"charge_cell_tables");upload(faceTables,in.faceTables,"charge_face_tables");
        upload(cellClass,in.cellClass,"charge_cell_class");upload(cellDofs,in.cellDofs,"charge_cell_dofs");
        upload(faceClass,in.faceClass,"charge_face_class");upload(faceDofs,in.faceDofs,"charge_face_dofs");
        upload(invalid,HostVector<int>{0},"charge_map_invalid_flag");
        B.rows=rows;B.offsets.allocate(rows+1);B.indices.allocate(nnz);B.values.allocate(nnz);
        flatChargeEntries<<<unsigned((rows+127)/128),128>>>(cells,faces,cm,fm,cw,fw,in.cellTables.size()/(cm*cw),in.faceTables.size()/(fm*fw),
            cellVertices.data,faceVertices.data,cellClass.data,cellDofs.data,faceClass.data,faceDofs.data,cellTables.data,faceTables.data,
            B.indices.data,B.values.data,invalid.data);check(cudaGetLastError());++launches;
        flatChargeOffsets<<<unsigned((rows+256)/256),256>>>(rows,cellRows,cw,fw,B.offsets.data);check(cudaGetLastError());++launches;
        int failed=0;check(cudaMemcpy(&failed,invalid.data,sizeof(int),cudaMemcpyDeviceToHost));
        transfers.push_back({"charge_map_flag","D2H",sizeof(int)});
        if(failed)throw std::runtime_error("Device HDiv charge map has a degenerate cell or face, invalid class or nonfinite value");
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
        scratchBytes=checkedProduct(flat.cells.size()+flat.faces.size()+in.cellTables.size()+in.faceTables.size(),sizeof(double))+
            checkedProduct(in.cellClass.size()+in.cellDofs.size()+in.faceClass.size()+in.faceDofs.size()+1,sizeof(int))+
            checkedProduct(4*nnz,sizeof(unsigned))+scratch;
    }
};
// Host bounds before any upload: the B rows must be the FlatInput charges in
// order, every index is range-checked and every table value is finite.
inline void validateFlatChargeMap(const radia_hdiv_python::FlatInput& flat,const radia_hdiv_python::ChargeMapInput& in){
    const auto cells=in.cellClass.size(),faces=in.faceClass.size(),cm=in.cellMonomials,fm=in.faceMonomials,cw=in.cellWidth,fw=in.faceWidth;
    if(!cells||cells>4096||!faces||!in.dofs||in.dofs>16384||!cm||cm>4||!fm||fm>6||!cw||cw>64||!fw||fw>cw)
        throw std::invalid_argument("Flat charge tables outside the bounded affine-TET scope");
    if(in.cellTables.empty()||in.cellTables.size()%(cm*cw)||in.cellTables.size()/(cm*cw)>24||
       in.faceTables.empty()||in.faceTables.size()%(fm*fw)||in.faceTables.size()/(fm*fw)>6)
        throw std::invalid_argument("Flat charge class tables outside the vertex-rank class bound");
    const auto rows=cells*cm+faces*fm;
    if(rows!=flat.host.size()||flat.kind.size()!=rows||flat.expo.size()!=3*rows||flat.cells.size()!=12*cells||flat.faces.size()!=9*faces||
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
