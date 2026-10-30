// Bounded analytic constant-charge block assembly. No centroid/cloud substitute.
#pragma once
#include "radia_field_exact.cuh"
#include "rad_charge_quadrature.h"
#include "hmatrix_host_fault.hpp"
namespace radia_hmatrix_cuda {
struct AnalyticSource {double vertices[12];std::size_t first,count;bool volume;};
struct QuadraturePoint {double point[3],weight;};
struct AnalyticBlock {HostVector<AnalyticSource> sources{HostFaultAllocator<AnalyticSource>("geometry.analytic-sources")};HostVector<QuadraturePoint> quadrature{HostFaultAllocator<QuadraturePoint>("geometry.analytic-quadrature")};};
template<class Allocator>
inline AnalyticBlock analyticSources(const std::vector<double,Allocator>& cells,const std::vector<double,Allocator>& faces,
                                    std::size_t sourceLimit=128){
    if(cells.size()%12||faces.size()%9||cells.size()/12+faces.size()/9==0||cells.size()/12+faces.size()/9>sourceLimit)
        throw std::invalid_argument("Analytic constant charge source count outside declared budget");
    AnalyticBlock block;HostVector<double> points{HostFaultAllocator<double>("geometry.duffy-points")},weights{HostFaultAllocator<double>("geometry.duffy-weights")};rad_gl4_duffy_tet(points,weights);
    for(int volume=1;volume>=0;--volume){
        const auto& geometry=volume?cells:faces;const std::size_t width=volume?12:9;
        for(std::size_t offset=0;offset<geometry.size();offset+=width){
            AnalyticSource source{};source.volume=volume!=0;source.first=block.quadrature.size();source.count=volume?64:7;
            std::copy(geometry.begin()+offset,geometry.begin()+offset+width,source.vertices);
            auto* v=source.vertices;double e1[3],e2[3],e3[3],cross[3];
            for(int k=0;k<3;++k){e1[k]=v[3+k]-v[k];e2[k]=v[6+k]-v[k];e3[k]=volume?v[9+k]-v[k]:0.;}
            cross[0]=e2[1]*e3[2]-e2[2]*e3[1];cross[1]=e2[2]*e3[0]-e2[0]*e3[2];cross[2]=e2[0]*e3[1]-e2[1]*e3[0];
            double measure;
            if(volume){const double volumeMeasure=std::fabs(e1[0]*cross[0]+e1[1]*cross[1]+e1[2]*cross[2])/6.;measure=6.*volumeMeasure;}
            else{cross[0]=e1[1]*e2[2]-e1[2]*e2[1];cross[1]=e1[2]*e2[0]-e1[0]*e2[2];cross[2]=e1[0]*e2[1]-e1[1]*e2[0];measure=.5*std::sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);}
            if(!std::isfinite(measure)||measure<=0)throw std::invalid_argument("Degenerate analytic charge geometry");
            for(std::size_t q=0;q<source.count;++q){
                QuadraturePoint point{};
                if(volume){
                    for(int k=0;k<3;++k)point.point[k]=v[k]+points[3*q]*(v[3+k]-v[k])+points[3*q+1]*(v[6+k]-v[k])+points[3*q+2]*(v[9+k]-v[k]);
                    point.weight=weights[q]*measure;
                }else{
                    for(int i=0;i<3;++i)for(int k=0;k<3;++k)point.point[k]+=RAD_DUN5[q][i]*v[3*i+k];
                    point.weight=RAD_DUN5[q][3]*measure;
                }
                block.quadrature.push_back(point);
            }
            block.sources.push_back(source);
        }
    }
    return block;
}
__device__ double analyticQuadDot(const AnalyticSource* sources,const QuadraturePoint* quadrature,
                                 std::size_t target,std::size_t source){
    const auto src=sources[source];const auto tgt=sources[target];double sum=0;
    for(std::size_t q=0;q<tgt.count;++q){
        const auto point=quadrature[tgt.first+q];double potential;
        if(src.volume)potential=rad_hdiv_cuda::PhiTet(reinterpret_cast<const double(*)[3]>(src.vertices),point.point);
        else potential=rad_hdiv_cuda::TriPotential(reinterpret_cast<const double(*)[3]>(src.vertices),point.point);
        sum+=point.weight*potential;
    }
    return sum*0.07957747154594766788;
}
__device__ double analyticSymmetricEntry(const AnalyticSource* sources,const QuadraturePoint* quadrature,
                                       std::size_t row,std::size_t column){
    // Canonical ordering preserves the dense CPU-compatible entry arithmetic,
    // including when ACA requests the same entry through its column direction.
    if(row>column){const auto swap=row;row=column;column=swap;}
    const double value=analyticQuadDot(sources,quadrature,row,column);
    return row==column?value:.5*(value+analyticQuadDot(sources,quadrature,column,row));
}
struct AnalyticEntryIndex {std::size_t row,column;};
__global__ void evaluateAnalyticEntries(const AnalyticSource* sources,const QuadraturePoint* quadrature,
        std::size_t sourceCount,const AnalyticEntryIndex* indices,std::size_t count,double* values,int* invalid){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=count)return;
    const auto index=indices[i];
    if(index.row>=sourceCount||index.column>=sourceCount){atomicExch(invalid,1);values[i]=0.;return;}
    const auto value=analyticSymmetricEntry(sources,quadrature,index.row,index.column);
    if(!isfinite(value))atomicExch(invalid,1);values[i]=value;
}
__global__ void assembleAnalyticConstant(const AnalyticSource* sources,const QuadraturePoint* quadrature,
        std::size_t n,double* matrix,int* invalid){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n*n)return;
    const auto row=i/n,col=i%n;if(row>col)return;double value=analyticSymmetricEntry(sources,quadrature,row,col);
    if(!isfinite(value))atomicExch(invalid,1);matrix[i]=value;if(row!=col)matrix[col*n+row]=value;
}
__global__ void analyticPhysicalScale(std::size_t n,const double* matrix,double* scale,int* invalid){
    auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n)return;
    double diagonal=matrix[i*n+i];if(!isfinite(diagonal)||diagonal<=0){atomicExch(invalid,1);return;}scale[i]=sqrt(diagonal);
}
__global__ void normalizeAnalytic(std::size_t n,double* matrix,const double* scale){
    auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n*n)matrix[i]*=(1./scale[i/n])*(1./scale[i%n]);
}
} // namespace radia_hmatrix_cuda
