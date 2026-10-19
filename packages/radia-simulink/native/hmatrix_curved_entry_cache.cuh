// Explicit indexed-entry source-potential reuse. No dense Gram allocation.
#pragma once
#include <map>
namespace radia_hmatrix_cuda {
__global__ void fillCachedDenseLeaf(const CurvedCharge* charges,const HostTargetPoint* points,
    std::size_t np,const HostTargetSpan* spans,const double* weights,std::size_t stride,const double* cache,
    const std::size_t* permutation,const double* scale,DeferredDenseSpan leaf,double* coefficients,int* invalid){
    const auto k=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(k>=leaf.rows*leaf.columns)return;
    const auto originalRow=permutation[leaf.row+k/leaf.columns],originalColumn=permutation[leaf.column+k%leaf.columns];
    auto row=originalRow,column=originalColumn;if(row>column){auto t=row;row=column;column=t;}
    double value=curvedCachedDot(charges,row,column,points,np,spans,cache,weights,stride);
    if(row!=column)value=.5*(value+curvedCachedDot(charges,column,row,points,np,spans,cache,weights,stride));
    value*=(1./scale[originalRow])*(1./scale[originalColumn]);
    if(!isfinite(value))atomicExch(invalid,1);coefficients[leaf.offset+k]=value;
}
__global__ void evaluateCachedCurvedOrdered(const CurvedCharge* charges,const HostTargetPoint* points,
    std::size_t np,const HostTargetSpan* spans,const double* weights,std::size_t stride,
    const double* cache,std::size_t sourceFirst,std::size_t sourceCount,
    const AnalyticEntryIndex* indices,std::size_t count,double* values,int* invalid){
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=count)return;
    const auto index=indices[i];
    if(index.column<sourceFirst||index.column-sourceFirst>=sourceCount){atomicExch(invalid,1);values[i]=0.;return;}
    const auto value=curvedCachedDot(charges,index.row,index.column-sourceFirst,points,np,spans,cache,weights,stride);
    if(!isfinite(value))atomicExch(invalid,1);values[i]=value;
}
class CurvedCachedEntryOracle {
    Buffer<CurvedCharge> charges;
    Buffer<HostTargetPoint> points;
    Buffer<HostTargetSpan> spans;
    Buffer<MappedHostPoint> mapped;
    Buffer<MomentSourceGroup> momentGroups;
    Buffer<std::size_t> groupIds;
    Buffer<double> weights,potentials,gl,gw,values;
    Buffer<AnalyticEntryIndex> indices;
    Buffer<int> invalid;
    HostVector<std::size_t> tileOffsets,tileTags,tileStamps;
    HostVector<AnalyticEntryIndex> fullIndices{HostFaultAllocator<AnalyticEntryIndex>("cache.full-indices")};
    HostVector<double> fullTerms{HostFaultAllocator<double>("cache.full-terms")};
    std::size_t capacity,n,np,stride,tileSources,slotElements,stamp=0;
    int device;
    std::mutex access;
    std::size_t loadTile(std::size_t tile,cudaStream_t stream){
        for(std::size_t slot=0;slot<tileTags.size();++slot)if(tileTags[slot]==tile){
            tileStamps[slot]=++stamp;++tileHits;return slot;
        }
        const auto slot=std::size_t(std::min_element(tileStamps.begin(),tileStamps.end())-tileStamps.begin());
        const auto sourceFirst=checkedProduct(tile,tileSources),sourcesHere=std::min(tileSources,n-sourceFirst);
        for(int kind=0;kind<2;++kind){const auto offset=2*tile+std::size_t(kind);
            const auto firstGroup=tileOffsets[offset],groupCount=tileOffsets[offset+1]-firstGroup;
            const auto work=checkedProduct(groupCount,np);
            constexpr std::size_t batchSize=std::size_t(1)<<22;
            for(std::size_t first=0;first<work;first+=batchSize){const auto count=std::min(batchSize,work-first);
                if(kind==0)cacheCurvedHostMomentPotentials<0><<<unsigned((count+127)/128),128,0,stream>>>(charges.data,n,points.data,np,spans.data,mapped.data,weights.data,stride,momentGroups.data,groupIds.data+firstGroup,groupCount,gl.data,gw.data,int(gl.size),first,count,sourceFirst,sourcesHere,potentials.data+slot*slotElements,invalid.data);
                else cacheCurvedHostMomentPotentials<1><<<unsigned((count+127)/128),128,0,stream>>>(charges.data,n,points.data,np,spans.data,mapped.data,weights.data,stride,momentGroups.data,groupIds.data+firstGroup,groupCount,gl.data,gw.data,int(gl.size),first,count,sourceFirst,sourcesHere,potentials.data+slot*slotElements,invalid.data);
                check(cudaGetLastError());
            }
        }
        tileTags[slot]=tile;tileStamps[slot]=++stamp;++tileLoads;return slot;
    }
public:
    std::size_t ownedBytes=0,tileLoads=0,tileHits=0,fullTileBatches=0,hostScratchBytes=0;
    CurvedCachedEntryOracle(const CurvedBlock& block,std::size_t batchCapacity,std::size_t cacheBudget=std::size_t(256)<<20,bool requireFull=false)
        :capacity(batchCapacity),n(block.charges.size()),np(0),stride(std::max(block.tet.size(),block.tri.size())){
        if(!n||n>8192||!capacity||capacity>EntryBatchCapacity||!stride||block.gl.empty()||block.gl.size()!=block.gw.size()||!cacheBudget||cacheBudget>(std::size_t(2)<<30))throw std::invalid_argument("Invalid cached curved oracle budget");
        check(cudaGetDevice(&device));
        HostVector<HostTargetPoint> hostPoints{HostFaultAllocator<HostTargetPoint>("cache.target-points")};
        HostVector<HostTargetSpan> hostSpans{HostFaultAllocator<HostTargetSpan>("cache.target-spans")};
        curvedHostTargets(block,hostPoints,hostSpans);np=hostPoints.size();
        const auto groups=curvedSourceMomentGroups(block,HostFaultAllocator<MomentSourceGroup>("cache.source-groups"));
        const auto sourceBytes=checkedProduct(np,sizeof(double));
        if(!sourceBytes||sourceBytes>cacheBudget)throw std::invalid_argument("One curved source cache column exceeds budget");
        const auto budgetSources=cacheBudget/sourceBytes;
        if(requireFull&&n>budgetSources)throw std::invalid_argument("Device dense fill requires a full source cache");
        if(n<=budgetSources)tileSources=n;
        else {if(budgetSources<2)throw std::invalid_argument("Two curved source cache tiles exceed budget");tileSources=budgetSources/2;}
        const auto tiles=1+(n-1)/tileSources,slots=std::min(std::size_t(2),tiles);
        slotElements=checkedProduct(tileSources,np);
        tileTags.assign(slots,std::numeric_limits<std::size_t>::max());tileStamps.assign(slots,0);
        if(tileSources==n){
            fullIndices.reserve(2*capacity);fullTerms.resize(2*capacity);
            hostScratchBytes=checkedProduct(fullIndices.capacity(),sizeof(AnalyticEntryIndex))+checkedProduct(fullTerms.capacity(),sizeof(double));
        }
        HostVector<std::size_t> ids{HostFaultAllocator<std::size_t>("cache.group-ids")};
        for(std::size_t tile=0;tile<tiles;++tile){const auto first=tile*tileSources,last=std::min(n,first+tileSources);
            for(int kind=0;kind<2;++kind){tileOffsets.push_back(ids.size());
                for(std::size_t g=0;g<groups.size();++g){if(block.charges[groups[g].member[0]].kind!=kind)continue;
                    bool active=false;for(std::size_t m=0;m<groups[g].count;++m)active|=groups[g].member[m]>=first&&groups[g].member[m]<last;
                    if(active)ids.push_back(g);
                }
            }
        }
        tileOffsets.push_back(ids.size());
        auto add=[&](std::size_t count,std::size_t width){const auto bytes=checkedProduct(count,width);
            if(bytes>std::numeric_limits<std::size_t>::max()-ownedBytes)throw std::overflow_error("Cached curved oracle allocation overflow");ownedBytes+=bytes;};
        add(n,sizeof(CurvedCharge));add(np,sizeof(HostTargetPoint)+sizeof(MappedHostPoint));add(hostSpans.size(),sizeof(HostTargetSpan));
        add(checkedProduct(n,stride),sizeof(double));add(groups.size(),sizeof(MomentSourceGroup));add(ids.size(),sizeof(std::size_t));
        add(checkedProduct(slots,slotElements),sizeof(double));add(block.gl.size()+block.gw.size(),sizeof(double));
        add(2*capacity,sizeof(AnalyticEntryIndex)+sizeof(double));add(1,sizeof(int));
        if(ownedBytes>(std::size_t(2)<<30))throw std::invalid_argument("Cached curved oracle workspace exceeds 2 GiB");
        std::size_t free,total;check(cudaMemGetInfo(&free,&total));const auto reserve=std::max(std::size_t(2)<<30,total/5);
        if(free<=reserve||ownedBytes>free-reserve)throw std::runtime_error("Cached curved oracle GPU reserve exceeded");
        charges.upload(block.charges);points.upload(hostPoints);spans.upload(hostSpans);momentGroups.upload(groups);groupIds.upload(ids);
        gl.upload(block.gl);gw.upload(block.gw);mapped.allocate(np);weights.allocate(checkedProduct(n,stride));
        potentials.allocate(checkedProduct(slots,slotElements));indices.allocate(2*capacity);values.allocate(2*capacity);invalid.allocate(1);
        mapCurvedHostPoints<<<unsigned((np+127)/128),128>>>(charges.data,points.data,np,mapped.data);check(cudaGetLastError());
        mapCurvedChargeWeights<<<unsigned((n*stride+127)/128),128>>>(charges.data,n,points.data,spans.data,stride,weights.data);check(cudaGetLastError());
        check(cudaStreamSynchronize(nullptr));
    }
    void fillDense(const Plan& plan,double* coefficients,const std::size_t* permutation,const double* scale,cudaStream_t stream){
        std::lock_guard<std::mutex> lock(access);int current;check(cudaGetDevice(&current));
        if(current!=device||tileSources!=n||plan.n!=n||plan.deferredDense.empty())throw std::invalid_argument("Invalid device dense fill context");
        check(cudaMemsetAsync(invalid.data,0,sizeof(int),stream));const auto slot=loadTile(0,stream);
        for(const auto& leaf:plan.deferredDense){
            const auto entries=checkedProduct(leaf.rows,leaf.columns);
            if(!entries||leaf.row>=n||leaf.column>=n||leaf.rows>n-leaf.row||leaf.columns>n-leaf.column||
                leaf.offset>plan.coefficientCount||entries>plan.coefficientCount-leaf.offset)throw std::invalid_argument("Invalid deferred coefficient range");
            fillCachedDenseLeaf<<<unsigned((entries+127)/128),128,0,stream>>>(charges.data,points.data,np,spans.data,weights.data,stride,potentials.data+slot*slotElements,permutation,scale,leaf,coefficients,invalid.data);
            check(cudaGetLastError());
        }
        int flag=0;check(cudaMemcpyAsync(&flag,invalid.data,sizeof(int),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
        if(flag)throw std::runtime_error("Device dense fill failed");
    }
    template<class Allocator>HostVector<double> evaluate(const std::vector<AnalyticEntryIndex,Allocator>& request,cudaStream_t stream){
        std::lock_guard<std::mutex> lock(access);int current;check(cudaGetDevice(&current));
        if(current!=device)throw std::runtime_error("Cached curved oracle belongs to another GPU");
        if(request.size()>capacity)throw std::invalid_argument("Cached curved entry batch exceeds workspace");
        if(tileSources==n){
            fullIndices.clear();
            for(const auto& index:request){auto row=index.row,col=index.column;
                if(row>=n||col>=n)throw std::invalid_argument("Cached curved entry index outside charge range");
                if(row>col)std::swap(row,col);
                fullIndices.push_back({row,col});if(row!=col)fullIndices.push_back({col,row});
            }
            HostVector<double> result(request.size(),0.,HostFaultAllocator<double>("cache.result"));if(request.empty())return result;
            int flag=0;check(cudaMemsetAsync(invalid.data,0,sizeof(int),stream));
            const auto slot=loadTile(0,stream);
            check(cudaMemcpyAsync(indices.data,fullIndices.data(),fullIndices.size()*sizeof(AnalyticEntryIndex),cudaMemcpyHostToDevice,stream));
            evaluateCachedCurvedOrdered<<<unsigned((fullIndices.size()+127)/128),128,0,stream>>>(charges.data,points.data,np,spans.data,weights.data,stride,potentials.data+slot*slotElements,0,n,indices.data,fullIndices.size(),values.data,invalid.data);
            check(cudaGetLastError());
            check(cudaMemcpyAsync(fullTerms.data(),values.data,fullIndices.size()*sizeof(double),cudaMemcpyDeviceToHost,stream));
            check(cudaMemcpyAsync(&flag,invalid.data,sizeof(int),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
            if(flag)throw std::runtime_error("Cached curved entry evaluation failed");
            std::size_t cursor=0;
            for(std::size_t i=0;i<request.size();++i){
                result[i]=request[i].row==request[i].column?fullTerms[cursor]:.5*(fullTerms[cursor]+fullTerms[cursor+1]);
                cursor+=request[i].row==request[i].column?1:2;
                if(!std::isfinite(result[i]))throw std::runtime_error("Cached curved symmetric entry overflow");
            }
            ++fullTileBatches;return result;
        }
        HostVector<AnalyticEntryIndex> ordered;HostVector<std::size_t> first(request.size()),second(request.size(),std::numeric_limits<std::size_t>::max());
        HostMap<std::size_t,HostVector<std::size_t>> byTile;
        for(std::size_t i=0;i<request.size();++i){auto row=request[i].row,col=request[i].column;
            if(row>=n||col>=n)throw std::invalid_argument("Cached curved entry index outside charge range");
            if(row>col)std::swap(row,col);
            first[i]=ordered.size();byTile[col/tileSources].push_back(ordered.size());ordered.push_back({row,col});
            if(row!=col){second[i]=ordered.size();byTile[row/tileSources].push_back(ordered.size());ordered.push_back({col,row});}
        }
        HostVector<double> result(request.size(),0.,HostFaultAllocator<double>("cache.result")),terms(ordered.size());if(request.empty())return result;
        int flag=0;check(cudaMemsetAsync(invalid.data,0,sizeof(int),stream));
        for(const auto& tile:byTile){const auto slot=loadTile(tile.first,stream),sourceFirst=tile.first*tileSources;
            HostVector<AnalyticEntryIndex> batch;for(auto index:tile.second)batch.push_back(ordered[index]);
            HostVector<double> output(batch.size());
            check(cudaMemcpyAsync(indices.data,batch.data(),batch.size()*sizeof(AnalyticEntryIndex),cudaMemcpyHostToDevice,stream));
            evaluateCachedCurvedOrdered<<<unsigned((batch.size()+127)/128),128,0,stream>>>(charges.data,points.data,np,spans.data,weights.data,stride,potentials.data+slot*slotElements,sourceFirst,std::min(tileSources,n-sourceFirst),indices.data,batch.size(),values.data,invalid.data);check(cudaGetLastError());
            check(cudaMemcpyAsync(output.data(),values.data,output.size()*sizeof(double),cudaMemcpyDeviceToHost,stream));
            // The final tile already waits for its output. Recover the shared
            // sticky error flag in that same ordered wait, before returning
            // any entries to ACA. Earlier tile errors remain set.
            if(tile.first==byTile.rbegin()->first)check(cudaMemcpyAsync(&flag,invalid.data,sizeof(int),cudaMemcpyDeviceToHost,stream));
            check(cudaStreamSynchronize(stream));
            for(std::size_t i=0;i<output.size();++i)terms[tile.second[i]]=output[i];
        }
        if(flag)throw std::runtime_error("Cached curved entry evaluation failed");
        for(std::size_t i=0;i<request.size();++i){result[i]=second[i]==std::numeric_limits<std::size_t>::max()?terms[first[i]]:.5*(terms[first[i]]+terms[second[i]]);
            if(!std::isfinite(result[i]))throw std::runtime_error("Cached curved symmetric entry overflow");}
        return result;
    }
};
}
