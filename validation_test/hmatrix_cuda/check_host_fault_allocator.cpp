/* Host container unwind contract; not an independent numerical oracle. */
#include "../../src/core/gpu/hdiv/hmatrix_host_fault.hpp"
#include "../../src/core/gpu/hdiv/hdiv_inverse_bh.hpp"
#include "../../src/core/gpu/hdiv/hdiv_resident_python_api.hpp"
#include <cassert>
#include <cmath>
#include "rad_charge_quadrature.h"
using namespace radia_hmatrix_cuda;
struct OwnedTest {
    HostVector<int> values;
    explicit OwnedTest(bool fail):values(8,3,HostFaultAllocator<int>("test.object-values")){
        if(fail)throw std::runtime_error("constructor failure");
    }
};
static void workload(){
    HostVector<double> qp{HostFaultAllocator<double>("test.quadrature-points")},qw{HostFaultAllocator<double>("test.quadrature-weights")};
    rad_gl4_duffy_tet(qp,qw);assert(qp.size()==192 && qw.size()==64);
    double sum=0.;for(double weight:qw)sum+=weight;assert(std::abs(sum-1./6.)<1.e-14);
    auto object=makeHostUnique<OwnedTest>("test.object",false);
    assert(object->values[7]==3);
    HostVector<double> a(17,2.,HostFaultAllocator<double>("test.values"));
    HostVector<bool> flags(100,false,HostFaultAllocator<bool>("test.bits"));
    HostVector<HostVector<int>> rows(3,HostVector<int>{HostFaultAllocator<int>("test.row")},HostFaultAllocator<HostVector<int>>("test.rows"));
    a.reserve(100);flags[99]=true;
    for(auto& r:rows){for(int i=0;i<12;++i)r.push_back(i);}
    auto copy=a;assert(copy[16]==2. && flags[99]);
    HostMap<int,HostVector<int>> tiles;
    for(int i=0;i<4;++i){tiles[i].push_back(i);tiles[i].reserve(20);}
    assert(tiles.at(3)[0]==3);
}
static void materialBoundaryWorkload(){
    radia_hdiv_bh::Law law{{0.,1.},{0.,0.,1.,0.},{0.,.5},1.,1.,1.,1.};
    radia_hdiv_bh::MaterialInput input{1,1,{0,1,1,1},{0},{1.},{1.}};
    auto copiedLaw=law;auto copiedInput=input;
    radia_hdiv_bh::Samples scalar;scalar.secant.resize(1);scalar.tangent.resize(1);scalar.energy.resize(1);
    radia_hdiv_bh::MaterialEvaluation output;output.gradient.resize(1);output.tangentAction.resize(1);
    assert(copiedLaw.knots[1]==1. && copiedInput.values[0]==1.);
}
static void newtonBoundaryWorkload(){
    radia_hdiv_python::NewtonResult result;result.x.resize(4);result.history.push_back({});
    auto saved=result;result.history.reserve(40);result.x.assign(saved.x.begin(),saved.x.end());
    assert(result.history.size()==1 && saved.x.size()==4);
}
static void adapterBoundaryWorkload(){
    radia_hdiv_python::LeafInput leaf{0,0,1,1,0,{2.},{}};
    radia_hdiv_python::MapInput map{1,1,{0,1},{0},{1.}};
    radia_hdiv_python::MaterialResult result{{3.},1,0.,8,0};
    HostVector<radia_hdiv_python::LeafInput> leaves;leaves.push_back(leaf);
    auto saved=map;HostVector<double> state{.2},direction{.1};
    assert(leaves[0].a1[0]==2. && saved.values[0]==1. && result.x[0]==3.);
}
static void sharedOwnerWorkload(){
    auto owner=std::allocate_shared<OwnedTest>(HostFaultAllocator<OwnedTest>("test.shared-owner"),false);
    std::weak_ptr<OwnedTest> weak=owner;
    assert(owner.use_count()==1);owner.reset();assert(weak.expired());weak.reset();
}
int main(){
    configureHostFault(0);sharedOwnerWorkload();const auto sharedCount=hostFaultState.attempts;
    assert(sharedCount==2 && hostFaultState.liveBytes==0);
    for(std::size_t at=1;at<=sharedCount;++at)for(int repeat=0;repeat<10;++repeat){
        configureHostFault(at);bool caught=false;
        try{sharedOwnerWorkload();}catch(const std::bad_alloc&){caught=true;}
        assert(caught && hostFaultState.attempts==at && hostFaultState.liveBytes==0);
    }

    configureHostFault(0);adapterBoundaryWorkload();const auto adapterCount=hostFaultState.attempts;
    assert(adapterCount>=10 && hostFaultState.liveBytes==0);
    for(std::size_t at=1;at<=adapterCount;++at)for(int repeat=0;repeat<10;++repeat){
        configureHostFault(at);bool caught=false;
        try{adapterBoundaryWorkload();}catch(const std::bad_alloc&){caught=true;}
        assert(caught && hostFaultState.attempts==at && hostFaultState.liveBytes==0);
    }

    configureHostFault(0);newtonBoundaryWorkload();const auto newtonCount=hostFaultState.attempts;
    assert(newtonCount>=5 && hostFaultState.liveBytes==0);
    for(std::size_t at=1;at<=newtonCount;++at)for(int repeat=0;repeat<10;++repeat){
        configureHostFault(at);bool caught=false;
        try{newtonBoundaryWorkload();}catch(const std::bad_alloc&){caught=true;}
        assert(caught && hostFaultState.attempts==at && hostFaultState.liveBytes==0);
    }

    configureHostFault(0);materialBoundaryWorkload();const auto boundaryCount=hostFaultState.attempts;
    assert(boundaryCount==19 && hostFaultState.liveBytes==0);
    for(std::size_t at=1;at<=boundaryCount;++at)for(int repeat=0;repeat<10;++repeat){
        configureHostFault(at);bool caught=false;
        try{materialBoundaryWorkload();}catch(const std::bad_alloc&){caught=true;}
        assert(caught && hostFaultState.attempts==at && hostFaultState.liveBytes==0);
    }

    auto registryWorkload=[](){
        using Map=HostMap<int,int>;
        auto registry=makeHostUnique<Map>("test.registry-object",std::less<int>{},Map::allocator_type("test.registry"));
        registry->emplace(1,2);
    };
    configureHostFault(0);registryWorkload();const auto registryCount=hostFaultState.attempts;
    assert(registryCount>=2 && hostFaultState.liveBytes==0);
    for(std::size_t at=1;at<=registryCount;++at){
        configureHostFault(at);bool caught=false;
        try{registryWorkload();}catch(const std::bad_alloc&){caught=true;}
        assert(caught && hostFaultState.attempts==at && hostFaultState.liveBytes==0);
    }
    configureHostFault(0);workload();
    const auto count=hostFaultState.attempts;
    assert(count>10 && hostFaultState.liveBytes==0);
    for(std::size_t at=1;at<=count;++at){
        for(int repeat=0;repeat<10;++repeat){
            configureHostFault(at);bool caught=false;
            try{workload();}catch(const std::bad_alloc&){caught=true;}
            assert(caught && hostFaultState.attempts==at);
            assert(hostFaultState.liveBytes==0);
        }
    }
    configureHostFault(0);workload();assert(hostFaultState.liveBytes==0);
    configureHostFaultSite(1,"test.selected");
    bool selectedCaught=false;
    try{
        HostVector<int> other(8,0,HostFaultAllocator<int>("test.other"));
        assert(hostFaultState.attempts==0 && hostFaultState.liveBytes==8*sizeof(int));
        HostVector<int> selected(4,0,HostFaultAllocator<int>("test.selected"));
    }catch(const std::bad_alloc&){selectedCaught=true;}
    assert(selectedCaught && hostFaultState.attempts==1 && hostFaultState.liveBytes==0);
    configureHostFault(0);assert(!hostFaultState.filter[0]);
    bool constructorCaught=false;
    try{auto object=makeHostUnique<OwnedTest>("test.object",true);}
    catch(const std::runtime_error&){constructorCaught=true;}
    assert(constructorCaught && hostFaultState.liveBytes==0);
    // Counter-only contract: actual CUDA resource construction is tested in MATLAB.
    for(std::size_t at=0;at<=2;++at)for(int repeat=0;repeat<10;++repeat){
        configureResourceFault(at);bool acquired=false,caught=false;
        try{resourceFaultPoint("test.create");resourceAcquired();acquired=true;resourceFaultPoint("test.after-create");}
        catch(const std::bad_alloc&){caught=true;}
        if(acquired)resourceReleased(true);
        assert(caught==(at!=0));assert(resourceFaultState.live==0 && resourceFaultState.freeErrors==0);
    }
    return 0;
}
