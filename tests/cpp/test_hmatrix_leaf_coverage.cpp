#include "../../src/core/gpu/hdiv/hmatrix_leaf_coverage.hpp"
#include <cassert>
#include <vector>
struct Block { std::size_t row,column,rows,columns; };
using Blocks=std::vector<Block>;
void reject(std::size_t n, const Blocks& blocks) {
    bool failed=false;
    try { radia_hmatrix_cuda::requireUpperLeafPartition(n,blocks); }
    catch(const std::invalid_argument&) { failed=true; }
    assert(failed);
}
int main() {
    const Blocks good{{0,0,2,2},{0,2,2,1},{2,2,1,1}};
    radia_hmatrix_cuda::requireUpperLeafPartition(3,good);
    auto duplicate=good;duplicate.push_back(good[1]);reject(3,duplicate);
    auto gap=good;gap.erase(gap.begin()+1);reject(3,gap);
    auto overlap=good;overlap.push_back({0,1,1,1});reject(3,overlap);
    reject(3,{{0,0,2,2},{0,2,3,1},{2,2,1,1}});
    reject(3,{{0,0,3,2}});reject(3,{{2,0,1,1}});
    reject(0,{});reject(4609,{});
    // Exercise bit-word boundaries and the admission cap without factor storage.
    radia_hmatrix_cuda::requireUpperLeafPartition(65,Blocks{{0,0,65,65}});
    radia_hmatrix_cuda::requireUpperLeafPartition(4608,Blocks{{0,0,4608,4608}});
    auto reordered=good;std::reverse(reordered.begin(),reordered.end());
    radia_hmatrix_cuda::requireUpperLeafPartition(3,reordered);
}
