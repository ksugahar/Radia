function result=checkFrameAdmission(mesh,cpu,large)
% Full HDiv numerical gate and conservative frame coefficient telemetry.
if nargin<3,large=false;end
validateattributes(large,{'logical'},{'scalar'});
checks=checkFinalResidual(mesh,cpu,large);
for mode={'baseline','candidate','repeated'}
    stats=checks.(mode{1}).operator_stats;
    assert(stats.assembly_frame_leaf_count>=stats.assembly_leaf_count && stats.assembly_frame_leaf_count>0);
    assert(stats.assembly_frame_upper_dense_coefficients>0);
    assert(stats.assembly_frame_storage_cap_coefficients>=stats.stored_coefficient_count);
    assert(8*stats.assembly_frame_storage_cap_coefficients<=2*2^30);
end
assert(radia_cuda_hmatrix('live_count')==0 && radia_cuda_projection('live_count')==0);
result=checks;result.frame_coefficient_bounds_cover_stored_values=true;
result.scope='Bounded full nonlinear HDiv solve and admitted C frame coefficient bounds; not complete memory-peak admission, larger-capacity or performance acceptance';
end
