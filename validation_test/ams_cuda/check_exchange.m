function result=check_exchange(path,largeCount)
% Cross-language checks use asymmetric shapes and zero-based int32 indices.
d=radia.gpu.readAMSExchange(path);
assert(isequal(d.matrix,[0 1 2 3 4;5 6 7 8 9;10 11 12 13 14]));
assert(isa(d.indices,'int32') && isequal(d.indices,int32([0;2;4])));
assert(islogical(d.snapshot.beta_zero) && ~d.snapshot.beta_zero);
assert(strcmp(d.snapshot.schema,'radia.ams-accelerator.v1'));
assert(isequal(d.snapshot.spaces{1}.levels{1}.l1,[1;2;3]));
assert(isempty(d.empty));
if largeCount>0
    assert(isa(d.large,'double') && numel(d.large)==largeCount);
    assert(all(d.large==1));
end
result=struct('accepted',true,'large_values',largeCount,'matrix_shape',size(d.matrix));
end
