function result = checkHDivHdf5(inputPath, outputPath)
%CHECKHDIVHDF5 Non-square CSR and metadata regression; no GPU computation.
[fixture, metadata]=radia.gpu.readHDivH5(inputPath);
assert(isequal(size(fixture.operator),[2 3]),'radia:HDF5Test','CSR shape changed');
assert(nnz(fixture.operator)==3,'radia:HDF5Test','CSR nonzero count changed');
assert(isequal(full(fixture.operator),[1 0 2;0 3 0]),'radia:HDF5Test','CSR payload changed');
assert(isrow(metadata.labels),'radia:HDF5Test','Labels must be canonical row list');
assert(isequal(metadata.labels,{'boundary',char([28201 24230])}), ...
    'radia:HDF5Test','UTF-8 labels changed');
dense=rmfield(fixture,'operator');
assert(iscolumn(fixture.complex_vector),'radia:HDF5Test','Rank-one arrays must read as columns');
dense.rhs=real(fixture.complex_vector(:));
dense.expected=dense.rhs;
radia.gpu.writeHDivH5(outputPath,dense,metadata);
[back, backMetadata]=radia.gpu.readHDivH5(outputPath);
assert(isequal(size(back.rhs),[2 1]) && isequal(size(back.expected),[2 1]), ...
    'radia:HDF5Test','RHS/expected columns changed');
assert(norm(back.rhs-back.expected)==0,'radia:HDF5Test','Vector comparison changed');
assert(isequaln(back,dense),'radia:HDF5Test','Dense or complex values changed');
assert(isequal(backMetadata,metadata),'radia:HDF5Test','Metadata changed');
result=struct('accepted',true,'csr_shape',size(fixture.operator), ...
    'csr_nnz',nnz(fixture.operator),'label_shape',size(backMetadata.labels), ...
    'scope','Typed HDF5 interchange only; no GPU or physical acceptance');
end
