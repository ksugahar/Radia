function grams=assembleHCurlCellCurlGrams(space,basis,options)
%ASSEMBLEHCURLCELLCURLGRAMS NGSolve-owned cell-local reduced curl Grams.
arguments
    space (1,1) radia.ngsolve.FESpace
    basis double {mustBeFinite}
    options.ElementIndices (:,1) {mustBeInteger,mustBeNonnegative}=int32.empty
end
grams=radia.internal.callMex('hcurl.topopt.cell_curl_grams', ...
    space.nativeHandle(),double(basis),int32(options.ElementIndices));
end
