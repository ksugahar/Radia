function result=assembleHCurlResistanceShapeTangents(space,basis,deformationModes,options)
%ASSEMBLEHCURLRESISTANCESHAPETANGENTS Reduced curl-curl matrix and Piola tangents.
arguments
    space (1,1) radia.ngsolve.FESpace
    basis double {mustBeFinite}
    deformationModes
    options.Conductivity (1,1) double {mustBePositive,mustBeFinite}=1
    options.ElementIndices (:,1) {mustBeInteger,mustBeNonnegative}=int32.empty
end
handles=localHandles(deformationModes);
result=radia.internal.callMex( ...
    'hcurl.topopt.resistance_shape_tangents',space.nativeHandle(), ...
    double(basis),handles,options.Conductivity,int32(options.ElementIndices));
end

function handles=localHandles(modes)
if isempty(modes), handles=zeros(0,1,'uint64'); return; end
if iscell(modes)
    handles=zeros(numel(modes),1,'uint64');
    for k=1:numel(modes)
        if ~isa(modes{k},'radia.ngsolve.GridFunction')
            error("radia:topopt:DeformationMode", ...
                "Every deformation mode must be a radia.ngsolve.GridFunction.");
        end
        handles(k)=modes{k}.nativeHandle();
    end
else
    if ~all(arrayfun(@(x)isa(x,'radia.ngsolve.GridFunction'),modes))
        error("radia:topopt:DeformationMode", ...
            "Every deformation mode must be a radia.ngsolve.GridFunction.");
    end
    handles=arrayfun(@(x)x.nativeHandle(),modes(:));
end
end
