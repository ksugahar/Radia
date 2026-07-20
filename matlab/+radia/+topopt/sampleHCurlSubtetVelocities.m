function velocity=sampleHCurlSubtetVelocities(referenceVertices,parentElements,deformationModes)
%SAMPLEHCURLSUBTETVELOCITIES Evaluate VectorH1 modes at analytic sub-TET vertices.
arguments
    referenceVertices double {mustBeFinite}
    parentElements (:,1) {mustBeInteger,mustBeNonnegative}
    deformationModes
end
if iscell(deformationModes)
    handles=zeros(numel(deformationModes),1,'uint64');
    for k=1:numel(deformationModes)
        if ~isa(deformationModes{k},'radia.ngsolve.GridFunction')
            error("radia:topopt:DeformationMode", ...
                "Every deformation mode must be a radia.ngsolve.GridFunction.");
        end
        handles(k)=deformationModes{k}.nativeHandle();
    end
else
    handles=zeros(numel(deformationModes),1,'uint64');
    for k=1:numel(deformationModes)
        if ~isa(deformationModes(k),'radia.ngsolve.GridFunction')
            error("radia:topopt:DeformationMode", ...
                "Every deformation mode must be a radia.ngsolve.GridFunction.");
        end
        handles(k)=deformationModes(k).nativeHandle();
    end
end
velocity=radia.internal.callMex( ...
    'hcurl.topopt.sample_subtet_velocities',double(referenceVertices), ...
    int32(parentElements),handles);
end
