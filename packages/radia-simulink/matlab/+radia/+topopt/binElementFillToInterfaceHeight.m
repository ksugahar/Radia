function result = binElementFillToInterfaceHeight( ...
        fillFraction, elementVolumes, firstCoordinate, secondCoordinate, ...
        firstEdges, secondEdges, binAreas)
%BINELEMENTFILLTOINTERFACEHEIGHT Conservatively map fill volume to height.
arguments
    fillFraction (:,1) double
    elementVolumes (:,1) double {mustBePositive}
    firstCoordinate (:,1) double
    secondCoordinate (:,1) double
    firstEdges (:,1) double
    secondEdges (:,1) double
    binAreas (:,:) double {mustBePositive}
end
count = numel(fillFraction);
if any([numel(elementVolumes),numel(firstCoordinate), ...
        numel(secondCoordinate)] ~= count)
    error("radia:topopt:InterfaceElementSize", ...
        "Fill, volume, and coordinates must have equal lengths.");
end
if numel(firstEdges) < 2 || numel(secondEdges) < 2 || ...
        any(diff(firstEdges) <= 0) || any(diff(secondEdges) <= 0)
    error("radia:topopt:InterfaceEdges", ...
        "Interface bin edges must be strictly increasing.");
end
shape = [numel(firstEdges)-1,numel(secondEdges)-1];
if ~isequal(size(binAreas),shape)
    error("radia:topopt:InterfaceAreaSize", ...
        "binAreas must contain one positive area per bin.");
end
[~,~,firstBin] = histcounts(firstCoordinate,firstEdges);
[~,~,secondBin] = histcounts(secondCoordinate,secondEdges);
if any(firstBin == 0) || any(secondBin == 0)
    error("radia:topopt:InterfaceCoordinateRange", ...
        "Every design element must lie inside the interface bins.");
end
signedVolume = accumarray([firstBin,secondBin], ...
    fillFraction .* elementVolumes,shape,@sum,0);
elementCount = accumarray([firstBin,secondBin],1,shape,@sum,0);
result = struct( ...
    first_edges=firstEdges,second_edges=secondEdges,bin_areas=binAreas, ...
    signed_volume=signedVolume,height_change=-signedVolume./binAreas, ...
    element_count=elementCount);
end
