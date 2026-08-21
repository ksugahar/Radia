function response = contractHDivElementFillResponse( ...
        responseRows, elementDofBlocks, designElements, fillPatterns)
%CONTRACTHDIVELEMENTFILLRESPONSE One material-fill column per HDiv element.
%   Local/global DOF ordering is supplied by NGSolve. This function does not
%   reconstruct basis functions, Piola maps, or orientation transformations.
arguments
    responseRows (:,:) double
    elementDofBlocks (1,:) cell
    designElements (:,1) double {mustBeInteger,mustBePositive}
    fillPatterns (1,:) cell
end
if numel(designElements) ~= numel(fillPatterns)
    error("radia:topopt:FillPatternCount", ...
        "fillPatterns must contain one entry per design element.");
end
response = zeros(size(responseRows,1),numel(designElements));
for column = 1:numel(designElements)
    element = designElements(column);
    if element > numel(elementDofBlocks)
        error("radia:topopt:DesignElementRange", ...
            "A design element is outside elementDofBlocks.");
    end
    dofs = elementDofBlocks{element};
    pattern = fillPatterns{column};
    if numel(dofs) ~= numel(pattern)
        error("radia:topopt:FillPatternSize", ...
            "A fill pattern does not match its element block.");
    end
    response(:,column) = responseRows(:,dofs) * pattern(:);
end
end
