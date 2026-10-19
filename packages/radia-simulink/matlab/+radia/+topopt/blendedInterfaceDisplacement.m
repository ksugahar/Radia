function displacement = blendedInterfaceDisplacement( ...
        normalCoordinate, interfaceCoordinate, lowerFixedCoordinate, ...
        upperFixedCoordinate, interfaceHeightChange)
%BLENDEDINTERFACEDISPLACEMENT Zero displacement at aperture and pole root.
arguments
    normalCoordinate double
    interfaceCoordinate double
    lowerFixedCoordinate double
    upperFixedCoordinate double
    interfaceHeightChange double
end
[normalCoordinate,interfaceCoordinate,lowerFixedCoordinate, ...
    upperFixedCoordinate,interfaceHeightChange] = compatibleArrays( ...
    normalCoordinate,interfaceCoordinate,lowerFixedCoordinate, ...
    upperFixedCoordinate,interfaceHeightChange);
if any(interfaceCoordinate <= lowerFixedCoordinate,"all") || ...
        any(upperFixedCoordinate <= interfaceCoordinate,"all")
    error("radia:topopt:InterfaceOrdering", ...
        "The interface must lie strictly between the fixed surfaces.");
end
weight = zeros(size(normalCoordinate));
below = normalCoordinate > lowerFixedCoordinate & ...
    normalCoordinate < interfaceCoordinate;
above = normalCoordinate >= interfaceCoordinate & ...
    normalCoordinate < upperFixedCoordinate;
weight(below) = (normalCoordinate(below)-lowerFixedCoordinate(below)) ./ ...
    (interfaceCoordinate(below)-lowerFixedCoordinate(below));
weight(above) = (upperFixedCoordinate(above)-normalCoordinate(above)) ./ ...
    (upperFixedCoordinate(above)-interfaceCoordinate(above));
displacement = interfaceHeightChange .* weight;
end

function varargout = compatibleArrays(varargin)
shape = [1,1];
for index = 1:nargin
    if ~isscalar(varargin{index})
        if isequal(shape,[1,1])
            shape = size(varargin{index});
        elseif ~isequal(size(varargin{index}),shape)
            error("radia:topopt:InterfaceBroadcast", ...
                "Nonscalar interface arrays must have equal sizes.");
        end
    end
end
varargout = varargin;
for index = 1:nargin
    if isscalar(varargout{index})
        varargout{index} = repmat(varargout{index},shape);
    end
end
end
