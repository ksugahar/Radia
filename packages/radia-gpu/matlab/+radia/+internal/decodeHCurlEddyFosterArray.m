function value = decodeHCurlEddyFosterArray(encoded, expectedShape)
%DECODEHCURLEDDYFOSTERARRAY Decode one row-major array of a Foster exchange.
%   The Python exporter writes {shape, values} in C order.  A matrix is the
%   MATLAB column-major reshape with the two dimensions swapped, then
%   transposed; a [force, state, port] tensor is rebuilt as
%   [port, state, force] and permuted back.

arguments
    encoded (1,1) struct
    expectedShape (1,:) double
end

if ~isfield(encoded, "shape") || ~isfield(encoded, "values") || ...
        ~isequal(double(encoded.shape(:)).', expectedShape)
    error("radia:simulink:HCurlFosterExchange", ...
        "exchange array shape does not match the declared dimensions.");
end
values = double(encoded.values(:));
if numel(values) ~= prod(expectedShape) || any(~isfinite(values))
    error("radia:simulink:HCurlFosterExchange", ...
        "exchange array values are incomplete or non-finite.");
end
switch numel(expectedShape)
    case 1
        value = values;
    case 2
        value = reshape(values, [expectedShape(2), expectedShape(1)]).';
    case 3
        raw = reshape(values, fliplr(expectedShape));
        value = permute(raw, [3, 2, 1]);
    otherwise
        error("radia:simulink:HCurlFosterExchange", ...
            "exchange arrays must have one to three dimensions.");
end
end
