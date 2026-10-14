function [U, singularValues, V, acaRank] = acaTsvd( ...
        matrix, modes, options)
%ACATSVD ACA+ compression followed by QR/TSVD recompression.
arguments
    matrix (:,:) double
    modes (1,1) double {mustBeInteger, mustBePositive}
    options.MaxRank (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    options.AcaTolerance (1,1) double {mustBePositive} = 1e-4
end
[U, singularValues, V, acaRank] = radia.internal.callMex( ...
    'stream.aca_tsvd', matrix, modes, options.MaxRank, ...
    options.AcaTolerance);
end
