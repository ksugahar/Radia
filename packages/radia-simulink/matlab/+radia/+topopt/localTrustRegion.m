function bounds=localTrustRegion(elementSizes,options)
%LOCALTRUSTREGION Per-cell normal-displacement limits from local mesh size.
arguments
 elementSizes (:,1) double {mustBePositive}
 options.Fraction (1,1) double {mustBePositive}=0.1
 options.Minimum double=double.empty
 options.Maximum double=double.empty
end
if options.Fraction>=1, error("radia:topopt:TrustFraction","Fraction must be below one."); end
bounds=options.Fraction*elementSizes;
if ~isempty(options.Minimum), bounds=max(bounds,options.Minimum); end
if ~isempty(options.Maximum), bounds=min(bounds,options.Maximum); end
end
