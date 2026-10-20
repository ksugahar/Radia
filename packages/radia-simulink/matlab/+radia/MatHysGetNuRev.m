function nuRev = MatHysGetNuRev(material)
%MATHYSGETNUREV Return the reversible reluctivity.

nuRev = radia.internal.callMex('radia.MatHysGetNuRev', double(material));
end
