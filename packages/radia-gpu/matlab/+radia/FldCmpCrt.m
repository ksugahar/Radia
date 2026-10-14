function FldCmpCrt(prcB, prcA, prcBInt, prcFrc, prcTrjCrd, prcTrjAng)
%FLDCMPCRT Set absolute field-computation accuracy criteria.

radia.internal.callMex('radia.FldCmpCrt', double(prcB), double(prcA), ...
    double(prcBInt), double(prcFrc), double(prcTrjCrd), double(prcTrjAng));
end
