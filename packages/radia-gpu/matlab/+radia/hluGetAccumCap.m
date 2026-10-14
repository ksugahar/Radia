function cap = hluGetAccumCap()
%HLUGETACCUMCAP Return the lazy low-rank update accumulator cap.
cap = radia.internal.callMex('hlu.get_accum_cap');
end
