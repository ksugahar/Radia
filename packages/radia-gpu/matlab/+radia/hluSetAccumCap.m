function hluSetAccumCap(cap)
%HLUSETACCUMCAP Set the lazy low-rank update accumulator cap.
radia.internal.callMex('hlu.set_accum_cap', double(cap));
end
