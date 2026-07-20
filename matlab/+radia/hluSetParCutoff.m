function hluSetParCutoff(cutoff)
%HLUSETPARCUTOFF Set the H-LU block-area threshold for parallel recursion.
radia.internal.callMex('hlu.set_par_cutoff', double(cutoff));
end
