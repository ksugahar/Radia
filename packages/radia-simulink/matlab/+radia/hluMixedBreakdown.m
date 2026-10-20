function breakdown = hluMixedBreakdown()
%HLUMIXEDBREAKDOWN Return H-LU mixed operand-kind counters.
breakdown = radia.internal.callMex('hlu.mixed_breakdown');
end
