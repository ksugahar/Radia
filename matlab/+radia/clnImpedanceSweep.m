function Z = clnImpedanceSweep(RDiag, LTridiag, frequencies)
%CLNIMPEDANCESWEEP Evaluate a CLN I impedance over a frequency vector.
Z = radia.internal.callMex( ...
    'cln.impedance_sweep', double(RDiag), double(LTridiag), double(frequencies));
end
