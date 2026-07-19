function Z = clnImpedance(RDiag, LTridiag, frequency)
%CLNIMPEDANCE Evaluate a CLN I impedance at one frequency in hertz.
Z = radia.internal.callMex( ...
    'cln.impedance', double(RDiag), double(LTridiag), double(frequency));
end
