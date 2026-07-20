function FldCmpPrc(options)
%FLDCMPPRC Set field-computation accuracy by option string.

radia.internal.callMex('radia.FldCmpPrc', char(string(options)));
end
