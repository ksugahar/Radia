function result = ObjCenFld(object, fieldType)
%OBJCENFLD Return an object's center and the field at that center.

if nargin < 2
    fieldType = "B";
end
result = radia.internal.callMex('radia.ObjCenFld', double(object), ...
    char(string(fieldType)));
end
