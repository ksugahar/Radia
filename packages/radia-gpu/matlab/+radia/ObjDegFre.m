function count = ObjDegFre(object)
%OBJDEGFRE Return the number of Radia interaction degrees of freedom.

count = radia.internal.callMex('radia.ObjDegFre', double(object));
end
