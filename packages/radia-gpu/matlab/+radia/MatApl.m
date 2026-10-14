function result = MatApl(object, material)
%MATAPL Apply a material to a Radia object.

result = radia.internal.callMex( ...
    'radia.MatApl', double(object), double(material));
end
