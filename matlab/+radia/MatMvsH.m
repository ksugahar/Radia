function magnetization = MatMvsH(material, component, hField)
%MATMVSH Evaluate magnetization from a material and H field.

magnetization = radia.internal.callMex('radia.MatMvsH', double(material), ...
    char(string(component)), double(hField));
end
