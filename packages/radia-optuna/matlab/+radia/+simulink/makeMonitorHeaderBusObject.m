function [busObject, value] = makeMonitorHeaderBusObject(options)
%MAKEMONITORHEADERBUSOBJECT Define the common Radia monitor header.
arguments
    options.AssignToBase (1,1) logical = true
end

value = struct( ...
    "schema_version", uint16(1), ...
    "status_code", 0, ...
    "time_s", 0, ...
    "revision", 0, ...
    "error_code", 0);
busObject = makeScalarBus(value, ...
    "Radia monitor lifecycle header v1");
if options.AssignToBase
    assignin("base", "RadiaMonitorHeaderV1", busObject);
end
end

function busObject = makeScalarBus(value, description)
names = string(fieldnames(value));
elements = repmat(Simulink.BusElement, numel(names), 1);
for index = 1:numel(names)
    fieldValue = value.(names(index));
    element = Simulink.BusElement;
    element.Name = char(names(index));
    element.Dimensions = 1;
    if islogical(fieldValue)
        element.DataType = "boolean";
    else
        element.DataType = class(fieldValue);
    end
    element.DimensionsMode = "Fixed";
    element.Complexity = "real";
    elements(index) = element;
end
busObject = Simulink.Bus;
busObject.Description = description;
busObject.Elements = elements;
end
