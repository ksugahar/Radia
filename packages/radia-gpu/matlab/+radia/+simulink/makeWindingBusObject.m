function busObject = makeWindingBusObject(runtime, options)
%MAKEWINDINGBUSOBJECT Create the fixed-width winding configuration Bus.

arguments
    runtime (1,1) struct
    options.Name (1,1) string = "RadiaWindingBus"
    options.AssignToBase (1,1) logical = true
end
busObject = makeNumericBus(runtime, ...
    "Fixed-width Radia winding, terminal, and .vol region contract");
if options.AssignToBase,assignin("base",char(options.Name),busObject);end
end

function bus = makeNumericBus(data,description)
names=string(fieldnames(data));
elements=repmat(Simulink.BusElement,numel(names),1);
for k=1:numel(names)
    value=data.(names(k));
    if ~(isnumeric(value)||islogical(value)) || ~isreal(value)
        error("radia:simulink:WindingBusField", ...
            "Runtime field '%s' must be real numeric or logical.",names(k));
    end
    element=Simulink.BusElement; element.Name=char(names(k));
    if isscalar(value),element.Dimensions=1;else,element.Dimensions=size(value);end
    if islogical(value),element.DataType="boolean";else,element.DataType=class(value);end
    element.DimensionsMode="Fixed"; element.Complexity="real";
    elements(k)=element;
end
bus=Simulink.Bus; bus.Description=description; bus.Elements=elements;
end
