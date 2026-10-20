function busObject = makeFieldStudyBusObject(runtime, options)
%MAKEFIELDSTUDYBUSOBJECT Define the fixed-width field-study Simulink.Bus.

arguments
    runtime (1,1) struct
    options.Name (1,1) string = "RadiaStudyBus"
    options.AssignToBase (1,1) logical = true
end
names=string(fieldnames(runtime));elements=repmat(Simulink.BusElement,numel(names),1);
for k=1:numel(names)
    value=runtime.(names(k));
    if ~(isnumeric(value)||islogical(value)) || ~isreal(value)
        error("radia:simulink:FieldStudyBusField", ...
            "Runtime field '%s' must be real numeric or logical.",names(k));
    end
    element=Simulink.BusElement;element.Name=char(names(k));
    if isscalar(value),element.Dimensions=1;else,element.Dimensions=size(value);end
    if islogical(value),element.DataType="boolean";else,element.DataType=class(value);end
    element.DimensionsMode="Fixed";element.Complexity="real";elements(k)=element;
end
busObject=Simulink.Bus;busObject.Description="Fixed-width Radia multiphysics study contract";
busObject.Elements=elements;
if options.AssignToBase,assignin("base",char(options.Name),busObject);end
end
