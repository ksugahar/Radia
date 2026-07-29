function buses = makeElectromechanicalBusObjects(options)
%MAKEELECTROMECHANICALBUSOBJECTS Define solver-independent dynamic I/O buses.

arguments
    options.MaxWindings (1,1) double {mustBeInteger,mustBePositive} = 16
    options.AssignToBase (1,1) logical = true
end
n=options.MaxWindings;
command=struct("schema_version",uint16(1),"winding_count",uint16(0), ...
    "excitation_mode",uint16(0), ... % 0=current, 1=voltage, 2=external circuit
    "terminal_voltage_V",zeros(n,1),"imposed_current_A",zeros(n,1), ...
    "rotor_angle_rad",0,"rotor_speed_rad_s",0,"load_torque_Nm",0, ...
    "translation_position_m",zeros(3,1), ...
    "translation_velocity_m_per_s",zeros(3,1),"load_force_N",zeros(3,1));
response=struct("schema_version",uint16(1),"winding_count",uint16(0), ...
    "terminal_current_A",zeros(n,1),"flux_linkage_Wb_turn",zeros(n,1), ...
    "back_emf_V",zeros(n,1),"electromagnetic_torque_Nm",0, ...
    "electromagnetic_force_N",zeros(3,1),"copper_loss_W",zeros(n,1), ...
    "iron_loss_W",0,"eddy_loss_W",0);
buses=struct("command",makeBus(command,"Radia electromechanical command"), ...
    "response",makeBus(response,"Radia electromechanical response"), ...
    "command_value",command,"response_value",response);
if options.AssignToBase
    assignin("base","RadiaMachineCommandBus",buses.command);
    assignin("base","RadiaMachineResponseBus",buses.response);
    assignin("base","radia_machine_command",command);
    assignin("base","radia_machine_response",response);
end
end

function bus=makeBus(data,description)
names=string(fieldnames(data)); elements=repmat(Simulink.BusElement,numel(names),1);
for k=1:numel(names)
    value=data.(names(k)); element=Simulink.BusElement; element.Name=char(names(k));
    if isscalar(value),element.Dimensions=1;else,element.Dimensions=size(value);end
    if islogical(value),element.DataType="boolean";else,element.DataType=class(value);end
    element.DimensionsMode="Fixed"; element.Complexity="real"; elements(k)=element;
end
bus=Simulink.Bus;bus.Description=description;bus.Elements=elements;
end
