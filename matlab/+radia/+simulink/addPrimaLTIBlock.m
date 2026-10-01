function block = addPrimaLTIBlock(blockPath, model, options)
%ADDPRIMALTIBLOCK Place a PRIMA port model in a Simulink LTI System block.
%   block = radia.simulink.addPrimaLTIBlock(blockPath, model) adds the
%   Control System Toolbox "LTI System" block at blockPath (the model must
%   be open) and points it at model.sys, stored in the model workspace.
%   model comes from radia.simulink.loadPrimaPortModel.
%
%   The LTI System block needs a proper model.  An inductive port is
%   improper in the impedance orientation (Z grows like s L); export it in
%   the admittance orientation instead.  An improper model is an error,
%   never silently differentiated.
%
%   Options:
%     VariableName  model-workspace variable holding the ss object
%                   (default "primaPortModel_<block name>")
%     InitialState  initial state vector of the explicit realization
%                   (default zeros)

arguments
    blockPath (1,1) string
    model (1,1) struct
    options.VariableName (1,1) string = ""
    options.InitialState double = []
end

if ~isfield(model, "sys") || ~isa(model.sys, "ss")
    error("radia:simulink:PrimaLTIBlock", ...
        "model must come from radia.simulink.loadPrimaPortModel.");
end
if ~isproper(model.sys)
    error("radia:simulink:PrimaImproper", ...
        "the %s model is improper; export it in the other orientation for the LTI System block.", ...
        model.orientation);
end
explicit = ss(model.sys, "explicit");
nState = size(explicit.A, 1);
x0 = options.InitialState;
if isempty(x0)
    x0 = zeros(nState, 1);
elseif ~isequal(size(x0(:)), [nState, 1]) || any(~isfinite(x0))
    error("radia:simulink:PrimaLTIBlock", ...
        "InitialState must be a finite vector of length %d.", nState);
end

parts = split(blockPath, "/");
if numel(parts) < 2 || ~bdIsLoaded(parts(1))
    error("radia:simulink:PrimaLTIBlock", ...
        "blockPath must lie in a loaded model: %s", blockPath);
end
name = options.VariableName;
if name == ""
    name = "primaPortModel_" + matlab.lang.makeValidName(parts(end));
end
if ~isvarname(name)
    error("radia:simulink:PrimaLTIBlock", "VariableName must be a MATLAB identifier.");
end
workspace = get_param(parts(1), "ModelWorkspace");
assignin(workspace, name, explicit);
assignin(workspace, name + "_x0", x0(:));

block = add_block("cstblocks/LTI System", blockPath);
set_param(block, "sys", name, "IC", name + "_x0");
end
