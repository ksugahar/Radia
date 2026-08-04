function [spec, notes] = normalizeIHGeometryRoles(spec)
%NORMALIZEIHGEOMETRYROLES Repair crossed .vol/.step/.sol geometry inputs.
%   [spec, notes] = normalizeIHGeometryRoles(spec) checks the geometry
%   path fields of an IHDesignSpec-shaped struct (wp_vol, coil_vol,
%   em_vol, peec_step, qsurf_sol) against the file extensions each slot
%   accepts.  These are the inputs users re-point most often, and the
%   extension identifies a file's role uniquely: when the filled slots
%   are crossed (e.g. the coil .step typed into wp_vol and the
%   workpiece .vol into peec_step) and exactly ONE arrangement of the
%   same values fits every slot, that arrangement is applied.  Each
%   repair is returned in notes (string column) and emitted as a
%   warning so it shows in run logs.  Without a unique arrangement this
%   errors immediately, naming every offending slot and its expected
%   extensions, instead of failing later inside a mesh or STEP reader.
%   Values are never moved into slots the user left empty -- that would
%   silently change the selected method.
%
%   MATLAB twin of radia.ih_design.IHDesignSpec.normalize_geometry_roles.

arguments
    spec (1,1) struct
end

slotNames = ["wp_vol", "coil_vol", "em_vol", "peec_step", "qsurf_sol"];
slotExtensions = { ...
    [".vol", ".vol.gz"], ...
    [".vol", ".vol.gz"], ...
    [".vol", ".vol.gz"], ...
    [".step", ".stp"], ...
    ".sol"};

notes = strings(0, 1);
filled = strings(0, 1);
filledExtensions = {};
filledValues = strings(0, 1);
for index = 1:numel(slotNames)
    name = slotNames(index);
    if ~isfield(spec, name)
        continue
    end
    value = strtrim(string(spec.(name)));
    if strlength(value) == 0
        continue
    end
    filled(end + 1, 1) = name; %#ok<AGROW>
    filledExtensions{end + 1, 1} = slotExtensions{index}; %#ok<AGROW>
    filledValues(end + 1, 1) = value; %#ok<AGROW>
end

wrong = find(arrayfun(@(k) ~fitsSlot(filledExtensions{k}, ...
    filledValues(k)), (1:numel(filled)).'));
if isempty(wrong)
    return
end

wrongValues = filledValues(wrong);
orders = perms(1:numel(wrong));
valid = strings(0, 1);
validOrders = {};
for row = 1:size(orders, 1)
    candidate = wrongValues(orders(row, :));
    if all(arrayfun(@(k) fitsSlot(filledExtensions{wrong(k)}, ...
            candidate(k)), (1:numel(wrong)).'))
        key = strjoin(candidate, "|");
        if ~any(valid == key)
            valid(end + 1, 1) = key; %#ok<AGROW>
            validOrders{end + 1, 1} = candidate; %#ok<AGROW>
        end
    end
end

if numel(validOrders) ~= 1
    lines = strings(numel(wrong), 1);
    for k = 1:numel(wrong)
        lines(k) = sprintf("  %s=%s expects %s", filled(wrong(k)), ...
            filledValues(wrong(k)), ...
            strjoin(string(filledExtensions{wrong(k)}), " / "));
    end
    hint = "";
    stepExtensions = [".step", ".stp"];
    if any(arrayfun(@(k) fitsSlot(stepExtensions, filledValues(k)), wrong))
        hint = newline + "  Hint: a .step coil belongs in peec_step " + ...
            "(PEEC methods); meshes are .vol.";
    end
    error("radia:simulink:IHGeometryRoles", ...
        "Geometry inputs do not match their slots and no unique " + ...
        "reassignment of the same values fits:%s%s", ...
        newline + strjoin(lines, newline), hint);
end

candidate = validOrders{1};
for k = 1:numel(wrong)
    slot = filled(wrong(k));
    newValue = candidate(k);
    if newValue ~= filledValues(wrong(k))
        note = sprintf("geometry input reassigned by extension: " + ...
            "%s <- '%s' (was '%s')", slot, newValue, ...
            filledValues(wrong(k)));
        notes(end + 1, 1) = note; %#ok<AGROW>
        warning("radia:simulink:IHGeometryRolesReassigned", "%s", note);
    end
    spec.(slot) = char(newValue);
end
end

function tf = fitsSlot(extensions, value)
tf = any(endsWith(lower(string(value)), lower(string(extensions))));
end
