function tests = test_ih_geometry_roles
%TEST_IH_GEOMETRY_ROLES Crossed .vol/.step/.sol inputs repair by extension.
%   MATLAB twin of tests/test_ih_geometry_roles.py: the geometry path
%   fields are the most re-pointed inputs, so an unambiguous crossing is
%   repaired deterministically and anything else fails loudly.
tests = functiontests(localfunctions);
end

function testSwappedWpVolAndPeecStepAreRepaired(testCase)
spec = struct("wp_vol", "coil.step", "peec_step", "wp.vol");
warned = warning("off", "radia:simulink:IHGeometryRolesReassigned");
restore = onCleanup(@() warning(warned));
[spec, notes] = radia.simulink.normalizeIHGeometryRoles(spec);
verifyEqual(testCase, string(spec.wp_vol), "wp.vol");
verifyEqual(testCase, string(spec.peec_step), "coil.step");
verifyEqual(testCase, numel(notes), 2);
verifyTrue(testCase, any(contains(notes, "wp_vol")));
end

function testStpAndCaseInsensitiveExtensions(testCase)
spec = struct("wp_vol", "coil.STP", "peec_step", "wp.VOL");
warned = warning("off", "radia:simulink:IHGeometryRolesReassigned");
restore = onCleanup(@() warning(warned));
spec = radia.simulink.normalizeIHGeometryRoles(spec);
verifyEqual(testCase, string(spec.wp_vol), "wp.VOL");
verifyEqual(testCase, string(spec.peec_step), "coil.STP");
end

function testMatchingInputsUntouchedWithoutNotes(testCase)
spec = struct("wp_vol", "wp.vol", "coil_vol", "coil.vol");
[normalized, notes] = radia.simulink.normalizeIHGeometryRoles(spec);
verifyEqual(testCase, normalized, spec);
verifyEmpty(testCase, notes);
end

function testStepInWpVolWithoutStepSlotErrors(testCase)
spec = struct("wp_vol", "coil.step", "coil_vol", "wp.vol");
verifyError(testCase, ...
    @() radia.simulink.normalizeIHGeometryRoles(spec), ...
    "radia:simulink:IHGeometryRoles");
end

function testUnknownExtensionErrors(testCase)
spec = struct("wp_vol", "wp.msh", "peec_step", "coil.step");
verifyError(testCase, ...
    @() radia.simulink.normalizeIHGeometryRoles(spec), ...
    "radia:simulink:IHGeometryRoles");
end

function testQsurfAndEmVolSwapRepaired(testCase)
spec = struct("qsurf_sol", "model_em.vol", "em_vol", "model_q.sol");
warned = warning("off", "radia:simulink:IHGeometryRolesReassigned");
restore = onCleanup(@() warning(warned));
[spec, notes] = radia.simulink.normalizeIHGeometryRoles(spec);
verifyEqual(testCase, string(spec.qsurf_sol), "model_q.sol");
verifyEqual(testCase, string(spec.em_vol), "model_em.vol");
verifyNotEmpty(testCase, notes);
end
