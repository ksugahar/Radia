function r=checkDeviceFaults(families,tracePath)
% Guards: every device allocation, stream wait, host vector, resource handle and mx allocation
% reached by the device commands, failed one at a time, raises radia:hdiv:device, publishes no
% handle, returns device and host memory to their level before the call, and leaves the
% prepared state usable: the same call then gives the same bytes. Only a stream wait in a
% destructor, which cannot raise, may pass silently; those are counted as absorbed_waits.
% The first 40 sites of a call are all failed and later ones stepped: see targets. Each destroy
% reads its command through mx, and a failed read must leave the handle usable; that is tested
% once per family, and the cleanup of a failed use is not itself failed. The solve, assembly and
% stats commands share one kernel across mesh families, so the first family sweeps them, each
% preconditioner path included; later families sweep what differs by family: creation, apply,
% mass and material evaluation. A failed linearization has overwritten device state, so the tangent
% of the earlier one must then be refused.
% Assumes: radia_hdiv_device_fault is the same kernel object as the shipped MEX, with the shim
% counting cudaMalloc, cudaFree and cudaStreamSynchronize. CUDA library workspaces are not counted.
% tracePath names the site under test before each injection, so a process abort can be located.
kinds={'device','sync','host','resource','mx'};
radia_hdiv_device_fault('fault_arm','none',0);
rows=struct('family',{},'command',{},'kind',{},'sites',{},'injected',{},'absorbed',{},'seconds',{});
for f=1:numel(families)
    family=families(f);
    op=radia_hdiv_device_fault('create_operator',family.charges);opOwner=onCleanup(@()radia_hdiv_device_fault('destroy',op));
    mass=radia_hdiv_device_fault('create_mass',family.evaluation,family.dofs);massOwner=onCleanup(@()radia_hdiv_device_fault('destroy',mass));
    material=radia_hdiv_device_fault('create_material',family.evaluation,family.dofs,family.law);materialOwner=onCleanup(@()radia_hdiv_device_fault('destroy',material));
    counters=radia_hdiv_device_fault('assembly',op).hmatrix;
    assert(strcmp(family.name,'hex-pair') || counters.lowrank_leaves>=1,'radia:hdiv:fault','%s reached no device ACA low-rank leaf',family.name);
    calls=familyCalls(family,op,mass,material);
    if f==1
        calls=[calls;{'create_analytic_gram',@analyticGram}];
    else
        calls=calls(~startsWith(calls(:,1),{'solve_','assembly','stats','inverse_bh_samples'}),:);
    end
    for c=1:size(calls,1)
        baseline=calls{c,2}();
        for k=1:numel(kinds)
            t=tic;[sites,absorbed]=sweep(kinds{k},[family.name ' ' calls{c,1}],calls{c,2},baseline,tracePath);
            rows(end+1)=struct('family',family.name,'command',calls{c,1},'kind',kinds{k},'sites',sites,'injected',numel(targets(sites)),'absorbed',absorbed,'seconds',toc(t)); %#ok<AGROW>
        end
    end
    checkFailedDestroy(op);
    checkFailedLinearization(material,1e3*sin((1:family.dofs)'),cos((1:family.dofs)'));
    clear materialOwner massOwner opOwner
    assert(radia_hdiv_device_fault('live_count')==0,'radia:hdiv:fault','%s left live handles',family.name);
end
final=radia_hdiv_device_fault('fault_report');
assert(final.device_live_bytes==0 && final.device_live_count==0 && final.free_errors==0 && ...
    final.resource_live==0 && final.resource_free_errors==0 && final.host_live_bytes==0, ...
    'radia:hdiv:fault','Resources remain after every handle is destroyed: %s',jsonencode(final));
sites=[rows.sites];
for k=1:numel(kinds)
    total=sum(sites(strcmp({rows.kind},kinds{k})));
    assert(total>0,'radia:hdiv:fault','No %s site was reached',kinds{k});
end
r=struct('accepted',true,'sites',sum(sites),'failure_events',sum([rows.injected]),'absorbed_waits',sum([rows.absorbed]),'rows',rows,'final',final, ...
    'scope','Failure paths of the device commands on the test build; CUDA library workspaces are not counted');
end

function calls=familyCalls(family,op,mass,material)
n=family.dofs;x=family.apply;
state=1e3*sin((1:n)');direction=cos((1:n)');
amplitudes=[1;1e3;1e5];
calls={'create_operator',@()createThen('create_operator',{family.charges},'apply',x); ...
    'create_mass',@()createThen('create_mass',{family.evaluation,n},'mass_apply',x); ...
    'create_material',@()createThen('create_material',{family.evaluation,n,family.law},'material_evaluate',state,direction); ...
    'apply',@()radia_hdiv_device_fault('apply',op,x); ...
    'assembly',@()counterFree(radia_hdiv_device_fault('assembly',op)); ...
    'stats',@()rmfield(radia_hdiv_device_fault('stats',op),'applications'); ...
    'solve_linear',@()radia_hdiv_device_fault('solve_linear',op,mass,family.rhs,family.inv_chi,family.tolerance,family.maximum_iterations,family.preconditioner); ...
    'solve_linear_cholesky',@()radia_hdiv_device_fault('solve_linear',op,mass,family.rhs,family.inv_chi,family.tolerance,family.maximum_iterations,'mass-cholesky'); ...
    'solve_nonlinear',@()radia_hdiv_device_fault('solve_nonlinear',op,material,family.rhs,family.initial,family.newton); ...
    'solve_nonlinear_chebyshev',@()radia_hdiv_device_fault('solve_nonlinear',op,material,family.rhs,family.initial,setfield(family.newton,'preconditioner','reference-chebyshev-graph')); ...
    'mass_apply',@()radia_hdiv_device_fault('mass_apply',mass,x); ...
    'material_evaluate',@()rmfield(radia_hdiv_device_fault('material_evaluate',material,state,direction),'linearization_epoch'); ...
    'material_linearize',@()rmfield(radia_hdiv_device_fault('material_linearize',material,state),'linearization_epoch'); ...
    'material_apply_tangent',@()tangent(material,state,direction); ...
    'inverse_bh_samples',@()radia_hdiv_device_fault('inverse_bh_samples',family.law,amplitudes)};
end

function y=createThen(create,inputs,use,varargin)
% The handle of a failed use is released here, as the MATLAB classes do.
h=radia_hdiv_device_fault(create,inputs{:});owner=onCleanup(@()release(h));
y=radia_hdiv_device_fault(use,h,varargin{:});
end

function release(h)
radia_hdiv_device_fault('fault_mx_pause',true);
radia_hdiv_device_fault('destroy',h);
radia_hdiv_device_fault('fault_mx_pause',false);
end

function checkFailedDestroy(h)
live=radia_hdiv_device_fault('live_count');
radia_hdiv_device_fault('fault_arm','mx',1);
failed=false;identifier='';
try
    radia_hdiv_device_fault('destroy',h);
catch e
    failed=true;identifier=e.identifier;
end
radia_hdiv_device_fault('fault_arm','none',0);
assert(failed && strcmp(identifier,'radia:hdiv:device'),'radia:hdiv:fault','A failed destroy gave [%s]',identifier);
assert(radia_hdiv_device_fault('live_count')==live,'radia:hdiv:fault','A failed destroy released the handle: %d live, expected %d',radia_hdiv_device_fault('live_count'),live);
end

function checkFailedLinearization(material,state,direction)
epoch=radia_hdiv_device_fault('material_linearize',material,state).linearization_epoch;
radia_hdiv_device_fault('material_apply_tangent',material,direction,epoch);
radia_hdiv_device_fault('fault_arm','sync',1);
failed='';
try
    radia_hdiv_device_fault('material_linearize',material,2*state);
catch e
    failed=e.identifier;
end
radia_hdiv_device_fault('fault_arm','none',0);
assert(strcmp(failed,'radia:hdiv:device'),'radia:hdiv:fault','A failed linearization gave [%s]',failed);
refused='';
try
    radia_hdiv_device_fault('material_apply_tangent',material,direction,epoch);
catch e
    refused=e.identifier;
end
assert(strcmp(refused,'radia:hdiv:device'),'radia:hdiv:fault','The tangent of the earlier linearization gave [%s] after a failed one',refused);
end

function y=tangent(material,state,direction)
% Each linearization starts a new epoch, so this call covers both commands.
e=radia_hdiv_device_fault('material_linearize',material,state).linearization_epoch;
y=radia_hdiv_device_fault('material_apply_tangent',material,direction,e);
end

function y=analyticGram()
first=[0 0 0;1 0 0;0 1 0;0 0 1];second=first+[3 0 0];
cells=reshape([first;second]',[],1);faces=reshape([first(1:3,:);second([1 3 4],:)]',[],1);
y=createThen('create_analytic_gram',{cells,faces},'apply',[1;-2;3;-4]);
end

function a=counterFree(a)
% Transfer counts grow with each apply; they are not the state under test.
a=rmfield(a,'transfer_ledger');
end

function [sites,absorbed]=sweep(kind,command,call,baseline,tracePath)
stateless=contains(command,' create_') || contains(command,' inverse_bh_samples');
slow=contains(command,' solve_');
radia_hdiv_device_fault('fault_arm',kind,0);call();
sites=reached(kind,radia_hdiv_device_fault('fault_report'));
radia_hdiv_device_fault('fault_arm','none',0);
before=radia_hdiv_device_fault('fault_report');handles=radia_hdiv_device_fault('live_count');
absorbed=0;
for at=targets(sites)
    where=sprintf('%s %s site %d of %d',command,kind,at,sites);
    trace=fopen(tracePath,'a');fprintf(trace,'%s\n',where);fclose(trace);
    radia_hdiv_device_fault('fault_arm',kind,at);
    failed=false;identifier='';
    result=[];
    try
        result=call();
    catch e
        failed=true;identifier=e.identifier;
    end
    report=radia_hdiv_device_fault('fault_report');
    radia_hdiv_device_fault('fault_arm','none',0);
    % A stream wait in a destructor cannot raise and the work has finished by then, so only such a
    % wait may pass silently, and the answer must still be the right one.
    if ~failed && strcmp(kind,'sync')
        assert(isequal(getByteStreamFromArray(result),getByteStreamFromArray(baseline)),'radia:hdiv:fault','%s: an ignored wait changed the answer',where);
        absorbed=absorbed+1;
    else
        assert(failed,'radia:hdiv:fault','%s: not refused',where);
        assert(strcmp(identifier,'radia:hdiv:device'),'radia:hdiv:fault','%s: identifier %s',where,identifier);
    end
    % Cleanup after a failed stream wait waits again, so only that count may pass the injected site.
    assert(reached(kind,report)==at || (strcmp(kind,'sync') && reached(kind,report)>at), ...
        'radia:hdiv:fault','%s: reached %d',where,reached(kind,report));
    assert(radia_hdiv_device_fault('live_count')==handles,'radia:hdiv:fault','%s: live handles %d, expected %d',where,radia_hdiv_device_fault('live_count'),handles);
    assert(report.device_live_bytes==before.device_live_bytes && report.free_errors==before.free_errors, ...
        'radia:hdiv:fault','%s: device %d bytes and %d free errors, expected %d and %d',where, ...
        report.device_live_bytes,report.free_errors,before.device_live_bytes,before.free_errors);
    assert(report.host_live_bytes==before.host_live_bytes,'radia:hdiv:fault','%s: host %d bytes, expected %d',where,report.host_live_bytes,before.host_live_bytes);
    assert(report.resource_live==before.resource_live && report.resource_free_errors==before.resource_free_errors, ...
        'radia:hdiv:fault','%s: resources %d and %d free errors, expected %d and %d',where, ...
        report.resource_live,report.resource_free_errors,before.resource_live,before.resource_free_errors);
    % A failed create leaves no prepared state, so one repeat at the end of its sweep is enough.
    % A solve is too slow to repeat at every site; every other stateful call repeats each time.
    if at==sites || (~stateless && (~slow || mod(at,8)==1)),same(call,baseline,where);end
end
end

function list=targets(sites)
% Every early site, where setup paths differ; later ones mostly repeat a loop, so they are stepped.
full=40;samples=30;
list=1:min(sites,full);
if sites>full,list=unique([list full+ceil((sites-full)/samples):ceil((sites-full)/samples):sites sites]);end
end

function same(call,baseline,where)
assert(isequal(getByteStreamFromArray(call()),getByteStreamFromArray(baseline)),'radia:hdiv:fault','%s: the next call differs',where);
end

function n=reached(kind,report)
n=report.([kind '_attempts']);
end
