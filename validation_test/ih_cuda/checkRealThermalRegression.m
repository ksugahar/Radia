function result=checkRealThermalRegression(fixturePath)
% Real thermal Krylov regression with a genuinely complex eddy system.
reference=jsondecode(fileread(fixturePath));cfg=makeResidentAlgebraConfig();
assert(radia_cuda_ih('live_count')==0);
cfg.bh_mode='linear';
h=radia_cuda_ih('create',cfg);owner=onCleanup(@()radia_cuda_ih('destroy',h));
temperature=zeros(6,3);heat=zeros(6,2);residuals=zeros(5,1);
for step=0:5
    [t,q,info]=radia_cuda_ih('outputs',h);
    temperature(step+1,:)=t;heat(step+1,:)=q;
    if step<5,radia_cuda_ih('update',h,1,300);end
end
electromagnetic=radia_cuda_ih('electromagnetic_state',h);
a=reference.matrix_real+1i*reference.matrix_imag;
b=reference.rhs_real+1i*reference.rhs_imag;
expected=reference.expected_electromagnetic_real+1i*reference.expected_electromagnetic_imag;
eddyError=norm(electromagnetic-expected)/norm(expected);
eddyResidual=norm(a*electromagnetic-b)/norm(b);
for step=1:5
    old=temperature(step,:)';increment=temperature(step+1,:)'-old;
    load=.01*(reference.heat_projection*heat(step+1,:)'-reference.stiffness*old ...
        +2*reference.convection*(300*ones(3,1)-old));
    residuals(step)=norm(reference.thermal_system*increment-load)/norm(load);
end
temperatureError=norm(temperature-reference.temperature,'fro')/norm(reference.temperature,'fro');
heatError=norm(heat(2:end,:)-reference.heat','fro')/norm(repmat(reference.heat',5,1),'fro');
assert(isreal(temperature) && all(isfinite(temperature),'all'));
assert(temperatureError<=reference.temperature_relative_limit && heatError<=reference.heat_relative_limit);
assert(eddyError<=1e-10 && eddyResidual<=reference.true_residual_limit);
assert(max(residuals)<=reference.true_residual_limit);
assert(info.eddy_true_residual<=reference.true_residual_limit && info.thermal_true_residual<=reference.true_residual_limit);
clear owner;assert(radia_cuda_ih('live_count')==0);
result=struct('accepted',true,'temperature_relative_error',temperatureError, ...
    'heat_relative_error',heatError,'complex_eddy_relative_error',eddyError, ...
    'independent_eddy_true_residual',eddyResidual,'independent_thermal_true_residuals',residuals, ...
    'complex_eddy_path_checked',true,'live_count',0,'scope','Synthetic fixed coupled algebra only');
end
