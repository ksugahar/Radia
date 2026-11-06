function input=readHDivParityInput(inputPath)
% The canonical device parity input written by hdiv_device_parity.py make-input.
nd=double(h5readatt(inputPath,'/','dofs'));
chargeFields={'cell_vertices','face_vertices','charge_host','charge_kind','charge_exponents','tet_rule', ...
    'tri_rule','tet_low_rule','tri_low_rule','cell_class','cell_dofs','cell_tables','cell_exponents', ...
    'face_class','face_dofs','face_tables','face_exponents'};
if strcmp(textAttribute(inputPath,'/charges','family'),'mixed')
    chargeFields=[chargeFields {'hex_rules','hex_bands','host_types','host_class','class_info','host_dofs','class_tables'}];
end
charges=struct();
% HDF5 row-major tables read as reversed MATLAB arrays: (:) is the native flat order.
for k=1:numel(chargeFields),charges.(chargeFields{k})=column(inputPath,['/charges/' chargeFields{k}]);end
charges.far_factor=double(h5readatt(inputPath,'/charges','far_factor'));
charges.compression=textAttribute(inputPath,'/charges','compression');
charges.family=textAttribute(inputPath,'/charges','family');
charges.dofs=nd;
charges.aca=struct('eps',double(h5readatt(inputPath,'/charges','aca_eps')),'eta',double(h5readatt(inputPath,'/charges','aca_eta')), ...
    'leaf_size',double(h5readatt(inputPath,'/charges','aca_leaf_size')),'max_rank',double(h5readatt(inputPath,'/charges','aca_max_rank')), ...
    'rank_cap',textAttribute(inputPath,'/charges','aca_rank_cap'));
evaluationFields={'rule','reference_shapes','vertices','element_class','element_dofs','element_points'};
groups=h5info(inputPath,'/evaluation').Groups;
for g=numel(groups):-1:1
    for k=1:numel(evaluationFields)
        evaluation(g,1).(evaluationFields{k})=column(inputPath,sprintf('/evaluation/group%d/%s',g,evaluationFields{k}));
    end
end
law=struct('knots',column(inputPath,'/law/knots'),'coefficients',h5read(inputPath,'/law/coefficients')', ...
    'integral_base',column(inputPath,'/law/integral_base'));
for name={'cap','barrier_slope','last_field','initial_reluctivity'}
    law.(name{1})=double(h5readatt(inputPath,'/law',name{1}));
end
newton=struct();
for name={'nonlinear_tolerance','linear_tolerance','maximum_iterations','maximum_linear_iterations','degree','lower_fraction'}
    newton.(name{1})=double(h5readatt(inputPath,'/newton',name{1}));
end
newton.preconditioner=textAttribute(inputPath,'/newton','preconditioner');
input=struct('name',textAttribute(inputPath,'/','case'),'dofs',nd,'charges',charges,'evaluation',evaluation,'law',law, ...
    'apply',column(inputPath,'/vectors/apply_input'),'rhs',column(inputPath,'/vectors/rhs'),'initial',column(inputPath,'/vectors/initial'), ...
    'inv_chi',double(h5readatt(inputPath,'/linear','inv_chi')),'tolerance',double(h5readatt(inputPath,'/linear','tolerance')), ...
    'maximum_iterations',double(h5readatt(inputPath,'/linear','maximum_iterations')), ...
    'preconditioner',textAttribute(inputPath,'/linear','preconditioner'),'newton',newton);
end

function value=column(path,name)
value=h5read(path,name);value=value(:);
end

function value=textAttribute(path,location,name)
value=char(string(h5readatt(path,location,name)));
end
