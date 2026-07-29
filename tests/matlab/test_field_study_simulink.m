classdef test_field_study_simulink < matlab.unittest.TestCase
    methods (Test)
        function compilesScalarPhysicsToFixedBusAndRequest(testCase)
            [meshFile,cleanup]=makeMesh(); %#ok<ASGLU>
            materials=struct("domain",radia.simulink.makeMaterialSpec( ...
                RelativePermittivity=2.5,Density_kg_per_m3=1000, ...
                SpecificHeat_J_per_kgK=500,ThermalConductivity_W_per_mK=4));
            materialContract=radia.simulink.compileMaterialDictionary(materials, ...
                MeshFile=meshFile);
            study=radia.simulink.makeFieldStudySpec(Physics="steady_heat", ...
                Formulation="planar",ElementFamily="P2",ModelDepth_m=0.2, ...
                DirichletValues=struct("left",300), ...
                RobinBoundaries=struct("right",struct( ...
                    "transfer_w_per_m2_k",10,"ambient_k",290)), ...
                VolumetricSources=struct("domain",100));
            contract=radia.simulink.compileFieldStudy(study,materialContract);

            testCase.verifyEqual(contract.runtime.physics_code,uint16(2));
            testCase.verifyEqual(contract.runtime.element_order,uint16(2));
            testCase.verifyEqual(contract.request.materials.domain.coefficient_si,4);
            testCase.verifyEqual(contract.request.dirichlet_values.left,[300 0]);
            testCase.verifyFalse(contract.runtime_policy.python_per_step);
        end

        function compilesLinearAndNonlinearHarmonicWinding(testCase)
            [meshFile,cleanup]=makeMesh(); %#ok<ASGLU>
            copper=radia.simulink.makeMaterialSpec(MuR=1, ...
                Conductivity_S_per_m=5.8e7);
            materialContract=radia.simulink.compileMaterialDictionary( ...
                struct("copper",copper),MeshFile=meshFile, ...
                RegionMaterials=struct("domain","copper"));
            windingContract=radia.simulink.compileWindingDictionary( ...
                struct("coil",radia.simulink.makeWindingSpec( ...
                    Regions="domain",Turns=20)),materialContract);
            study=radia.simulink.makeFieldStudySpec(Physics="harmonic_eddy", ...
                Frequency_Hz=1000,DirichletValues=struct("outer",0), ...
                BranchCurrent_A=2+3i);
            contract=radia.simulink.compileFieldStudy(study,materialContract, ...
                WindingContract=windingContract);

            testCase.verifyEqual(contract.request.formulation,"planar");
            testCase.verifyEqual(contract.request.branches(1).turns,20);
            testCase.verifyEqual(contract.request.branch_current_a,[2 3]);
            testCase.verifyEqual( ...
                contract.request.materials.domain.conductivity_s_per_m,5.8e7);

            output=string(tempname("C:\temp"))+".json";
            cleanupOutput=onCleanup(@()deleteIfPresent(output));
            radia.simulink.writeFieldStudyRequest(contract,output);
            text=fileread(output);
            testCase.verifyMatches(text,'"branches"\s*:\s*\[');
            testCase.verifyMatches(text,'"branch_current_a"\s*:\s*\[\s*\[');

            study.dirichlet_values.outer=1;
            testCase.verifyError(@()radia.simulink.compileFieldStudy( ...
                study,materialContract,WindingContract=windingContract), ...
                "radia:simulink:FieldStudyHarmonicBoundary");

            nonlinear=radia.simulink.makeMaterialSpec(Conductivity_S_per_m=3, ...
                BH_B_T=[0;0.2;0.8;1.5],BH_H_A_per_m=[0;80;600;20000]);
            nonlinearContract=radia.simulink.compileMaterialDictionary( ...
                struct("steel",nonlinear),MeshFile=meshFile, ...
                RegionMaterials=struct("domain","steel"));
            nonlinearWinding=radia.simulink.compileWindingDictionary( ...
                struct("coil",radia.simulink.makeWindingSpec( ...
                    Regions="domain",Turns=20)),nonlinearContract);
            nonlinearStudy=radia.simulink.makeFieldStudySpec( ...
                Physics="harmonic_eddy",Frequency_Hz=1000, ...
                DirichletValues=struct("outer",0),BranchCurrent_A=2+3i);
            nonlinearCompiled=radia.simulink.compileFieldStudy( ...
                nonlinearStudy,nonlinearContract,WindingContract=nonlinearWinding);
            testCase.verifyEqual(numel( ...
                nonlinearCompiled.request.materials.domain.bh_curve),4);
            testCase.verifyEqual( ...
                nonlinearCompiled.request.materials.domain.bh_curve(3).b_t,0.8);
        end

        function compilesTransientHeatAndElectrostaticSystem(testCase)
            [meshFile,cleanup]=makeMesh(); %#ok<ASGLU>
            materials=struct("domain",radia.simulink.makeMaterialSpec( ...
                RelativePermittivity=2.5,Density_kg_per_m3=1000, ...
                SpecificHeat_J_per_kgK=500,ThermalConductivity_W_per_mK=4));
            materialContract=radia.simulink.compileMaterialDictionary(materials, ...
                MeshFile=meshFile);
            thermal=radia.simulink.makeFieldStudySpec(Physics="transient_heat", ...
                ElementFamily="Q1",Time_s=[0;0.1;0.2],Theta=1, ...
                InitialTemperature_K=300,VolumetricSources=struct("domain",8), ...
                VolumetricSourceHistory_W_per_m3=struct("domain",[8;8;8]));
            thermalContract=radia.simulink.compileFieldStudy(thermal,materialContract);

            testCase.verifyEqual(thermalContract.runtime.physics_code,uint16(5));
            testCase.verifyEqual(thermalContract.runtime.time_sample_count,uint16(3));
            testCase.verifyEqual(thermalContract.request.materials.domain. ...
                volumetric_heat_capacity_j_per_m3_k,5e5);
            testCase.verifyEqual(thermalContract.request.time_s,[0 0.1 0.2]);

            electrostatic=radia.simulink.makeFieldStudySpec( ...
                Physics="electrostatic_system",Conductors=["left";"right";"outer"], ...
                AppliedVoltages_V=[1;-1;0],ForceBoundaries=["left";"right"]);
            electrostaticContract=radia.simulink.compileFieldStudy( ...
                electrostatic,materialContract);
            testCase.verifyEqual(electrostaticContract.runtime.physics_code,uint16(4));
            testCase.verifyEqual(electrostaticContract.runtime.conductor_count,uint16(3));
            testCase.verifyEqual(string(electrostaticContract.request.conductors), ...
                ["left";"right";"outer"]);
            testCase.verifyEqual(electrostaticContract.request.applied_voltages_v,[1 -1 0]);
        end

        function writesRequestWithNoComplexJsonValues(testCase)
            [meshFile,cleanup]=makeMesh(); %#ok<ASGLU>
            materialContract=radia.simulink.compileMaterialDictionary( ...
                struct("domain",radia.simulink.makeMaterialSpec( ...
                    Conductivity_S_per_m=1)),MeshFile=meshFile);
            study=radia.simulink.makeFieldStudySpec(Physics="current_flow", ...
                Frequency_Hz=50,DirichletValues=struct("left",1+2i,"right",0));
            contract=radia.simulink.compileFieldStudy(study,materialContract);
            output=string(tempname("C:\temp"))+".json";
            cleanupOutput=onCleanup(@()deleteIfPresent(output));
            radia.simulink.writeFieldStudyRequest(contract,output);
            decoded=jsondecode(fileread(output));
            testCase.verifyEqual(decoded.dirichlet_values.left,[1;2]);

            dcStudy=radia.simulink.makeFieldStudySpec(Physics="current_flow", ...
                DirichletValues=struct("left",0,"right",1));
            dcContract=radia.simulink.compileFieldStudy(dcStudy,materialContract);
            testCase.verifyFalse(isfield(dcContract.request.materials.domain, ...
                "relative_permittivity"));
        end

        function collapsesRepeatedNetgenBoundaryIdsByExactName(testCase)
            path=string(tempname("C:\temp"))+".vol";
            cleanup=onCleanup(@()deleteIfPresent(path));
            text="mesh3d"+newline+"dimension"+newline+"2"+newline+ ...
                "materials"+newline+"1"+newline+"1 domain"+newline+ ...
                "bcnames"+newline+"4"+newline+"1 outer"+newline+ ...
                "2 conductor"+newline+"3 outer"+newline+"4 outer"+newline+ ...
                "endmesh"+newline;
            file=fopen(path,"w","n","UTF-8");fprintf(file,"%s",text);fclose(file);
            inventory=radia.simulink.inspectVolBoundaries(path);

            testCase.verifyEqual(inventory.boundary_names,["outer";"conductor"]);
            testCase.verifyEqual(inventory.boundary_id_groups{1},uint32([1;3;4]));
            testCase.verifyEqual(inventory.boundary_count,uint16(2));
            testCase.verifyEqual(inventory.raw_boundary_id_count,uint16(4));
            testCase.verifyEqual(inventory.raw_boundary_ids,uint32([1;2;3;4]));
            testCase.verifyEqual(inventory.raw_boundary_names, ...
                ["outer";"conductor";"outer";"outer"]);

            materialContract=radia.simulink.compileMaterialDictionary( ...
                struct("domain",radia.simulink.makeMaterialSpec( ...
                    RelativePermittivity=2.5)),MeshFile=path);
            study=radia.simulink.makeFieldStudySpec(Physics="electrostatic", ...
                DirichletValues=struct("outer",0,"conductor",1));
            compiled=radia.simulink.compileFieldStudy(study,materialContract);
            testCase.verifyEqual(compiled.runtime.boundary_count,uint16(4));
            testCase.verifyEqual(compiled.runtime.boundary_id(1:4),uint32([1;2;3;4]));
            testCase.verifyEqual(compiled.runtime.boundary_kind(1:4),uint16(ones(4,1)));
            testCase.verifyEqual(compiled.runtime.boundary_value_real(1:4),[0;1;0;0]);
        end

        function savedLibraryExposesApplicationAndTypedStudyBus(testCase)
            root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
            library=fullfile(root,"matlab","radia_simulink_library.slx");
            name="radia_simulink_library";
            if bdIsLoaded(name),close_system(name,0);end
            load_system(library);
            cleanup=onCleanup(@()closeIfLoaded(name));

            application=name+"/Applications/Field Study";
            coupling=name+"/Coupling/Field Study";
            testCase.verifyEqual(string(get_param(application,"FunctionName")), ...
                "radia_application_sfun");
            testCase.verifySubstring(string(get_param(application,"Parameters")), ...
                "'field'");
            testCase.verifyEqual(string(get_param(coupling+"/Compiled Study Bus", ...
                "OutDataTypeStr")),"Bus: RadiaStudyBus");
            testCase.verifyEqual(string(get_param(coupling+"/study","Port")),"1");
        end
    end
end

function [path,cleanup]=makeMesh()
path=string(tempname("C:\temp"))+".vol";
text="mesh3d"+newline+"dimension"+newline+"2"+newline+ ...
    "materials"+newline+"1"+newline+"1 domain"+newline+ ...
    "bcnames"+newline+"3"+newline+"1 left"+newline+ ...
    "2 right"+newline+"3 outer"+newline+"endmesh"+newline;
file=fopen(path,"w","n","UTF-8");fprintf(file,"%s",text);fclose(file);
cleanup=onCleanup(@()deleteIfPresent(path));
end

function deleteIfPresent(path)
if isfile(path),delete(path);end
end

function closeIfLoaded(name)
if bdIsLoaded(name),close_system(name,0);end
end
