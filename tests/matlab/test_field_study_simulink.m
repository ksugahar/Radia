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

        function compilesHarmonicWindingAndRejectsNonlinearMaterial(testCase)
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
