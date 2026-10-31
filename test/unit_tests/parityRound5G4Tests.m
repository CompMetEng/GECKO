%run this test case with the command
%results = runtests('parityRound5G4Tests.m')

function tests = parityRound5G4Tests
    tests = functiontests(localfunctions);
end

function testBestParsimoniousSkipsNaNRmse_tc0600(testCase)
    % A NaN rmse never sets the best rmse and is never chosen.
    points = struct('threshold',{0,3,6}, 'nChanged',{4000,3000,1000}, ...
        'rmse',{NaN,0.905,1.4});
    point = bestParsimonious(points);
    verifyEqual(testCase, point.nChanged, 3000)
end

function testFillCarbonNumWithoutFormulasWarnsAndWeightsOne_tc0601(testCase)
    model = exchangeToyModel();
    verifyWarning(testCase, @() fillCarbonNum(model), 'GECKO:fillCarbonNum:noFormulas')
    model.metFormulas = {'';''};
    verifyWarning(testCase, @() fillCarbonNum(model), 'GECKO:fillCarbonNum:noFormulas')
    warning('off','GECKO:fillCarbonNum:noFormulas')
    cleanup = onCleanup(@() warning('on','GECKO:fillCarbonNum:noFormulas'));
    withoutField = rmfield(model,'metFormulas');
    out = fillCarbonNum(withoutField);
    verifyEqual(testCase, out.excarbon(strcmp(out.rxns,'EX_a')), 1)
    verifyEqual(testCase, out.excarbon(strcmp(out.rxns,'EX_b')), 1)
end

function testFillCarbonNumUsesFormulaCarbons_tc0602(testCase)
    model = exchangeToyModel();
    model.metFormulas = {'C6H12O6';'O2'};
    out = fillCarbonNum(model);
    verifyEqual(testCase, out.excarbon(strcmp(out.rxns,'EX_a')), 6)
end

function testSaveEcModelWritesDefaultDescriptionOnlyWhenNone_tc0603(testCase)
    adapter = ModelAdapterManager.getAdapter(fullfile(findGECKOroot,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    project = tempname;
    testCase.addTeardown(@() rmdir(project, 's'));
    copyfile(fullfile(adapter.params.path, 'data'), fullfile(project, 'data'));
    adapter.params.path = project;
    mkdir(fullfile(project, 'models'));
    model = getECfromGEM(makeEcModel(getGeckoTestModel(), false, adapter));
    model.annotation.description = '';
    saveEcModel(model, 'default.yml', adapter);
    text = fileread(fullfile(project, 'models', 'default.yml'));
    verifyNotEmpty(testCase, strfind(text, ['- description: Enzyme-constrained model of ' model.id]))
    model.annotation.description = 'my own description';
    saveEcModel(model, 'own.yml', adapter);
    text = fileread(fullfile(project, 'models', 'own.yml'));
    verifyNotEmpty(testCase, strfind(text, '- description: my own description'))
    verifyEmpty(testCase, strfind(text, 'Enzyme-constrained model of'))
end

function model = exchangeToyModel()
    model = struct();
    model.id = 'toy';
    model.rxns = {'EX_a';'EX_b'};
    model.rxnNames = {'exchange a';'exchange b'};
    model.mets = {'a';'b'};
    model.metNames = {'a';'b'};
    model.S = sparse([-1 0;0 -1]);
    model.lb = [-1000;-1000];
    model.ub = [1000;1000];
    model.rev = [1;1];
    model.c = [0;0];
    model.b = [0;0];
    model.comps = {'c'};
    model.compNames = {'cytosol'};
    model.metComps = [1;1];
end
