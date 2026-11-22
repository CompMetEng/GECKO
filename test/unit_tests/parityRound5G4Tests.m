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

function testLoadEvotuneDataReadsTheExperimentFiles_tc0601(testCase)
    adapter = ModelAdapterManager.getAdapter(fullfile(findGECKOroot,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    project = tempname;
    mkdir(fullfile(project, 'data'));
    testCase.addTeardown(@() rmdir(project, 's'));
    adapter.params.path = project;
    tab = char(9);
    writeLines(fullfile(project, 'data', 'fluxData.tsv'), ...
        {['Condition' tab 'Ptot' tab 'grRate' tab 'glucose (EX_glc)' tab 'weight'], ...
        ['c1' tab '0.5' tab '0.4' tab '-10' tab '0'], ...
        ['c2' tab '0.5' tab '0.3' tab '-8' tab '1']});
    writeLines(fullfile(project, 'data', 'maxGrowth.tsv'), ...
        {['Condition' tab 'Ptot' tab 'grRate' tab 'glucose (EX_glc)'], ...
        ['glucose' tab 'NaN' tab '0.41' tab '-1000']});
    writeLines(fullfile(project, 'data', 'zeroFlux.tsv'), {'Rxns', 'EX_a', 'EX_b'});
    experimentData = loadEvotuneData('modelAdapter', adapter);
    verifyEqual(testCase, experimentData.fluxData.conds, {'c1'; 'c2'})
    verifyEqual(testCase, experimentData.fluxData.weight, [0; 1])
    verifyEqual(testCase, experimentData.maxGrowth.exchRxnIDs, {'EX_glc'})
    verifyEqual(testCase, experimentData.zeroFlux, {'EX_a'; 'EX_b'})
end

function testLoadEvotuneDataWithoutFilesIsEmpty_tc0602(testCase)
    adapter = ModelAdapterManager.getAdapter(fullfile(findGECKOroot,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    experimentData = loadEvotuneData('modelAdapter', adapter);
    verifyEmpty(testCase, experimentData.fluxData)
    verifyEmpty(testCase, experimentData.maxGrowth)
    verifyEmpty(testCase, experimentData.zeroFlux)
end

function testLoadEvotuneDataRejectsTheLegacyWeightColumn_tc0604(testCase)
    adapter = ModelAdapterManager.getAdapter(fullfile(findGECKOroot,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    project = tempname;
    mkdir(fullfile(project, 'data'));
    testCase.addTeardown(@() rmdir(project, 's'));
    adapter.params.path = project;
    tab = char(9);
    writeLines(fullfile(project, 'data', 'fluxData.tsv'), ...
        {['Condition' tab 'Ptot' tab 'grRate' tab 'glucose (EX_glc)' tab 'evotuneRMSEweight'], ...
        ['c1' tab '0.5' tab '0.4' tab '-10' tab '1']});
    verifyError(testCase, @() loadEvotuneData('modelAdapter', adapter), 'loadEvotuneData:legacyWeightColumn')
end

function testEvotuneScoreMaximisesGrowthWhateverTheModelObjective_tc0605(testCase)
    % The score uses the growth reaction of the adapter as objective: a model whose
    % objective is another reaction scores like one whose objective is growth, and
    % no objective warning is issued.
    adapter = ModelAdapterManager.getAdapter(fullfile(findGECKOroot,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    adapter.params.bioRxn = 'biomass';
    project = tempname;
    mkdir(fullfile(project, 'data'));
    testCase.addTeardown(@() rmdir(project, 's'));
    adapter.params.path = project;
    tab = char(9);
    writeLines(fullfile(project, 'data', 'maxGrowth.tsv'), ...
        {['Condition' tab 'Ptot' tab 'grRate' tab 'glucose (EX_glc)'], ...
        ['glucose' tab 'NaN' tab '0.41' tab '-10']});
    experimentData = loadEvotuneData('modelAdapter', adapter);
    model = struct();
    model.id = 'toy';
    model.rxns = {'EX_glc';'biomass';'EX_x'};
    model.rxnNames = model.rxns;
    model.mets = {'glc_e'};
    model.metNames = {'glucose'};
    model.metFormulas = {'C6H12O6'};
    model.metComps = 1;
    model.comps = {'c'};
    model.compNames = {'cytosol'};
    model.S = sparse([-1 -1 -1]);
    model.lb = [-1000; 0; 0];
    model.ub = [0; 1000; 1000];
    model.rev = [1; 0; 0];
    model.b = 0;
    model.genes = {};
    model.grRules = repmat({''}, 3, 1);
    model.rxnGeneMat = sparse(3, 0);
    modelGrowth = model;
    modelGrowth.c = [0; 1; 0];
    modelOther = model;
    modelOther.c = [0; 0; 1];
    lastwarn('');
    scoreGrowth = evotuneScore(modelGrowth, experimentData, 'modelAdapter', adapter);
    scoreOther = evotuneScore(modelOther, experimentData, 'modelAdapter', adapter);
    [~, warningId] = lastwarn;
    verifyEqual(testCase, scoreOther, scoreGrowth)
    verifyNotEqual(testCase, warningId, 'simulateConditions:objectiveNotGrowth')
    % the score of the model that is passed in does not change its objective
    verifyEqual(testCase, modelOther.c, [0; 0; 1])
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
    model.annotation.description = '';
    model.annotation.note = 'loaded description';
    saveEcModel(model, 'note.yml', adapter);
    text = fileread(fullfile(project, 'models', 'note.yml'));
    verifyNotEmpty(testCase, strfind(text, '- note: loaded description'))
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

function file = writeLines(file, lines)
fid = fopen(file, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end
