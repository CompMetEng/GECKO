%run this test case with the command
%results = runtests('kcatContractTests.m')

function tests = kcatContractTests
    tests = functiontests(localfunctions);
end

function testCollectKcatsCountsEachSourceFileOnce_tc0300(testCase)
    % database_top and database_bottom both read brenda.tsv, which has an
    % empty log10_sd (NaN) column: the file is one list, so R3 has the two
    % values of the file as candidates, not four.
    [model, adapter, project] = kcatProject();
    adapter.params.kcat_assignment = struct('isozymes', false, 'standard', false, 'stage', ...
        struct('action', 'select', 'sources', {{'database_top', 'database_bottom'}}));
    [baseRxns, enz] = getEcRxnKeys(model);
    i = find(strcmp(model.ec.rxns, 'R3'));
    writeLines(fullfile(project, 'data', 'kcats', 'brenda.tsv'), {'# source: brenda', ...
        sprintf('rxn\tenzymes\tsubstrate\tkcat\tlog10_sd\tsource\tdetail'), ...
        sprintf('%s\t%s\tm1\t5\t\tbrenda\twc=0; origin=1', baseRxns{i}, enz{i}), ...
        sprintf('%s\t%s\tm1\t20\t\tbrenda\twc=0; origin=1', baseRxns{i}, enz{i})});
    decisions = collectKcats(model, 'modelAdapter', adapter, 'writeReport', false);
    row = strcmp(decisions.ecRxn, 'R3');
    verifyEqual(testCase, decisions.kcat(row), 20)
    verifyEqual(testCase, decisions.candidates(row), 2)
    verifyEqual(testCase, decisions.spreadLog10(row), log10(4), 'AbsTol', 1e-12)
end

function testMergeKcatsMatchesOkpTaggedSources_tc0301(testCase)
    % The OKP- tag is dropped for matching: OKP-BRENDA is an exact database
    % hit and OKP-CataPro matches the priority entry catapro; the rows keep
    % their labels and win over DLKcat.
    list = struct('source', 'x', 'rxns', {{'A'; 'A'; 'B'; 'B'}}, 'kcats', [3; 9; 4; 9], ...
        'kcatSource', {{'OKP-BRENDA'; 'DLKcat'; 'OKP-CataPro'; 'DLKcat'}});
    merged = mergeKcats({list}, {'database_exact', 'catapro', 'dlkcat'});
    verifyEqual(testCase, merged.rxns, {'A'; 'B'})
    verifyEqual(testCase, merged.kcats, [3; 4])
    verifyEqual(testCase, merged.kcatSource, {'OKP-BRENDA'; 'OKP-CataPro'})
end

function [model, adapter, project] = kcatProject()
    adapter = ModelAdapterManager.getAdapter(fullfile(findGECKOroot, 'test', 'unit_tests', 'ecTestGEM', 'TestGEMAdapter.m'));
    project = tempname;
    copyfile(fullfile(adapter.params.path, 'data'), fullfile(project, 'data'));
    delete(fullfile(project, 'data', 'customKcats.tsv'));
    mkdir(fullfile(project, 'data', 'kcats'));
    adapter.params.path = project;
    model = makeEcModel(getGeckoTestModel(), false, adapter);
    model = getECfromGEM(model);
end

function writeLines(file, lines)
    fid = fopen(file, 'w');
    fprintf(fid, '%s\n', lines{:});
    fclose(fid);
end
