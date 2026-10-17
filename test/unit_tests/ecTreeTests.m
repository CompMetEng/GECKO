%run this test case with the command
%results = runtests('ecTreeTests.m')

function tests = ecTreeTests
    tests = functiontests(localfunctions);
end

function testEcTreeWorkedExample_tc0084(testCase)
    % log10 kcats: global median 2; class 1.1 median 1 (5 values), 1.1.1
    % median 1.5 (4 values). geckopy's test_ec_tree_worked_example has the
    % same numbers.
    m = table_();
    [v, n, depth, ec] = ecTreeEstimate({'1.1.1.1'}, m, 5);
    % 1.1: 0.5*1 + 0.5*2 = 1.5; 1.1.1: 1.5; 1.1.1.1: 3/8*1 + 5/8*1.5.
    verifyEqual(testCase, v, 1.3125, 'AbsTol', 1e-12)
    verifyEqual(testCase, [n, depth], [3, 4])
    verifyEqual(testCase, ec, {'1.1.1.1'})
    [v, n, depth] = ecTreeEstimate({'1.1.1.9'; '1.1.1.-'}, m, 5);
    verifyEqual(testCase, v, [1.5; 1.5], 'AbsTol', 1e-12)
    verifyEqual(testCase, [n, depth], [0 3; 0 3])
    [~, ~, ~, ec] = ecTreeEstimate({'1.1.1.2; 1.1.1.1'; '1.1.1.9;1.1.1.8'}, m, 5);
    verifyEqual(testCase, ec, {'1.1.1.1'; '1.1.1.9'})
    [v, ~, depth] = ecTreeEstimate({'9.9.9.9'; ''}, m, 5);
    verifyEqual(testCase, depth, [0; 0])
    verifyEqual(testCase, v, [2; 2])
    verifyEqual(testCase, ecTreeEstimate({'1.1.1.1'}, m, 10), 403/273, 'AbsTol', 1e-12)
    verifyError(testCase, @() ecTreeEstimate({'1.1.1.1'}, m, 0), 'ecTreeEstimate:k')
end

function testBundledMeasurements_tc0085(testCase)
    m = loadBrendaMeasurements('folder', fetchBrendaAsset());
    verifyTrue(testCase, startsWith(m.release, 'BRENDA release '))
    verifyGreaterThan(testCase, numel(m.ec), 50000)
    verifyTrue(testCase, all(m.kcat > 0))
    verifyFalse(testCase, any(contains(m.ec, '-')))
    % A file that is not gzipped in the release, and a second call that
    % reuses the downloads.
    folder = fetchBrendaAsset({'kcatCombineParams.tsv'});
    verifyTrue(testCase, isfile(fullfile(folder, 'kcatCombineParams.tsv')))
    verifyEqual(testCase, fetchBrendaAsset(), folder)
    verifyError(testCase, @() fetchBrendaAsset({'other.tsv'}), 'fetchBrendaAsset:file')
end

function testKcatsFromECTreeWritesFile_tc0086(testCase)
    [model, adapter, project] = buildProject();
    m = measurementsFor(model);
    kcatList = kcatsFromECTree(model, 'measurements', m, 'modelAdapter', adapter);
    withEc = model.ec.rxns(~cellfun(@isempty, model.ec.eccodes));
    verifyNotEmpty(testCase, withEc)
    verifyEqual(testCase, sort(kcatList.rxns), sort(withEc))
    [~, i] = ismember(kcatList.rxns, model.ec.rxns);
    verifyEqual(testCase, log10(kcatList.kcats), ecTreeEstimate(model.ec.eccodes(i), m, 5), 'AbsTol', 1e-12)
    [~, loaded] = loadKcatSource(model, fullfile(project, 'data', 'kcats', 'ec_tree.tsv'));
    records = loaded.records;
    meta = loaded.header;
    verifyEqual(testCase, meta.parameters, '{"k":5}')
    verifyEqual(testCase, meta.data, 'BRENDA release test')
    verifyEqual(testCase, unique(records.source), {'ec_tree'})
    verifyTrue(testCase, all(~cellfun(@isempty, regexp(records.detail, '^n=\d+; depth=[234]$', 'once'))))
end

function testBrendaEcTreePreset_tc0087(testCase)
    % The default preset: custom.tsv, then BRENDA at origin 1, then the EC
    % tree, then DLKcat.
    [model, adapter, project] = buildProject();
    verifyFalse(testCase, isfield(adapter.params, 'kcat_assignment'))
    writeMeasurements(fullfile(project, 'data', 'measurements.tsv'), measurementsFor(model));
    report = updateKcatSources(model, 'dryRun', true, 'modelAdapter', adapter);
    verifyEqual(testCase, report.source', {'custom','brenda','ec_tree','dlkcat'})
    verifyEqual(testCase, report.status{3}, 'missing')
    fuzzy = kcatsFromBrenda(model, 'aggregation', 'median', 'modelAdapter', adapter);
    kcatsFromECTree(model, 'modelAdapter', adapter);
    [decisions, rep] = collectKcats(model, 'modelAdapter', adapter, 'writeReport', false);
    verifyEqual(testCase, rep.preset, 'brenda_ec_tree')
    origin1 = fuzzy.rxns(fuzzy.origin == 1 & fuzzy.wildcardLvl == 0);
    stage1 = strcmp(decisions.stage, '1');
    fromBrenda = stage1 & startsWith(decisions.source, 'brenda');
    verifyTrue(testCase, all(ismember(decisions.ecRxn(fromBrenda), origin1)))
    verifyTrue(testCase, all(ismember(decisions.source(stage1 & ~fromBrenda), {'custom','ec_tree','DLKcat'})))
    verifyTrue(testCase, any(strcmp(decisions.source, 'ec_tree')))
    % Made from the adapter's measurement table with k = 5, as the preset
    % asks: current. Another k makes it outdated.
    report = updateKcatSources(model, 'dryRun', true, 'modelAdapter', adapter);
    verifyEqual(testCase, report.status{3}, 'current')
    kcatsFromECTree(model, 'k', 2, 'modelAdapter', adapter);
    report = updateKcatSources(model, 'dryRun', true, 'modelAdapter', adapter);
    verifyEqual(testCase, report.status{3}, 'outdated')
end

function writeMeasurements(file, m)
    fid = fopen(file, 'w');
    fprintf(fid, '# %s\nec\tkcat\n', m.release);
    for i = 1:numel(m.ec)
        fprintf(fid, '%s\t%.17g\n', m.ec{i}, m.kcat(i));
    end
    fclose(fid);
end

function m = table_()
    ec = {'1.1.1.1'; '1.1.1.1'; '1.1.1.1'; '1.1.1.2'; '1.1.2.1'; '2.7.1.1'; '2.7.1.1'};
    m = struct('ec', {ec}, 'log10kcat', [0; 1; 2; 3; -1; 2; 2], 'release', 'BRENDA release test');
    m.kcat = 10 .^ m.log10kcat;
end

function m = measurementsFor(model)
    % Two measurements for each EC number of the model.
    ecs = unique(strtrim(split(strjoin(model.ec.eccodes(~cellfun(@isempty, model.ec.eccodes)), ';'), ';')));
    ecs = ecs(~cellfun(@isempty, ecs));
    v = repelem((1:numel(ecs))', 2);
    m = struct('ec', {repelem(ecs, 2)}, 'log10kcat', v, 'kcat', 10 .^ v, 'release', 'BRENDA release test');
end

function [base, adapter, project] = buildProject()
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    project = tempname;
    copyfile(fullfile(adapter.params.path, 'data'), fullfile(project, 'data'));
    adapter.params.path = project;
    model = getGeckoTestModel();
    base = makeEcModel(model, false, adapter);
    base = getECfromGEM(base);
end
