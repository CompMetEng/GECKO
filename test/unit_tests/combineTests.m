%run this test case with the command
%results = runtests('combineTests.m')

function tests = combineTests
    tests = functiontests(localfunctions);
end

function testOkpPresetCombines_tc0088(testCase)
    % CataPro 100 /s (similarity 30: bin <40, bias 0.2, variance 2) and
    % UniKP 10 /s (similarity 50: bin 40-60, bias -0.3, variance 1):
    % (1.8 * 0.5 + 1.3 * 1) / 1.5; sd sqrt(1/1.5 + 0.25) with the okp scale
    % 1. geckopy's test_okp_presets has the same numbers.
    [model, adapter, rxns] = buildProject('okp');
    report = updateKcatSources(model, 'dryRun', true, 'modelAdapter', adapter);
    verifyEqual(testCase, report.source', {'custom', 'okp_catapro', 'okp_unikp'})
    saveKcatSource(model, okpList(rxns(1:3), 100, 30, 'CataPro'), 'okp_catapro', 'modelAdapter', adapter);
    saveKcatSource(model, okpList(rxns(1), 10, 50, 'UniKP'), 'okp_unikp', 'modelAdapter', adapter);
    decisions = collectKcats(model, 'modelAdapter', adapter, 'writeReport', false);
    i = find(strcmp(decisions.ecRxn, rxns{1}));
    verifyEqual(testCase, decisions.kcat(i), 10 ^ ((1.8 * 0.5 + 1.3) / 1.5), 'RelTol', 1e-12)
    verifyEqual(testCase, decisions.source{i}, sprintf('combined (sd=%.2f, n=2)', sqrt(1/1.5 + 0.25)))
    verifyEqual(testCase, decisions.stage{i}, '2')
    verifyEqual(testCase, decisions.candidates(i), 2)
    verifyEqual(testCase, decisions.sdLog10(i), sqrt(1/1.5 + 0.25), 'AbsTol', 1e-12)
    j = find(strcmp(decisions.ecRxn, rxns{2}));
    verifyEqual(testCase, decisions.kcat(j), 10 ^ 1.8, 'RelTol', 1e-12)
end

function testEcTreeOkpPreset_tc0089(testCase)
    % EC tree 10 /s (n=4, depth 4: bin ec4:3-5, bias 0.1, variance 1) and
    % CataPro 100 /s (bin <40): (0.9 + 0.9) / 1.5 = 1.2; the ec_tree+okp
    % scale is 2.
    [model, adapter, rxns] = buildProject('ec_tree_okp');
    saveKcatSource(model, okpList(rxns(1), 100, 30, 'CataPro'), 'okp_catapro', 'modelAdapter', adapter);
    tree = struct('source', 'ec_tree', 'rxns', {rxns(1)}, 'kcats', 10, 'detail', {{'n=4; depth=4'}});
    saveKcatSource(model, tree, 'ec_tree', 'modelAdapter', adapter);
    decisions = collectKcats(model, 'modelAdapter', adapter, 'writeReport', false);
    i = find(strcmp(decisions.ecRxn, rxns{1}));
    verifyEqual(testCase, decisions.kcat(i), 10 ^ 1.2, 'RelTol', 1e-12)
    verifyEqual(testCase, decisions.sdLog10(i), sqrt(1/1.5 + 0.25) * 2, 'AbsTol', 1e-12)
    verifyEqual(testCase, decisions.stage{i}, '2')
end

function testExactHomologDominates_tc0098(testCase)
    % EC tree 10 /s (ec4:3-5: bias 0.1, variance 1) and a homolog at 100%
    % identity measured at 100 /s (exact: bias 0, variance 0.01):
    % (0.9 * 1 + 2 * 100) / 101. geckopy's test_exact_homolog_dominates
    % has the same numbers.
    [model, adapter, rxns] = buildProject('ec_tree_homologs');
    tree = struct('source', 'ec_tree', 'rxns', {rxns(1)}, 'kcats', 10, 'detail', {{'n=4; depth=4'}});
    saveKcatSource(model, tree, 'ec_tree', 'modelAdapter', adapter);
    hom = struct('source', 'homologs', 'rxns', {rxns(1)}, 'kcats', 100, 'substrates', {{'m1'}}, ...
        'detail', {{'identity=1.000; n=1'}});
    saveKcatSource(model, hom, 'homologs', 'modelAdapter', adapter);
    decisions = collectKcats(model, 'modelAdapter', adapter, 'writeReport', false);
    i = find(strcmp(decisions.ecRxn, rxns{1}));
    verifyEqual(testCase, decisions.kcat(i), 10 ^ ((0.9 + 200) / 101), 'RelTol', 1e-12)
end

function testCombineErrors_tc0090(testCase)
    % An OKP method without combination parameters, and an unknown method.
    [model, adapter, rxns] = buildProject('okp');
    adapter.params.kcat_sources.okp.methods = {'EITLEM'};
    saveKcatSource(model, okpList(rxns(1), 100, 30, 'EITLEM'), 'okp_eitlem', 'modelAdapter', adapter);
    verifyError(testCase, @() collectKcats(model, 'modelAdapter', adapter, 'writeReport', false), ...
        'collectKcats:combine')
    adapter.params.kcat_sources.okp.methods = {'CataPro', 'NoSuchMethod'};
    verifyError(testCase, @() kcatAssignmentConfig(adapter.params), 'kcatAssignment:config')
end

function list = okpList(rxns, kcat, similarity, method)
    n = numel(rxns);
    list = struct('source', 'OpenKineticsPredictor', 'rxns', {rxns(:)}, 'kcats', repmat(kcat, n, 1), ...
        'kcatSource', {repmat({['OKP-' method]}, n, 1)}, ...
        'detail', {repmat({sprintf('similarity=%g', similarity)}, n, 1)});
end

function [model, adapter, rxns] = buildProject(preset)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    project = tempname;
    copyfile(fullfile(adapter.params.path, 'data'), fullfile(project, 'data'));
    delete(fullfile(project, 'data', 'customKcats.tsv'));
    adapter.params.path = project;
    adapter.params.kcat_assignment = struct('preset', preset);
    adapter.params.kcat_sources = struct('okp', struct('methods', {{'CataPro', 'UniKP'}}));
    % TestGEMAdapter's BRENDA folder is the project's data folder.
    fid = fopen(fullfile(project, 'data', 'kcatCombineParams.tsv'), 'w');
    fprintf(fid, ['# Combination parameters for tests\n# measurement_sd: 0.5\n# min_variance: 0.01\n' ...
        '# interval_scale ec_tree+okp: 2\n# interval_scale okp: 1\nsource\tbin\tbias\tvariance\tn\n' ...
        'ec_tree\tec4:3-5\t0.1\t1\t50\nec_tree\tec3\t0\t4\t50\nokp_catapro\t<40\t0.2\t2\t50\n' ...
        'okp_catapro\t>80\t0\t0.5\t50\nokp_unikp\t40-60\t-0.3\t1\t50\n' ...
        'homologs\t90\t0\t1\t50\nhomologs\texact\t0\t0.01\t0\n']);
    fclose(fid);
    model = makeEcModel(getGeckoTestModel(), false, adapter);
    model = getECfromGEM(model);
    rxns = model.ec.rxns(any(model.ec.rxnEnzMat, 2));
end
