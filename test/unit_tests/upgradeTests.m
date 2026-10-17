%run this test case with the command
%results = runtests('upgradeTests.m')

function tests = upgradeTests
    tests = functiontests(localfunctions);
end

function testConvertCustomKcats_tc0091(testCase)
    % The ecTestGEM customKcats.tsv: P3 200 and P1 + P2 100 (proteins only),
    % then R2_REV and R5 50 (reactions only), which replaces P1 + P2 on
    % R2_REV, as applyCustomKcats does. One row per reaction, keyed by
    % conventional reaction and enzyme set.
    [model, adapter, project] = buildProject();
    [file, report] = convertCustomKcats(model, 'modelAdapter', adapter);
    verifyEqual(testCase, file, fullfile(project, 'data', 'kcats', 'custom.tsv'))
    text = strsplit(strtrim(fileread(file)), newline);
    verifyTrue(testCase, any(strcmp(text, '# kind: curated')))
    verifyTrue(testCase, any(strcmp(text, '# source: custom')))
    [~, loaded] = loadKcatSource(model, file);
    rows = strcat(loaded.records.rxn, '|', loaded.records.enzymes, '|', ...
        cellstr(string(loaded.records.kcat)));
    verifyTrue(testCase, any(strcmp(rows, 'R2_REV|P1 + P2|50')))
    verifyTrue(testCase, any(strcmp(rows, 'R2|P1 + P2|100')))
    verifyEqual(testCase, report.rows, numel(rows))
    verifyError(testCase, @() convertCustomKcats(model, 'modelAdapter', adapter), 'convertCustomKcats:exists')
end

function testConvertedCustomKcatsGiveSameKcats_tc0092(testCase)
    % The ecTestGEM customKcats.tsv: assignKcats with data/customKcats.tsv
    % and with the converted data/kcats/custom.tsv give the same kcats.
    [model, adapter, project] = buildProject();
    adapter.params.kcat_assignment = struct('preset', 'gecko_median');
    kcatsFromBrenda(model, 'modelAdapter', adapter);
    a = collectKcats(model, 'modelAdapter', adapter, 'writeReport', false);
    convertCustomKcats(model, 'modelAdapter', adapter);
    b = collectKcats(model, 'modelAdapter', adapter, 'writeReport', false);
    verifyEqual(testCase, b.kcat, a.kcat)
    verifyEqual(testCase, strcmp(b.source, 'custom'), strcmp(a.source, 'custom'))
    verifyTrue(testCase, any(strcmp(b.source, 'custom')))
    verifyTrue(testCase, isfile(fullfile(project, 'data', 'kcats', 'custom.tsv')))
end

function testUpgradeEcModel_tc0093(testCase)
    % A model with GECKO 3 protein reactions is turned forward; a GECKO 4
    % model is left as it is.
    [model] = buildProject();
    old = gecko3(model);
    verifyNotEqual(testCase, old.S, model.S)
    [up, changed] = upgradeEcModel(old);
    verifyEqual(testCase, up.S, model.S)
    verifyEqual(testCase, [up.lb up.ub], [model.lb model.ub])
    verifyEqual(testCase, up.c, model.c)
    verifyEqual(testCase, numel(changed), sum(startsWith(model.rxns, 'usage_prot_')) + 1)
    [again, changed] = upgradeEcModel(up);
    verifyEmpty(testCase, changed)
    verifyEqual(testCase, again.S, model.S)
end

function testUpgradeGECKOproject_tc0094(testCase)
    % A GECKO 3 project: adapter class, customKcats.tsv, a saved GECKO 3
    % ecModel and a script calling a renamed function.
    [model, ~, project] = buildProject();
    copyfile(fullfile(findGECKOroot, 'test', 'unit_tests', 'ecTestGEM', 'TestGEMAdapter.m'), project);
    mkdir(fullfile(project, 'models'));
    copyfile(fullfile(findGECKOroot, 'test', 'unit_tests', 'ecTestGEM', 'models', 'testModel.xml'), ...
        fullfile(project, 'models'));
    writeYAMLmodel(gecko3(model), fullfile(project, 'models', 'ecModel.yml'));
    writeLines(fullfile(project, 'run.m'), {'u = enzymeUsage(ecModel);', 'x = getEnzymeUsage(ecModel);'});
    % DLKcat output for the test model, as a stored prediction.
    writeLines(fullfile(project, 'data', 'DLKcat.tsv'), {sprintf('rxns\tgenes\tsubstrates\tsmiles\tsequence\tkcat'), ...
        sprintf('R2_EXP_1\tG1\tm1\tC\tMRAL\t12.5'), sprintf('R3\tG4\tm1\tC\tMDFM\t3.25')});

    dry = upgradeGECKOproject(project);
    verifyEqual(testCase, dry.step(1:5)', {'adapter', 'custom kcats', 'predictions', 'saved ecModel', 'project format'})
    verifyEqual(testCase, dry.status(1:5)', {'to do', 'skipped', 'skipped', 'to do', 'to do'})
    verifyFalse(testCase, isfile(fullfile(project, 'model_adapter.toml')))
    scripts = dry.note(strcmp(dry.step, 'script'));
    verifyEqual(testCase, scripts, {'run.m:1: enzymeUsage -> getEnzymeUsage'})

    done = upgradeGECKOproject(project, 'dryRun', false, 'ecModel', model);
    verifyEqual(testCase, done.status(strcmp(done.step, 'predictions')), {'done'})
    verifyTrue(testCase, isfile(fullfile(project, 'data', 'kcats', 'dlkcat.tsv')))
    verifyTrue(testCase, all(strcmp(done.status(strcmp(done.step, 'custom kcats')), 'done')))
    toml = fileread(fullfile(project, 'model_adapter.toml'));
    verifyTrue(testCase, startsWith(toml, 'project_format = 4'))
    verifyTrue(testCase, isfile(fullfile(project, 'data', 'kcats', 'custom.tsv')))
    legacy = fullfile(project, 'legacy', char(datetime('now', 'Format', 'yyyy-MM-dd')));
    verifyTrue(testCase, isfile(fullfile(legacy, 'data', 'customKcats.tsv')))
    verifyTrue(testCase, isfile(fullfile(legacy, 'models', 'ecModel.yml')))
    saved = readYAMLmodel(fullfile(project, 'models', 'ecModel.yml'));
    [~, changed] = upgradeEcModel(saved);
    verifyEmpty(testCase, changed)
    params = loadAdapterParams(fullfile(project, 'model_adapter.toml'));
    verifyEqual(testCase, params.projectFormat, 4)

    again = upgradeGECKOproject(project);
    verifyTrue(testCase, all(ismember(again.status, {'by hand'})))
end

function old = gecko3(model)
    % The model as GECKO 3 defined it: usage and pool exchange reversed.
    r = startsWith(model.rxns, 'usage_prot_') | strcmp(model.rxns, 'prot_pool_exchange');
    old = model;
    old.S(:, r) = -model.S(:, r);
    old.lb(r) = -model.ub(r);
    old.ub(r) = -model.lb(r);
    old.c(r) = -model.c(r);
end

function [model, adapter, project] = buildProject()
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    project = tempname;
    copyfile(fullfile(adapter.params.path, 'data'), fullfile(project, 'data'));
    adapter.params.path = project;
    model = makeEcModel(getGeckoTestModel(), false, adapter);
    model = getECfromGEM(model);
end

function writeLines(file, lines)
    fid = fopen(file, 'w');
    fprintf(fid, '%s\n', lines{:});
    fclose(fid);
end
