%run this test case with the command
%results = runtests('kcatEvidenceTests.m')

function tests = kcatEvidenceTests
    tests = functiontests(localfunctions);
end

function testEvidenceListsEverySource_tc0748(testCase)
    [model, adapter] = buildProject();
    rows = find(~cellfun(@isempty, model.ec.eccodes) & any(model.ec.rxnEnzMat, 2), 3);
    sel = model.ec.rxns(rows);
    a = struct('rxns', {sel}, 'kcats', [1; 2; 4]);
    saveKcatSource(model, a, 'srcA', 'modelAdapter', adapter);
    b = struct('rxns', {sel(1)}, 'kcats', 8);
    saveKcatSource(model, b, 'srcB', 'modelAdapter', adapter);
    [evidence, summary, idx] = kcatEvidence(model, 'rxns', sel, 'sources', {'srcA', 'srcB'}, ...
        'modelAdapter', adapter);
    verifyEqual(testCase, idx(:), rows(:))
    verifyEqual(testCase, summary.ecRxn, sel(:))
    verifyEqual(testCase, summary.src_srcA, [1; 2; 4])
    verifyEqual(testCase, summary.src_srcB(1), 8)
    verifyEqual(testCase, summary.Properties.VariableNames, {'ecRxn','kcat','tuned','fold','prior', ...
        'origin','enzymes','eccodes','src_srcA','src_srcB'})
    [~, summaryOrdered] = kcatEvidence(model, 'rxns', sel, 'sources', {'srcB', 'srcA'}, 'modelAdapter', adapter);
    verifyEqual(testCase, summaryOrdered.Properties.VariableNames(end-1:end), {'src_srcA','src_srcB'})
    verifyEqual(testCase, evidence.Properties.VariableNames, {'ecRxn','kcat','source','substrate','detail','assigned'})
    verifyTrue(testCase, all(isnan(summary.src_srcB(2:3))))
    verifyEqual(testCase, height(evidence), 4)
    verifyEqual(testCase, summary.kcat, model.ec.kcat(rows(:)))
    % A source that is not listed is not read.
    [~, summary] = kcatEvidence(model, 'rxns', rows, 'sources', {'srcB'}, 'modelAdapter', adapter);
    verifyFalse(testCase, ismember('src_srcA', summary.Properties.VariableNames))
    verifyError(testCase, @() kcatEvidence(model, 'rxns', {'nope'}, 'modelAdapter', adapter), 'kcatEvidence:rxns')
    verifyError(testCase, @() kcatEvidence(model, 'compute', {'x'}, 'modelAdapter', adapter), 'kcatEvidence:compute')
end

function testEvidenceOfTunedKcat_tc0749(testCase)
    [model, adapter] = buildProject();
    i = find(any(model.ec.rxnEnzMat, 2), 1);
    prior = 0.5;
    model.ec.kcat(i) = 5;
    model.ec.source{i} = 'custom | tuned=sensitivity (prior=0.5)';
    a = struct('rxns', {model.ec.rxns(i)}, 'kcats', prior);
    saveKcatSource(model, a, 'srcA', 'modelAdapter', adapter);
    tunedKcats = struct('rxns', {model.ec.rxns(i)});
    if model.ec.geckoLight
        tunedKcats.rxns = extractAfter(model.ec.rxns(i), 4);
    end
    [evidence, summary] = kcatEvidence(model, 'tunedKcats', tunedKcats, 'sources', {'srcA'}, ...
        'modelAdapter', adapter);
    verifyEqual(testCase, summary.ecRxn, model.ec.rxns(i))
    verifyTrue(testCase, summary.tuned)
    verifyEqual(testCase, [summary.kcat, summary.prior, summary.fold], [5, 0.5, 10], 'AbsTol', 1e-12)
    verifyEqual(testCase, summary.origin, {'custom'})
    verifyTrue(testCase, evidence.assigned)
    % The prior in ec.source has six significant digits: a source kcat that
    % equals the assigned one up to that rounding is still the assigned one.
    model.ec.kcat(i) = 123.456789;
    model.ec.source{i} = 'custom | tuned=sensitivity (prior=12.3457)';
    saveKcatSource(model, struct('rxns', {model.ec.rxns(i)}, 'kcats', 12.3456789), 'srcA', 'modelAdapter', adapter);
    evidence = kcatEvidence(model, 'rxns', model.ec.rxns(i), 'sources', {'srcA'}, 'modelAdapter', adapter);
    verifyTrue(testCase, evidence.assigned)
    % A tuning status without a prior is tuned, and has no fold.
    model.ec.source{i} = 'custom | tuned=cmaes';
    [~, summary] = kcatEvidence(model, 'rxns', i, 'sources', {}, 'modelAdapter', adapter);
    verifyTrue(testCase, summary.tuned && isnan(summary.prior) && isnan(summary.fold))
    model.ec.kcat(i) = 5;
    model.ec.source{i} = 'custom | tuned=sensitivity (prior=0.5)';
    file = [tempname '.tsv'];
    kcatEvidence(model, 'tunedKcats', tunedKcats, 'sources', {'srcA'}, 'saveAs', file, 'modelAdapter', adapter);
    saved = readtable(file, 'FileType', 'text', 'Delimiter', '\t');
    verifyEqual(testCase, saved.fold, 10, 'AbsTol', 1e-9)
end

function testEvidenceByEnzymeAndEcTree_tc0750(testCase)
    [model, adapter] = buildProject();
    e = find(any(model.ec.rxnEnzMat ~= 0, 1), 1);
    expected = find(model.ec.rxnEnzMat(:, e) ~= 0);
    m = measurementsFor(model);
    [evidence, summary, idx] = kcatEvidence(model, 'enzymes', model.ec.enzymes(e), 'sources', {}, ...
        'compute', {'ec_tree'}, 'measurements', m, 'modelAdapter', adapter);
    verifyEqual(testCase, idx(:), expected(:))
    withEc = ~cellfun(@isempty, model.ec.eccodes(expected));
    verifyEqual(testCase, sum(~isnan(summary.src_ec_tree)), sum(withEc))
    verifyEqual(testCase, unique(evidence.source), {'ec_tree'})
    verifyError(testCase, @() kcatEvidence(model, 'enzymes', {'X00000'}, 'modelAdapter', adapter), 'kcatEvidence:enzymes')
    % A computed EC-tree estimate replaces a file of that name.
    e = find(~cellfun(@isempty, model.ec.eccodes), 1);
    saveKcatSource(model, struct('rxns', {model.ec.rxns(e)}, 'kcats', 1e9), 'ec_tree', 'modelAdapter', adapter);
    [~, summary] = kcatEvidence(model, 'rxns', e, 'compute', {'ec_tree'}, 'measurements', m, 'modelAdapter', adapter);
    verifyEqual(testCase, summary.src_ec_tree, 10 .^ ecTreeEstimate(model.ec.eccodes(e), m, 5), 'RelTol', 1e-9)
end

function testEmptySelectionSelectsNothing_tc0755(testCase)
    [model, adapter] = buildProject();
    for args = {{'rxns', {}}, {'rxns', []}, {'enzymes', {}}, {'tunedKcats', []}, {'tunedKcats', struct('rxns', {{}})}}
        [~, summary, rows] = kcatEvidence(model, args{1}{:}, 'sources', {}, 'modelAdapter', adapter);
        verifyEmpty(testCase, rows)
        verifyEqual(testCase, height(summary), 0)
    end
    [~, ~, rows] = kcatEvidence(model, 'sources', {}, 'modelAdapter', adapter);
    verifyEqual(testCase, numel(rows), numel(model.ec.rxns))
    verifyError(testCase, @() kcatEvidence(model, 'rxns', 10^6, 'modelAdapter', adapter), 'kcatEvidence:rxns')
    verifyError(testCase, @() kcatEvidence(model, 'rxns', true(2, 1), 'modelAdapter', adapter), 'kcatEvidence:rxns')
    verifyError(testCase, @() kcatsForReactions(model, {}, 'modelAdapter', adapter), 'kcatsForReactions:none')
    verifyError(testCase, @() kcatsForReactions(model, 1, 'measurement', 1, 'modelAdapter', adapter), 'kcatsForReactions:arguments')
end

function testConventionalNameSelectsTheIsozymes_tc0756(testCase)
    [model, adapter] = buildProject(true);
    base = getEcRxnKeys(model);
    [names, ~, g] = unique(base);
    [~, most] = max(accumarray(g, 1));
    name = names{most};
    expected = find(strcmp(base, name));
    verifyGreaterThan(testCase, numel(expected), 1)
    [~, ~, rows] = kcatEvidence(model, 'rxns', name, 'sources', {}, 'modelAdapter', adapter);
    verifyEqual(testCase, rows(:), expected(:))
    [~, ~, rows] = kcatEvidence(model, 'tunedKcats', struct('rxns', {{name}}), 'sources', {}, 'modelAdapter', adapter);
    verifyEqual(testCase, rows(:), expected(:))
end

function testKcatsForReactionsWithoutOkp_tc0751(testCase)
    [model, adapter] = buildProject();
    i = find(~cellfun(@isempty, model.ec.eccodes) & any(model.ec.rxnEnzMat, 2), 1);
    [~, summary] = kcatsForReactions(model, i, 'measurements', measurementsFor(model), 'modelAdapter', adapter);
    verifyEqual(testCase, summary.ecRxn, model.ec.rxns(i))
    verifyFalse(testCase, isnan(summary.src_ec_tree))
    verifyError(testCase, @() kcatsForReactions(model, i, 'okp', 'NoSuchMethod', 'modelAdapter', adapter), ...
        'kcatsForReactions:okp')
    % An empty okp asks for no OKP method: no job is submitted and no file is written.
    [~, summary] = kcatsForReactions(model, i, 'okp', [], 'measurements', measurementsFor(model), 'modelAdapter', adapter);
    verifyEqual(testCase, height(summary), 1)
end

function testOkpJobTagAndNames_tc0752(testCase)
    [model, adapter, project] = buildProject();
    verifyError(testCase, @() kcatsFromOKP(model, 'method', 'KinForm-H', 'ecRxns', {'nope'}, ...
        'modelAdapter', adapter), 'kcatsFromOKP:ecRxns')
    % A tagged job has its own files: nothing is stored for it, whatever is stored for the method.
    writeSharedJob(project);
    verifyError(testCase, @() kcatsFromOKP(model, 'method', 'KinForm-H', 'jobTag', 'query', ...
        'useStored', true, 'modelAdapter', adapter), 'kcatsFromOKP:noStored')
    try
        kcatsFromOKP(model, 'method', 'KinForm-H', 'jobTag', 'query', 'useStored', true, 'modelAdapter', adapter);
    catch err
        verifyTrue(testCase, contains(err.message, 'OKP_job_kinform_h__query.txt'))
    end
    verifyError(testCase, @() kcatsFromOKP(model, 'method', 'KinForm-H', 'jobTag', '!!', 'useStored', true, ...
        'modelAdapter', adapter), 'kcatsFromOKP:jobTag')
    verifyError(testCase, @() kcatsFromOKP(model, 'method', 'KinForm-H', 'ecRxns', 10^6, ...
        'modelAdapter', adapter), 'kcatsFromOKP:ecRxns')
    verifyError(testCase, @() kcatsFromOKP(model, 'method', 'KinForm-H', 'ecRxns', {}, ...
        'modelAdapter', adapter), 'kcatsFromOKP:ecRxns')
    verifyError(testCase, @() kcatsFromOKP(model, 'method', 'KinForm-H', 'ecRxns', true(2, 1), ...
        'modelAdapter', adapter), 'kcatsFromOKP:ecRxns')
end

function writeSharedJob(project)
    fid = fopen(fullfile(project, 'data', 'OKP_job_kinform_h.txt'), 'w');
    fprintf(fid, 'jobId: abc\nmethod: KinForm-H\nresultOf: abc\n');
    fclose(fid);
    fid = fopen(fullfile(project, 'data', 'OKP_output_kinform_h.csv'), 'w');
    fprintf(fid, 'x\n');
    fclose(fid);
end

function m = measurementsFor(model)
    % Two measurements for each EC number of the model.
    ecs = unique(strtrim(split(strjoin(model.ec.eccodes(~cellfun(@isempty, model.ec.eccodes)), ';'), ';')));
    ecs = ecs(~cellfun(@isempty, ecs));
    v = repelem((1:numel(ecs))', 2);
    m = struct('ec', {repelem(ecs, 2)}, 'log10kcat', v, 'kcat', 10 .^ v, 'release', 'BRENDA release test');
end

function [base, adapter, project] = buildProject(light)
    if nargin < 1
        light = false;
    end
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    project = tempname;
    copyfile(fullfile(adapter.params.path, 'data'), fullfile(project, 'data'));
    adapter.params.path = project;
    model = getGeckoTestModel();
    base = makeEcModel(model, light, adapter);
    base = getECfromGEM(base);
end
