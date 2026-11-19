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
end

function testKcatsForReactionsWithoutOkp_tc0751(testCase)
    [model, adapter] = buildProject();
    i = find(~cellfun(@isempty, model.ec.eccodes) & any(model.ec.rxnEnzMat, 2), 1);
    [~, summary] = kcatsForReactions(model, i, 'measurements', measurementsFor(model), 'modelAdapter', adapter);
    verifyEqual(testCase, summary.ecRxn, model.ec.rxns(i))
    verifyFalse(testCase, isnan(summary.src_ec_tree))
    verifyError(testCase, @() kcatsForReactions(model, i, 'okp', 'NoSuchMethod', 'modelAdapter', adapter), ...
        'kcatsForReactions:okp')
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
        verifyTrue(testCase, contains(err.message, 'OKP_job_kinform_h_query.txt'))
    end
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
