%run this test case with the command
%results = runtests('homologsTests.m')

function tests = homologsTests
    tests = functiontests(localfunctions);
end

% P4 (reaction R3, substrate m1) has an exact copy h1, measured on m1 at 10
% and 40 /s, a copy h2 with 60 of 300 residues changed (80% identity),
% measured at 100 /s, and an unrelated h3 at 1000 /s. geckopy's
% test_homologs uses the same fixture (ecTestGEM/data/homologFixture).

function testKcatsFromHomologs_tc0095(testCase)
    [model, adapter, project] = buildProject();
    kcatList = kcatsFromHomologs(model, 'modelAdapter', adapter);
    verifyEqual(testCase, kcatList.rxns, {'R3'})
    % h1 is the same protein (100% identity): its median alone, 20 /s.
    verifyEqual(testCase, kcatList.kcats, 20, 'RelTol', 1e-12)
    verifyEqual(testCase, kcatList.detail, {'identity=1.000; n=1'})
    verifyEqual(testCase, kcatList.substrates, {'m1'})
    [~, loaded] = loadKcatSource(model, fullfile(project, 'data', 'kcats', 'homologs.tsv'));
    verifyEqual(testCase, loaded.header.parameters, '{"coverage":0.8,"k":10,"min_identity":0.4}')
    verifyEqual(testCase, loaded.header.data, 'BRENDA release test')
    verifyEqual(testCase, loaded.records.rxn, {'R3'})
    % Without h1, the homologs below 95% identity are averaged by identity:
    % h2 (80%) alone.
    seqs = fullfile(project, 'data', 'sequences.fasta');
    lines = strsplit(strtrim(fileread(seqs)), newline);
    keep = ~ismember(1:numel(lines), [1 2]); % the >h1 header and its sequence
    fid = fopen(seqs, 'w'); fprintf(fid, '%s\n', lines{keep}); fclose(fid);
    rest = kcatsFromHomologs(model, 'save', false, 'modelAdapter', adapter);
    verifyEqual(testCase, rest.kcats, 100, 'RelTol', 1e-12)
    verifyEqual(testCase, rest.detail, {'identity=0.800; n=1'})
    % Above 90% identity nothing is left.
    high = kcatsFromHomologs(model, 'minIdentity', 0.9, 'save', false, 'modelAdapter', adapter);
    verifyEmpty(testCase, high.rxns)
end

function testSubstrateMatching_tc0096(testCase)
    % ChEBI and MetaNetX annotations first, then the name; geckopy's
    % test_metabolite_inchikeys has the same case.
    [model, adapter, project] = buildProject();
    fid = fopen(fullfile(project, 'data', 'substrates.tsv'), 'w');
    fprintf(fid, 'inchikey\tkind\tvalue\nAAAAAAAAAAAAAA-BBBBBBBBSA\tchebi\t12345\n');
    fclose(fid);
    i = find(strcmp(model.mets, 'm1c'));
    model.metNames{i} = 'not m1';
    model.metMiriams{i} = struct('name', {{'chebi'}}, 'value', {{'CHEBI:12345'}});
    kcatList = kcatsFromHomologs(model, 'save', false, 'modelAdapter', adapter);
    verifyEqual(testCase, kcatList.rxns, {'R3'})
    pairs = extractEnzymeSubstratePairs(model, 'modelAdapter', adapter);
    verifyTrue(testCase, any(model.ec.rxns(pairs.ecRxn) == "R3" & pairs.met == i))
end

function testHomologPresets_tc0097(testCase)
    [~, adapter] = buildProject();
    adapter.params.kcat_assignment = struct('preset', 'ec_tree_homologs');
    config = kcatAssignmentConfig(adapter.params);
    verifyEqual(testCase, {config.stages.action}, {'select', 'combine', 'combine'})
    verifyEqual(testCase, config.stages(2).files, {'ec_tree', 'homologs'})
    adapter.params.kcat_assignment = struct('preset', 'ec_tree_okp');
    config = kcatAssignmentConfig(adapter.params);
    verifyEqual(testCase, config.stages(2).files(1:2), {'ec_tree', 'homologs'})
end

function [model, adapter, project] = buildProject()
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    project = tempname;
    copyfile(fullfile(adapter.params.path, 'data'), fullfile(project, 'data'));
    % TestGEMAdapter's BRENDA folder is the project's data folder.
    fixture = fullfile(geckoPath,'test','unit_tests','ecTestGEM','data','homologFixture');
    for f = {'measurements.tsv', 'sequences.fasta', 'substrates.tsv'}
        copyfile(fullfile(fixture, f{1}), fullfile(project, 'data', f{1}));
    end
    adapter.params.path = project;
    model = makeEcModel(getGeckoTestModel(), false, adapter);
    model = getECfromGEM(model);
    lines = strsplit(strtrim(fileread(fullfile(fixture, 'enzymes.fasta'))), newline);
    for k = 1:2:numel(lines)
        hit = strcmp(model.ec.enzymes, lines{k}(2:end));
        model.ec.sequence(hit) = lines(k+1);
    end
end
