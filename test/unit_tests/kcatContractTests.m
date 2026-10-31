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

function testLoadKcatSourceQuotedRows_tc0302(testCase)
    % A row of empty quoted fields is blank and a quoted comment is a
    % comment; a field that opens a quote and does not close it is an error.
    [model, ~, project] = kcatProject();
    file = fullfile(project, 'odd.tsv');
    writeLines(file, {'# source: odd', sprintf('rxn\tenzymes\tkcat\tsource'), ...
        sprintf('R3\tP4\t25\todd'), '"# a quoted comment"', sprintf('""\t""\t""\t""'), ...
        sprintf('R5\tP5\t30\todd')});
    [list, report] = loadKcatSource(model, file);
    verifyEqual(testCase, numel(report.records.kcat), 2)
    verifyEmpty(testCase, report.invalid)
    writeLines(file, {'# source: odd', sprintf('rxn\tenzymes\tkcat\tsource'), ...
        sprintf('R3\tP4\t25\todd'), '"# a comment that was quoted'});
    verifyError(testCase, @() loadKcatSource(model, file), 'loadKcatSource:quote')
    writeLines(file, {'# source: odd', sprintf('rxn\tenzymes\tkcat\tsource'), ...
        sprintf('"R3\tP4\t25\todd')});
    verifyError(testCase, @() loadKcatSource(model, file), 'loadKcatSource:quote')
end

function testReadDLKcatOutputKcatFields_tc0303(testCase)
    % A kcat field that is not a finite number ('inf', '1.2.3', 'NA', 'None')
    % drops the row before the reactions are checked; a file without a
    % number raises.
    [model, adapter, project] = kcatProject();
    file = fullfile(project, 'DLKcat_out.tsv');
    head = sprintf('rxns\tgenes\tsubstrates\tsmiles\tsequence\tkcat');
    row = @(r, s, k) sprintf('%s\tG1\t%s\tsmiles\tseq\t%s', r, s, k);
    read = @() readDLKcatOutput(model, 'outFile', file, 'modelAdapter', adapter);

    writeLines(file, {head, row('R2_EXP_1', 'm1', '12.5'), row('R2_EXP_2', 'm1', 'inf'), ...
        row('R2_EXP_2', 'm1', '1.2.3')});
    kcatList = read();
    verifyEqual(testCase, kcatList.rxns, {'R2_EXP_1'})
    verifyEqual(testCase, kcatList.kcats, 12.5)

    writeLines(file, {head, row('R2_EXP_1', 'm1', 'None'), row('R2_EXP_2', 'm1', 'None')});
    verifyError(testCase, read, 'readDLKcatOutput:noKcats')

    writeLines(file, {head, row('R2_EXP_1', 'm1', '12.5'), row('NOTAREACTION_EXP_1', 'm1', 'NA')});
    kcatList = read();
    verifyEqual(testCase, kcatList.rxns, {'R2_EXP_1'})

    writeLines(file, {head, row('R2_EXP_1', 'm1', '12.5'), row('NOTAREACTION_EXP_1', 'm1', '7')});
    verifyError(testCase, read, 'readDLKcatOutput:unknownReaction')
end

function testReadDLKcatOutputUnknownSubstrate_tc0304(testCase)
    % A row with a substrate that is not in the model is left out with a
    % warning; strict raises.
    [model, adapter, project] = kcatProject();
    file = fullfile(project, 'DLKcat_out.tsv');
    writeLines(file, {sprintf('rxns\tgenes\tsubstrates\tsmiles\tsequence\tkcat'), ...
        sprintf('R2_EXP_1\tG1\tm1\tsmiles\tseq\t12.5'), sprintf('R2_EXP_2\tG1\tNOTAREALMETABOLITE\tsmiles\tseq\t7')});
    verifyWarning(testCase, @() readDLKcatOutput(model, 'outFile', file, 'modelAdapter', adapter), ...
        'readDLKcatOutput:unknownSubstrate')
    state = warning('off', 'readDLKcatOutput:unknownSubstrate');
    restore = onCleanup(@() warning(state));
    kcatList = readDLKcatOutput(model, 'outFile', file, 'modelAdapter', adapter);
    verifyEqual(testCase, kcatList.rxns, {'R2_EXP_1'})
    verifyError(testCase, @() readDLKcatOutput(model, 'outFile', file, 'modelAdapter', adapter, ...
        'strict', true), 'readDLKcatOutput:unknownSubstrate')
end

function testWriteDLKcatInputEmptySelection_tc0305(testCase)
    % A selection without a reaction writes an empty file and returns no
    % rows; the rows of a selection follow model.ec.rxns.
    [model, adapter, project] = kcatProject();
    file = fullfile(project, 'DLKcat_in.tsv');
    none = false(numel(model.ec.rxns), 1);
    evalc(['written = writeDLKcatInput(model, ''ecRxns'', none, ''modelAdapter'', adapter, ' ...
        '''onlyWithSmiles'', false, ''filename'', file, ''overwrite'', true);']);
    verifyEqual(testCase, size(written), [6 0])
    info = dir(file);
    verifyEqual(testCase, info.bytes, 0)
    some = ismember(model.ec.rxns, {'R5', 'R3'});
    evalc(['written = writeDLKcatInput(model, ''ecRxns'', some, ''modelAdapter'', adapter, ' ...
        '''onlyWithSmiles'', false, ''filename'', file, ''overwrite'', true);']);
    order = find(some);
    verifyEqual(testCase, unique(written(1, :), 'stable'), model.ec.rxns(order)')
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
