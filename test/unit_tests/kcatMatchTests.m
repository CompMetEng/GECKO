%run this test case with the command
%results = runtests('kcatMatchTests.m')

function tests = kcatMatchTests
    tests = functiontests(localfunctions);
end

function testLoadBRENDAdataDropsKcatFarBelowItsEc_tc0220(testCase)
    % A kcat more than four orders of magnitude below the median of its EC
    % is dropped by default and kept when the filter is switched off.
    adapter = brendaAdapter(testCase, ...
        [kcatRows('1.1.1.1', {'o1','o2','o3','o4','o5','o6'}, 10, 10), ...
         kcatRows('1.1.1.1', {'o7'}, 1e-4, 1e-4)]);
    KCATcell = loadBRENDAdata('modelAdapter', adapter);
    verifyEqual(testCase, KCATcell{4}, repmat(10,6,1))
    KCATcell = loadBRENDAdata('modelAdapter', adapter, 'filterOutlierKcats', false);
    verifyEqual(testCase, sort(KCATcell{4}), sort([repmat(10,6,1); 1e-4]))
end

function testLoadBRENDAdataReferenceUsesMaxAndMedian_tc0221(testCase)
    % The reference of an EC is built from the max and the median values,
    % whichever aggregation is read: a median that is far below its EC is an
    % outlier even when the max of the same row is not.
    adapter = brendaAdapter(testCase, ...
        [kcatRows('2.1.1.1', {'o1','o2','o3','o4'}, 100, 100), ...
         kcatRows('2.1.1.1', {'o5'}, 100, 0.001)]);
    KCATcell = loadBRENDAdata('modelAdapter', adapter, 'kcatAggregation', 'median');
    verifyEqual(testCase, KCATcell{4}, repmat(100,4,1))
    KCATcell = loadBRENDAdata('modelAdapter', adapter, 'kcatAggregation', 'max');
    verifyEqual(testCase, KCATcell{4}, repmat(100,5,1))
end

function testLoadBRENDAdataReferenceNeedsFiveValues_tc0222(testCase)
    % The values of the max and the median of a row both count towards the
    % reference of its EC. With fewer than five of them at every EC level
    % there is no reference and no row is dropped; with five or more at a
    % shorter EC prefix, that prefix is the reference.
    adapter = brendaAdapter(testCase, ...
        [kcatRows('3.1.1.1', {'o1'}, 10, 10), kcatRows('3.1.1.1', {'o2'}, 1e-9, 1e-9)]);
    KCATcell = loadBRENDAdata('modelAdapter', adapter);
    verifyEqual(testCase, numel(KCATcell{4}), 2)
    adapter = brendaAdapter(testCase, ...
        [kcatRows('4.1.1.2', {'o1','o2','o3'}, 10, 10), kcatRows('4.1.1.1', {'o4'}, 1e-9, 1e-9)]);
    KCATcell = loadBRENDAdata('modelAdapter', adapter);
    verifyEqual(testCase, KCATcell{4}, repmat(10,3,1))
end

function testFuzzyKcatMatchingFullWildcard_tc0223(testCase)
    % An EC number without any match at the wildcard levels 1 to 3 is matched
    % against every EC number in BRENDA, at wildcard level 4.
    [ecModel, fuzzyAdapter] = fuzzyModel(testCase, {'5.5.5.5\tm1\ttestus testus\t42\t42\t1\t*\n'});
    ecModel.ec.eccodes{strcmp(ecModel.ec.rxns, 'R2_EXP_1')} = '9.9.9.9';
    kcatList = fuzzyKcatMatching(ecModel, ismember(ecModel.ec.rxns,{'R2_EXP_1'}), fuzzyAdapter);
    verifyEqual(testCase, kcatList.kcats, 42)
    verifyEqual(testCase, kcatList.wildcardLvl, 4)
end

function testFuzzyKcatMatchingTieReportsTheSelectedEc_tc0224(testCase)
    % Two EC numbers of a reaction give the same kcat, one through a wildcard
    % and one exactly: the wildcard level and origin are those of the exact one.
    [ecModel, fuzzyAdapter] = fuzzyModel(testCase, ...
        {'1.1.1.1\tm1\ttestus testus\t7\t7\t1\t*\n', ...
         '1.1.1.9\tm1\ttestus testus\t7\t7\t1\t*\n'});
    ecModel.ec.eccodes{strcmp(ecModel.ec.rxns, 'R2_EXP_1')} = '1.1.1.5;1.1.1.1';
    kcatList = fuzzyKcatMatching(ecModel, ismember(ecModel.ec.rxns,{'R2_EXP_1'}), fuzzyAdapter);
    verifyEqual(testCase, kcatList.kcats, 7)
    verifyEqual(testCase, kcatList.wildcardLvl, 0)
    verifyEqual(testCase, kcatList.origin, 1)
end

function testFuzzyKcatMatchingSubsetKeepsOrderAndRepeats_tc0225(testCase)
    % Indices select the reactions in the order given, a repeated index
    % gives a repeated reaction, and a logical vector gives model order.
    [ecModel, fuzzyAdapter] = fuzzyModel(testCase, {'1.1.1.1\tm1\ttestus testus\t7\t7\t1\t*\n'});
    ids = {'R3','R2_EXP_1','R2_EXP_1'};
    idx = cellfun(@(r) find(strcmp(ecModel.ec.rxns, r)), ids);
    kcatList = fuzzyKcatMatching(ecModel, 'ecRxns', idx, 'modelAdapter', fuzzyAdapter);
    verifyEqual(testCase, kcatList.rxns(:)', ids)
    verifyEqual(testCase, numel(kcatList.kcats), 3)
    kcatList = fuzzyKcatMatching(ecModel, 'ecRxns', ismember(ecModel.ec.rxns, ids), ...
        'modelAdapter', fuzzyAdapter);
    verifyEqual(testCase, kcatList.rxns, ecModel.ec.rxns(ismember(ecModel.ec.rxns, ids)))
end

function testFuzzyKcatMatchingEmptySubsetSelectsNone_tc0226(testCase)
    % [] selects all reactions; any other empty selection selects none.
    [ecModel, fuzzyAdapter] = fuzzyModel(testCase, {'1.1.1.1\tm1\ttestus testus\t7\t7\t1\t*\n'});
    kcatList = fuzzyKcatMatching(ecModel, 'ecRxns', [], 'modelAdapter', fuzzyAdapter);
    verifyEqual(testCase, numel(kcatList.rxns), numel(ecModel.ec.rxns))
    kcatList = fuzzyKcatMatching(ecModel, 'ecRxns', zeros(0,1), 'modelAdapter', fuzzyAdapter);
    verifyEmpty(testCase, kcatList.rxns)
    verifyEmpty(testCase, kcatList.kcats)
    kcatList = fuzzyKcatMatching(ecModel, 'ecRxns', false(numel(ecModel.ec.rxns),1), ...
        'modelAdapter', fuzzyAdapter);
    verifyEmpty(testCase, kcatList.rxns)
    adapter = getTestAdapter();
    verifyError(testCase, @() kcatsFromBrenda(ecModel, 'ecRxns', zeros(0,1), 'save', false, ...
        'modelAdapter', adapter), 'kcatsFromBrenda:emptySubset')
end

function testFuzzyKcatMatchingForceWildcardEscalatesFromTheRight_tc0231(testCase)
    % forceWClvl replaces the rightmost level that is not a wildcard yet, for
    % every EC number of a reaction, also a provisional one: level 4 gives
    % -.-.-.- and level 1 gives a wildcard at the last level.
    [ecModel, fuzzyAdapter] = fuzzyModel(testCase, {'1.1.1.1\tm1\ttestus testus\t7\t7\t1\t*\n'});
    ecModel.ec.eccodes{strcmp(ecModel.ec.rxns, 'R2_EXP_1')} = '1.1.1.1;3.5.1.n3';
    sel = ismember(ecModel.ec.rxns,{'R2_EXP_1'});
    kcatList = fuzzyKcatMatching(ecModel, 'ecRxns', sel, 'modelAdapter', fuzzyAdapter, 'forceWClvl', 1);
    verifyEqual(testCase, kcatList.eccodes, {'1.1.1.-;3.5.1.-'})
    kcatList = fuzzyKcatMatching(ecModel, 'ecRxns', sel, 'modelAdapter', fuzzyAdapter, 'forceWClvl', 4);
    verifyEqual(testCase, kcatList.eccodes, {'-.-.-.-;-.-.-.-'})
end

function testApplyCustomKcatsSkipsABlankKcat_tc0227(testCase)
    % An entry without a kcat leaves the reaction as it is.
    [ecModel, adapter] = customModel();
    isR3 = strcmp(ecModel.ec.rxns,'R3');
    isR2 = strcmp(ecModel.ec.rxns,'R2_EXP_1');
    custom = customStruct({''; ''}, [NaN; 50], {'R3'; 'R2'}, {'';''});
    ecModel = applyCustomKcats(ecModel, custom, adapter);
    verifyEqual(testCase, ecModel.ec.kcat(isR3), 10)
    verifyFalse(testCase, strcmp(ecModel.ec.source{isR3}, 'custom'))
    verifyEqual(testCase, ecModel.ec.kcat(isR2), 50)
    verifyEqual(testCase, ecModel.ec.source{isR2}, 'custom')
end

function testApplyCustomKcatsEmptyNoteIsNotAppended_tc0228(testCase)
    % An empty note does not add a separator to the note of the reaction.
    [ecModel, adapter] = customModel();
    custom = customStruct({''; ''}, [50; 60], {'R3'; 'R3'}, {'first'; ''});
    ecModel = applyCustomKcats(ecModel, custom, adapter);
    verifyEqual(testCase, ecModel.ec.notes{strcmp(ecModel.ec.rxns,'R3')}, 'first')
    custom = customStruct({''; ''}, [50; 60], {'R3'; 'R3'}, {''; 'second'});
    ecModel = applyCustomKcats(customModelAgain(), custom, adapter);
    verifyEqual(testCase, ecModel.ec.notes{strcmp(ecModel.ec.rxns,'R3')}, 'second')
end

function testApplyCustomKcatsMatchesIdsExactly_tc0229(testCase)
    % Protein and reaction ids are case sensitive.
    [ecModel, adapter] = customModel();
    before = ecModel.ec.kcat;
    custom = customStruct({'p3'}, 123, {''}, {''});
    updated = applyCustomKcats(ecModel, custom, adapter);
    verifyEqual(testCase, updated.ec.kcat, before)
    custom = customStruct({'P3'}, 123, {'r3'}, {''});
    updated = applyCustomKcats(ecModel, custom, adapter);
    verifyEqual(testCase, updated.ec.kcat, before)
    custom = customStruct({'P3'}, 123, {''}, {''});
    updated = applyCustomKcats(ecModel, custom, adapter);
    verifyTrue(testCase, any(updated.ec.kcat == 123))
end

function testApplyCustomKcatsLightModel_tc0230(testCase)
    % A light ecModel is updated by the proteins of its rows, and all three
    % outputs are returned.
    adapter = getTestAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), true, adapter);
    ecModel = getECfromGEM(ecModel);
    custom = customStruct({'P3'}, 123, {''}, {'light'});
    [updated, rxnUpdated, notMatch] = applyCustomKcats(ecModel, custom, adapter);
    verifyNotEmpty(testCase, rxnUpdated)
    verifyTrue(testCase, all(updated.ec.kcat(ismember(updated.ec.rxns, rxnUpdated)) == 123))
    verifyTrue(testCase, istable(notMatch))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adapter = getTestAdapter()
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
end

function rows = kcatRows(ec, organisms, kcatMax, kcatMedian)
    rows = cell(1, numel(organisms));
    for i = 1:numel(organisms)
        rows{i} = sprintf('%s\tm1\t%s\t%g\t%g\t1\t*\n', ec, organisms{i}, kcatMax, kcatMedian);
    end
end

function writeBrenda(brendaDir, kcatLines)
    fid = fopen(fullfile(brendaDir, 'kcat.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - kcat in 1/s\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tkcat_max\tkcat_median\tn\treferences\n');
    for i = 1:numel(kcatLines)
        fprintf(fid, kcatLines{i});
    end
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'mw.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - molecular weight in g/mol\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tmw\tn\treferences\n');
    fprintf(fid, '0.0.0.0\t*\tplaceholder\t1\t1\t*\n');
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'sa.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - specific activity in umol/min/mg\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tsa_max\tsa_median\tn\treferences\n');
    fprintf(fid, '0.0.0.0\t*\tplaceholder\t1\t1\t1\t*\n');
    fclose(fid);
end

function adapter = brendaAdapter(testCase, kcatLines)
    scratchDir = tempname;
    brendaDir = fullfile(scratchDir, 'data');
    mkdir(brendaDir);
    testCase.addTeardown(@() rmdir(scratchDir, 's'));
    writeBrenda(brendaDir, kcatLines);
    adapter = getTestAdapter();
    adapter.params.path = scratchDir;
end

function [ecModel, fuzzyAdapter] = fuzzyModel(testCase, kcatLines)
    fuzzyAdapter = brendaAdapter(testCase, kcatLines);
    geckoPath = findGECKOroot;
    copyfile(fullfile(geckoPath,'test','unit_tests','ecTestGEM','data','PhylDist.mat'), ...
        fullfile(fuzzyAdapter.params.path, 'data', 'PhylDist.mat'));
    adapter = getTestAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    ecModel = getECfromGEM(ecModel);
end

function [ecModel, adapter] = customModel()
    adapter = getTestAdapter();
    ecModel = customModelAgain();
end

function ecModel = customModelAgain()
    adapter = getTestAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(:) = 10;
end

function custom = customStruct(proteins, kcats, rxns, notes)
    custom.proteins = proteins;
    custom.kcat     = kcats;
    custom.rxns     = rxns;
    custom.notes    = notes;
end
