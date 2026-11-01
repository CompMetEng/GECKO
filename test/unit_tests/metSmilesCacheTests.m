%run this test case with the command
%results = runtests('metSmilesCacheTests.m')

function tests = metSmilesCacheTests
    tests = functiontests(localfunctions);
end

function testFindMetSmilesUsesLastRowOfDuplicateName_tc0407(testCase)
    adapter = cacheAdapter(testCase, sprintf('a\tFIRST\nb\tCCO\na\tLAST\n'));
    model = struct('metNames', {{'a';'b';'a'}});
    [~, out] = evalc('findMetSmiles(model, ''modelAdapter'', adapter, ''verbose'', false)');
    verifyEqual(testCase, out.metSmiles, {'LAST';'CCO';'LAST'})
end

function testFindMetSmilesSkipsEmptyNameWithWarning_tc0408(testCase)
    adapter = cacheAdapter(testCase, sprintf('a\tCCO\n'));
    model = struct('metNames', {{'a';'';'prot_P1'}});
    verifyWarning(testCase, @() quietFind(model, adapter), 'findMetSmiles:emptyName')
    warning('off', 'findMetSmiles:emptyName')
    cleanup = onCleanup(@() warning('on', 'findMetSmiles:emptyName')); %#ok<NASGU>
    [~, out] = evalc('findMetSmiles(model, ''modelAdapter'', adapter, ''verbose'', false)');
    verifyEqual(testCase, out.metSmiles, {'CCO';'';''})
end

function quietFind(model, adapter)
    evalc('findMetSmiles(model, ''modelAdapter'', adapter, ''verbose'', false)');
end

function adapter = cacheAdapter(testCase, cacheText)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    folder = tempname();
    mkdir(fullfile(folder, 'data'));
    testCase.addTeardown(@() rmdir(folder, 's'))
    fid = fopen(fullfile(folder, 'data', 'smilesDB.tsv'), 'w');
    fprintf(fid, '%s', cacheText);
    fclose(fid);
    adapter.params.path = folder;
end
