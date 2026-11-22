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
    % The empty name is not queried, so nothing is added to the cache.
    verifyEqual(testCase, fileread(fullfile(adapter.params.path, 'data', 'smilesDB.tsv')), sprintf('a\tCCO\n'))
end

function testFindMetSmilesCachesNoMatchButNotFailures_tc0409(testCase)
    adapter = cacheAdapter(testCase, sprintf('a\tCCO\n'));
    fakeWebread(testCase, ...
        {'nomatch','MATLAB:webservices:HTTP404StatusCodeError'; ...
         'badname','MATLAB:webservices:HTTP400StatusCodeError'; ...
         'broken','MATLAB:webservices:HTTP500StatusCodeError'; ...
         'ok',''});
    model = struct('metNames', {{'a';'nomatch';'badname';'broken';'ok';'prot_P1'}});
    warnState = warning('off', 'findMetSmiles:requestFailed');
    [~, out] = evalc('findMetSmiles(model, ''modelAdapter'', adapter, ''verbose'', false)');
    warning(warnState)
    verifyEqual(testCase, out.metSmiles, {'CCO';'';'';'';'C1CC1';''})
    % A clean no match is cached with an empty SMILES; the failed request and
    % the skipped prot_ name are not.
    cache = fileread(fullfile(adapter.params.path, 'data', 'smilesDB.tsv'));
    verifyEqual(testCase, sort(strsplit(strtrim(cache), newline)), ...
        sort({sprintf('a\tCCO'), sprintf('nomatch\t'), sprintf('badname\t'), sprintf('ok\tC1CC1')}))
    verifyWarning(testCase, @() quietFind(model, adapter), 'findMetSmiles:requestFailed')
end

function fakeWebread(testCase, answers)
    % Put a webread and a pause on the path that answer from a table: the
    % last path element of the URL selects the name, a non-empty identifier
    % raises that error, otherwise a SMILES is returned.
    folder = tempname();
    mkdir(folder);
    testCase.addTeardown(@() rmdir(folder, 's'))
    fid = fopen(fullfile(folder, 'webread.m'), 'w');
    fprintf(fid, 'function out = webread(url, varargin)\n');
    fprintf(fid, 'tok = regexp(url, ''name/(.*)/property'', ''tokens'', ''once'');\n');
    fprintf(fid, 'switch tok{1}\n');
    for i = 1:size(answers,1)
        if isempty(answers{i,2})
            fprintf(fid, 'case ''%s''\n out = sprintf(''C1CC1\\nOTHER\\n'');\n', answers{i,1});
        else
            fprintf(fid, 'case ''%s''\n error(''%s'', ''fake'');\n', answers{i,1}, answers{i,2});
        end
    end
    fprintf(fid, 'end\nend\n');
    fclose(fid);
    fid = fopen(fullfile(folder, 'pause.m'), 'w');
    fprintf(fid, 'function pause(varargin)\nend\n');
    fclose(fid);
    shadow = warning('off', 'MATLAB:dispatcher:nameConflict');
    addpath(folder);
    testCase.addTeardown(@() rmpath(folder))
    testCase.addTeardown(@() warning(shadow))
    rehash
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
