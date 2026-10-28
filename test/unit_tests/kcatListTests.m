%run this test case with the command
%results = runtests('kcatListTests.m')

function tests = kcatListTests
    tests = functiontests(localfunctions);
end

function testAssignKcatValuesSourceAfterZeroRows_tc0104(testCase)
    % Rows with a zero kcat are dropped together with their source.
    adapter = getTestAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    kcatList = struct();
    kcatList.source = 'test';
    kcatList.rxns = ecModel.ec.rxns(1:2);
    kcatList.kcats = [0; 5];
    kcatList.kcatSource = {'a'; 'b'};
    ecModel = assignKcatValues(ecModel, kcatList);
    verifyEqual(testCase, ecModel.ec.kcat(2), 5)
    verifyEqual(testCase, ecModel.ec.source{2}, 'b')
    verifyEqual(testCase, ecModel.ec.kcat(1), 0)
end

function testLoadBRENDAdataLowercaseEc_tc0105(testCase)
    % An EC code with a lowercase letter (provisional EC numbers end in n)
    % is found in the molecular weight table.
    adapter = getTestAdapter();
    scratchDir = tempname;
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    mkdir(fullfile(scratchDir,'data'));
    writeTable(fullfile(scratchDir,'data','kcat.tsv'), ...
        'ec_code\tsubstrate\torganism\tkcat_max\tkcat_median\tn\treferences\n', ...
        '3.5.1.n3\tm1\ttestorg\t99\t5\t3\tPMID:1\n');
    writeTable(fullfile(scratchDir,'data','mw.tsv'), ...
        'ec_code\tsubstrate\torganism\tmw\tn\treferences\n', ...
        '3.5.1.n3\t*\ttestorg\t60000\t1\t*\n');
    writeTable(fullfile(scratchDir,'data','sa.tsv'), ...
        'ec_code\tsubstrate\torganism\tsa_max\tsa_median\tn\treferences\n', ...
        '3.5.1.n3\t*\ttestorg\t120\t7\t2\tPMID:2\n');
    brendaAdapter = adapter;
    brendaAdapter.params.path = scratchDir;
    [~, SAcell] = loadBRENDAdata('modelAdapter', brendaAdapter);
    verifyEqual(testCase, SAcell{1}, {'3.5.1.n3'})
    verifyEqual(testCase, SAcell{3}, (120/60) * (60000/1000), 'AbsTol', 1e-9)
end

function adapter = getTestAdapter()
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
end

function writeTable(file, header, row)
    fid = fopen(file, 'w');
    fprintf(fid, '# BRENDA release test-fixture\n');
    fprintf(fid, header);
    fprintf(fid, row);
    fclose(fid);
end
