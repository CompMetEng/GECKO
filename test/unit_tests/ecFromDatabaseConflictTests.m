%run this test case with the command
%results = runtests('ecFromDatabaseConflictTests.m')

function tests = ecFromDatabaseConflictTests
    tests = functiontests(localfunctions);
end

function testGetECfromDatabaseDisplayWarnsAndUsesFirst_tc0404(testCase)
    [ecModel, adapter] = conflictModel(testCase);
    verifyWarning(testCase, @() getECfromDatabase(ecModel, [], 'display', adapter), 'getECfromDatabase:conflicts')
    ecModel2 = getECfromDatabase(ecModel, [], 'ignore', adapter);
    verifyEqual(testCase, ecModel2.ec.eccodes, {'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3'})
    warning('off', 'getECfromDatabase:conflicts')
    cleanup = onCleanup(@() warning('on', 'getECfromDatabase:conflicts')); %#ok<NASGU>
    displayed = getECfromDatabase(ecModel, [], 'display', adapter);
    verifyEqual(testCase, displayed.ec.eccodes, ecModel2.ec.eccodes)
end

function testGetECfromDatabaseIgnoresConflictOutsideSelection_tc0405(testCase)
    [ecModel, adapter] = conflictModel(testCase);
    selected = strcmp(ecModel.ec.rxns, 'R2_EXP_1');
    lastwarn('')
    ecModel2 = getECfromDatabase(ecModel, selected, 'display', adapter);
    verifyEmpty(testCase, lastwarn())
    verifyEqual(testCase, ecModel2.ec.eccodes(selected), {'1.1.1.1'})
    verifyEqual(testCase, ecModel2.ec.eccodes(~selected), repmat({''}, nnz(~selected), 1))
end

function [ecModel, adapter] = conflictModel(testCase)
    % ecTestGEM with a second UniProt entry for gene G4 and another EC number.
    geckoPath = findGECKOroot;
    source = fullfile(geckoPath,'test','unit_tests','ecTestGEM');
    adapter = ModelAdapterManager.getAdapter(fullfile(source, 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    folder = tempname();
    mkdir(fullfile(folder, 'data'));
    testCase.addTeardown(@() rmdir(folder, 's'))
    copyfile(fullfile(source, 'data', 'uniprot.tsv'), fullfile(folder, 'data', 'uniprot.tsv'))
    copyfile(fullfile(source, 'data', 'kegg.tsv'), fullfile(folder, 'data', 'kegg.tsv'))
    fid = fopen(fullfile(folder, 'data', 'uniprot.tsv'), 'a');
    fprintf(fid, '\nP4b\tG4\t1.1.2.2\t40001\tMDFM\n');
    fclose(fid);
    adapter.params.path = folder;
end
