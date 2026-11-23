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

function testGetECfromDatabaseEmptySelectionChangesNothing_tc0406(testCase)
    % An empty logical selection searches nothing and returns the model as it
    % is; [] selects every ec reaction.
    [ecModel, adapter] = conflictModel(testCase);
    ecModel.ec.eccodes = repmat({'9.9.9.9'}, numel(ecModel.ec.rxns), 1);
    lastwarn('')
    same = getECfromDatabase(ecModel, false(0,1), 'display', adapter);
    verifyEqual(testCase, same, ecModel)
    verifyEmpty(testCase, lastwarn())
    warning('off', 'getECfromDatabase:conflicts')
    cleanup = onCleanup(@() warning('on', 'getECfromDatabase:conflicts')); %#ok<NASGU>
    all = getECfromDatabase(ecModel, [], 'ignore', adapter);
    verifyEqual(testCase, all.ec.eccodes{1}, '1.1.1.1')
end

function testGetECfromDatabaseSelectionMatchesFullRun_tc0733(testCase)
    % Selecting entries gives the same eccodes for them as a run over all
    % entries; an all-false selection of full length changes nothing.
    [ecModel, adapter] = conflictModel(testCase);
    warning('off', 'getECfromDatabase:conflicts')
    cleanup = onCleanup(@() warning('on', 'getECfromDatabase:conflicts')); %#ok<NASGU>
    full = getECfromDatabase(ecModel, [], 'ignore', adapter);
    n = numel(ecModel.ec.rxns);
    selected = false(n,1);
    selected(1:2:n) = true;
    part = getECfromDatabase(ecModel, selected, 'ignore', adapter);
    verifyEqual(testCase, part.ec.eccodes(selected), full.ec.eccodes(selected))
    verifyEqual(testCase, part.ec.eccodes(~selected), repmat({''}, nnz(~selected), 1))
    none = getECfromDatabase(ecModel, false(n,1), 'ignore', adapter);
    verifyEqual(testCase, none.ec.eccodes, repmat({''}, n, 1))
end

function testGetECfromDatabaseRejectsUnsupportedAction_tc0409(testCase)
    [ecModel, adapter] = conflictModel(testCase);
    verifyError(testCase, @() getECfromDatabase(ecModel, [], 'add', adapter), 'getECfromDatabase:unsupportedAction')
    verifyError(testCase, @() getECfromDatabase(ecModel, [], 'other', adapter), 'getECfromDatabase:unsupportedAction')
end

function testLoadDatabasesSplitsGeneCells_tc0748(testCase)
    % A gene cell with several names is one row per name, whatever separates
    % them; an empty cell stays one row, and a repeated ID across file rows is
    % still a duplicate.
    [~, adapter] = conflictModel(testCase);
    file = fullfile(adapter.params.path, 'data', 'uniprot.tsv');
    before = loadDatabases('uniprot', adapter);
    fid = fopen(file, 'a');
    fprintf(fid, 'P7\tG7  G8\tx\t70000\tMKAL\nP8\t\ty\t80000\tMNTD\n');
    fprintf(fid, 'P9\tG9; G10;G11\tz\t90000\tMSYN\nP10\t ; \tw\t10000\tMAAA\n');
    fclose(fid);
    db = loadDatabases('uniprot', adapter);
    i = strcmp(db.uniprot.ID, 'P7');
    verifyEqual(testCase, db.uniprot.genes(i), {'G7';'G8'})
    verifyEqual(testCase, db.uniprot.eccodes(i), {'x';'x'})
    verifyEqual(testCase, db.uniprot.MW(i), [70000;70000])
    verifyEqual(testCase, db.uniprot.seq(i), {'MKAL';'MKAL'})
    verifyEqual(testCase, db.uniprot.genes(strcmp(db.uniprot.ID, 'P9')), {'G9';'G10';'G11'})
    verifyEqual(testCase, db.uniprot.genes(strcmp(db.uniprot.ID, 'P8')), {''})
    verifyEqual(testCase, db.uniprot.genes(strcmp(db.uniprot.ID, 'P10')), {''})
    verifyEqual(testCase, numel(db.uniprot.ID), numel(before.uniprot.ID) + 2 + 3 + 1 + 1)
    fid = fopen(file, 'a');
    fprintf(fid, 'P7\tG12\tz\t90000\tMSYN\n');
    fclose(fid);
    verifyError(testCase, @() loadDatabases('uniprot', adapter), ?MException)
end

function testGetECfromDatabaseMatchesSecondGeneName_tc0749(testCase)
    % A model gene that is the second name in a gene cell finds its UniProt
    % entry; before, only a gene equal to the whole cell matched.
    [ecModel, adapter] = conflictModel(testCase);
    file = fullfile(adapter.params.path, 'data', 'uniprot.tsv');
    txt = fileread(file);
    txt = strrep(txt, sprintf('P5\tG5\t\t'), sprintf('P5\tXX; G5\t3.3.3.3\t'));
    verifyNotEqual(testCase, txt, fileread(file))
    fid = fopen(file, 'w');
    fwrite(fid, txt);
    fclose(fid);
    warning('off', 'getECfromDatabase:conflicts')
    cleanup = onCleanup(@() warning('on', 'getECfromDatabase:conflicts')); %#ok<NASGU>
    r = find(strcmp(ecModel.ec.rxns, 'R5'));
    verifyNumElements(testCase, r, 1)
    model = getECfromDatabase(ecModel, [], 'ignore', adapter);
    verifyEqual(testCase, model.ec.eccodes{r}, '3.3.3.3')
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
