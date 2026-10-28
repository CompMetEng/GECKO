%run this test case with the command
%results = runtests('protDataTests.m')

function tests = protDataTests
    tests = functiontests(localfunctions);
end

function testLoadProtDataRsdFilterDropsWholeRow_tc0101(testCase)
    % A protein that exceeds maxRSD is discarded in every replicate, not
    % only in the first one.
    adapter = getTestAdapter();
    file = writeProteomics({'X1','X2','X3'}, [1 1 10; 5 5 5; 6 6 6]);
    protData = loadProtData(3, 'protDataFile', file, 'modelAdapter', adapter, ...
        'maxRSD', 1, 'addStdevs', 0, 'cutLowest', 0);
    verifyEqual(testCase, sort(protData.uniprotIDs), {'X2';'X3'})
    verifyEqual(testCase, sort(protData.abundances), [5;6])
end

function testLoadProtDataCutLowestIsUsed_tc0102(testCase)
    % The lowest cutLowest percent of the values of a condition is discarded.
    adapter = getTestAdapter();
    ids = strcat('P', string(1:40));
    values = repmat((1:40)', 1, 3);
    file = writeProteomics(cellstr(ids), values);
    p0  = loadProtData(3, 'protDataFile', file, 'modelAdapter', adapter, 'addStdevs', 0, 'cutLowest', 0);
    p5  = loadProtData(3, 'protDataFile', file, 'modelAdapter', adapter, 'addStdevs', 0);
    p20 = loadProtData(3, 'protDataFile', file, 'modelAdapter', adapter, 'addStdevs', 0, 'cutLowest', 20);
    verifyEqual(testCase, numel(p0.uniprotIDs), 40)
    verifyEqual(testCase, numel(p5.uniprotIDs), 38)
    verifyEqual(testCase, numel(p20.uniprotIDs), 32)
    verifyEqual(testCase, min(p20.abundances), 9)
end

function testApplyComplexDataStoichiometryByProteinId_tc0103(testCase)
    % The complex lists its proteins in another order than model.ec.enzymes;
    % each protein receives its own stoichiometry.
    adapter = getTestAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    complexInfo = struct('complexID', 'R2Compl', 'name', 'R2Compl', ...
        'specie', 'Testus testus; 000000', 'geneName', {{'G2';'G1'}}, ...
        'protID', {{'P2';'P1'}}, 'stochiometry', [2;1]);
    ecModel = applyComplexData(ecModel, complexInfo, adapter, false);
    rxn = find(strcmp(ecModel.ec.rxns, 'R2_EXP_1'));
    p1 = strcmp(ecModel.ec.enzymes, 'P1');
    p2 = strcmp(ecModel.ec.enzymes, 'P2');
    verifyEqual(testCase, full(ecModel.ec.rxnEnzMat(rxn, p1)), 1)
    verifyEqual(testCase, full(ecModel.ec.rxnEnzMat(rxn, p2)), 2)
end

function adapter = getTestAdapter()
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
end

function file = writeProteomics(ids, values)
    file = [tempname '.tsv'];
    fid = fopen(file, 'w');
    fprintf(fid, 'uniprot\tr1\tr2\tr3\n');
    for i = 1:numel(ids)
        fprintf(fid, '%s\t%g\t%g\t%g\n', ids{i}, values(i,:));
    end
    fclose(fid);
end
