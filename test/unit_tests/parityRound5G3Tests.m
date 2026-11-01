%run this test case with the command
%results = runtests('parityRound5G3Tests.m')

function tests = parityRound5G3Tests
    tests = functiontests(localfunctions);
end

function testAddNewRxnsToECWithoutNewEnzymes_tc0500(testCase)
    % The genes of the new reaction already have enzymes: the reaction is
    % added and no enzyme is.
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    nEnzymes = numel(ecModel.ec.enzymes);
    nEcRxns = numel(ecModel.ec.rxns);
    newRxns.rxns      = {'RNEW'};
    newRxns.rxnNames  = {'RNEW'};
    newRxns.equations = {'m1[c] => e2[e]'};
    newRxns.grRules   = {'G1'};
    [ecModel, rxnsAdded, enzAdded] = addNewRxnsToEC(ecModel, newRxns, [], adapter);
    verifyEqual(testCase, rxnsAdded, {'RNEW'})
    verifyEmpty(testCase, enzAdded)
    verifyEqual(testCase, numel(ecModel.ec.enzymes), nEnzymes)
    verifyEqual(testCase, numel(ecModel.ec.rxns), nEcRxns + 1)
    [~, rxnIdx] = ismember('RNEW', ecModel.ec.rxns);
    [~, enzIdx] = ismember('P1', ecModel.ec.enzymes);
    verifyEqual(testCase, full(ecModel.ec.rxnEnzMat(rxnIdx, :)), double((1:nEnzymes) == enzIdx))
end

function testApplyComplexDataPartialZeroStoichiometry_tc0501(testCase)
    % A subunit with stoichiometry 0 counts as a subunit with stoichiometry 1.
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    complexInfo = jsondecode(['[{"complexID":"R2Compl","name":"R2Compl","specie":"x",' ...
        '"geneName":["G1","G2"],"protID":["P1","P2"],"stochiometry":[2,0]}]']);
    [ecModel, foundComplex] = applyComplexData(ecModel, 'complexInfo', complexInfo, ...
        'modelAdapter', adapter, 'verbose', false);
    [~, p12] = ismember({'P1','P2'}, ecModel.ec.enzymes);
    verifyEqual(testCase, full(ecModel.ec.rxnEnzMat(1, p12)), [2 1])
    verifyEqual(testCase, height(foundComplex), 2)
end

function testCalculateFfactorEmptyEnzymeList_tc0502(testCase)
    % An empty cell array means no enzymes, [] means every enzyme of the model.
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    protData.uniprotIDs = {'P1';'P2';'X1'};
    protData.abundances = [1;1;2];
    verifyEqual(testCase, calculateFfactor(ecModel, 'protData', protData, ...
        'enzymes', {}, 'modelAdapter', adapter), 0)
    verifyEqual(testCase, calculateFfactor(ecModel, 'protData', protData, ...
        'modelAdapter', adapter), 0.5)
    verifyEqual(testCase, calculateFfactor(ecModel, 'protData', protData, ...
        'enzymes', {'P1'}, 'modelAdapter', adapter), 0.25)
end

function testConstrainEnzConcsNegativeConcentrationRaises_tc0503(testCase)
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    ecModel.ec.concs(:) = NaN;
    ecModel.ec.concs(1) = -0.5;
    verifyError(testCase, @() constrainEnzConcs(ecModel), 'constrainEnzConcs:negativeConcentration')
    % removing the constraints does not read the concentrations
    cleared = constrainEnzConcs(ecModel, 'removeConstraints', true);
    verifyEqual(testCase, cleared.ub(startsWith(cleared.rxns, 'usage_prot_')), ...
        1000 * ones(numel(ecModel.ec.enzymes), 1))
end

function testFillEnzConcsNegativeAbundanceRaises_tc0504(testCase)
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    protData.uniprotIDs = {'P1';'X1'};
    protData.abundances = [-0.5;-1];
    verifyError(testCase, @() fillEnzConcs(ecModel, protData), 'fillEnzConcs:negativeAbundance')
    % a negative abundance of a protein that is not an enzyme of the model is not read
    protData.uniprotIDs = {'P1';'X1'};
    protData.abundances = [0.5;-1];
    ecModel = fillEnzConcs(ecModel, protData);
    verifyEqual(testCase, ecModel.ec.concs(1), 0.5)
end

function testConstrainFluxDataInvalidArguments_tc0505(testCase)
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    fluxData.conds = {'c1'};
    fluxData.grRate = 0.1;
    fluxData.exchFluxes = 1;
    fluxData.exchRxnIDs = {'S1'};
    fluxData.exchMets = {'m'};
    constrained = constrainFluxData(ecModel, 'fluxData', fluxData, 'modelAdapter', adapter);
    verifyEqual(testCase, constrained.ub(strcmp(constrained.rxns, 'S1')), 1)
    verifyError(testCase, @() constrainFluxData(ecModel, 'fluxData', fluxData, ...
        'maxMinGrowth', 'median', 'modelAdapter', adapter), 'constrainFluxData:invalidMaxMinGrowth')
    fluxData.exchRxnIDs = {'NOT_A_RXN'};
    verifyError(testCase, @() constrainFluxData(ecModel, 'fluxData', fluxData, ...
        'modelAdapter', adapter), 'constrainFluxData:missingReaction')
    fluxData.exchRxnIDs = {'S1'};
    verifyError(testCase, @() constrainFluxData(ecModel, 'fluxData', fluxData, ...
        'condition', 'nope', 'modelAdapter', adapter), 'constrainFluxData:unknownCondition')
end

function testLoadFluxDataMalformedFiles_tc0506(testCase)
    adapter = getAdapter();
    folder = tempname;
    mkdir(folder)
    cleanup = onCleanup(@() rmdir(folder, 's'));
    tab = char(9);

    verifyError(testCase, @() loadFluxData(fullfile(folder, 'missing.tsv'), adapter), ...
        'loadFluxData:fileNotFound')

    fixedOnly = writeLines(fullfile(folder, 'fixed.tsv'), ...
        {['cond' tab 'Ptot' tab 'grRate'], ['c1' tab '0.5' tab '0.4']});
    verifyError(testCase, @() loadFluxData(fixedOnly, adapter), 'loadFluxData:noExchangeColumns')

    ragged = writeLines(fullfile(folder, 'ragged.tsv'), ...
        {['cond' tab 'Ptot' tab 'grRate' tab 'Glc (EX_a)' tab 'O2 (EX_b)'], ...
        ['c1' tab '0.5' tab '0.4' tab '-1' tab '2'], ...
        ['c2' tab '0.5' tab '0.4' tab '-1']});
    verifyError(testCase, @() loadFluxData(ragged, adapter), 'loadFluxData:raggedRow')

    nested = writeLines(fullfile(folder, 'nested.tsv'), ...
        {['cond' tab 'Ptot' tab 'grRate' tab 'X (EX_a(e))' tab 'O2 (EX_b)'], ...
        ['c1' tab '0.5' tab '0.4' tab '-1' tab '2']});
    fluxData = loadFluxData(nested, adapter);
    verifyEqual(testCase, fluxData.exchMets, {'X', 'O2'})
    verifyEqual(testCase, fluxData.exchRxnIDs, {'EX_a(e)', 'EX_b'})
end

function testGetEnzymeUsageLeavesOutMissingUsageReaction_tc0507(testCase)
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    ecModel = removeReactions(ecModel, 'usage_prot_P3');
    fluxes = zeros(numel(ecModel.rxns), 1);
    verifyWarning(testCase, @() getEnzymeUsage(ecModel, fluxes), 'getEnzymeUsage:missingUsageReaction')
    warning('off', 'getEnzymeUsage:missingUsageReaction')
    restore = onCleanup(@() warning('on', 'getEnzymeUsage:missingUsageReaction'));
    usageData = getEnzymeUsage(ecModel, fluxes);
    verifyEqual(testCase, usageData.protID, setdiff(ecModel.ec.enzymes, {'P3'}, 'stable'))
    verifyEqual(testCase, numel(usageData.absUsage), numel(usageData.protID))
end

function testReportEnzymeUsageWithoutZeroUsage_tc0508(testCase)
    % getEnzymeUsage with zero=false filters the enzymes; the report reads
    % the usage of the filtered enzymes.
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(:) = 100;
    ecModel.ec.source(:) = {'manual'};
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
    ecModel.ub(strcmp(ecModel.rxns, 'usage_prot_P5')) = 10;
    fluxes = zeros(numel(ecModel.rxns), 1);
    fluxes(strcmp(ecModel.rxns, 'usage_prot_P5')) = 10;
    fluxes(strcmp(ecModel.rxns, 'R5')) = 10;
    usageData = getEnzymeUsage(ecModel, fluxes, 'zero', false);
    verifyEqual(testCase, usageData.protID, {'P5'})
    report = reportEnzymeUsage(ecModel, usageData);
    verifyEqual(testCase, report.highCapUsage.protID, {'P5'})
    verifyEqual(testCase, report.highCapUsage.absUsage, 10)
    verifyEqual(testCase, report.highCapUsage.capUsage, 1)
    verifyEqual(testCase, report.topAbsUsage.protID, {'P5'})
end

function testReportEnzymeUsageNegativeTopAbsUsageRaises_tc0509(testCase)
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(:) = 100;
    ecModel.ec.source(:) = {'manual'};
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
    usageData = getEnzymeUsage(ecModel, zeros(numel(ecModel.rxns), 1));
    verifyError(testCase, @() reportEnzymeUsage(ecModel, usageData, 'topAbsUsage', -2), ...
        'reportEnzymeUsage:invalidTopAbsUsage')
end

function testTruncateValuesInfinite_tc0510(testCase)
    values = {Inf, -Inf, 1234565, 0.5};
    truncated = truncateValues(values, 1:4);
    verifyEqual(testCase, truncated{1}, Inf)
    verifyEqual(testCase, truncated{2}, -Inf)
    verifyEqual(testCase, truncated{3}, 1234570)
    verifyEqual(testCase, truncated{4}, 0.5)
end

function testRelaxProteomicsGreedyEmptyEnzymeSet_tc0511(testCase)
    % An empty cell array makes no enzyme eligible; [] makes every
    % constrained enzyme eligible.
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(:) = 10;
    ecModel.ec.source(:) = {'manual'};
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel, 0.5, 0.5, 0.5, adapter);
    ecModel.ec.concs(:) = NaN;
    ecModel.ec.concs(strcmp(ecModel.ec.enzymes, 'P5')) = 10;
    ecModel = constrainEnzConcs(ecModel);

    warning('off', 'all')
    restore = onCleanup(@() warning('on', 'all'));
    none = relaxProteomicsGreedy(ecModel, 'minimalGrowth', 50, 'enzymeSet', {});
    verifyFalse(testCase, none.converged)
    verifyEmpty(testCase, fieldnames(none.relaxed))
    all = relaxProteomicsGreedy(ecModel, 'minimalGrowth', 50, 'enzymeSet', []);
    verifyTrue(testCase, all.converged)
    verifyEqual(testCase, fieldnames(all.relaxed), {'P5'})
end

function testApplyKcatConstraintsLightUpdateIds_tc0512(testCase)
    % On a light model the update list takes cobra reaction ids and ec
    % reaction ids (with the ###_ prefix); an id of neither kind raises.
    adapter = getAdapter();
    lecModel = makeEcModel(getGeckoTestModel(), true, adapter);
    lecModel.ec.kcat(:) = 10;
    poolRow = strcmp(lecModel.mets, 'prot_pool');
    all = applyKcatConstraints(lecModel);
    stale = all;
    stale.S(poolRow, strcmp(stale.rxns, 'R3')) = 0;
    byCobra = applyKcatConstraints(stale, 'updateRxns', {'R3'});
    verifyEqual(testCase, byCobra.S, all.S)
    ecId = lecModel.ec.rxns(endsWith(lecModel.ec.rxns, '_R3'));
    byEc = applyKcatConstraints(stale, 'updateRxns', ecId);
    verifyEqual(testCase, byEc.S, all.S)
    verifyError(testCase, @() applyKcatConstraints(lecModel, 'updateRxns', {'NOT_A_RXN'}), ...
        'applyKcatConstraints:unknownReaction')
    lecModel.rxns{end+1} = lecModel.ec.rxns{1};
    verifyError(testCase, @() applyKcatConstraints(lecModel, 'updateRxns', lecModel.ec.rxns(1)), ...
        'applyKcatConstraints:ambiguousReaction')
end

function adapter = getAdapter()
geckoPath = findGECKOroot;
adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
end

function file = writeLines(file, lines)
fid = fopen(file, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end
