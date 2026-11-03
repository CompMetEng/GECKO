%run this test case with the command
%results = runtests('isozymeTests.m')

function tests = isozymeTests
    tests = functiontests(localfunctions);
end

function testGetIsozymeCoverageReportsPartialGroups_tc0287(testCase)
    % W, X and Y have a measured and an unmeasured entry; Z has one entry only.
    % P4 also catalyses Z, so zeroing it would limit Z as well.
    model = isozymeModel([3; NaN; 2; NaN; NaN]);
    [coverage, zeroable] = getIsozymeCoverage(model);
    verifyEqual(testCase, coverage.rxn, {'W';'X';'Y'})
    verifyEqual(testCase, coverage.nEntries, [2;2;2])
    verifyEqual(testCase, coverage.nMeasured, [1;1;1])
    verifyEqual(testCase, coverage.nUnmeasured, [1;1;1])
    verifyEqual(testCase, coverage.measuredEnzymes, {{'P1'};{'P1'};{'P3'}})
    verifyEqual(testCase, coverage.unmeasuredEnzymes, {{'P5'};{'P2'};{'P4'}})
    verifyEqual(testCase, coverage.sharedEnzymes, {cell(0,1);cell(0,1);{'P4'}})
    verifyEqual(testCase, zeroable, {'P2';'P5'})
end

function testGetIsozymeCoverageCountsAMeasuredSubunitAsMeasuredEntry_tc0288(testCase)
    % The entry V has the subunits P1 (measured) and P2 (not measured): it is a
    % measured entry, and P2 is not an unmeasured isozyme.
    model = isozymeModel([3; NaN; NaN; NaN; NaN]);
    model.ec.rxns      = {'V';'V_EXP_1'};
    model.ec.rxnEnzMat = sparse([1 1 0 0 0; 0 0 0 1 1]);
    [coverage, zeroable] = getIsozymeCoverage(model);
    verifyEqual(testCase, coverage.rxn, {'V'})
    verifyEqual(testCase, coverage.measuredEnzymes, {{'P1'}})
    verifyEqual(testCase, coverage.unmeasuredEnzymes, {{'P4';'P5'}})
    verifyEqual(testCase, zeroable, {'P4';'P5'})
end

function testGetIsozymeCoverageIgnoresMeasuredZero_tc0289(testCase)
    % A concentration of 0 is neither a measured sibling nor an unmeasured entry
    model = isozymeModel([0; NaN; 2; NaN; NaN]);
    [coverage, zeroable] = getIsozymeCoverage(model);
    verifyEqual(testCase, coverage.rxn, {'Y'})
    verifyEmpty(testCase, zeroable)
end

function testGetIsozymeCoverageEmptyWithoutMeasurements_tc0290(testCase)
    model = isozymeModel(nan(5, 1));
    [coverage, zeroable] = getIsozymeCoverage(model);
    verifyEqual(testCase, height(coverage), 0)
    verifyEmpty(testCase, zeroable)
    model.ec.concs = [];
    verifyEqual(testCase, height(getIsozymeCoverage(model)), 0)
end

function testGetIsozymeCoverageRejectsLightModel_tc0291(testCase)
    model = isozymeModel([3; NaN; 2; NaN; NaN]);
    model.ec.geckoLight = true;
    verifyError(testCase, @() getIsozymeCoverage(model), 'getIsozymeCoverage:geckoLight')
end

function testZeroUnmeasuredIsozymesSetsUsageBoundToZero_tc0292(testCase)
    model = isozymeModel([3; NaN; 2; NaN; NaN]);
    concs = model.ec.concs;
    [model, changed] = zeroUnmeasuredIsozymes(model);
    verifyEqual(testCase, changed.uniprot, {'P2';'P5'})
    verifyEqual(testCase, changed.newUpperBound, [0; 0])
    verifyEqual(testCase, usageUb(model, {'P1','P2','P3','P4','P5'}), [1000 0 1000 1000 0])
    verifyEqual(testCase, model.ec.concs, concs)
end

function testZeroUnmeasuredIsozymesDetectionLimit_tc0293(testCase)
    model = isozymeModel([3; NaN; 2; NaN; NaN]);
    % The default limit is the smallest concentration above 0
    modelDefault = zeroUnmeasuredIsozymes(model, 'method', 'detectionLimit');
    verifyEqual(testCase, usageUb(modelDefault, {'P2','P5'}), [2 2])
    modelLimit = zeroUnmeasuredIsozymes(model, 'method', 'detectionLimit', 'limit', 0.5);
    verifyEqual(testCase, usageUb(modelLimit, {'P2','P5','P4'}), [0.5 0.5 1000])
    % An upper bound below the limit is not raised
    model.ub(strcmp(model.rxns, 'usage_prot_P2')) = 0.1;
    modelLow = zeroUnmeasuredIsozymes(model, 'method', 'detectionLimit', 'limit', 0.5);
    verifyEqual(testCase, usageUb(modelLow, {'P2','P5'}), [0.1 0.5])
end

function testZeroUnmeasuredIsozymesNoChangeWithoutPartialGroup_tc0294(testCase)
    model = isozymeModel([0; NaN; 2; NaN; NaN]);
    ub = model.ub;
    [model, changed] = zeroUnmeasuredIsozymes(model);
    verifyEqual(testCase, model.ub, ub)
    verifyEqual(testCase, height(changed), 0)
end

function testZeroUnmeasuredIsozymesRejectsUnknownMethod_tc0295(testCase)
    model = isozymeModel([3; NaN; 2; NaN; NaN]);
    verifyError(testCase, @() zeroUnmeasuredIsozymes(model, 'method', 'half'), ...
        'zeroUnmeasuredIsozymes:method')
    verifyError(testCase, @() zeroUnmeasuredIsozymes(model, 'method', 'detectionLimit', 'limit', -1), ...
        'zeroUnmeasuredIsozymes:limit')
    noConcs = isozymeModel(nan(5, 1));
    verifyError(testCase, @() zeroUnmeasuredIsozymes(noConcs, 'method', 'detectionLimit'), ...
        'zeroUnmeasuredIsozymes:limit')
end

function testFlexibilizeEnzConcsLeavesMeasuredZeroAtZero_tc0296(testCase)
    % E1 is measured at 1 and limits growth; E2 is measured at 0 on a parallel
    % route. E1 is raised, E2 stays at 0 and the function warns about it.
    adapter = getAdapter();
    model = zeroToy([1 0]);
    verifyWarning(testCase, @() flexibilizeEnzConcs(model, 'expGrowth', 2, 'modelAdapter', adapter, 'verbose', false), ...
        'flexibilizeEnzConcs:measuredZero')
    warning('off', 'flexibilizeEnzConcs:measuredZero')
    restore = onCleanup(@() warning('on', 'flexibilizeEnzConcs:measuredZero'));
    [~, model, flexEnz] = evalc("flexibilizeEnzConcs(model, 'expGrowth', 2, 'modelAdapter', adapter, 'verbose', false)");
    verifyEqual(testCase, flexEnz.uniprotIDs, {'E1'})
    verifyEqual(testCase, usageUb(model, {'E2'}), 0)
    sol = solveLP(model);
    verifyEqual(testCase, sol.f, 2, 'RelTol', 5e-3)
end

function testFlexibilizeEnzConcsCannotRaiseMeasuredZero_tc0297(testCase)
    % The only route is measured at 0: no enzyme is limiting, the bound stays 0.
    % The test pins this behaviour, which a change that raises a bound of 0 breaks.
    adapter = getAdapter();
    model = zeroToy([NaN 0]);
    model.ub(strcmp(model.rxns, 'R_AB1')) = 0;
    warning('off', 'flexibilizeEnzConcs:measuredZero')
    restore = onCleanup(@() warning('on', 'flexibilizeEnzConcs:measuredZero'));
    [~, model, flexEnz] = evalc("flexibilizeEnzConcs(model, 'expGrowth', 2, 'modelAdapter', adapter, 'verbose', false)");
    verifyEmpty(testCase, flexEnz.uniprotIDs)
    verifyEqual(testCase, usageUb(model, {'E2'}), 0)
    verifyEqual(testCase, solveLP(model).f, 0, 'AbsTol', 1e-9)
end


function testGetIsozymeCoverageKeepsEnzymeOfMeasuredEntry_tc0314(testCase)
    % P2 is not measured, and it is in the measured entry X (with P1) as well
    % as in the unmeasured entry X_EXP_1, so it is shared and not zeroable.
    model = isozymeModel([3; NaN; 2; NaN; NaN]);
    model.ec.rxnEnzMat(1, :) = [1 1 0 0 0];
    [coverage, zeroable] = getIsozymeCoverage(model);
    verifyEqual(testCase, coverage.sharedEnzymes{strcmp(coverage.rxn, 'X')}, {'P2'})
    verifyEqual(testCase, zeroable, {'P5'})
end

function testLevelComplexConcsKeepsMeasuredZero_tc0315(testCase)
    % B is measured at 0 in a complex with A and C: it stays 0 and is not
    % filled with the level of the complex
    model.ec.enzymes   = {'A';'B';'C'};
    model.ec.rxns      = {'R1'};
    model.ec.genes     = {'g_A';'g_B';'g_C'};
    model.ec.mw        = [100; 100; 100];
    model.ec.concs     = [2; 0; 4];
    model.ec.rxnEnzMat = sparse([1 1 1]);
    [model, result] = levelComplexConcs(model);
    verifyEqual(testCase, model.ec.concs, [3; 0; 3], 'AbsTol', 1e-12)
    verifyEmpty(testCase, result.added)
end


function testLevelComplexConcsSkipsComplexWithOnlyZeroSubunits_tc0316(testCase)
    model.ec.enzymes   = {'A';'B';'C'};
    model.ec.rxns      = {'R1'};
    model.ec.genes     = {'g_A';'g_B';'g_C'};
    model.ec.mw        = [100; 100; 100];
    model.ec.concs     = [0; NaN; NaN];
    model.ec.rxnEnzMat = sparse([1 1 1]);
    [model, result] = levelComplexConcs(model);
    verifyEqual(testCase, model.ec.concs, [0; NaN; NaN])
    verifyEqual(testCase, numel(result.skipped), 1)
    verifyEmpty(testCase, result.added)
end

function testFlexibilizePerComplexKeepsMeasuredZeroOutOfTheReport_tc0317(testCase)
    % X1 has a concentration of 0 but an open usage bound, X2 limits the
    % complex. X2 is raised and carries the count; X1 is not reported.
    adapter = getAdapter();
    model = complexZeroToy([0 1], [1000 1]);
    warning('off', 'flexibilizeEnzConcs:measuredZero')
    restore = onCleanup(@() warning('on', 'flexibilizeEnzConcs:measuredZero'));
    [~, model, flexEnz] = evalc("flexibilizeEnzConcs(model, 'expGrowth', 2, 'modelAdapter', adapter, 'verbose', false, 'perComplex', true)");
    verifyEqual(testCase, flexEnz.uniprotIDs, {'X2'})
    verifyEqual(testCase, flexEnz.frequence, 1)
    verifyFalse(testCase, any(isnan(flexEnz.ratioIncr)))
    verifyEqual(testCase, usageUb(model, {'X1'}), 1000)
end


function adapter = getAdapter()
    adapter = ModelAdapterManager.getAdapter(fullfile(findGECKOroot,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
end

function ub = usageUb(model, enzymes)
    [~, idx] = ismember(strcat('usage_prot_', enzymes), model.rxns);
    ub = model.ub(idx)';
end

function model = isozymeModel(concs)
    % Five enzymes with usage reactions (upper bound 1000) and the entries
    % X, X_EXP_1: P1, P2;  Y, Y_EXP_1: P3, P4;  Z: P4;  W, W_EXP_1: P5, P1
    enzymes = {'P1';'P2';'P3';'P4';'P5'};
    model.rxns = strcat('usage_prot_', enzymes);
    model.ub   = 1000 * ones(5, 1);
    model.ec.geckoLight = false;
    model.ec.enzymes   = enzymes;
    model.ec.concs     = concs(:);
    model.ec.rxns      = {'X';'X_EXP_1';'Y';'Y_EXP_1';'Z';'W';'W_EXP_1'};
    model.ec.rxnEnzMat = sparse([1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; ...
                                 0 0 0 1 0; 0 0 0 0 1; 1 0 0 0 0]);
end

function model = zeroToy(concs)
    % A -> B by E1 or by E2 (R_AB1, R_AB2); growth (R4) equals the flux.
    mets = {'A_e';'A_c';'B_c';'prot_pool';'prot_E1';'prot_E2'};
    spec = {'EX_A', {'A_e',-1}, -1000, 0; ...
            'TR_A', {'A_e',-1,'A_c',1}, 0, 1000; ...
            'R_AB1', {'A_c',-1,'B_c',1,'prot_E1',-1}, 0, 1000; ...
            'R_AB2', {'A_c',-1,'B_c',1,'prot_E2',-1}, 0, 1000; ...
            'R4', {'B_c',-1}, 0, 1000; ...
            'prot_pool_exchange', {'prot_pool',1}, 0, 1000; ...
            'usage_prot_E1', {'prot_pool',-1,'prot_E1',1}, 0, 1000; ...
            'usage_prot_E2', {'prot_pool',-1,'prot_E2',1}, 0, 1000};
    S = zeros(numel(mets), size(spec, 1));
    for r = 1:size(spec, 1)
        entries = spec{r, 2};
        for n = 1:2:numel(entries)
            S(strcmp(mets, entries{n}), r) = entries{n+1};
        end
    end
    model.id       = 'toy';
    model.mets     = mets;
    model.metNames = mets;
    model.rxns     = spec(:, 1);
    model.rxnNames = spec(:, 1);
    model.S        = sparse(S);
    model.lb       = cell2mat(spec(:, 3));
    model.ub       = cell2mat(spec(:, 4));
    model.c        = double(strcmp(model.rxns, 'R4'));
    model.b        = zeros(numel(mets), 1);
    model.rev      = double(model.lb < 0);
    model.ec.enzymes   = {'E1';'E2'};
    model.ec.rxns      = {'R_AB1';'R_AB2'};
    model.ec.genes     = {'g_E1';'g_E2'};
    model.ec.mw        = [100; 100];
    model.ec.concs     = concs(:);
    model.ec.rxnEnzMat = sparse([1 0; 0 1]);
    model = constrainEnzConcs(model);
end

function model = complexZeroToy(concs, ubs)
    % A -> B is catalysed by the complex X1 + X2; growth (R4) equals the flux.
    mets = {'A_e';'A_c';'B_c';'prot_pool';'prot_X1';'prot_X2'};
    spec = {'EX_A', {'A_e',-1}, -1000, 0; ...
            'TR_A', {'A_e',-1,'A_c',1}, 0, 1000; ...
            'R_AB', {'A_c',-1,'B_c',1,'prot_X1',-1,'prot_X2',-1}, 0, 1000; ...
            'R4', {'B_c',-1}, 0, 1000; ...
            'prot_pool_exchange', {'prot_pool',1}, 0, 1000; ...
            'usage_prot_X1', {'prot_pool',-1,'prot_X1',1}, 0, ubs(1); ...
            'usage_prot_X2', {'prot_pool',-1,'prot_X2',1}, 0, ubs(2)};
    S = zeros(numel(mets), size(spec, 1));
    for r = 1:size(spec, 1)
        entries = spec{r, 2};
        for n = 1:2:numel(entries)
            S(strcmp(mets, entries{n}), r) = entries{n+1};
        end
    end
    model.id       = 'toy';
    model.mets     = mets;
    model.metNames = mets;
    model.rxns     = spec(:, 1);
    model.rxnNames = spec(:, 1);
    model.S        = sparse(S);
    model.lb       = cell2mat(spec(:, 3));
    model.ub       = cell2mat(spec(:, 4));
    model.c        = double(strcmp(model.rxns, 'R4'));
    model.b        = zeros(numel(mets), 1);
    model.rev      = double(model.lb < 0);
    model.ec.enzymes   = {'X1';'X2'};
    model.ec.rxns      = {'R_AB'};
    model.ec.genes     = {'g_X1';'g_X2'};
    model.ec.mw        = [100; 100];
    model.ec.concs     = concs(:);
    model.ec.rxnEnzMat = sparse([1 1]);
end
