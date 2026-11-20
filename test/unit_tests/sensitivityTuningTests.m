%run this test case with the command
%results = runtests('sensitivityTuningTests.m')

function tests = sensitivityTuningTests
    tests = functiontests(localfunctions);
end

function testGrowthTargetReportsStatus(testCase)
    [ecModel, adapter] = buildTuningModel(testCase);
    [~, tunedKcats] = sensitivityTuning(ecModel, 'modelAdapter', adapter, ...
        'desiredGrowthRate', 120, 'verbose', false);
    verifyEqual(testCase, tunedKcats.status, 'reached')
    verifyEqual(testCase, tunedKcats.rxns, {'R5'})
end

function testConditionRaisesTheKcatThatLimitsIt(testCase)
    [ecModel, adapter] = buildTuningModel(testCase);
    fluxData = buildFluxData(120, -120);
    [tunedModel, tunedKcats] = sensitivityTuning(ecModel, 'modelAdapter', adapter, ...
        'fluxData', fluxData, 'verbose', false);
    verifyEqual(testCase, tunedKcats.status, 'reached')
    verifyEqual(testCase, tunedKcats.rxns, {'R5'})
    verifyEqual(testCase, tunedKcats.oldKcat, 10)
    verifyEqual(testCase, tunedKcats.newKcat, 100)
    %the model keeps the tuned kcat and none of the condition's bounds
    verifyEqual(testCase, tunedModel.ec.kcat(strcmp(tunedModel.ec.rxns, 'R5')), 100)
    verifyEqual(testCase, tunedModel.lb(strcmp(tunedModel.rxns, 'S1')), ecModel.lb(strcmp(ecModel.rxns, 'S1')))
    verifyEqual(testCase, tunedModel.lb(strcmp(tunedModel.rxns, 'R5')), ecModel.lb(strcmp(ecModel.rxns, 'R5')))
end

function testConditionByNameAndGrowthOverride(testCase)
    [ecModel, adapter] = buildTuningModel(testCase);
    fluxData = buildFluxData(1, -120);
    fluxData.conds = {'chemostat'};
    [~, tunedKcats] = sensitivityTuning(ecModel, 'modelAdapter', adapter, ...
        'fluxData', fluxData, 'condition', 'chemostat', 'desiredGrowthRate', 120, ...
        'verbose', false);
    verifyEqual(testCase, tunedKcats.status, 'reached')
    verifyEqual(testCase, tunedKcats.rxns, {'R5'})
end

function testConditionInfeasibleInTheGemIsMarkedAndLeftAlone(testCase)
    %growth 120 needs 120 of e1, and the measured uptake is 50
    [ecModel, adapter] = buildTuningModel(testCase);
    fluxData = buildFluxData(120, -50);
    output = evalc(['[tunedModel, tunedKcats] = sensitivityTuning(ecModel, ' ...
        '''modelAdapter'', adapter, ''fluxData'', fluxData, ''verbose'', false);']);
    verifySubstring(testCase, output, 'infeasible without any enzyme limit')
    verifyEqual(testCase, tunedKcats.status, 'infeasibleInGem')
    verifyEmpty(testCase, tunedKcats.rxns)
    verifyEqual(testCase, tunedModel.ec.kcat, ecModel.ec.kcat)
end

function testPerFluxToleranceWidensTheBand(testCase)
    %uptake 80 does not allow growth 120 within 10%, but does within 50%
    [ecModel, adapter] = buildTuningModel(testCase);
    fluxData = buildFluxData(120, -80);
    [~, tight] = sensitivityTuning(ecModel, 'modelAdapter', adapter, ...
        'fluxData', fluxData, 'verbose', false);
    verifyEqual(testCase, tight.status, 'infeasibleInGem')
    [~, loose] = sensitivityTuning(ecModel, 'modelAdapter', adapter, ...
        'fluxData', fluxData, 'fluxTolerances', {'S1', 0.5}, 'verbose', false);
    verifyEqual(testCase, loose.status, 'reached')
end

function testInfeasibleConditionsListsTheNames(testCase)
    [ecModel, adapter] = buildTuningModel(testCase);
    fluxData = buildFluxData(120, -120);
    fluxData.conds = {'ok'};
    fluxData.grRate = 120;
    fluxData.exchFluxes = [-120];
    both = fluxData;
    both.conds = {'ok'; 'tooLow'};
    both.Ptot = [NaN; NaN];
    both.grRate = [120; 120];
    both.exchFluxes = [-120; -50];
    names = infeasibleConditions(ecModel, both, 'modelAdapter', adapter);
    verifyEqual(testCase, names, {'tooLow'})
end

function testDiagnoseNamesTheMeasurementThatHasToMove(testCase)
    %growth 120 needs 120 of e1, and the measured uptake is 50: growth is the
    %cheaper value to move, relative to its size
    [ecModel, adapter] = buildTuningModel(testCase);
    conflicts = diagnoseCondition(ecModel, buildFluxData(120, -50), 'modelAdapter', adapter);
    verifyEqual(testCase, conflicts.reaction{1}, 'R5')
    verifyGreaterThan(testCase, conflicts.stretch(1), 0)
    verifyEqual(testCase, conflicts.measured(1), 120)
    verifyEqual(testCase, height(conflicts), 2)
end

function testInfeasibleConditionCarriesTheConflicts(testCase)
    [ecModel, adapter] = buildTuningModel(testCase);
    evalc(['[~, tunedKcats] = sensitivityTuning(ecModel, ''modelAdapter'', adapter, ' ...
        '''fluxData'', buildFluxData(120, -50), ''verbose'', false);']);
    verifyEqual(testCase, tunedKcats.status, 'infeasibleInGem')
    verifyEqual(testCase, tunedKcats.conflicts.reaction{1}, 'R5')
end

function testReachedConditionHasNoConflicts(testCase)
    [ecModel, adapter] = buildTuningModel(testCase);
    evalc(['[~, tunedKcats] = sensitivityTuning(ecModel, ''modelAdapter'', adapter, ' ...
        '''fluxData'', buildFluxData(120, -120), ''verbose'', false);']);
    verifyEmpty(testCase, tunedKcats.conflicts)
end

function testSingleEnzymeLimitIsUnfixable(testCase)
    %growth 120 is feasible without enzyme limits, but the usage of P5 is
    %capped far below what it needs, so lifting the pool does not help
    [ecModel, adapter] = buildTuningModel(testCase);
    ecModel.ub(strcmp(ecModel.rxns, 'usage_prot_P5')) = 1;
    fluxData = buildFluxData(120, -120);
    evalc(['[~, tunedKcats] = sensitivityTuning(ecModel, ''modelAdapter'', adapter, ' ...
        '''fluxData'', fluxData, ''verbose'', false);']);
    verifyEqual(testCase, tunedKcats.status, 'unfixable')
end

function testUnknownConditionAndExchangeRaise(testCase)
    [ecModel, adapter] = buildTuningModel(testCase);
    fluxData = buildFluxData(120, -120);
    verifyError(testCase, @() sensitivityTuning(ecModel, 'modelAdapter', adapter, ...
        'fluxData', fluxData, 'condition', 'missing', 'verbose', false), ...
        'sensitivityTuning:unknownCondition')
    fluxData.exchRxnIDs = {'noSuchRxn'};
    verifyError(testCase, @() sensitivityTuning(ecModel, 'modelAdapter', adapter, ...
        'fluxData', fluxData, 'verbose', false), 'sensitivityTuning:missingExchange')
end


function [ecModel, adapter] = buildTuningModel(testCase)
    %Test model with kcat 10 for all enzymes and a protein pool of 0.5. Growth
    %is the flux of R5, which is 90 at most before any tuning.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    oldBioRxn = adapter.params.bioRxn;
    adapter.params.bioRxn = 'R5';
    testCase.addTeardown(@() setBioRxn(adapter, oldBioRxn));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(:) = 10;
    ecModel.ec.source(:) = {'manual'};
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel, 0.5, 0.5, 0.5, adapter);
end

function setBioRxn(adapter, bioRxn)
    adapter.params.bioRxn = bioRxn;
end

function fluxData = buildFluxData(growth, uptake)
    fluxData = struct();
    fluxData.conds      = {'c1'};
    fluxData.Ptot       = NaN;
    fluxData.grRate     = growth;
    fluxData.exchFluxes = uptake;
    fluxData.exchMets   = {'e1'};
    fluxData.exchRxnIDs = {'S1'};
end
