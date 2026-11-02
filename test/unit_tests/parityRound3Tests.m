%run this test case with the command
%results = runtests('parityRound3Tests.m')

function tests = parityRound3Tests
    tests = functiontests(localfunctions);
end

function testApplyKcatConstraintsAllZeroClearsFullModel_tc0250(testCase)
    % With every kcat 0 the enzyme constraints of a full ecModel are removed
    % and no new ones are written. A light ecModel is returned unchanged.
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    ecModel.ec.kcat(:) = 10;
    ecModel = applyKcatConstraints(ecModel);
    protMets = startsWith(ecModel.mets, 'prot_') & ~strcmp(ecModel.mets, 'prot_pool');
    metabolic = ~startsWith(ecModel.rxns, 'usage_prot_') & ~strcmp(ecModel.rxns, 'prot_pool_exchange');
    verifyGreaterThan(testCase, nnz(ecModel.S(protMets, metabolic)), 0)
    ecModel.ec.kcat(:) = 0;
    cleared = applyKcatConstraints(ecModel);
    verifyEqual(testCase, nnz(cleared.S(protMets, metabolic)), 0)
    % the usage reactions and the pool exchange keep their fixed entries
    verifyEqual(testCase, cleared.S(:, ~metabolic), ecModel.S(:, ~metabolic))
end

function testApplyKcatConstraintsAllZeroClearsLightModel_tc0250light(testCase)
    % With every kcat 0 the pool coefficients of a light ecModel are cleared,
    % as the enzyme constraints of a full ecModel are.
    adapter = getAdapter();
    lecModel = makeEcModel(getGeckoTestModel(), true, adapter);
    lecModel.ec.kcat(:) = 10;
    lecModel = applyKcatConstraints(lecModel);
    poolRow = strcmp(lecModel.mets, 'prot_pool');
    verifyGreaterThan(testCase, nnz(lecModel.S(poolRow, :)), 1)
    lecModel.ec.kcat(:) = 0;
    cleared = applyKcatConstraints(lecModel);
    poolRow = strcmp(cleared.mets, 'prot_pool');
    metabolic = ~strcmp(cleared.rxns, 'prot_pool_exchange');
    verifyEqual(testCase, nnz(cleared.S(poolRow, metabolic)), 0)
    verifyEqual(testCase, cleared.S(:, ~metabolic), lecModel.S(:, ~metabolic))
end

function testSensitivityTuningStopsOnGrowthBelowTolerance_tc0251(testCase)
    % The solver reports a growth of 1 in the first iteration and 1 + 5e-10 in
    % the second: the gain is below 1e-9, so the loop stops after the second
    % iteration (3 solves with the one before the loop). An exact-equality
    % loop would go on until two growths are identical (4 solves).
    [ecModel, tuneAdapter] = tuningModel();
    cleanup = shadowSolveLP('growth', [1 1 1+5e-10]); %#ok<NASGU>
    tunedModel = sensitivityTuning(ecModel, 'desiredGrowthRate', 50, ...
        'modelAdapter', tuneAdapter, 'protToIgnore', {'P4'}, 'verbose', false);
    verifyEqual(testCase, getappdata(0, 'gecko_test_solvelp_count'), 3)
    r5 = tunedModel.ec.kcat(strcmp(tunedModel.ec.rxns, 'R5'));
    verifyGreaterThan(testCase, r5, 200)
end

function testSensitivityTuningEndsWhenTheLPBecomesInfeasible_tc0254(testCase)
    % The third solve (the second iteration) has no solution: RAVEN returns
    % an empty objective value. The loop ends with a warning and returns
    % what was tuned in the first iteration.
    [ecModel, tuneAdapter] = tuningModel();
    cleanup = shadowSolveLP('infeasibleFrom', 3); %#ok<NASGU>
    out = evalc(['tunedModel = sensitivityTuning(ecModel, ''desiredGrowthRate'', 50, ' ...
        '''modelAdapter'', tuneAdapter, ''protToIgnore'', {''P4''}, ''verbose'', false);']);
    verifyTrue(testCase, contains(out, 'no feasible solution'))
    verifyEqual(testCase, getappdata(0, 'gecko_test_solvelp_count'), 3)
    r5 = tunedModel.ec.kcat(strcmp(tunedModel.ec.rxns, 'R5'));
    verifyGreaterThan(testCase, r5, 200)
end

function testSensitivityTuningStopsWhenNothingCarriesFlux_tc0253(testCase)
    % R5, the growth reaction, is blocked: growth is 0 from the start and no
    % usage reaction carries flux. The loop ends with a warning instead of
    % tuning an arbitrary reaction or failing, and no kcat changes.
    adapter = getAdapter();
    model = loadConventionalGEM('modelAdapter', adapter);
    for id = {'R2', 'R4', 'R5'}
        idx = strcmp(model.rxns, id{1});
        model.lb(idx) = 0; model.ub(idx) = 0;
    end
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R3')) = 0.05;
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
    tuneAdapter = adapter;
    tuneAdapter.params.bioRxn = 'R5';
    [tunedModel, tunedKcats] = sensitivityTuning(ecModel, 'desiredGrowthRate', 50, ...
        'modelAdapter', tuneAdapter, 'verbose', false);
    verifyEqual(testCase, tunedModel.ec.kcat, ecModel.ec.kcat)
    verifyEmpty(testCase, tunedKcats.rxns)
end

function testGetKcatAcrossIsozymesUsesPositiveSiblingsOnly_tc0252(testCase)
    % Only isozymes with a kcat above 0 enter the mean: a negative or NaN
    % sibling does not fill a kcat of 0.
    adapter = getAdapter();
    ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
    ecModel.ec.kcat(:) = 0;
    ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R3')) = 8;
    % R2: a zero next to a negative sibling is not filled
    ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R2_EXP_2')) = -3;
    % R2_REV: a zero next to a NaN sibling is not filled
    ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R2_REV_EXP_2')) = NaN;
    out = getKcatAcrossIsozymes(ecModel);
    verifyEqual(testCase, out.ec.kcat(strcmp(out.ec.rxns, 'R2_EXP_1')), 0)
    verifyEqual(testCase, out.ec.kcat(strcmp(out.ec.rxns, 'R2_REV_EXP_1')), 0)
    % a positive sibling fills the zero (control: the same before and after)
    ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R2_EXP_2')) = 6;
    out = getKcatAcrossIsozymes(ecModel);
    j = strcmp(out.ec.rxns, 'R2_EXP_1');
    verifyEqual(testCase, out.ec.kcat(j), 6)
    verifyEqual(testCase, out.ec.source{j}, 'isozymes')
    % a NaN sibling next to a positive one is left out of the mean: three
    % isozymes of R2 with kcats 0, NaN and 6 give 6 to the first
    ecModel.ec.kcat(:) = 0;
    ecModel.ec.rxns{strcmp(ecModel.ec.rxns, 'R3')} = 'R2_EXP_3';
    ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R2_EXP_2')) = NaN;
    ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R2_EXP_3')) = 6;
    out = getKcatAcrossIsozymes(ecModel);
    j = strcmp(out.ec.rxns, 'R2_EXP_1');
    verifyEqual(testCase, out.ec.kcat(j), 6)
    verifyEqual(testCase, out.ec.source{j}, 'isozymes')
    verifyTrue(testCase, isnan(out.ec.kcat(strcmp(out.ec.rxns, 'R2_EXP_2'))))
end

function testSigmaFitterReturnsStartOfGrowthPlateau_tc0255(testCase)
% The growth target is above what the uptake-limited model reaches, so growth
% is flat from sigma 0.06 on. GLPK returns the plateau growth with rounding
% noise (49.999999999999986 at sigma 0.06, 50 at sigma 1). The fitted sigma is
% the start of the plateau, 0.06, whatever the noise, and the same on Gurobi.
previous = getpref('RAVEN', 'solver', 'glpk');
restore = onCleanup(@() setRavenSolver(previous)); %#ok<NASGU>
setRavenSolver('glpk');
sigmaGlpk = plateauSigma();
verifyEqual(testCase, sigmaGlpk, 0.06, 'AbsTol', 1e-12)
hasGurobi = false;
if exist('gurobi', 'file')
    try
        setRavenSolver('gurobi');
        hasGurobi = true;
    catch
    end
end
assumeTrue(testCase, hasGurobi, 'Gurobi is not available')
verifyEqual(testCase, plateauSigma(), sigmaGlpk, 'AbsTol', 1e-12)
end

function testGetStandardKcatSubSystemNamesIgnoreCase_tc0256(testCase)
% R3 is in 'Alpha', R5 in 'alpha', the entries of R2 in 'Beta'. The names
% 'Alpha' and 'alpha' are one subSystem, so R6 ('ALPHA') gets the median of
% R3 and R5 (5) and not the model-wide median (100).
[ecModel, adapter] = standardKcatModel(struct('R5', 'alpha', 'R2', 'Beta', 'R6', 'ALPHA'));
ecModel = getStandardKcat(ecModel, 'modelAdapter', adapter, 'threshold', 2);
verifyEqual(testCase, ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R6')), 5)
end

function testGetStandardKcatReactionWithoutSubSystem_tc0257(testCase)
% R7 has an empty subSystem while the other reactions have one: it gets the
% model-wide median kcat (100), R6 the one of its subSystem 'Alpha' (5).
[ecModel, adapter] = standardKcatModel(struct('R2', 'Beta', 'R7', ''));
ecModel = getStandardKcat(ecModel, 'modelAdapter', adapter, 'threshold', 2);
verifyEqual(testCase, ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R7')), 100)
verifyEqual(testCase, ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R6')), 5)
end

function testGetSubsetEcModelKeepsEverythingForTheFullModel_tc0258(testCase)
% The context model is the conventional model, so no gene is removed and the
% ecModel is returned as it is, for the full and the light ecModel. The
% reaction identifiers of the light ecModel have a prefix; they are matched
% by the reaction they stand for, so they are not all removed.
adapter = getAdapter();
model = loadConventionalGEM('modelAdapter', adapter);
for light = [false true]
    bigEcModel = makeEcModel(model, light, adapter);
    smallEcModel = getSubsetEcModel(bigEcModel, model);
    verifyEqual(testCase, smallEcModel, bigEcModel)
    verifyNotEmpty(testCase, smallEcModel.ec.rxns)
end
end

function testGetSubsetEcModelLightDropsEcEntriesOfRemovedRxns_tc0259(testCase)
% The context model lacks R2: the light ecModel loses the entry 001_R2 of
% its ec structure and keeps the others.
adapter = getAdapter();
model = loadConventionalGEM('modelAdapter', adapter);
small = removeReactions(model, {'R2'}, 'removeUnusedMets', true, 'removeUnusedGenes', true);
bigEcModel = makeEcModel(model, true, adapter);
smallEcModel = getSubsetEcModel(bigEcModel, small);
verifyFalse(testCase, ismember('001_R2', smallEcModel.ec.rxns))
verifyTrue(testCase, all(ismember({'001_R3', '001_R5'}, smallEcModel.ec.rxns)))
verifyEqual(testCase, size(smallEcModel.ec.rxnEnzMat, 1), numel(smallEcModel.ec.rxns))
end

function [ecModel, adapter] = standardKcatModel(subSystems)
% ecTestGEM with R6 and R7 added (no genes) and the subSystems given as a
% struct of reaction id to name (an empty name is no subSystem); reactions
% not listed are in 'Alpha'. The kcats
% are R3 2, R5 8, the entries of R2 100.
adapter = getAdapter();
model = loadConventionalGEM('modelAdapter', adapter);
rxnsToAdd.rxns      = {'R6'; 'R7'};
rxnsToAdd.equations = {'m1[c] => m2[c]'; 'm1[c] => m2[c]'};
rxnsToAdd.grRules   = {''; ''};
model = addRxns(model, rxnsToAdd, 3);
model.subSystems = cell(numel(model.rxns), 1);
for i = 1:numel(model.rxns)
    name = 'Alpha';
    if isfield(subSystems, model.rxns{i})
        name = subSystems.(model.rxns{i});
    end
    if isempty(name)
        model.subSystems{i} = {};
    else
        model.subSystems{i} = {name};
    end
end
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);
ecModel.ec.kcat(:) = 100;
ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R3')) = 2;
ecModel.ec.kcat(strcmp(ecModel.ec.rxns, 'R5')) = 8;
ecModel.ec.source(:) = {'manual'};
end

function sigma = plateauSigma()
adapter = getAdapter();
model = loadConventionalGEM('modelAdapter', adapter);
model.lb(strcmp(model.rxns, 'S1')) = -50;
for id = {'R2', 'R4'}
    idx = strcmp(model.rxns, id{1});
    model.lb(idx) = 0; model.ub(idx) = 0;
end
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);
kcats = struct('R3', 5, 'R5', 200);
for id = fieldnames(kcats)'
    j = strcmp(ecModel.ec.rxns, id{1});
    ecModel.ec.kcat(j) = kcats.(id{1});
    ecModel.ec.source{j} = 'manual';
end
ecModel.c = double(strcmp(ecModel.rxns, 'R5'));
ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
[~, sigma] = sigmaFitter(ecModel, 'growthRate', 100, 'Ptot', 0.5, 'f', 4, ...
    'makePlot', false, 'modelAdapter', adapter);
end

function [ecModel, tuneAdapter] = tuningModel()
% R3 limits growth, its enzyme P4 is ignored by the tests, so tuning raises R5.
adapter = getAdapter();
model = loadConventionalGEM('modelAdapter', adapter);
for id = {'R2', 'R4'}
    idx = strcmp(model.rxns, id{1});
    model.lb(idx) = 0; model.ub(idx) = 0;
end
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);
kcats = struct('R3', 0.05, 'R5', 200);
for id = fieldnames(kcats)'
    j = strcmp(ecModel.ec.rxns, id{1});
    ecModel.ec.kcat(j) = kcats.(id{1});
    ecModel.ec.source{j} = 'manual';
end
ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
tuneAdapter = adapter;
tuneAdapter.params.bioRxn = 'R5';
end

function cleanup = shadowSolveLP(mode, values)
% Puts a solveLP on the path that calls RAVEN's and then reports the growth
% of solve number n as values(min(n, end)) ('growth'), or no solution from
% solve number values on ('infeasibleFrom'). The number of solves is kept in
% the root appdata gecko_test_solvelp_count. The path is restored when the
% returned object is cleared.
setappdata(0, 'gecko_test_solvelp_real', @solveLP);
setappdata(0, 'gecko_test_solvelp_count', 0);
setappdata(0, 'gecko_test_solvelp_mode', mode);
setappdata(0, 'gecko_test_solvelp_values', values);
stubDir = tempname;
mkdir(stubDir);
fid = fopen(fullfile(stubDir, 'solveLP.m'), 'w');
fprintf(fid, '%s\n', ...
    'function [res,hs] = solveLP(varargin)', ...
    'n = getappdata(0,''gecko_test_solvelp_count'') + 1;', ...
    'setappdata(0,''gecko_test_solvelp_count'',n);', ...
    'real = getappdata(0,''gecko_test_solvelp_real'');', ...
    'mode = getappdata(0,''gecko_test_solvelp_mode'');', ...
    'values = getappdata(0,''gecko_test_solvelp_values'');', ...
    '[res,hs] = real(varargin{:});', ...
    'switch mode', ...
    '    case ''growth''', ...
    '        res.f = values(min(n,numel(values)));', ...
    '    case ''infeasibleFrom''', ...
    '        if n >= values', ...
    '            res.f = [];', ...
    '            res.x = [];', ...
    '        end', ...
    'end', ...
    'end');
fclose(fid);
addpath(stubDir, '-begin');
cleanup = onCleanup(@() removeShadow(stubDir));
end

function removeShadow(stubDir)
rmpath(stubDir);
rmdir(stubDir, 's');
for name = {'real', 'mode', 'values'}
    if isappdata(0, ['gecko_test_solvelp_' name{1}])
        rmappdata(0, ['gecko_test_solvelp_' name{1}]);
    end
end
end

function adapter = getAdapter()
geckoPath = findGECKOroot;
adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
end
