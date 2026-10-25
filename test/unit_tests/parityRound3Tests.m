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

function testApplyKcatConstraintsAllZeroLightIsUnchanged_tc0250light(testCase)
    % A light ecModel with every kcat 0 is returned unchanged (control: the
    % behaviour is the same before and after the full-model change).
    adapter = getAdapter();
    lecModel = makeEcModel(getGeckoTestModel(), true, adapter);
    lecModel.ec.kcat(:) = 10;
    lecModel = applyKcatConstraints(lecModel);
    poolRow = strcmp(lecModel.mets, 'prot_pool');
    verifyGreaterThan(testCase, nnz(lecModel.S(poolRow, :)), 1)
    lecModel.ec.kcat(:) = 0;
    same = applyKcatConstraints(lecModel);
    verifyEqual(testCase, same.S, lecModel.S)
end

function testSensitivityTuningStopsOnGrowthBelowTolerance_tc0251(testCase)
    % R3 limits growth and its enzyme P4 is ignored, so tuning only raises
    % R5. Growth stops rising once R3 alone limits it: the loop ends on the
    % first iteration whose growth gain is below 1e-9, instead of raising
    % the kcat of R5 until growth is exactly unchanged.
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
    tunedModel = sensitivityTuning(ecModel, 'desiredGrowthRate', 50, ...
        'modelAdapter', tuneAdapter, 'protToIgnore', {'P4'}, 'verbose', false);
    r5 = tunedModel.ec.kcat(strcmp(tunedModel.ec.rxns, 'R5'));
    verifyGreaterThan(testCase, r5, 200)
    verifyLessThan(testCase, r5, 1e12)
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
end

function adapter = getAdapter()
geckoPath = findGECKOroot;
adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
end
