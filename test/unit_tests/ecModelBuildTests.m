%run this test case with the command
%results = runtests('ecModelBuildTests.m')

function tests = ecModelBuildTests
    tests = functiontests(localfunctions);
end

function testLightNestedGeneRuleExpandsToIsozymes_tc0210(testCase)
    % 'G1 and (G2 or G3)' is the two isozymes {G1,G2} and {G1,G3}; the
    % brackets must not be dropped, which would make the second one {G3}.
    adapter = getAdapter();
    model = getGeckoTestModel();
    model = changeGrRules(model, {'R2'}, {'G1 and (G2 or G3)'}, true);
    ecModel = makeEcModel(model, true, adapter);
    verifyEqual(testCase, sort(coupledGenes(ecModel, '001_R2')), {'G1';'G2'})
    verifyEqual(testCase, sort(coupledGenes(ecModel, '002_R2')), {'G1';'G3'})
    verifyEqual(testCase, sum(strcmp(ecModel.ec.rxns, '003_R2')), 0)
end

function testAddNewRxnsToECNestedRuleAndGeneIdWithOr_tc0211(testCase)
    % The gene ids of the rule are whole words, so 'orfX' is a gene and not
    % the operator 'or' followed by 'fX', and the nested rule is split into
    % one reaction per isozyme: 'orfX and GNEW2' and 'orfX and GNEW3'.
    adapter = getAdapter();
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);

    newRxns.rxns      = {'RNEWA'};
    newRxns.rxnNames  = {'RNEWA'};
    newRxns.equations = {'m1[c] => e2[e]'};
    newRxns.grRules   = {'orfX and (GNEW2 or GNEW3)'};
    newEnzymes.enzymes = {'ENEW1';'ENEW2';'ENEW3'};
    newEnzymes.genes   = {'orfX';'GNEW2';'GNEW3'};
    newEnzymes.mw      = [15000;16000;17000];

    ecModel = addNewRxnsToEC(ecModel, newRxns, newEnzymes, adapter);
    verifyEqual(testCase, ecModel.grRules(strcmp(ecModel.rxns, 'RNEWA_EXP_1')), {'orfX and GNEW2'})
    verifyEqual(testCase, ecModel.grRules(strcmp(ecModel.rxns, 'RNEWA_EXP_2')), {'orfX and GNEW3'})
end

function testMakeEcModelCouplesTheKnownSubunitsOfAComplex_tc0212(testCase)
    % GX is in neither database: the reaction is still coupled to the enzyme
    % of G4, as the databases know it.
    adapter = getAdapter();
    model = getGeckoTestModel();
    model = changeGrRules(model, {'R3'}, {'G4 and GX'}, true);
    ecModel = makeEcModel(model, false, adapter);
    verifyEqual(testCase, coupledGenes(ecModel, 'R3'), {'G4'})
end

function testMapRxnsToConvStripsOnlyTheSuffixes_tc0213(testCase)
    % R_REVX is a conventional reaction whose id contains _REV; its reverse
    % variant R_REVX_REV carries flux that maps back to it with the opposite sign.
    adapter = getAdapter();
    model = getGeckoTestModel();
    model.rxns{strcmp(model.rxns, 'R1')} = 'R_REVX';
    ecModel = makeEcModel(model, false, adapter);
    fluxVect = zeros(numel(ecModel.rxns), 1);
    fluxVect(strcmp(ecModel.rxns, 'R_REVX')) = 5;
    fluxVect(strcmp(ecModel.rxns, 'R_REVX_REV')) = 2;
    mappedFlux = mapRxnsToConv(ecModel, model, fluxVect);
    verifyEqual(testCase, mappedFlux(strcmp(model.rxns, 'R_REVX')), 3)
end

function testEcFVAStripsOnlyTheSuffixes_tc0214(testCase)
    % The same reaction id in ecFVA, which groups the variants of each
    % conventional reaction by id.
    adapter = getAdapter();
    model = getGeckoTestModel();
    model.rxns{strcmp(model.rxns, 'R1')} = 'R_REVX';
    ecModel = makeEcModel(model, false, adapter);
    [minFlux, maxFlux] = ecFVA(ecModel, model, 'runParallel', false);
    verifySize(testCase, minFlux, [numel(model.rxns) 1])
    verifySize(testCase, maxFlux, [numel(model.rxns) 1])
end

function testApplyKcatConstraintsLightSkipsIsozymeWithoutEnzymes_tc0215(testCase)
    % R2 has two isozymes; the second one has no enzyme but a kcat. It cannot
    % be the cheapest isozyme: the cost of R2 is the one of the first.
    adapter = getAdapter();
    model = getGeckoTestModel();
    lecModel = makeEcModel(model, true, adapter);
    lecModel.ec.kcat(:) = 10;
    first  = strcmp(lecModel.ec.rxns, '001_R2');
    second = strcmp(lecModel.ec.rxns, '002_R2');
    lecModel.ec.rxnEnzMat(second,:) = 0;
    lecModel = applyKcatConstraints(lecModel);
    expected = -(lecModel.ec.rxnEnzMat(first,:) * lecModel.ec.mw) / 10 / 3600;
    actual = full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'), strcmp(lecModel.rxns, 'R2')));
    verifyGreaterThan(testCase, abs(expected), 0)
    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-10)
end


function adapter = getAdapter()
geckoPath = findGECKOroot;
adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
end

function genes = coupledGenes(ecModel, ecRxn)
% The genes of the enzymes that ec reaction ecRxn is coupled to.
row = strcmp(ecModel.ec.rxns, ecRxn);
genes = ecModel.ec.genes(logical(ecModel.ec.rxnEnzMat(row,:)));
end
