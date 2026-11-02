%run this test case with the command
%results = runtests('levelComplexTests.m')

function tests = levelComplexTests
    tests = functiontests(localfunctions);
end

function testComplexFamiliesGroupIsozymesAndFlagSharedSubunits_tc0262(testCase)
    % B is in an isozyme variant with C, so A, B and C are one complex; S is
    % in two complexes of equal weight and belongs to none.
    model = ecOnly(nan(9, 1), toyRows());
    families = getComplexFamilies(model);
    verifyEqual(testCase, sortedComplexes(families), {{'A','B','C'}, {'D','F','S'}, {'E','G','S'}})
    verifyEqual(testCase, families.shared, {'S'})
    [~, b] = ismember('B', model.ec.enzymes);
    [~, s] = ismember('S', model.ec.enzymes);
    [~, h] = ismember('H', model.ec.enzymes);
    verifyEqual(testCase, families.complexes{families.owner(b)}.enzymes(families.complexes{families.owner(b)}.owned), {'A';'B';'C'})
    verifyEqual(testCase, families.owner(s), 0)
    verifyEqual(testCase, families.owner(h), 0)
end

function testComplexFamiliesEnzymeThatWorksAloneIsShared_tc0263(testCase)
    rows = toyRows();
    rows(end+1, :) = {'R5', struct('A', 1)};
    families = getComplexFamilies(ecOnly(nan(9, 1), rows));
    verifyTrue(testCase, ismember('A', families.shared))
end

function testComplexFamiliesCarrierDoesNotBridgeComplexes_tc0264(testCase)
    % B, C, D take carrier S or carrier H (two variants); S also pairs with
    % E, like cytochrome c with a dehydrogenase. A single shared protein does
    % not merge the two.
    rows = {'R1_EXP_1', struct('B',1,'C',1,'D',1,'S',1); ...
            'R1_EXP_2', struct('B',1,'C',1,'D',1,'H',1); ...
            'R2',       struct('E',1,'S',1)};
    families = getComplexFamilies(ecOnly(nan(9, 1), rows));
    verifyEqual(testCase, sortedComplexes(families), {{'B','C','D','H','S'}, {'E','S'}})
    verifyTrue(testCase, ismember('S', families.shared))
end

function testComplexFamiliesSubunitInTwoComplexesGoesToTheHeavier_tc0265(testCase)
    % A, B, C, D is a complex on three reactions; A, B with E, F, G, H on
    % one (an assembly step): A and B belong to the first complex.
    big = struct('A',1,'B',1,'C',1,'D',1);
    rows = {'R1', big; 'R2', big; 'R3', big; ...
            'R4', struct('A',1,'B',1,'E',1,'F',1,'G',1,'H',1)};
    model = ecOnly(nan(9, 1), rows);
    families = getComplexFamilies(model);
    verifyEmpty(testCase, families.shared)
    verifyEqual(testCase, sort(ownedBy(families, model, 'A')), {'A';'B';'C';'D'})
    verifyEqual(testCase, sort(ownedBy(families, model, 'E')), {'E';'F';'G';'H'})
    % with equal weights the two subunits are shared
    rows = {'R1', big; 'R2', struct('A',1,'B',1,'E',1,'F',1,'G',1,'H',1)};
    families = getComplexFamilies(ecOnly(nan(9, 1), rows));
    verifyTrue(testCase, all(ismember({'A','B'}, families.shared)))
end

function testComplexFamiliesWarnWhenVariantsDisagreeOnCopies_tc0266(testCase)
    rows = {'R1', struct('A',1,'B',2); 'R2', struct('A',1,'B',3,'C',1)};
    model = ecOnly(nan(9, 1), rows);
    verifyWarning(testCase, @() getComplexFamilies(model), 'getComplexFamilies:copiesDisagree')
    warning('off', 'getComplexFamilies:copiesDisagree');
    families = getComplexFamilies(model);
    warning('on', 'getComplexFamilies:copiesDisagree');
    [~, b] = ismember('B', families.complexes{1}.enzymes);
    verifyEqual(testCase, families.complexes{1}.copies(b), 2)
end

function testLevelComplexConcsMedianScalesByCopies_tc0267(testCase)
    % Molar amount per copy: A 1, B 8 / 2 = 4, C 2; the median 2 gives
    % A = 2 * 1 copy, B = 2 * 2 copies, C = 2 * 1 copy.
    concs = nan(9, 1);
    concs(1) = 100; concs(2) = 1600; concs(3) = 100;
    model = ecOnly(concs, toyRows());
    [model, result] = levelComplexConcs(model);
    verifyEqual(testCase, model.ec.concs(1:3), [200; 800; 100], 'RelTol', 1e-12)
    verifyTrue(testCase, isnan(model.ec.concs(9)))
    verifyEqual(testCase, sort(result.adjusted), {'A';'B'})
    verifyEmpty(testCase, result.added)
end

function testLevelComplexConcsMeanIsAvailable_tc0268(testCase)
    % Molar amount per copy: A 1, B 4 / 2 = 2, C 6; the mean is 3.
    concs = nan(9, 1);
    concs(1) = 100; concs(2) = 800; concs(3) = 300;
    model = levelComplexConcs(ecOnly(concs, toyRows()), 'statistic', 'mean');
    verifyEqual(testCase, model.ec.concs(1:3), [300; 1200; 150], 'RelTol', 1e-12)
    model = levelComplexConcs(ecOnly(concs, toyRows()), 'median');
    verifyEqual(testCase, model.ec.concs(1:3), [200; 800; 100], 'RelTol', 1e-12)
end

function testLevelComplexConcsSharedSubunitStaysAndMissingOneIsFilled_tc0269(testCase)
    concs = nan(9, 1);
    concs(4) = 100;      % D
    concs(6) = 900;      % S, shared between two complexes
    [model, result] = levelComplexConcs(ecOnly(concs, toyRows()));
    verifyEqual(testCase, model.ec.concs([4 5 6]), [100; 100; 900])
    verifyEqual(testCase, result.shared, {'S'})
    verifyEqual(testCase, result.added, {'F'})
    % filling can be switched off
    model = levelComplexConcs(ecOnly(concs, toyRows()), 'fillMissing', false);
    verifyTrue(testCase, isnan(model.ec.concs(5)))
    verifyEqual(testCase, model.ec.concs(6), 900)
end

function testLevelComplexConcsFillsSubunitOfIsozymeVariant_tc0270(testCase)
    % C is in the second variant and is not measured: it gets the
    % median per copy (A 1, B 4 / 2 = 2: 1.5) at its molecular weight.
    concs = nan(9, 1);
    concs(1) = 100; concs(2) = 800;
    [model, result] = levelComplexConcs(ecOnly(concs, toyRows()));
    verifyEqual(testCase, model.ec.concs(1:3), [150; 600; 75], 'RelTol', 1e-12)
    verifyEqual(testCase, result.added, {'C'})
end

function testLevelComplexConcsReportsComplexWithoutMeasurement_tc0271(testCase)
    concs = nan(9, 1);
    concs(1) = 100;      % complex A, B, C is measured
    [model, result] = levelComplexConcs(ecOnly(concs, toyRows()));
    verifyEqual(testCase, result.skipped, {{'D','F','S'}; {'E','G','S'}})
    verifyTrue(testCase, all(isnan(model.ec.concs(4:9))))
end

function testLevelComplexConcsRejectsUnknownStatistic_tc0272(testCase)
    model = ecOnly(nan(9, 1), toyRows());
    verifyError(testCase, @() levelComplexConcs(model, 'statistic', 'max'), 'levelComplexConcs:statistic')
end

function testFlexibilizePerComplexRaisesSubunitsByTheSameFactor_tc0273(testCase)
    adapter = getAdapter();
    % Without perComplex the limiting subunit alone is raised
    [model, flexEnz] = flexibilizeEnzConcs(complexToy([1 3], [1 3]), 2, [], [], adapter, false);
    verifyEqual(testCase, flexEnz.uniprotIDs, {'X1'})
    verifyEqual(testCase, usageUb(model, 'X2'), 3)
    % With perComplex both are raised by one factor
    [model, flexEnz] = flexibilizeEnzConcs(complexToy([1 3], [1 3]), 'expGrowth', 2, ...
        'modelAdapter', adapter, 'verbose', false, 'perComplex', true);
    verifyEqual(testCase, sort(flexEnz.uniprotIDs), {'X1';'X2'})
    [~, order] = ismember({'X1';'X2'}, flexEnz.uniprotIDs);
    verifyEqual(testCase, flexEnz.flexConcs(order), [1.995; 5.985], 'RelTol', 1e-6)
    verifyEqual(testCase, usageUb(model, 'X1') / 1, usageUb(model, 'X2') / 3, 'RelTol', 1e-9)
    sol = solveLP(model);
    verifyEqual(testCase, sol.f, 2, 'RelTol', 5e-3)
    % a cap is never lowered below its measurement
    verifyGreaterThanOrEqual(testCase, usageUb(model, 'X1'), 1)
    verifyGreaterThanOrEqual(testCase, usageUb(model, 'X2'), 3)
    verifyEqual(testCase, flexEnz.ratioIncr, flexEnz.flexConcs ./ flexEnz.oldConcs)
end

function testFlexibilizePerComplexKeepsLevelledComplexLevelled_tc0274(testCase)
    % Both subunits capped at 1 bind together: no single subunit shows a gain
    adapter = getAdapter();
    [model, flexEnz] = flexibilizeEnzConcs(complexToy([1 1], [1 1]), 'expGrowth', 2, ...
        'modelAdapter', adapter, 'verbose', false, 'perComplex', true);
    verifyEqual(testCase, usageUb(model, 'X1'), usageUb(model, 'X2'), 'RelTol', 1e-9)
    verifyEqual(testCase, sort(flexEnz.uniprotIDs), {'X1';'X2'})
    sol = solveLP(model);
    verifyEqual(testCase, sol.f, 2, 'RelTol', 5e-3)
end

function testFlexibilizePerComplexFalseIsTheDefault_tc0275(testCase)
    adapter = getAdapter();
    [m0, f0] = flexibilizeEnzConcs(complexToy([1 3], [1 3]), 'expGrowth', 2, 'modelAdapter', adapter, 'verbose', false);
    [m1, f1] = flexibilizeEnzConcs(complexToy([1 3], [1 3]), 'expGrowth', 2, 'modelAdapter', adapter, 'verbose', false, 'perComplex', false);
    verifyEqual(testCase, m1, m0)
    verifyEqual(testCase, f1, f0)
    verifyEqual(testCase, f0.uniprotIDs, {'X1'})
end

function testFlexibilizePerComplexDoesNotPickComplexForPoolPrice_tc0276(testCase)
    % A -> B by the complex X1, X2, X3 and B -> P by the single enzyme E,
    % with a protein pool that binds. At flux 1 every subunit uses 80% of a
    % 1.25 cap, E uses all of its cap of 1. Of the enzymes, E alone limits.
    adapter = getAdapter();
    [model, flexEnz] = flexibilizeEnzConcs(seriesToy(), 'expGrowth', 1.1, 'modelAdapter', adapter, ...
        'verbose', false, 'perComplex', true);
    verifyEqual(testCase, flexEnz.uniprotIDs(:), {'prot_pool';'E'})
    verifyEqual(testCase, usageUb(model, 'E'), 1.09725, 'RelTol', 1e-4)
    for e = {'X1','X2','X3'}
        verifyEqual(testCase, usageUb(model, e{1}), 1.25, 'RelTol', 1e-9)
    end
    sol = solveLP(model);
    verifyGreaterThanOrEqual(testCase, sol.f, 1.1 * 0.9975 - 1e-9)
end

% ---- helpers

function adapter = getAdapter()
    adapter = ModelAdapterManager.getAdapter(fullfile(findGECKOroot,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
end

function rows = toyRows()
    rows = {'R1_EXP_1', struct('A',1,'B',2); ...
            'R1_EXP_2', struct('A',1,'B',2,'C',1); ...  % isozyme variant: same complex
            'R2',       struct('D',1,'F',1,'S',1); ...
            'R3',       struct('E',1,'G',1,'S',1); ...  % S is shared with the complex of R2
            'R4',       struct('H',1)};                 % a single enzyme: not a complex
end

function model = ecOnly(concs, rows)
    % A model that carries nothing but the ec data: rows is {rxn, struct(enzyme, copies)}
    enzymes = {'A';'B';'C';'D';'F';'S';'E';'G';'H'};
    mw = [100; 200; 50; 100; 100; 100; 100; 100; 100];
    mat = zeros(size(rows, 1), numel(enzymes));
    for r = 1:size(rows, 1)
        names = fieldnames(rows{r, 2});
        for n = 1:numel(names)
            mat(r, strcmp(enzymes, names{n})) = rows{r, 2}.(names{n});
        end
    end
    model.ec.rxns      = rows(:, 1);
    model.ec.genes     = strcat('g_', enzymes);
    model.ec.enzymes   = enzymes;
    model.ec.mw        = mw;
    model.ec.concs     = concs(:);
    model.ec.rxnEnzMat = sparse(mat);
end

function sets = sortedComplexes(families)
    sets = cellfun(@(c) sort(c.enzymes(:))', families.complexes, 'UniformOutput', false);
    [~, order] = sort(cellfun(@(c) strjoin(c, ''), sets, 'UniformOutput', false));
    sets = sets(order)';
end

function names = ownedBy(families, model, enzyme)
    [~, e] = ismember(enzyme, model.ec.enzymes);
    c = families.complexes{families.owner(e)};
    names = c.enzymes(c.owned);
end

function ub = usageUb(model, enzyme)
    ub = model.ub(strcmp(model.rxns, ['usage_prot_' enzyme]));
end

function model = complexToy(concs, ubs)
    % A -> B is catalysed by the complex X1 + X2; growth (R4) equals the flux.
    model = buildLP({'A_e';'A_c';'B_c';'prot_pool';'prot_X1';'prot_X2'}, ...
        {'EX_A', {'A_e',-1}, -1000, 0; ...
         'TR_A', {'A_e',-1,'A_c',1}, 0, 1000; ...
         'R_AB', {'A_c',-1,'B_c',1,'prot_X1',-1,'prot_X2',-1}, 0, 1000; ...
         'R4', {'B_c',-1}, 0, 1000; ...
         'prot_pool_exchange', {'prot_pool',1}, 0, 1000; ...
         'usage_prot_X1', {'prot_pool',-1,'prot_X1',1}, 0, ubs(1); ...
         'usage_prot_X2', {'prot_pool',-1,'prot_X2',1}, 0, ubs(2)});
    model.ec.enzymes   = {'X1';'X2'};
    model.ec.rxns      = {'R_AB'};
    model.ec.genes     = {'g_X1';'g_X2'};
    model.ec.mw        = [100; 100];
    model.ec.concs     = concs(:);
    model.ec.rxnEnzMat = sparse([1 1]);
end

function model = seriesToy()
    % A -> B by the complex X1, X2, X3, B -> P by the enzyme E; the pool
    % holds exactly 4 mg.
    model = buildLP({'A';'B';'P';'prot_pool';'prot_X1';'prot_X2';'prot_X3';'prot_E'}, ...
        {'EX_A', {'A',1}, 0, 1000; ...
         'R1', {'A',-1,'B',1,'prot_X1',-1,'prot_X2',-1,'prot_X3',-1}, 0, 1000; ...
         'R2', {'B',-1,'P',1,'prot_E',-1}, 0, 1000; ...
         'R4', {'P',-1}, 0, 1000; ...
         'prot_pool_exchange', {'prot_pool',1}, 0, 4; ...
         'usage_prot_X1', {'prot_pool',-1,'prot_X1',1}, 0, 1.25; ...
         'usage_prot_X2', {'prot_pool',-1,'prot_X2',1}, 0, 1.25; ...
         'usage_prot_X3', {'prot_pool',-1,'prot_X3',1}, 0, 1.25; ...
         'usage_prot_E',  {'prot_pool',-1,'prot_E',1}, 0, 1});
    model.ec.enzymes   = {'X1';'X2';'X3';'E'};
    model.ec.rxns      = {'R1';'R2'};
    model.ec.genes     = {'g1';'g2';'g3';'gE'};
    model.ec.mw        = [100; 100; 100; 100];
    model.ec.concs     = [1.25; 1.25; 1.25; 1];
    model.ec.rxnEnzMat = sparse([1 1 1 0; 0 0 0 1]);
end

function model = buildLP(mets, rxnSpec)
    % rxnSpec rows: id, {metabolite, coefficient, ...}, lb, ub. The reaction
    % named R4 is the objective, which is the biomass reaction of the adapter.
    nR = size(rxnSpec, 1);
    S = zeros(numel(mets), nR);
    for r = 1:nR
        entries = rxnSpec{r, 2};
        for n = 1:2:numel(entries)
            S(strcmp(mets, entries{n}), r) = entries{n+1};
        end
    end
    model.id      = 'toy';
    model.mets    = mets;
    model.metNames = mets;
    model.rxns    = rxnSpec(:, 1);
    model.rxnNames = rxnSpec(:, 1);
    model.S       = sparse(S);
    model.lb      = cell2mat(rxnSpec(:, 3));
    model.ub      = cell2mat(rxnSpec(:, 4));
    model.c       = double(strcmp(model.rxns, 'R4'));
    model.b       = zeros(numel(mets), 1);
    model.rev     = double(model.lb < 0);
end
