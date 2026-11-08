%run this test case with the command
%results = runtests('kcatSourceFormatTests.m')

function tests = kcatSourceFormatTests
    tests = functiontests(localfunctions);
end

function testFormatKcatSourceExamples_tc0701(testCase)
    verifyEqual(testCase, formatKcatSource('custom'), 'custom')
    verifyEqual(testCase, formatKcatSource('standard'), 'standard')
    verifyEqual(testCase, formatKcatSource('brenda', struct('wc', 0, 'origin', 1)), 'brenda (wc=0, origin=1)')
    verifyEqual(testCase, formatKcatSource('okp_catapro', struct('similarity', 42)), 'okp_catapro (similarity=42)')
    verifyEqual(testCase, formatKcatSource('combined', struct('sd', sprintf('%.2f', 0.61), 'n', 6)), 'combined (sd=0.61, n=6)')
    verifyEqual(testCase, formatKcatSource('homologs', struct('identity', 0.953)), 'homologs (identity=0.953)')
    verifyEqual(testCase, formatKcatSource('brenda', struct('wc', 0, 'origin', 2), {'cmaes'}, 6.55), ...
        'brenda (wc=0, origin=2) | tuned=cmaes (prior=6.55)')
    verifyEqual(testCase, formatKcatSource('okp_catapro', struct('similarity', 42), {'sensitivity', 'cmaes'}, 5.48), ...
        'okp_catapro (similarity=42) | tuned=sensitivity+cmaes (prior=5.48)')
end

function testFormatKcatSourceNormalisesOriginAndKeys_tc0702(testCase)
    verifyEqual(testCase, formatKcatSource('OKP-CataPro'), 'okp_catapro')
    verifyEqual(testCase, formatKcatSource('  DLKcat '), 'dlkcat')
    verifyEqual(testCase, formatKcatSource('a - b'), 'a_b')
    verifyEqual(testCase, formatKcatSource('_x_'), 'x')
    verifyEqual(testCase, formatKcatSource('Brenda', struct('WC', 0, 'Origin', 1)), 'brenda (wc=0, origin=1)')
    % an empty origin is the empty text, unless there is more to write
    verifyEqual(testCase, formatKcatSource(''), '')
    verifyEqual(testCase, formatKcatSource('', struct('wc', 0)), 'unknown (wc=0)')
    verifyEqual(testCase, formatKcatSource('', [], {'cmaes'}, 3), 'unknown | tuned=cmaes (prior=3)')
end

function testFormatKcatSourceDetailValuesAndStatus_tc0703(testCase)
    % numbers use %.6g, text is verbatim, NaN and [] leave the pair out and
    % the brackets with it
    verifyEqual(testCase, formatKcatSource('x', struct('v', 1234567)), 'x (v=1.23457e+06)')
    verifyEqual(testCase, formatKcatSource('x', struct('v', 0.1 + 0.2)), 'x (v=0.3)')
    verifyEqual(testCase, formatKcatSource('x', struct('v', '0.60')), 'x (v=0.60)')
    verifyEqual(testCase, formatKcatSource('x', struct('a', NaN)), 'x')
    verifyEqual(testCase, formatKcatSource('x', struct('a', [])), 'x')
    verifyEqual(testCase, formatKcatSource('x', struct('a', NaN, 'b', 2)), 'x (b=2)')
    % the prior is written when there is a tuning status
    verifyEqual(testCase, formatKcatSource('x', [], {}, 5), 'x')
    verifyEqual(testCase, formatKcatSource('x', [], {'cmaes'}, NaN), 'x | tuned=cmaes')
    verifyEqual(testCase, formatKcatSource('x', [], {'cmaes'}, []), 'x | tuned=cmaes')
    verifyEqual(testCase, formatKcatSource('x', [], 'cmaes', 2), 'x | tuned=cmaes (prior=2)')
end

function testParseKcatSourceExamples_tc0704(testCase)
    info = parseKcatSource('brenda (wc=0, origin=1)');
    verifyEqual(testCase, info.origin, 'brenda')
    verifyEqual(testCase, info.detail, struct('wc', 0, 'origin', 1))
    verifyEqual(testCase, info.tuned, {})
    verifyTrue(testCase, isnan(info.prior))

    info = parseKcatSource('combined (sd=0.61, n=6)');
    verifyEqual(testCase, info.origin, 'combined')
    verifyEqual(testCase, info.detail, struct('sd', 0.61, 'n', 6))

    info = parseKcatSource('okp_catapro (similarity=42) | tuned=sensitivity+cmaes (prior=5.48)');
    verifyEqual(testCase, info.origin, 'okp_catapro')
    verifyEqual(testCase, info.detail, struct('similarity', 42))
    verifyEqual(testCase, info.tuned, {'sensitivity', 'cmaes'})
    verifyEqual(testCase, info.prior, 5.48)

    info = parseKcatSource('standard | tuned=cmaes');
    verifyEqual(testCase, info.origin, 'standard')
    verifyEqual(testCase, info.tuned, {'cmaes'})
    verifyTrue(testCase, isnan(info.prior))

    % a value that is not a number stays text; whitespace around is stripped
    info = parseKcatSource('  x ( a = abc ,b=2,junk, c=NaN )  ');
    verifyEqual(testCase, info.origin, 'x')
    verifyEqual(testCase, info.detail, struct('a', 'abc', 'b', 2, 'c', 'NaN'))
end

function testParseKcatSourceLegacyAndNonText_tc0705(testCase)
    verifyEqual(testCase, parseKcatSource('brenda').origin, 'brenda')
    % an unclosed bracket is part of the origin, with no detail
    info = parseKcatSource('brenda (wc=1');
    verifyEqual(testCase, info.origin, 'brenda (wc=1')
    verifyEqual(testCase, info.detail, struct())
    % the origin keeps the case it was written in
    info = parseKcatSource('OKP-CataPro');
    verifyEqual(testCase, info.origin, 'OKP-CataPro')
    verifyEqual(testCase, info.detail, struct())
    verifyEqual(testCase, info.tuned, {})
    verifyEqual(testCase, parseKcatSource('sensitivityTuning').origin, 'sensitivityTuning')
    verifyEqual(testCase, parseKcatSource('other (x)').origin, 'other')
    verifyEqual(testCase, parseKcatSource('other (x)').detail, struct())
    verifyEqual(testCase, parseKcatSource('no closing (a=1').origin, 'no closing (a=1')
    for bad = {'', NaN, [], 3}
        info = parseKcatSource(bad{1});
        verifyEqual(testCase, info.origin, '')
        verifyEqual(testCase, info.detail, struct())
        verifyEqual(testCase, info.tuned, {})
        verifyTrue(testCase, isnan(info.prior))
    end
    verifyEqual(testCase, parseKcatSource("custom").origin, 'custom')
    verifyEqual(testCase, parseKcatSource(string(missing)).origin, '')
    % a cell array gives a struct array of its shape
    info = parseKcatSource({'brenda (wc=0)', ''; 'custom', NaN});
    verifyEqual(testCase, size(info), [2 2])
    verifyEqual(testCase, {info.origin}, {'brenda', 'custom', '', ''})
    verifyEqual(testCase, info(1).detail, struct('wc', 0))
    verifyEqual(testCase, size(parseKcatSource({})), [0 0])
end

function testFormatOfParseGivesTheText_tc0706(testCase)
    texts = {'brenda (wc=0, origin=1)', 'custom', 'okp_catapro (similarity=42)', 'combined (sd=0.61, n=6)', ...
        'brenda (wc=0, origin=2) | tuned=cmaes (prior=6.55)', ...
        'okp_catapro (similarity=42) | tuned=sensitivity+cmaes (prior=5.48)', ...
        'standard', 'manual', 'isozymes', 'homologs (identity=0.953)', 'standard | tuned=cmaes'};
    for i = 1:numel(texts)
        info = parseKcatSource(texts{i});
        verifyEqual(testCase, formatKcatSource(info.origin, info.detail, info.tuned, info.prior), texts{i})
    end
    % numbers are written again with %.6g
    info = parseKcatSource('combined (sd=0.60, n=6)');
    verifyEqual(testCase, formatKcatSource(info.origin, info.detail, info.tuned, info.prior), 'combined (sd=0.6, n=6)')
end

function testAssignKcatValuesWritesTheSourceText_tc0707(testCase)
    % The source text of a match: the label in lower case snake_case, the
    % wildcard level with the origin when both are known, the level alone
    % when the origin is not, and no detail when the level is not known.
    model = ecFixture();
    rxns = model.ec.rxns(1:5);
    kcatList = struct('source', 'x', 'rxns', {rxns}, 'kcats', [5; 6; 7; 8; 9], ...
        'kcatSource', {{'BRENDA'; 'DLKcat'; 'OKP-CataPro'; 'brenda'; 'Brenda'}}, ...
        'wildcardLvl', [0; NaN; NaN; 2; NaN], 'origin', [1; NaN; 3; NaN; 4]);
    model = assignKcatValues(model, kcatList);
    verifyEqual(testCase, model.ec.source(1:5), {'brenda (wc=0, origin=1)'; 'dlkcat'; 'okp_catapro'; 'brenda (wc=2)'; 'brenda'})
end

function testSetKcatForReactionsWritesManualWithoutStatus_tc0708(testCase)
    model = ecFixture();
    model.ec.source(:) = {'brenda (wc=0, origin=1) | tuned=cmaes (prior=2)'};
    j = strcmp(model.ec.rxns, 'R3');
    out = setKcatForReactions(model, 'R3', 7);
    verifyEqual(testCase, out.ec.source{j}, 'manual')
    verifyEqual(testCase, out.ec.kcat(j), 7)
    verifyEqual(testCase, out.ec.source(~j), model.ec.source(~j))
end

function testGetStandardKcatKeepsATunedStandardRow_tc0709(testCase)
    % R3 holds a standard kcat on a real enzyme: it is set to zero and its
    % source cleared. R5 holds a tuned standard kcat: it is kept as it is.
    % The row added for the reaction without GPR (R1) reads 'standard'.
    [model, adapter] = ecFixture();
    model.ec.kcat(:) = 0;
    model.ec.kcat(strcmp(model.ec.rxns, 'R3')) = 50;
    model.ec.kcat(strcmp(model.ec.rxns, 'R5')) = 999;
    model.ec.source(:) = {'manual'};
    model.ec.source{strcmp(model.ec.rxns, 'R3')} = 'standard';
    model.ec.source{strcmp(model.ec.rxns, 'R5')} = 'standard | tuned=cmaes (prior=3)';
    evalc('out = getStandardKcat(model, ''modelAdapter'', adapter, ''threshold'', 1, ''fillZeroKcat'', false);');
    verifyEqual(testCase, out.ec.source{strcmp(out.ec.rxns, 'R3')}, '')
    verifyEqual(testCase, out.ec.kcat(strcmp(out.ec.rxns, 'R3')), 0)
    verifyEqual(testCase, out.ec.source{strcmp(out.ec.rxns, 'R5')}, 'standard | tuned=cmaes (prior=3)')
    verifyEqual(testCase, out.ec.kcat(strcmp(out.ec.rxns, 'R5')), 999)
    verifyEqual(testCase, out.ec.source{strcmp(out.ec.rxns, 'R1')}, 'standard')
end

function testRemoveStandardKcatKeepsATunedStandardRow_tc0710(testCase)
    % The standard rows are cleared whatever their letter case; the one with
    % a tuning status and the others are kept.
    model = ecFixture();
    model.ec.kcat(:) = 10;
    model.ec.source(:) = {'brenda'};
    model.ec.source{1} = 'standard';
    model.ec.source{2} = 'Standard';
    model.ec.source{3} = 'standard | tuned=sensitivity (prior=5)';
    out = removeStandardKcat(model);
    verifyEqual(testCase, out.ec.source(1:2), {''; ''})
    verifyEqual(testCase, out.ec.kcat(1:2), [0; 0])
    verifyEqual(testCase, out.ec.source(3:end), model.ec.source(3:end))
    verifyEqual(testCase, out.ec.kcat(3:end), model.ec.kcat(3:end))
end

function testReviewKcatAssignmentRecognisesOriginsOfLabels_tc0711(testCase)
    % Every kcat is outside the magnitude bounds. The standard rows, in any
    % letter case and with a tuning status, are left out; the custom row is
    % reported.
    findings = reviewKcatAssignment({'A'; 'B'; 'C'}, [1e-5; 1e-5; 1e-5], ...
        {'Custom (x=1)'; 'Standard'; 'standard | tuned=cmaes (prior=1)'}, {'1.1.1.1'; '1.1.1.2'; '1.1.1.3'});
    verifyEqual(testCase, findings.rxnId, {'A'})
    verifyEqual(testCase, findings.checks, {'magnitude'})
    findings = reviewKcatAssignment({'A'; 'B'}, [1e-5; 1e-5], {'standard (x=1)'; 'brenda'}, {'1.1.1.1'; '1.1.1.2'}, ...
        'include', {'standard-fallback'});
    verifyEqual(testCase, findings.rxnId, {'A'; 'B'})
    verifyEqual(testCase, findings.checks, {'magnitude,standard-fallback'; 'magnitude'})
end

function testResolveTunableContextGroupsATunedSourceByItsOrigin_tc0712(testCase)
    model = ecFixture();
    model.ec.kcat(1:4) = 10;
    model.ec.source(1:4) = {'brenda | tuned=cmaes (prior=2)'; 'Combined (sd=0.10, n=2) | tuned=sensitivity (prior=1)'; ...
        'DLKcat'; 'other | tuned=cmaes (prior=2)'};
    params = struct('kcatSources', {{'brenda', 'combined', 'dlkcat'}}, 'sigma0logSource', [0.2; 0.3; 0.4], ...
        'sigma0logDefault', 0.5, 'tieIsozymes', false);
    mask = false(numel(model.ec.rxns), 1);
    mask(1:4) = true;
    [~, ~, ~, groups, sigma0log] = resolveTunableContext(model, params, mask);
    verifyEqual(testCase, groups, {'brenda'; 'combined'; 'dlkcat'; 'unlabelled'})
    verifyEqual(testCase, sigma0log, [0.2; 0.3; 0.4; 0.5])
end

function testSensitivityTuningExtendsTheStatusOnce_tc0713(testCase)
    % Every kcat the tuning changes gets the origin it had, the method
    % sensitivity and the kcat before the first tuning as prior. A second
    % tuning adds no second method and keeps the prior. The notes are left
    % as they were.
    [ecModel, tuneAdapter] = tuningModel();
    ecModel.ec.source(:) = {'brenda (wc=0, origin=1)'};
    cleanup = shadowSolveLP([1 1 50]); %#ok<NASGU>
    first = sensitivityTuning(ecModel, 'desiredGrowthRate', 50, ...
        'modelAdapter', tuneAdapter, 'protToIgnore', {'P4'}, 'verbose', false);
    changed = find(first.ec.kcat ~= ecModel.ec.kcat);
    verifyTrue(testCase, any(strcmp(first.ec.rxns(changed), 'R5')))
    iR5 = strcmp(ecModel.ec.rxns, 'R5');
    verifyEqual(testCase, ecModel.ec.kcat(iR5), 200)
    verifyEqual(testCase, first.ec.source{iR5}, 'brenda (wc=0, origin=1) | tuned=sensitivity (prior=200)')
    for i = changed'
        verifyEqual(testCase, first.ec.source{i}, ...
            sprintf('brenda (wc=0, origin=1) | tuned=sensitivity (prior=%.6g)', ecModel.ec.kcat(i)))
    end
    verifyEqual(testCase, first.ec.source(first.ec.kcat == ecModel.ec.kcat), ecModel.ec.source(first.ec.kcat == ecModel.ec.kcat))
    verifyEqual(testCase, first.ec.notes, ecModel.ec.notes)

    clear cleanup
    cleanup2 = shadowSolveLP([1 1 50]); %#ok<NASGU>
    second = sensitivityTuning(first, 'desiredGrowthRate', 50, ...
        'modelAdapter', tuneAdapter, 'protToIgnore', {'P4'}, 'verbose', false);
    changedFirst = find(first.ec.kcat ~= ecModel.ec.kcat);
    changedSecond = find(second.ec.kcat ~= first.ec.kcat);
    both = intersect(changedFirst(:), changedSecond(:));
    verifyNotEmpty(testCase, both)
    for k = 1:numel(both)
        i = both(k);
        verifyEqual(testCase, second.ec.source{i}, ...
            sprintf('brenda (wc=0, origin=1) | tuned=sensitivity (prior=%.6g)', ecModel.ec.kcat(i)))
    end
    onlySecond = setdiff(changedSecond(:), changedFirst(:));
    for k = 1:numel(onlySecond)
        i = onlySecond(k);
        verifyEqual(testCase, second.ec.source{i}, ...
            sprintf('brenda (wc=0, origin=1) | tuned=sensitivity (prior=%.6g)', first.ec.kcat(i)))
    end
end

function testSensitivityTuningWritesUnknownForAnEmptyOrigin_tc0714(testCase)
    [ecModel, tuneAdapter] = tuningModel();
    ecModel.ec.source(:) = {''};
    cleanup = shadowSolveLP([1 1 50]); %#ok<NASGU>
    out = sensitivityTuning(ecModel, 'desiredGrowthRate', 50, ...
        'modelAdapter', tuneAdapter, 'protToIgnore', {'P4'}, 'verbose', false);
    changed = find(out.ec.kcat ~= ecModel.ec.kcat);
    verifyNotEmpty(testCase, changed)
    for i = changed'
        verifyEqual(testCase, out.ec.source{i}, sprintf('unknown | tuned=sensitivity (prior=%.6g)', ecModel.ec.kcat(i)))
    end
end

function [ecModel, adapter] = ecFixture()
% The ecTestGEM ecModel with the enzyme data from the model.
adapter = getAdapter();
ecModel = makeEcModel(getGeckoTestModel(), false, adapter);
ecModel = getECfromGEM(ecModel);
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

function cleanup = shadowSolveLP(values)
% Puts a solveLP on the path that calls RAVEN's and then reports the growth
% of solve number n as values(min(n, end)). The path is restored when the
% returned object is cleared.
setappdata(0, 'gecko_test_solvelp_real', @solveLP);
setappdata(0, 'gecko_test_solvelp_count', 0);
setappdata(0, 'gecko_test_solvelp_values', values);
stubDir = tempname;
mkdir(stubDir);
fid = fopen(fullfile(stubDir, 'solveLP.m'), 'w');
fprintf(fid, '%s\n', ...
    'function [res,hs] = solveLP(varargin)', ...
    'n = getappdata(0,''gecko_test_solvelp_count'') + 1;', ...
    'setappdata(0,''gecko_test_solvelp_count'',n);', ...
    'real = getappdata(0,''gecko_test_solvelp_real'');', ...
    'values = getappdata(0,''gecko_test_solvelp_values'');', ...
    '[res,hs] = real(varargin{:});', ...
    'res.f = values(min(n,numel(values)));', ...
    'end');
fclose(fid);
addpath(stubDir, '-begin');
cleanup = onCleanup(@() removeShadow(stubDir));
end

function removeShadow(stubDir)
rmpath(stubDir);
rmdir(stubDir, 's');
for name = {'real', 'values'}
    if isappdata(0, ['gecko_test_solvelp_' name{1}])
        rmappdata(0, ['gecko_test_solvelp_' name{1}]);
    end
end
if isappdata(0, 'gecko_test_solvelp_count')
    rmappdata(0, 'gecko_test_solvelp_count');
end
end

function adapter = getAdapter()
geckoPath = findGECKOroot;
adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
end

function testParseKcatSourceReadsANumberOnlyInPlainDecimalForm_tc0715(testCase)
    info = parseKcatSource('x (a=1_0, b=inf, c=1e3, d=.5, e=-2., f=0x10)');
    verifyEqual(testCase, info.detail, struct('a', '1_0', 'b', 'inf', 'c', 1000, 'd', 0.5, 'e', -2, 'f', '0x10'))
end

function testParseKcatSourceIgnoresAPairWithoutALetterKey_tc0716(testCase)
    info = parseKcatSource('x (1a=1, a b=2, ok=3, _c=4)');
    verifyEqual(testCase, info.detail, struct('ok', 3))
end

function testParseKcatSourceNeedsOneSpaceBeforeThePrior_tc0717(testCase)
    info = parseKcatSource('x | tuned=cmaes(prior=3)');
    verifyEqual(testCase, info.tuned, {'cmaes(prior=3)'})
    verifyTrue(testCase, isnan(info.prior))
end

function testFormatKcatSourceLeavesOutUnusablePairsAndPriors_tc0718(testCase)
    verifyEqual(testCase, formatKcatSource('x', struct('ok', 3, 'a1', 2)), 'x (ok=3, a1=2)')
    verifyEqual(testCase, formatKcatSource('x', struct('a', Inf, 'b', 1)), 'x (b=1)')
    verifyEqual(testCase, formatKcatSource('x', struct('a', -Inf)), 'x')
    verifyEqual(testCase, formatKcatSource('x', [], {'cmaes'}, Inf), 'x | tuned=cmaes')
    verifyEqual(testCase, formatKcatSource('x', [], 'cmaes', 2), 'x | tuned=cmaes (prior=2)')
end

function testParseKcatSourceNeedsExactlyOneSpaceBeforeThePrior_tc0719(testCase)
    info = parseKcatSource('x | tuned=cmaes  (prior=3)');
    verifyEqual(testCase, info.tuned, {'cmaes  (prior=3)'})
    verifyTrue(testCase, isnan(info.prior))
end

function testParseKcatSourceReadsAClosingBracketInThePriorAsAMethod_tc0720(testCase)
    info = parseKcatSource('x | tuned=cmaes (prior=1))');
    verifyEqual(testCase, info.tuned, {'cmaes (prior=1))'})
    verifyTrue(testCase, isnan(info.prior))
end

function testParseKcatSourceTrimsTheHeadBeforeTheTuningStatus_tc0721(testCase)
    info = parseKcatSource('x (a=1)  | tuned=cmaes');
    verifyEqual(testCase, info.origin, 'x')
    verifyEqual(testCase, info.detail, struct('a', 1))
    verifyEqual(testCase, info.tuned, {'cmaes'})
end

function testParseKcatSourceTakesTheFirstPriorMarker_tc0722(testCase)
    info = parseKcatSource('x | tuned=cmaes (prior=3) (prior=4)');
    verifyEqual(testCase, info.tuned, {'cmaes (prior=3) (prior=4)'})
    verifyTrue(testCase, isnan(info.prior))
end

