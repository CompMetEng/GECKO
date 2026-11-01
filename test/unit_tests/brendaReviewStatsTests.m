%run this test case with the command
%results = runtests('brendaReviewStatsTests.m')

function tests = brendaReviewStatsTests
    tests = functiontests(localfunctions);
end

function testEcStatsFromBrendaMedianMaximumAndCase_tc0400(testCase)
    % Per EC number: the rows and maximum of the max view, the median of the
    % median view, the values of both views, and the specific-activity rows.
    % The EC keys keep their case.
    kcatMax = {{'1.1.1.1';'1.1.1.1';'3.5.1.n3';'3.5.1.N5'}, {'a';'b';'c';'d'}, {'o';'o';'o';'o'}, [10;4;40;7]};
    kcatMed = {{'1.1.1.1';'1.1.1.1';'1.1.1.1';'3.5.1.n3'}, {'a';'b';'e';'c'}, {'o';'o';'o';'o'}, [2;3;6;4]};
    sa = {{'4.1.1.1';'1.1.1.1'}, {'o';'o'}, [1;2], [1;1]};
    stats = ecStatsFromBrenda(kcatMax, sa, kcatMed);
    verifyEqual(testCase, sort(stats.keys()), {'1.1.1.1','3.5.1.N5','3.5.1.n3','4.1.1.1'})
    s = stats('1.1.1.1');
    verifyEqual(testCase, [s.nKcat s.kcatMax s.kcatMedian s.nSa], [2 10 3 1])
    verifyEqual(testCase, s.values(:)', [2 3 4 6 10])
    s = stats('3.5.1.n3');
    verifyEqual(testCase, [s.nKcat s.kcatMax s.kcatMedian s.nSa], [1 40 4 0])
    verifyEqual(testCase, s.values(:)', [4 40])
    s = stats('3.5.1.N5');
    verifyEqual(testCase, [s.nKcat s.kcatMax s.nSa], [1 7 0])
    verifyTrue(testCase, isnan(s.kcatMedian))
    s = stats('4.1.1.1');
    verifyEqual(testCase, [s.nKcat s.nSa], [0 1])
    verifyTrue(testCase, isnan(s.kcatMax) && isnan(s.kcatMedian))
    verifyEmpty(testCase, s.values)
end

function testEcStatsFromBrendaWithoutMedianViewWarns_tc0401(testCase)
    kcatMax = {{'1.1.1.1';'1.1.1.1'}, {'a';'b'}, {'o';'o'}, [10;4]};
    verifyWarning(testCase, @() ecStatsFromBrenda(kcatMax, {[],[],[],[]}), 'ecStatsFromBrenda:noMedianView')
    warning('off', 'ecStatsFromBrenda:noMedianView')
    cleanup = onCleanup(@() warning('on', 'ecStatsFromBrenda:noMedianView')); %#ok<NASGU>
    stats = ecStatsFromBrenda(kcatMax, {[],[],[],[]});
    s = stats('1.1.1.1');
    verifyTrue(testCase, isnan(s.kcatMedian))
    verifyEqual(testCase, [s.nKcat s.kcatMax], [2 10])
    verifyEqual(testCase, s.values(:)', [4 10])
end

function testReviewFlagsEcMaximumWithLoaderStatistics_tc0402(testCase)
    % A kcat equal to the maximum of its EC, ten times the median, is flagged
    % 'ec-maximum', also for an EC with a lowercase letter.
    kcatMax = {{'1.1.1.1';'3.5.1.n3'}, {'a';'c'}, {'o';'o'}, [10;40]};
    kcatMed = {{'1.1.1.1';'3.5.1.n3'}, {'a';'c'}, {'o';'o'}, [2;4]};
    stats = ecStatsFromBrenda(kcatMax, {[],[],[],[]}, kcatMed);
    f = reviewKcatAssignment({'R_A';'R_F'}, [10;40], {'brenda';'brenda'}, {'1.1.1.1';'3.5.1.n3'}, 'ecStats', stats);
    verifyEqual(testCase, f.rxnId', {'R_A','R_F'})
    verifyTrue(testCase, all(contains(f.checks, 'ec-maximum')))
end
