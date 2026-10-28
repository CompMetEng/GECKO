%run this test case with the command
%results = runtests('parityDecisionTests.m')

function tests = parityDecisionTests
    tests = functiontests(localfunctions);
end

function testMergeKcatsDefaultTopOriginLimitIs4_tc0240(testCase)
    % A BRENDA match of origin 5 (specific activity, correct organism) ranks
    % after DLKcat by default and above it with topOriginLimit 6.
    [fuzzy, dlkcat] = originFiveLists();
    tiers = {'database_top', 'dlkcat', 'database_bottom'};
    merged = mergeKcats({fuzzy, dlkcat}, tiers);
    verifyEqual(testCase, merged.kcats(strcmp(merged.rxns, 'R1')), 50)
    merged = mergeKcats({fuzzy, dlkcat}, tiers, 'topOriginLimit', 4);
    verifyEqual(testCase, merged.kcats(strcmp(merged.rxns, 'R1')), 50)
    merged = mergeKcats({fuzzy, dlkcat}, tiers, 'topOriginLimit', 6);
    verifyEqual(testCase, merged.kcats(strcmp(merged.rxns, 'R1')), 5)
end

function testMergeDLKcatAndFuzzyKcatsDefaultTopOriginLimitIs4_tc0241(testCase)
    [fuzzy, dlkcat] = originFiveLists();
    merged = mergeDLKcatAndFuzzyKcats(dlkcat, fuzzy);
    verifyEqual(testCase, merged.kcats(strcmp(merged.rxns, 'R1')), 50)
    merged = mergeDLKcatAndFuzzyKcats(dlkcat, fuzzy, 'topOriginLimit', 6);
    verifyEqual(testCase, merged.kcats(strcmp(merged.rxns, 'R1')), 5)
end

function testFindECInDBSingleGeneDropsWildcardParent_tc0242(testCase)
    % One gene whose protein carries '1.1.1.4 1.1.1.- 1.1.1.303' keeps the
    % two specific codes and loses the wildcard of the same prefix.
    geneHashMap = containers.Map({'G1'}, {1});
    geneIndex = {1};
    DBecNum = {'1.1.1.4 1.1.1.- 1.1.1.303'};
    DBMW = 10000;
    EC = findECInDB({'G1'}, DBecNum, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '1.1.1.4;1.1.1.303')
    % A wildcard without a specific code of its prefix stays
    DBecNum = {'1.1.1.- 2.2.2.2'};
    EC = findECInDB({'G1'}, DBecNum, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '1.1.1.-;2.2.2.2')
end

function testFindECInDBSingleGeneShortCodesAndTrailingSpace_tc0243(testCase)
    % A code with fewer than four levels is neither the wildcard parent nor
    % the child of a four-level code, and a trailing space leaves no empty code.
    geneHashMap = containers.Map({'G1'}, {1});
    geneIndex = {1};
    DBMW = 10000;
    EC = findECInDB({'G1'}, {'1.1.1 1.1.1.4'}, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '1.1.1;1.1.1.4')
    EC = findECInDB({'G1'}, {'1.1.1.4 '}, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '1.1.1.4')
end

function testFindECInDBComplexDropsWildcardParent_tc0244(testCase)
    % The shared EC set of an enzyme complex loses a wildcard code next to a
    % specific code of its prefix, as the set of a single gene does.
    DBMW = [10000; 20000; 30000];
    geneHashMap = containers.Map({'G1','G2','G3'}, {1,2,3});
    geneIndex = {1; 2; 3};
    EC = findECInDB({'G1','G2'}, {'1.1.1.4 1.1.1.-'; '1.1.1.4 1.1.1.-'; ''}, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '1.1.1.4')
    % A wildcard that only the union carries is dropped too
    EC = findECInDB({'G1','G2','G3'}, {'1.1.1.4'; '1.1.1.- 2.2.2.2'; '5.5.5.5'}, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '1.1.1.4;2.2.2.2;5.5.5.5')
    % A wildcard shared by every gene stays, as nothing more specific exists
    EC = findECInDB({'G1','G2'}, {'1.1.1.-'; '1.1.1.-'; ''}, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '1.1.1.-')
end

function testFindECInDBGeneWithoutECIsLeftOut_tc0245(testCase)
    % A gene without an EC number does not change the codes of the others:
    % the complex is treated as the genes that have codes, in either order.
    DBMW = [10000; 20000];
    geneHashMap = containers.Map({'G1','G2'}, {1,2});
    geneIndex = {1; 2};
    EC = findECInDB({'G1','G2'}, {'1.1.1.4 1.1.1.-'; ''}, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '1.1.1.4')
    EC = findECInDB({'G1','G2'}, {''; '1.1.1.4 1.1.1.-'}, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '1.1.1.4')
    % No gene has a code
    EC = findECInDB({'G1','G2'}, {''; ''}, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '')
    % A gene missing from the database is left out as well
    geneHashMap = containers.Map({'G2'}, {2});
    geneIndex = {[]; 2};
    EC = findECInDB({'G1','G2'}, {''; '1.1.1.4 1.1.1.-'}, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '1.1.1.4')
end

function [fuzzy, dlkcat] = originFiveLists()
    fuzzy.rxns        = {'R1'};
    fuzzy.kcats       = 5;
    fuzzy.genes       = {{}};
    fuzzy.substrates  = {{}};
    fuzzy.eccodes     = {'1.1.1.1'};
    fuzzy.wildcardLvl = 0;
    fuzzy.origin      = 5;
    fuzzy.source      = 'brenda';
    dlkcat.rxns       = {'R1'};
    dlkcat.kcats      = 50;
    dlkcat.genes      = {'G1'};
    dlkcat.substrates = {'m1'};
    dlkcat.source     = 'DLKcat';
end
