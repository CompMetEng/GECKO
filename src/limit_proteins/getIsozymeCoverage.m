function [coverage, zeroable] = getIsozymeCoverage(model)
% getIsozymeCoverage  Report isozyme groups with partial proteomics coverage.
%
% The reactions in model.ec.rxns that share an identifier once the suffix
% _EXP_<n> is removed are the isozymes of one reaction of the conventional
% GEM. An entry is measured when at least one of its enzymes has a
% concentration above 0 in model.ec.concs, and unmeasured when none of its
% enzymes has a concentration (all NaN). A group has partial coverage when
% it has both a measured and an unmeasured entry. An unmeasured entry next to
% a measured one can be a protein below the detection limit as well as a
% protein that is absent, so the report makes no statement on which. The kcat
% of an entry is not considered.
%
% Parameters
% ----------
% model : struct
%     a full ecModel (with ecModel.ec structure), not
%     GECKO light, with concentrations in model.ec.concs (see fillEnzConcs).
%
% Returns
% -------
% coverage : table
%     one row per group with partial coverage, sorted by reaction, with the
%     variables:
%
%     - rxn : the reaction identifier without the _EXP_<n> suffix.
%     - nEntries : number of entries of the group in model.ec.rxns.
%     - nMeasured : number of measured entries.
%     - nUnmeasured : number of unmeasured entries.
%     - measuredEnzymes : cell array with the enzymes of the measured entries
%       that have a concentration above 0.
%     - unmeasuredEnzymes : cell array with the enzymes of the unmeasured
%       entries.
%     - sharedEnzymes : cell array with the unmeasured enzymes that also
%       catalyse an entry outside the unmeasured entries of groups with
%       partial coverage. zeroUnmeasuredIsozymes leaves them unchanged.
% zeroable : cell array
%     the enzymes that zeroUnmeasuredIsozymes acts on: the enzymes in
%     unmeasuredEnzymes that are not in sharedEnzymes of any group.
%
% Examples
% --------
%     coverage = getIsozymeCoverage(model);
%     [coverage, zeroable] = getIsozymeCoverage(model);
%
% See also
% --------
% zeroUnmeasuredIsozymes, getKcatAcrossIsozymes, fillEnzConcs

if isfield(model, 'ec') && isfield(model.ec, 'geckoLight') && model.ec.geckoLight
    error('getIsozymeCoverage:geckoLight', ...
        'This function does not work on GECKO light models.')
end
if ~isfield(model, 'ec') || ~isfield(model.ec, 'rxnEnzMat') || isempty(model.ec.rxnEnzMat)
    error('getIsozymeCoverage:noEnzymeMatrix', ...
        'The model has no ec.rxnEnzMat. Run makeEcModel first.')
end

enzymes = model.ec.enzymes(:);
nEnz    = numel(enzymes);
concs   = nan(nEnz, 1);
if isfield(model.ec, 'concs')
    nConcs  = min(numel(model.ec.concs), nEnz);
    concs(1:nConcs) = model.ec.concs(1:nConcs);
end
matrix  = logical(model.ec.rxnEnzMat);

entryMeasured   = (matrix * double(concs > 0)) > 0;
entryUnmeasured = (matrix * double(~isnan(concs))) == 0;
baseRxn         = regexprep(model.ec.rxns, '_EXP_\d+', '');
[groups, ~, groupIdx] = unique(baseRxn);

partial = false(numel(groups), 1);
for g = 1:numel(groups)
    inGroup = groupIdx == g;
    partial(g) = any(entryMeasured(inGroup)) && any(entryUnmeasured(inGroup));
end

% An unmeasured entry of a group with partial coverage
inPartial = entryUnmeasured & partial(groupIdx);
% An enzyme is zeroable when it is unmeasured and every entry it catalyses
% is such an entry
catalyses = matrix';
nEntries  = full(sum(catalyses, 2));
nInPartial = full(catalyses * double(inPartial));
zeroableMask = isnan(concs) & nEntries > 0 & nInPartial == nEntries;
zeroable = enzymes(zeroableMask);

rows = find(partial);
n = numel(rows);
rxn = groups(rows);
nEntriesGroup = zeros(n, 1);
nMeasured = zeros(n, 1);
nUnmeasured = zeros(n, 1);
measuredEnzymes = cell(n, 1);
unmeasuredEnzymes = cell(n, 1);
sharedEnzymes = cell(n, 1);
for k = 1:n
    inGroup = groupIdx == rows(k);
    nEntriesGroup(k) = sum(inGroup);
    nMeasured(k)     = sum(inGroup & entryMeasured);
    nUnmeasured(k)   = sum(inGroup & entryUnmeasured);
    measuredMask   = any(matrix(inGroup & entryMeasured, :), 1)' & concs > 0;
    unmeasuredMask = any(matrix(inGroup & entryUnmeasured, :), 1)' & isnan(concs);
    measuredEnzymes{k}   = enzymes(measuredMask);
    unmeasuredEnzymes{k} = enzymes(unmeasuredMask);
    sharedEnzymes{k}     = enzymes(unmeasuredMask & ~zeroableMask);
end
coverage = table(rxn, nEntriesGroup, nMeasured, nUnmeasured, measuredEnzymes, ...
    unmeasuredEnzymes, sharedEnzymes, ...
    'VariableNames', {'rxn','nEntries','nMeasured','nUnmeasured', ...
    'measuredEnzymes','unmeasuredEnzymes','sharedEnzymes'});
end
