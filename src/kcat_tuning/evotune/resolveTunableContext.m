function [tunableIdx,ecRxnIdsTunable,kcat0,groups,sigma0log,tieMap] = resolveTunableContext(ecModel,params,tunableMask)
% resolveTunableContext  The tunable kcat subset plus everything derived
% from it: source-group trust tiers and isozyme tie groups.
%
% Shared by screenKcatLeverage and cmaesKcatTuning so a screen and the
% tuning run it feeds are always looking at exactly the same parameters.
%
% A kcat belongs to the group of the first params.kcatSources entry that
% matches its source, and is labelled with the name of that group
% (params.kcatGroups, one per entry). Without groups the error is raised:
% No evotune source groups: set [evotune.source_groups] in model_adapter.toml. An entry matches the origin of the kcat's ec.source
% or its whole ec.source (case-insensitive, see parseKcatSource), so
% 'brenda (wc=0, origin=1)' is in the group of 'brenda',
% 'combined (sd=0.10, n=2)' in the group of 'combined' and
% 'brenda (wc=0, origin=1) | tuned=cmaes (prior=2)' in the group of 'brenda',
% with that group's sigma0logSource. An entry that ends in * is a prefix
% pattern on the origin ('okp_*' covers okp_eitlem_kinetics and
% okp_kinform_h). Exact entries are tried first, over all entries in order;
% then the patterns, the longest prefix first, and in the order of
% params.kcatSources among equal prefixes. Kcats that match no entry get
% sigma0logDefault; when more than 20 % of the tunable kcats match none,
% a warning (resolveTunableContext:unlabelled) lists up to five origins,
% the most common first and in alphabetical order among equal counts.
%
% Parameters
% ----------
% ecModel : struct
%     ecModel with a populated ec.kcat.
% params : struct
%     modelAdapter.params.evotune (kcatSources, kcatGroups, sigma0logSource,
%     sigma0logDefault, tieIsozymes).
% tunableMask : logical
%     same length as ecModel.ec.rxns; [] to consider every kcat > 0.
%
% Returns
% -------
% tunableIdx : double
%     indices into ecModel.ec.rxns/ec.kcat of the tunable subset.
% ecRxnIdsTunable : cell
% kcat0 : double
% groups : cell
%     source-group name per tunable row ('unlabelled' if unmatched).
% sigma0log : double
%     per-row prior std dev in log-space.
% tieMap : double
%     isozymeTieMap(ecRxnIdsTunable,kcat0,'sources',groups) if
%     params.tieIsozymes, otherwise every row is its own free parameter.
%
% See also
% --------
% screenKcatLeverage, cmaesKcatTuning, isozymeTieMap
isTunable = ecModel.ec.kcat > 0;
if ~isempty(tunableMask)
    isTunable = isTunable & tunableMask(:);
end
tunableIdx = find(isTunable);
if isempty(tunableIdx)
    if isempty(tunableMask)
        error('No tunable kcats: ecModel.ec.kcat is all <= 0.')
    else
        error('No tunable kcats: ecModel.ec.kcat is all <= 0, or tunableMask excludes every one with a kcat.')
    end
end

ecRxnIdsTunable = ecModel.ec.rxns(tunableIdx);
kcat0           = ecModel.ec.kcat(tunableIdx);
sources         = ecModel.ec.source(tunableIdx);
sourceInfo      = parseKcatSource(sources);
sourceNames     = {sourceInfo.origin}';

groups    = repmat({'unlabelled'}, numel(tunableIdx), 1);
sigma0log = params.sigma0logDefault * ones(numel(tunableIdx),1);
assigned  = false(numel(tunableIdx),1);
entries   = params.kcatSources(:)';
if ~isfield(params, 'kcatGroups') || isempty(params.kcatGroups)
    error('GECKO:resolveTunableContext', ...
        'No evotune source groups: set [evotune.source_groups] in model_adapter.toml')
end
entryGroup = params.kcatGroups(:)';
isPattern = endsWith(entries, '*');
patternLength = cellfun(@numel, entries);
[~, order] = sort(patternLength(isPattern), 'descend'); % stable: ties keep the group order
patternIdx = find(isPattern);
for i = [find(~isPattern), patternIdx(order)]
    if isPattern(i)
        idx = ~assigned & startsWith(sourceNames, entries{i}(1:end-1), 'IgnoreCase', true);
    else
        idx = ~assigned & (strcmpi(sourceNames, entries{i}) | strcmpi(sources, entries{i}));
    end
    assigned(idx)  = true;
    groups(idx)    = entryGroup(i);
    sigma0log(idx) = params.sigma0logSource(i);
end
if mean(~assigned) > 0.2
    [names, ~, k] = unique(sourceNames(~assigned));
    [~, byCount] = sort(accumarray(k, 1), 'descend');
    names = names(byCount(1:min(5, end)));
    warning('resolveTunableContext:unlabelled', ...
        '%d of %d tunable kcats match no evotune source group and get sigma0_log_default (%g); origins: %s', ...
        sum(~assigned), numel(assigned), params.sigma0logDefault, strjoin(names, ', '))
end

if isfield(params,'tieIsozymes') && params.tieIsozymes
    tieMap = isozymeTieMap(ecRxnIdsTunable, kcat0, 'sources', groups);
else
    tieMap = (1:numel(kcat0))';
end
end
