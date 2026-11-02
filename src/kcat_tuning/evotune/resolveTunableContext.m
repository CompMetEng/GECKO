function [tunableIdx,ecRxnIdsTunable,kcat0,groups,sigma0log,tieMap] = resolveTunableContext(ecModel,params,tunableMask)
% resolveTunableContext  The tunable kcat subset plus everything derived
% from it: source-group trust tiers and isozyme tie groups.
%
% Shared by screenKcatLeverage and cmaesKcatTuning so a screen and the
% tuning run it feeds are always looking at exactly the same parameters.
%
% A kcat belongs to the group whose params.kcatSources entry equals its
% ec.source, or the part of ec.source before ' (' (case-insensitive), so
% 'brenda (wc=0, origin=1)' is in the group of 'brenda' and
% 'combined (sd=0.10, n=2)' in the group of 'combined', with that group's
% sigma0logSource.
%
% Parameters
% ----------
% ecModel : struct
%     ecModel with a populated ec.kcat.
% params : struct
%     modelAdapter.params.evotune (kcatSources, sigma0logSource,
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
% A source with a detail, such as 'combined (sd=0.10, n=2)', is also matched by
% its name before ' ('.
sourceNames     = strtrim(regexprep(sources, ' \(.*$', ''));

groups    = repmat({'unlabelled'}, numel(tunableIdx), 1);
sigma0log = params.sigma0logDefault * ones(numel(tunableIdx),1);
for i = 1:numel(params.kcatSources)
    idx = strcmpi(sources, params.kcatSources{i}) | strcmpi(sourceNames, params.kcatSources{i});
    groups(idx)    = params.kcatSources(i);
    sigma0log(idx) = params.sigma0logSource(i);
end

if isfield(params,'tieIsozymes') && params.tieIsozymes
    tieMap = isozymeTieMap(ecRxnIdsTunable, kcat0, 'sources', groups);
else
    tieMap = (1:numel(kcat0))';
end
end
