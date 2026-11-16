%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kcatAggregation = resolveKcatAggregation(kcatAggregation, modelAdapter)
% resolveKcatAggregation  Settle which BRENDA aggregate to read.
%
% The BRENDA snapshot records both a maximum and a median per
% (EC, substrate, organism) triple. Which one a run uses can be given per
% call. Without one, it is [kcat_sources.brenda] aggregation of the kcat
% assignment configuration (kcatAssignmentConfig): the preset's, which is
% 'median' for the default preset brenda_ec_tree and 'max' for
% gecko_classic, or the one that model_adapter.toml sets. A configuration
% in which no stage sets an aggregation for brenda gives 'max'.
%
% Input Arguments
% ---------------
% kcatAggregation : char
%     the per-call choice, or empty to fall back to the adapter.
% modelAdapter : ModelAdapter
%     a loaded model adapter.
%
% Returns
% -------
% kcatAggregation : char
%     'max' or 'median'.
%
% See also
% --------
% loadBRENDAdata, fuzzyKcatMatching

if isempty(kcatAggregation)
    config = kcatAssignmentConfig(modelAdapter.getParameters());
    kcatAggregation = 'max';
    if isfield(config.sources, 'brenda') && isfield(config.sources.brenda, 'aggregation')
        kcatAggregation = config.sources.brenda.aggregation;
    end
end

if ~(ischar(kcatAggregation) || isstring(kcatAggregation)) || ...
        ~any(strcmpi(kcatAggregation,{'max','median'}))
    error('resolveKcatAggregation:invalidKcatAggregation', ...
        'aggregation must be max or median')
end
kcatAggregation = lower(char(kcatAggregation));
end
