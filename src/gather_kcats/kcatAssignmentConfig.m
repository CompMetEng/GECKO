function config = kcatAssignmentConfig(params)
% kcatAssignmentConfig  The kcat assignment settings of a model adapter.
%
% Reads [kcat_assignment] and [kcat_sources.*] of model_adapter.toml (or
% params.kcat_assignment and params.kcat_sources of an adapter class), and
% expands a preset into its stages, as collectKcats and updateKcatSources
% use them. Without a preset or stages the preset is brenda_ec_tree.
% geckopy's kcat_assignment_config gives the same settings.
%
% Parameters
% ----------
% params : struct
%     the parameters of a model adapter (modelAdapter.params).
%
% Returns
% -------
% config : struct
%     - preset: the preset name, or '' when stages are given;
%     - isozymes, standard: logical;
%     - stages: struct array with action ('select' or 'combine'), sources
%       (cell, as given), files (the source files they read, with okp
%       replaced by the files of its methods), criteria, topOriginLimit,
%       bottomOriginLimit and wildcardLimit;
%     - okpMethods: the OKP methods of the okp source group ([kcat_sources.okp]
%       methods; default all methods with combination parameters);
%     - sources: per source name its options (from [kcat_sources] and the
%       preset).
%
% Examples
% --------
%     config = kcatAssignmentConfig(modelAdapter.params);
%     config.stages(1).sources   % the order of the sources
%
% See also
% --------
% collectKcats, assignKcats, updateKcatSources

a = struct();
if isfield(params, 'kcat_assignment') && isstruct(params.kcat_assignment)
    a = params.kcat_assignment;
end
config.preset = '';
config.isozymes = getOr(a, 'isozymes', true);
config.standard = getOr(a, 'standard', true);
config.sources = struct();
hasPreset = isfield(a, 'preset') && ~isempty(a.preset);
hasStages = isfield(a, 'stage') && ~isempty(a.stage);
if hasPreset && hasStages
    error('kcatAssignment:config', ['[kcat_assignment] gives both a preset and stages; ' ...
        'give one of the two.'])
end
if hasStages
    config.stages = arrayfun(@readStage, a.stage(:)');
else
    preset = 'brenda_ec_tree';
    if hasPreset
        preset = char(a.preset);
    end
    [config.stages, config.sources] = presetStages(preset);
    config.preset = preset;
end
if isfield(params, 'kcat_sources') && isstruct(params.kcat_sources)
    names = fieldnames(params.kcat_sources);
    for k = 1:numel(names)
        config.sources.(names{k}) = params.kcat_sources.(names{k});
    end
end
% The okp group: the OKP methods of [kcat_sources.okp] methods.
known = {'CataPro','CatPred','CatRange','DLKcat','EITLEM','KinForm-H','KinForm-L', ...
    'OmniESI','OmniESI-O2DENet','TurNuP','UniKP'};
config.okpMethods = known(~strcmp(known, 'CatRange')); % CatRange has no combination parameters
if isfield(config.sources, 'okp') && isfield(config.sources.okp, 'methods') && ~isempty(config.sources.okp.methods)
    config.okpMethods = cellstr(config.sources.okp.methods);
    config.okpMethods = config.okpMethods(:)';
end
unknown = setdiff(config.okpMethods, known);
if ~isempty(unknown)
    error('kcatAssignment:config', 'Unknown OKP methods in [kcat_sources.okp] methods: %s; methods: %s.', ...
        strjoin(unknown, ', '), strjoin(known, ', '))
end
okpFiles = cellfun(@(m) ['okp_' lower(regexprep(regexprep(m, '[^A-Za-z0-9]+', '_'), '^_+|_+$', ''))], ...
    config.okpMethods, 'UniformOutput', false);
for s = 1:numel(config.stages)
    files = {};
    for t = config.stages(s).sources
        if strcmp(t{1}, 'okp')
            files = [files, okpFiles]; %#ok<AGROW>
        else
            files{end+1} = t{1}; %#ok<AGROW>
        end
    end
    [~, first] = unique(files, 'stable');
    config.stages(s).files = files(sort(first));
end
end

function [stages, sources] = presetStages(preset)
% The presets of kcat assignment; geckopy defines the same.
stage = @(action, srcs, top) struct('action', action, 'sources', {srcs}, 'criteria', 'max', ...
    'topOriginLimit', top, 'bottomOriginLimit', 6, 'wildcardLimit', 3);
order = {'custom', 'database_top', 'dlkcat', 'database_bottom'};
own = {'custom', 'database_top'}; % with topOriginLimit 1: the organism's own BRENDA value
sources = struct();
switch preset
    case 'brenda_ec_tree'
        stages = [stage('select', [own, {'ec_tree'}], 1), stage('combine', {'okp'}, 4)];
        sources.brenda = struct('aggregation', 'median');
        sources.ec_tree = struct('k', 5);
    case 'ec_tree_homologs'
        stages = [stage('select', own, 1), stage('combine', {'ec_tree', 'homologs'}, 4), ...
            stage('combine', {'okp'}, 4)];
        sources.brenda = struct('aggregation', 'median');
        sources.ec_tree = struct('k', 5);
    case 'ec_tree_okp'
        stages = [stage('select', own, 1), stage('combine', {'ec_tree', 'homologs', 'okp'}, 4), ...
            stage('combine', {'okp'}, 4)];
        sources.brenda = struct('aggregation', 'median');
        sources.ec_tree = struct('k', 5);
    case 'okp'
        stages = [stage('select', {'custom'}, 4), stage('combine', {'okp'}, 4)];
    case 'gecko_median'
        stages = stage('select', order, 4);
        sources.brenda = struct('aggregation', 'median');
    case 'gecko_classic'
        stages = stage('select', order, 6);
        sources.brenda = struct('aggregation', 'max');
    otherwise
        error('kcatAssignment:config', 'Unknown kcat assignment preset "%s". Presets: %s.', ...
            preset, 'brenda_ec_tree, ec_tree_homologs, ec_tree_okp, okp, gecko_median, gecko_classic')
end
end

function stage = readStage(s)
stage.action = getOr(s, 'action', 'select');
if ~any(strcmp(stage.action, {'select', 'combine'}))
    error('kcatAssignment:config', ['Unknown stage action "%s"; this GECKO version knows ' ...
        '"select" and "combine".'], stage.action)
end
sources = getOr(s, 'sources', {});
if ischar(sources)
    sources = {sources};
end
if isempty(sources)
    error('kcatAssignment:config', 'A [[kcat_assignment.stage]] needs sources.')
end
stage.sources = sources(:)';
stage.criteria = validatestring(getOr(s, 'criteria', 'max'), {'max','min','median','mean'});
stage.topOriginLimit = getOr(s, 'top_origin_limit', 4);
stage.bottomOriginLimit = getOr(s, 'bottom_origin_limit', 6);
stage.wildcardLimit = getOr(s, 'wildcard_limit', 3);
end

function v = getOr(s, field, default)
v = default;
if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
    v = s.(field);
end
end
